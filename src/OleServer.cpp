#include "OleServer.h"
#include "MainWindow.h"

#include <QApplication>
#include <QMimeData>
#include <QTimer>
#include <atomic>
#include <cstring>
#include <string>

namespace ole {

const CLSID kClsid = {0xf6858801, 0x14d9, 0x488c, {0xb7, 0x9a, 0x8c, 0x7d, 0x07, 0x45, 0x28, 0x62}};
static const wchar_t* kStream = L"Penzene";  // the pages, .penz JSON
static bool serving = false;                 // penzene.exe -Embedding: quit once the drawing is done with

bool isRegistered() {
    HKEY key = nullptr;
    if (RegOpenKeyExW(HKEY_CLASSES_ROOT, L"CLSID\\{F6858801-14D9-488C-B79A-8C7D07452862}\\LocalServer32", 0, KEY_READ,
                      &key) != ERROR_SUCCESS)
        return false;
    RegCloseKey(key);
    return true;
}

static bool writePages(IStorage* storage, const std::vector<Sheet>& sheets) {
    const QByteArray json = sheetsToJson(sheets);
    IStream* stream = nullptr;
    if (FAILED(storage->CreateStream(kStream, STGM_CREATE | STGM_WRITE | STGM_SHARE_EXCLUSIVE, 0, 0, &stream)))
        return false;
    ULONG written = 0;
    const bool ok = SUCCEEDED(stream->Write(json.constData(), ULONG(json.size()), &written)) && written == ULONG(json.size());
    stream->Release();
    return ok;
}

static std::vector<Sheet> readPages(IStorage* storage) {
    IStream* stream = nullptr;  // read whole and closed: the storage stays usable for Save
    if (FAILED(storage->OpenStream(kStream, nullptr, STGM_READ | STGM_SHARE_EXCLUSIVE, 0, &stream))) return {};
    STATSTG stat{};
    QByteArray json;
    if (SUCCEEDED(stream->Stat(&stat, STATFLAG_NONAME)) && stat.cbSize.QuadPart < (1ull << 30)) {
        json.resize(qsizetype(stat.cbSize.QuadPart));
        ULONG read = 0;
        if (FAILED(stream->Read(json.data(), ULONG(json.size()), &read))) read = 0;
        json.truncate(read);
    }
    stream->Release();
    return sheetsFromJson(json);
}

IStorage* embedSource(const std::vector<Sheet>& sheets) {
    if (sheets.empty()) return nullptr;
    ILockBytes* bytes = nullptr;
    if (FAILED(CreateILockBytesOnHGlobal(nullptr, TRUE, &bytes))) return nullptr;
    IStorage* storage = nullptr;
    const HRESULT hr = StgCreateDocfileOnILockBytes(bytes, STGM_CREATE | STGM_READWRITE | STGM_SHARE_EXCLUSIVE, 0, &storage);
    bytes->Release();  // the storage keeps it
    if (FAILED(hr)) return nullptr;
    if (FAILED(WriteClassStg(storage, kClsid)) || !writePages(storage, sheets) || FAILED(storage->Commit(STGC_DEFAULT))) {
        storage->Release();
        return nullptr;
    }
    return storage;
}

// Page 1 as Copy draws it; the caller deletes it. Null for an empty drawing.
static HENHMETAFILE picture(const MainWindow& window) {
    const QByteArray emf = window.embeddedPicture();
    return emf.isEmpty() ? nullptr : SetEnhMetaFileBits(UINT(emf.size()), reinterpret_cast<const BYTE*>(emf.constData()));
}

// The same picture as an old-style metafile, which some containers cache instead.
static HGLOBAL metafilePict(HENHMETAFILE emf) {
    ENHMETAHEADER header{};
    GetEnhMetaFileHeader(emf, sizeof header, &header);
    HDC screen = GetDC(nullptr);
    std::vector<BYTE> bits(GetWinMetaFileBits(emf, 0, nullptr, MM_ANISOTROPIC, screen));
    const UINT size = bits.empty() ? 0 : GetWinMetaFileBits(emf, UINT(bits.size()), bits.data(), MM_ANISOTROPIC, screen);
    ReleaseDC(nullptr, screen);
    HMETAFILE wmf = size ? SetMetaFileBitsEx(size, bits.data()) : nullptr;
    HGLOBAL global = wmf ? GlobalAlloc(GMEM_MOVEABLE, sizeof(METAFILEPICT)) : nullptr;
    if (!global) {
        if (wmf) DeleteMetaFile(wmf);
        return nullptr;
    }
    auto* pict = static_cast<METAFILEPICT*>(GlobalLock(global));
    *pict = {MM_ANISOTROPIC, header.rclFrame.right - header.rclFrame.left, header.rclFrame.bottom - header.rclFrame.top, wmf};
    GlobalUnlock(global);
    return global;
}

namespace {

class Object : public IOleObject, public IDataObject, public IPersistStorage {
public:
    static inline Object* live = nullptr;  // the one a -Embedding process serves

    explicit Object(MainWindow& window) : window_(window) {
        CreateOleAdviseHolder(&oleAdvise_);
        CreateDataAdviseHolder(&dataAdvise_);
        live = this;
    }

    // The window was closed (its changes already saved or discarded), or the container closed us.
    void closed() {
        if (closed_) return;
        closed_ = true;
        window_.hide();
        if (site_) site_->OnShowWindow(FALSE);  // the container stops hatching the object
        if (oleAdvise_) oleAdvise_->SendOnClose();
        if (serving) QCoreApplication::quit();
    }

    // IUnknown
    STDMETHODIMP QueryInterface(REFIID iid, void** out) override {
        if (iid == IID_IUnknown || iid == IID_IOleObject) *out = static_cast<IOleObject*>(this);
        else if (iid == IID_IDataObject) *out = static_cast<IDataObject*>(this);
        else if (iid == IID_IPersistStorage || iid == IID_IPersist) *out = static_cast<IPersistStorage*>(this);
        else return *out = nullptr, E_NOINTERFACE;
        AddRef();
        return S_OK;
    }
    STDMETHODIMP_(ULONG) AddRef() override { return ++refs_; }
    STDMETHODIMP_(ULONG) Release() override {
        const ULONG left = --refs_;
        if (left) return left;
        if (live == this) live = nullptr;
        delete this;
        if (serving) QCoreApplication::quit();  // the container is done with the drawing
        return 0;
    }

    // IOleObject: open editing only, in Penzene's window; the registry answers the rest (OLE_S_USEREG).
    STDMETHODIMP SetClientSite(IOleClientSite* site) override {
        if (site) site->AddRef();
        if (site_) site_->Release();
        site_ = site;
        return S_OK;
    }
    STDMETHODIMP GetClientSite(IOleClientSite** site) override {
        if ((*site = site_)) site_->AddRef();
        return S_OK;
    }
    STDMETHODIMP SetHostNames(LPCOLESTR app, LPCOLESTR document) override {
        window_.setEmbeddedIn(QString::fromWCharArray(document && *document ? document : app ? app : L""));
        return S_OK;
    }
    STDMETHODIMP Close(DWORD option) override {
        // Prompting or not, an edit is kept: the container asks when its document closes.
        if (option != OLECLOSE_NOSAVE && !window_.isClean()) window_.save();
        closed();
        return S_OK;
    }
    STDMETHODIMP SetMoniker(DWORD, IMoniker*) override { return E_NOTIMPL; }
    STDMETHODIMP GetMoniker(DWORD, DWORD, IMoniker**) override { return E_NOTIMPL; }
    STDMETHODIMP InitFromData(IDataObject*, BOOL, DWORD) override { return E_NOTIMPL; }
    STDMETHODIMP GetClipboardData(DWORD, IDataObject**) override { return E_NOTIMPL; }
    STDMETHODIMP DoVerb(LONG verb, LPMSG, IOleClientSite*, LONG, HWND, LPCRECT) override {
        if (verb == OLEIVERB_HIDE) {
            window_.hide();
            if (site_) site_->OnShowWindow(FALSE);
            return S_OK;
        }
        if (verb == OLEIVERB_DISCARDUNDOSTATE) return S_OK;
        window_.show();  // Edit, Open, Show, and in-place requests too: this window is the editor
        window_.raise();
        window_.activateWindow();
        if (site_) site_->ShowObject(), site_->OnShowWindow(TRUE);
        return S_OK;
    }
    STDMETHODIMP EnumVerbs(IEnumOLEVERB**) override { return OLE_S_USEREG; }
    STDMETHODIMP Update() override { return S_OK; }
    STDMETHODIMP IsUpToDate() override { return S_OK; }
    STDMETHODIMP GetUserClassID(CLSID* clsid) override { return *clsid = kClsid, S_OK; }
    STDMETHODIMP GetUserType(DWORD, LPOLESTR*) override { return OLE_S_USEREG; }
    STDMETHODIMP SetExtent(DWORD, SIZEL*) override { return E_FAIL; }  // sized by its drawing
    STDMETHODIMP GetExtent(DWORD aspect, SIZEL* size) override {
        if (aspect != DVASPECT_CONTENT) return E_INVALIDARG;
        HENHMETAFILE emf = picture(window_);
        if (!emf) return E_FAIL;
        ENHMETAHEADER header{};
        GetEnhMetaFileHeader(emf, sizeof header, &header);
        DeleteEnhMetaFile(emf);
        *size = {header.rclFrame.right - header.rclFrame.left, header.rclFrame.bottom - header.rclFrame.top};  // 0.01 mm
        return S_OK;
    }
    STDMETHODIMP Advise(IAdviseSink* sink, DWORD* connection) override { return oleAdvise_->Advise(sink, connection); }
    STDMETHODIMP Unadvise(DWORD connection) override { return oleAdvise_->Unadvise(connection); }
    STDMETHODIMP EnumAdvise(IEnumSTATDATA** out) override { return oleAdvise_->EnumAdvise(out); }
    STDMETHODIMP GetMiscStatus(DWORD, DWORD*) override { return OLE_S_USEREG; }
    STDMETHODIMP SetColorScheme(LOGPALETTE*) override { return E_NOTIMPL; }

    // IDataObject: the picture, for the container's cache.
    STDMETHODIMP GetData(FORMATETC* format, STGMEDIUM* medium) override {
        if (const HRESULT hr = QueryGetData(format); hr != S_OK) return hr;
        HENHMETAFILE emf = picture(window_);
        if (!emf) return E_FAIL;
        medium->pUnkForRelease = nullptr;
        if (format->cfFormat == CF_ENHMETAFILE) {
            medium->tymed = TYMED_ENHMF;
            medium->hEnhMetaFile = emf;
            return S_OK;
        }
        medium->tymed = TYMED_MFPICT;
        medium->hMetaFilePict = metafilePict(emf);
        DeleteEnhMetaFile(emf);
        return medium->hMetaFilePict ? S_OK : E_FAIL;
    }
    STDMETHODIMP GetDataHere(FORMATETC*, STGMEDIUM*) override { return E_NOTIMPL; }
    STDMETHODIMP QueryGetData(FORMATETC* format) override {
        if (format->dwAspect != DVASPECT_CONTENT) return DV_E_DVASPECT;
        const bool emf = format->cfFormat == CF_ENHMETAFILE && (format->tymed & TYMED_ENHMF);
        const bool wmf = format->cfFormat == CF_METAFILEPICT && (format->tymed & TYMED_MFPICT);
        return emf || wmf ? S_OK : DV_E_FORMATETC;
    }
    STDMETHODIMP GetCanonicalFormatEtc(FORMATETC*, FORMATETC* out) override {
        out->ptd = nullptr;
        return DATA_S_SAMEFORMATETC;
    }
    STDMETHODIMP SetData(FORMATETC*, STGMEDIUM*, BOOL) override { return E_NOTIMPL; }
    STDMETHODIMP EnumFormatEtc(DWORD, IEnumFORMATETC**) override { return OLE_S_USEREG; }
    STDMETHODIMP DAdvise(FORMATETC* format, DWORD flags, IAdviseSink* sink, DWORD* connection) override {
        return dataAdvise_->Advise(static_cast<IDataObject*>(this), format, flags, sink, connection);
    }
    STDMETHODIMP DUnadvise(DWORD connection) override { return dataAdvise_->Unadvise(connection); }
    STDMETHODIMP EnumDAdvise(IEnumSTATDATA** out) override { return dataAdvise_->EnumAdvise(out); }

    // IPersistStorage: the pages, as one .penz stream; read whole on Load, so no storage is kept.
    STDMETHODIMP GetClassID(CLSID* clsid) override { return *clsid = kClsid, S_OK; }
    STDMETHODIMP IsDirty() override { return window_.isClean() ? S_FALSE : S_OK; }
    STDMETHODIMP InitNew(IStorage*) override {  // Insert → Object: an empty drawing
        edit({});
        return S_OK;
    }
    STDMETHODIMP Load(IStorage* storage) override {
        const std::vector<Sheet> sheets = readPages(storage);
        if (sheets.empty()) return STG_E_READFAULT;
        edit(sheets);
        return S_OK;
    }
    STDMETHODIMP Save(IStorage* storage, BOOL) override {
        return writePages(storage, window_.sheets()) ? S_OK : STG_E_WRITEFAULT;
    }
    STDMETHODIMP SaveCompleted(IStorage*) override {
        if (oleAdvise_) oleAdvise_->SendOnSave();
        return S_OK;
    }
    STDMETHODIMP HandsOffStorage() override { return S_OK; }

private:
    ~Object() {
        if (site_) site_->Release();
        if (oleAdvise_) oleAdvise_->Release();
        if (dataAdvise_) dataAdvise_->Release();
    }

    // Save in the window: a new picture for the container's cache, then the container saves the object.
    void edit(const std::vector<Sheet>& sheets) {
        window_.editEmbedded(sheets, [this] {
            if (!site_) return false;
            if (dataAdvise_) dataAdvise_->SendOnDataChange(static_cast<IDataObject*>(this), 0, 0);
            return SUCCEEDED(site_->SaveObject());
        });
    }

    MainWindow& window_;
    std::atomic<ULONG> refs_ = 1;
    IOleClientSite* site_ = nullptr;
    IOleAdviseHolder* oleAdvise_ = nullptr;
    IDataAdviseHolder* dataAdvise_ = nullptr;
    bool closed_ = false;
};

// Single use: each embedded drawing gets its own Penzene process, and window.
class Factory : public IClassFactory {
public:
    explicit Factory(MainWindow& window) : window_(window) {}
    STDMETHODIMP QueryInterface(REFIID iid, void** out) override {
        if (iid != IID_IUnknown && iid != IID_IClassFactory) return *out = nullptr, E_NOINTERFACE;
        *out = static_cast<IClassFactory*>(this);
        return S_OK;
    }
    STDMETHODIMP_(ULONG) AddRef() override { return 2; }  // lives as long as serve()
    STDMETHODIMP_(ULONG) Release() override { return 1; }
    STDMETHODIMP CreateInstance(IUnknown* outer, REFIID iid, void** out) override {
        *out = nullptr;
        if (outer) return CLASS_E_NOAGGREGATION;
        if (Object::live) return E_FAIL;  // one drawing per window
        IUnknown* object = newObject(window_);
        const HRESULT hr = object->QueryInterface(iid, out);
        object->Release();
        return hr;
    }
    STDMETHODIMP LockServer(BOOL) override { return S_OK; }

private:
    MainWindow& window_;
};

}  // namespace

IUnknown* newObject(MainWindow& window) { return static_cast<IOleObject*>(new Object(window)); }

int serve(MainWindow& window) {
    serving = true;
    QApplication::setQuitOnLastWindowClosed(false);  // quit once the container knows it's closed
    QObject::connect(qApp, &QGuiApplication::lastWindowClosed, [] {
        if (Object::live) Object::live->closed();
        else QCoreApplication::quit();
    });
    OleInitialize(nullptr);
    Factory factory(window);
    DWORD cookie = 0;
    if (FAILED(CoRegisterClassObject(kClsid, &factory, CLSCTX_LOCAL_SERVER, REGCLS_SINGLEUSE, &cookie))) return 1;
    QTimer::singleShot(120000, [] {  // started by hand, or the container gave up: don't linger hidden
        if (!Object::live) QCoreApplication::quit();
    });
    const int result = QApplication::exec();
    if (Object::live) CoDisconnectObject(static_cast<IOleObject*>(Object::live), 0);
    CoRevokeClassObject(cookie);
    OleUninitialize();
    return result;
}

}  // namespace ole

static const char* kPenzMime = "application/x-penzene";  // MainWindow's copy of the drawing

static CLIPFORMAT clipboardFormat(const wchar_t* name) { return CLIPFORMAT(RegisterClipboardFormatW(name)); }

bool EmbedClipboard::canConvertFromMime(const FORMATETC& format, const QMimeData* mime) const {
    if (!mime->hasFormat(kPenzMime) || !ole::isRegistered()) return false;
    return (format.cfFormat == clipboardFormat(L"Embed Source") && (format.tymed & TYMED_ISTORAGE)) ||
           (format.cfFormat == clipboardFormat(L"Object Descriptor") && (format.tymed & TYMED_HGLOBAL));
}

bool EmbedClipboard::convertFromMime(const FORMATETC& format, const QMimeData* mime, STGMEDIUM* medium) const {
    if (!canConvertFromMime(format, mime)) return false;
    medium->pUnkForRelease = nullptr;
    if (format.cfFormat == clipboardFormat(L"Embed Source")) {
        medium->pstg = ole::embedSource(sheetsFromJson(mime->data(kPenzMime)));
        medium->tymed = TYMED_ISTORAGE;
        return medium->pstg;
    }
    // Object Descriptor: the class and its name, which Paste Special lists ("Penzene Drawing Object").
    const std::wstring type = L"Penzene Drawing", source = L"Penzene";
    const DWORD typeBytes = DWORD((type.size() + 1) * sizeof(wchar_t)), sourceBytes = DWORD((source.size() + 1) * sizeof(wchar_t));
    const DWORD size = sizeof(OBJECTDESCRIPTOR) + typeBytes + sourceBytes;
    HGLOBAL global = GlobalAlloc(GMEM_MOVEABLE | GMEM_ZEROINIT, size);
    if (!global) return false;
    auto* d = static_cast<OBJECTDESCRIPTOR*>(GlobalLock(global));
    d->cbSize = size;
    d->clsid = ole::kClsid;
    d->dwDrawAspect = DVASPECT_CONTENT;
    d->dwFullUserTypeName = sizeof(OBJECTDESCRIPTOR);
    d->dwSrcOfCopy = sizeof(OBJECTDESCRIPTOR) + typeBytes;
    std::memcpy(reinterpret_cast<char*>(d) + d->dwFullUserTypeName, type.c_str(), typeBytes);
    std::memcpy(reinterpret_cast<char*>(d) + d->dwSrcOfCopy, source.c_str(), sourceBytes);
    GlobalUnlock(global);
    medium->tymed = TYMED_HGLOBAL;
    medium->hGlobal = global;
    return true;
}

QList<FORMATETC> EmbedClipboard::formatsForMime(const QString& type, const QMimeData*) const {
    if (type != kPenzMime || !ole::isRegistered()) return {};
    return {FORMATETC{clipboardFormat(L"Embed Source"), nullptr, DVASPECT_CONTENT, -1, TYMED_ISTORAGE},
            FORMATETC{clipboardFormat(L"Object Descriptor"), nullptr, DVASPECT_CONTENT, -1, TYMED_HGLOBAL}};
}
