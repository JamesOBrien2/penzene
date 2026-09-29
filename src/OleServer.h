#pragma once
// Windows only (#229): Penzene as an OLE embedding server. Word and PowerPoint keep a
// "Penzene Drawing Object" (the .penz pages and an EMF picture of page 1, which OLE's default
// handler draws without Penzene); double-clicking it starts penzene.exe -Embedding, which
// opens the drawing in a window, and Save puts it back. The installer registers kClsid
// (cmake/penzene.iss). Editing is in Penzene's own window, as ChemDraw does, not in place.
#include "Document.h"

#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <ole2.h>

#include <QVariant>
#include <QWindowsMimeConverter>
#include <vector>

class MainWindow;

namespace ole {
extern const CLSID kClsid;  // {F6858801-14D9-488C-B79A-8C7D07452862}, "Penzene.Drawing"
bool isRegistered();        // by the installer, so Office can open what it pastes
// The storage an embedded drawing lives in: the class, and the pages as .penz JSON. The
// caller releases it. Copy offers it as "Embed Source".
IStorage* embedSource(const std::vector<Sheet>& sheets);
// An embedded drawing, edited in window: IOleObject, IDataObject and IPersistStorage.
IUnknown* newObject(MainWindow& window);
// penzene.exe -Embedding: serves one drawing through window until it's closed or released.
int serve(MainWindow& window);
}  // namespace ole

// Copy's drawing (application/x-penzene) as "Embed Source" and "Object Descriptor", so Word and
// PowerPoint paste an object that opens in Penzene. Offered only once Penzene is registered.
struct EmbedClipboard : QWindowsMimeConverter {
    bool canConvertFromMime(const FORMATETC& format, const QMimeData* mime) const override;
    bool convertFromMime(const FORMATETC& format, const QMimeData* mime, STGMEDIUM* medium) const override;
    QList<FORMATETC> formatsForMime(const QString& type, const QMimeData*) const override;
    bool canConvertToMime(const QString&, IDataObject*) const override { return false; }
    QVariant convertToMime(const QString&, IDataObject*, QMetaType) const override { return {}; }
    QString mimeForFormat(const FORMATETC&) const override { return {}; }
};
