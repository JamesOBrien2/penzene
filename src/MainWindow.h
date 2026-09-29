#pragma once
#include <QMainWindow>
#include <functional>
#include <string>
#include <vector>

#ifdef Q_OS_MACOS
#include <QUtiMimeConverter>
// ChemDraw's pasteboard type (binary CDX) as chemical/x-cdx, com.adobe.pdf as application/pdf,
// and public.png as image/png (the bytes unchanged, so an exported drawing stays in them).
struct ChemDrawPasteboard : QUtiMimeConverter {
    ChemDrawPasteboard();
    QString mimeForUti(const QString& uti) const override;
    QString utiForMime(const QString& mime) const override;
    QVariant convertToMime(const QString&, const QList<QByteArray>& data, const QString&) const override;
    QList<QByteArray> convertFromMime(const QString&, const QVariant& data, const QString&) const override;
};
#endif
#ifdef Q_OS_WIN
#include <QVariant>
#include <QWindowsMimeConverter>
// The copied drawing's Enhanced Metafile (image/x-emf) as CF_ENHMETAFILE, which Office pastes as a vector picture.
struct EmfClipboard : QWindowsMimeConverter {
    bool canConvertFromMime(const FORMATETC& format, const QMimeData* mime) const override;
    bool convertFromMime(const FORMATETC& format, const QMimeData* mime, STGMEDIUM* medium) const override;
    QList<FORMATETC> formatsForMime(const QString& type, const QMimeData*) const override;
    bool canConvertToMime(const QString&, IDataObject*) const override { return false; }
    QVariant convertToMime(const QString&, IDataObject*, QMetaType) const override { return {}; }
    QString mimeForFormat(const FORMATETC&) const override { return {}; }
};
#endif

class Canvas;
class QLabel;
class QPrinter;
class QUndoStack;
#include "Document.h"

// The drawing at its export size, centred on the page; shrunk to fit if it's bigger.
bool printDocument(QPrinter& printer, const Document& doc);

class QColor;
class QPixmap;
class QStyle;
// The platform style, with message box icons that show in every theme: a dark glyph (macOS's
// question mark is black) is drawn in the palette's text colour (#418). Install it before the palette.
QStyle* themedStyle();
QPixmap inkGlyph(const QPixmap& icon, const QColor& ink);  // `icon` in `ink` if it's dark; a coloured or light one unchanged

// RDKit's formula as HTML: counts subscripted, a trailing charge ("-2") superscripted as "2−".
QString formulaHtml(const std::string& formula);

class MainWindow : public QMainWindow {
    Q_OBJECT
public:
    MainWindow();
    ~MainWindow() override;
    bool openFile(const QString& path);
    // Crash recovery: offer the autosaved document left behind by a crash, if any.
    void offerRecovery();
    // A newer release on GitHub? By hand (Help menu) it always answers; quietly (the weekly
    // check, off unless turned on in Preferences) it only speaks up if there is one.
    void checkForUpdates(bool quietly);
    void maybeCheckForUpdates();  // at startup: weekly, if turned on
    // The first start of a newer version shows its CHANGELOG.md section (not on a fresh install).
    void maybeShowWhatsNew();
    void showWhatsNew();
    void autosave();  // writes unsaved changes to autosavePath() (every minute)
    static QString autosavePath();
    QStringList recentFiles() const;
    void showPreferences();  // theme, default style, export resolution and background
    // Interface languages: those with a translation built in (:/i18n/penzene_<code>.qm), and
    // installing one (empty: the system's language; "en": none, the source strings).
    static QStringList languages();
    static bool installTranslations(const QString& language);
    QWidget* checkStructure();  // lists problems; clicking one selects its atoms
    std::vector<Sheet> sheets() const;  // every page, the one on the canvas as it is now
    bool isClean() const;               // no unsaved change on any page
    bool save();
#ifdef Q_OS_WIN
    // A drawing embedded in Word or PowerPoint (#229, OleServer.cpp): its pages replace the
    // window's, and Save calls save, which puts them back in the document, instead of a file.
    void editEmbedded(const std::vector<Sheet>& sheets, std::function<bool()> save);
    void setEmbeddedIn(const QString& document);  // named in the title
    QByteArray embeddedPicture() const;  // page 1 as an EMF, rendered as Copy renders it
#endif

protected:
    void closeEvent(QCloseEvent* e) override;
    void resizeEvent(QResizeEvent* e) override;
    void dragEnterEvent(QDragEnterEvent* e) override;  // a file from Finder or Explorer
    void dropEvent(QDropEvent* e) override;            // opens, as File → Open does
    bool eventFilter(QObject* watched, QEvent* e) override;

private:
    void buildTools();
    void buildMenus();
    void buildWelcome();
    void paintExamples();
    bool saveTo(const QString& path, bool v3000 = false);
    bool saveAs();
    bool maybeSave();
    // With the preference on: lists doc's problems (Check Structure's) and asks whether to go on.
    bool confirmStructure(const Document& doc, const QString& title, const QString& proceed);
    void exportImage();
    void exportDescriptors();
    void importSmiles();
    void importSequence();
    void importName();
    void print();
    bool copy();  // false when nothing was copied (empty, or the structure warning cancelled)
    void paste();
    void updateTitle();
    void updateInfo();
    void updateProfile();  // the Properties panel, while it's visible
    void updateMassSpec();  // the Mass Spec panel, likewise
    void updateNmr();       // and the NMR panel
    void applyTheme(const QString& name);
    void remember(const QString& path);  // most recent first, at most 10

    QUndoStack* undo_;  // the page on the canvas's history
    Canvas* canvas_;
    // The file's pages, tabs along the bottom. The canvas holds pages_[page_]'s drawing (its
    // entry here is stale while it's shown); every page keeps its own undo history.
    struct PageState {
        QString name;
        Document doc;
        QUndoStack* undo;
    };
    std::vector<PageState> pages_;
    int page_ = 0;
    bool pagesEdited_ = false;  // pages added, deleted, renamed or moved since the last save
    bool penz1_ = false;  // Save As chose "Penzene 1": .penz saves as version 1 until the next New or Open
    class QUndoGroup* undoGroup_;
    class QTabBar* pageTabs_;
    void setPages(const std::vector<Sheet>& sheets);  // New, Open: every page replaced, no history
    void showPage(int i);
    void renumberPages();  // compound numbers run on from page to page (#569)
    void addPage();
    void renamePage(int i);
    void deletePage(int i);
    void moveSelectionToPage(int i);
    QString path_;
#ifdef Q_OS_WIN
    std::function<bool()> embeddedSave_;
    QString embeddedIn_;
#endif
    QLabel* info_;
    class QFrame* welcome_ = nullptr;
    class QAction* shortcutsAction_ = nullptr;
    std::vector<std::pair<class QToolButton*, Document>> welcomeExamples_;
    class QActionGroup* themeGroup_ = nullptr;
    class QDockWidget* profileDock_;
    class QDockWidget* templateDock_;
    class QTreeWidget* templates_;
    void fillTemplates();
    void insertTemplate(class QTreeWidgetItem* item);
    void saveTemplate();
    QLabel* profile_;
    QString profileText_;  // plain-text copy of the panel, for the Copy button
    class QDockWidget* massDock_;
    class QComboBox* ion_;
    class SpectrumView* spectrum_;
    class QLabel* eiIons_;  // under [M]: the EI ions worth checking
    class QDockWidget* nmrDock_;
    class QComboBox* nucleus_;
    class NmrView* nmr_;
    std::vector<std::pair<QAction*, std::function<QIcon()>>> icons_;
    std::vector<class QFrame*> flyouts_;  // the tool rail's group flyouts
};
