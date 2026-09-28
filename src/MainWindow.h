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

protected:
    void closeEvent(QCloseEvent* e) override;
    void dragEnterEvent(QDragEnterEvent* e) override;  // a file from Finder or Explorer
    void dropEvent(QDropEvent* e) override;            // opens, as File → Open does
    bool eventFilter(QObject* watched, QEvent* e) override;

private:
    void buildTools();
    void buildMenus();
    void buildWelcome();
    void paintExamples();
    bool saveTo(const QString& path, bool v3000 = false);
    bool save();
    bool saveAs();
    bool maybeSave();
    // With the preference on: lists doc's problems (Check Structure's) and asks whether to go on.
    bool confirmStructure(const Document& doc, const QString& title, const QString& proceed);
    void exportImage();
    void exportDescriptors();
    void importSmiles();
    void importName();
    void print();
    void copy();
    void paste();
    void updateTitle();
    void updateInfo();
    void updateProfile();  // the Properties panel, while it's visible
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
    class QUndoGroup* undoGroup_;
    class QTabBar* pageTabs_;
    void setPages(const std::vector<Sheet>& sheets);  // New, Open: every page replaced, no history
    std::vector<Sheet> sheets() const;
    void showPage(int i);
    void addPage();
    void renamePage(int i);
    void deletePage(int i);
    void moveSelectionToPage(int i);
    bool isClean() const;  // no unsaved change on any page
    QString path_;
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
    std::vector<std::pair<QAction*, std::function<QIcon()>>> icons_;
    std::vector<class QFrame*> flyouts_;  // the tool rail's group flyouts
};
