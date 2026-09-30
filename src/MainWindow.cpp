#include "MainWindow.h"
#include "Canvas.h"
#include "Chem.h"
#include "Edit.h"
#include "Online.h"
#include "Templates.h"
#include "WhatsNew.h"
#ifdef Q_OS_WIN
#include "OleServer.h"
#endif

#include <QAccessibleWidget>
#include <QActionGroup>
#include <QButtonGroup>
#include <QTextBrowser>
#include <QNetworkReply>
#include <QNetworkAccessManager>
#include <QDesktopServices>
#include <QDateTime>
#include <QCheckBox>
#include <QHash>
#include <QTreeWidget>
#include <QPainter>
#include <QPrintDialog>
#include <QPrinter>
#include <QApplication>
#include <QClipboard>
#include <QCloseEvent>
#include <QColorDialog>
#include <QFile>
#include <QFileOpenEvent>
#include <QLockFile>
#include <QSaveFile>
#include <QTabBar>
#include <QUndoGroup>
#include <QComboBox>
#include <QFileDialog>
#include <QSpinBox>
#include <QFormLayout>
#include <QVBoxLayout>
#include <QListWidget>
#include <QDialogButtonBox>
#include <QDialog>
#include <QDoubleSpinBox>
#include <QHBoxLayout>
#include <QPushButton>
#include <QDockWidget>
#include <QFileInfo>
#include <QInputDialog>
#include <memory>
#include <QSet>
#include <QWidgetAction>
#include <QToolButton>
#include <QMenu>
#include <QGridLayout>
#include <QFrame>
#include <QDir>
#include <QLibraryInfo>
#include <QTranslator>
#include <QSettings>
#include <QSignalBlocker>
#include <QMouseEvent>
#include <algorithm>
#include <QStandardPaths>
#include <QTimer>
#include <QProxyStyle>
#include <QStyle>
#include <QStyleHints>
#include <QtMath>
#include <functional>
#include <QLabel>
#include <QRegularExpression>
#include <QMenuBar>
#include <QMessageBox>
#include <QMimeData>
#include <QDropEvent>
#include <QStatusBar>
#include <QStackedWidget>
#include <QToolBar>
#include <QUndoStack>
#include <array>

static const char* kMolMime = "chemical/x-mdl-molfile";
static const char* kPenzMime = "application/x-penzene";  // full fidelity: arrows and text too
static const char* kWinPngMime = "application/x-qt-windows-mime;value=\"PNG\"";  // the clipboard format Office reads
#ifdef Q_OS_WIN
static const char* kEmfMime = "image/x-emf";  // offered to Office as CF_ENHMETAFILE by EmfClipboard
#endif
static const char* kDocsUrl = "https://penzene.readthedocs.io/";
static const QSize kExampleIcon(168, 84);

// The Mass Spec panel's stick spectrum: m/z along the bottom, the main peaks labelled.
class SpectrumView : public QWidget {
public:
    std::vector<chem::Peak> peaks;
    SpectrumView() { setMinimumSize(260, 200); }

protected:
    void paintEvent(QPaintEvent*) override {
        QPainter p(this);
        p.setRenderHint(QPainter::Antialiasing);
        const QColor ink = palette().color(QPalette::WindowText);
        if (peaks.empty()) {
            p.setPen(ink);
            p.drawText(rect(), Qt::AlignCenter | Qt::TextWordWrap,
                       QCoreApplication::translate("MainWindow", "Draw or select a valid structure."));
            return;
        }
        const QFontMetrics fm(font());
        const QRectF plot = QRectF(rect()).adjusted(8, fm.height() + 6, -8, -fm.height() - 8);
        const double lo = peaks.front().mz - 1, hi = peaks.back().mz + 1;
        auto x = [&](double mz) { return plot.left() + (mz - lo) / (hi - lo) * plot.width(); };
        p.setPen(QPen(ink, 1));
        p.drawLine(plot.bottomLeft(), plot.bottomRight());
        for (int m = int(std::ceil(lo)); m <= hi; ++m) p.drawLine(QPointF(x(m), plot.bottom()), QPointF(x(m), plot.bottom() + 3));
        p.setPen(QPen(palette().color(QPalette::Highlight), 2, Qt::SolidLine, Qt::FlatCap));
        for (const auto& k : peaks)
            p.drawLine(QPointF(x(k.mz), plot.bottom()), QPointF(x(k.mz), plot.bottom() - k.intensity / 100 * plot.height()));
        // Labels on the tallest stick of each nominal mass, if it's 5% or more.
        p.setPen(ink);
        for (const auto& k : peaks) {
            const bool main = k.intensity >= 5 && std::none_of(peaks.begin(), peaks.end(), [&](const chem::Peak& o) {
                                  return std::abs(o.mz - k.mz) < 0.5 && o.intensity > k.intensity;
                              });
            if (!main) continue;
            const QString text = QString::number(k.mz, 'f', 4);
            const double y = plot.bottom() - k.intensity / 100 * plot.height() - 4;
            p.drawText(QPointF(x(k.mz) - fm.horizontalAdvance(text) / 2.0, y), text);
        }
        p.drawText(QRectF(plot.left(), plot.bottom() + 4, plot.width(), fm.height() + 4), Qt::AlignRight, "m/z");
    }
};

// The NMR panel's predicted spectrum: one stick per set of equivalent atoms, as tall as the atoms (or H)
// it stands for, ppm falling left to right. The stick under the pointer, or stepped to with Left/Right,
// is labelled and lights its atoms on the canvas.
class NmrView : public QWidget {
public:
    std::vector<chem::NmrStick> sticks;
    bool proton = false;
    Document doc;  // the sticks' own: names their atoms, whatever the canvas holds by now
    std::vector<int> atoms;  // the molecule they're for
    Document molecule;       // and it alone, drawn as the legend
    std::function<void(const std::vector<int>&)> light;
    int current = -1;
    NmrView() {
        setMinimumSize(260, 200);
        setMouseTracking(true);
        setFocusPolicy(Qt::StrongFocus);
    }
    QString describe(int k) const {
        const auto& s = sticks[k];
        QStringList atoms;
        for (int a : s.atoms) atoms << QString::fromStdString(chem::symbol(doc.atoms[a].z)) + QString::number(a + 1);
        return QCoreApplication::translate("MainWindow", "δ %1%2, %3 %4: %5")
            .arg(s.weak ? "~" : "")
            .arg(s.ppm, 0, 'f', proton ? 2 : 1)
            .arg(s.count)
            .arg(proton ? "H, " + s.multiplicity() : "C", atoms.join(", "));
    }
    QString reading() const {  // for screen readers: the stick in hand, else all of them
        QStringList all;
        for (int k = 0; k < int(sticks.size()); ++k)
            if (current < 0 || k == current) all << describe(k);
        return all.join("; ");
    }
    void setCurrent(int k) {
        if (k == current) return;
        current = k;
        setAccessibleDescription(reading());
        light(k >= 0 ? sticks[k].atoms : std::vector<int>{});
        update();
    }

protected:
    QRectF plot() const {
        const QFontMetrics fm(font());
        return QRectF(rect()).adjusted(8, fm.height() + 6, -8, -fm.height() - 8);
    }
    std::pair<double, double> range() const {  // ppm at the left and right edges
        double left = proton ? 12 : 220, right = 0;
        for (const auto& s : sticks) left = std::max(left, s.ppm + 1), right = std::min(right, s.ppm - 1);
        return {left, right};
    }
    double x(double ppm) const {
        const auto [l, r] = range();
        return plot().left() + (l - ppm) / (l - r) * plot().width();
    }
    void paintEvent(QPaintEvent*) override {
        QPainter p(this);
        p.setRenderHint(QPainter::Antialiasing);
        const QColor ink = palette().color(QPalette::WindowText);
        if (sticks.empty()) {
            p.setPen(ink);
            p.drawText(rect(), Qt::AlignCenter | Qt::TextWordWrap,
                       QCoreApplication::translate("MainWindow", "Draw or select a structure with predicted shifts."));
            return;
        }
        const QFontMetrics fm(font());
        const QRectF plot = this->plot();
        const auto [l, r] = range();
        p.setPen(QPen(ink, 1));
        p.drawLine(plot.bottomLeft(), plot.bottomRight());
        const int step = proton ? 1 : 20;
        for (int v = int(std::floor(l / step)) * step; v >= r; v -= step) {
            p.drawLine(QPointF(x(v), plot.bottom()), QPointF(x(v), plot.bottom() + 3));
            const QString t = QString::number(v);
            if (v % (2 * step) == 0) p.drawText(QPointF(x(v) - fm.horizontalAdvance(t) / 2.0, plot.bottom() + 4 + fm.ascent()), t);
        }
        int tallest = 1;
        for (const auto& s : sticks) tallest = std::max(tallest, s.count);
        // Each stick's lines as (x, height): 1H multiplets as n + 1 lines, Pascal's triangle tall, drawn wider
        // than a real J so they can be read.
        std::vector<std::vector<QPointF>> lines;
        std::vector<QPointF> all;
        for (const auto& s : sticks) {
            const int n = proton ? std::min(s.coupled, 6) : 0;
            std::vector<double> row{1};
            for (int k = 0; k < n; ++k) {
                row.push_back(0);
                for (int j = k + 1; j > 0; --j) row[j] += row[j - 1];
            }
            const double peak = *std::max_element(row.begin(), row.end());
            lines.emplace_back();
            for (int j = 0; j <= n; ++j)
                lines.back().push_back({x(s.ppm) + (j - n / 2.0) * 3, row[j] / peak * s.count / tallest}), all.push_back(lines.back().back());
        }
        // The molecule as the legend, as large as fits a corner, the sticks shrunk to keep clear of it and their label.
        const QRectF bounds = documentBounds(molecule);
        const double fit = std::min({0.45 * plot.width() / std::max(bounds.width(), 1.0), 0.5 * plot.height() / std::max(bounds.height(), 1.0),
                                     22 / kBondLength});  // a bond no longer than 22 px
        const Legend legend = placeLegend(plot, bounds.size() * fit, all, fm.height() + 6);
        auto top = [&](int k) { return plot.bottom() - double(sticks[k].count) / tallest * legend.scale * plot.height(); };
        for (size_t k = 0; k < sticks.size(); ++k) {
            const QColor c = int(k) == current ? ink : palette().color(QPalette::Highlight);
            p.setPen(QPen(c, int(k) == current ? 3 : 2, Qt::SolidLine, Qt::FlatCap));
            for (QPointF l : lines[k]) p.drawLine(QPointF(l.x(), plot.bottom()), QPointF(l.x(), plot.bottom() - l.y() * legend.scale * plot.height()));
        }
        p.save();
        p.translate(legend.rect.topLeft());
        p.scale(fit, fit);
        p.translate(-bounds.topLeft());
        paintDocument(p, molecule, {ink});
        if (current >= 0) {  // the stick in hand's atoms, ringed
            p.setPen(QPen(palette().color(QPalette::Highlight), 1.5 / fit));
            p.setBrush(Qt::NoBrush);
            for (int a : sticks[current].atoms) p.drawEllipse(doc.atoms[a].pos, 7.0, 7.0);
        }
        p.restore();
        p.setPen(ink);
        if (current >= 0) {
            const QString t = describe(current).section(':', 0, 0);
            const double left = std::clamp(x(sticks[current].ppm) - fm.horizontalAdvance(t) / 2.0, 0.0, width() - fm.horizontalAdvance(t) - 0.0);
            p.drawText(QPointF(left, std::max(top(current) - 4, double(fm.ascent()))), t);
        }
        p.drawText(QRectF(plot.left(), 0, plot.width(), fm.height() + 4), Qt::AlignRight, "δ / ppm");
    }
    void mouseMoveEvent(QMouseEvent* e) override {
        int best = -1;
        double nearest = 6;  // px
        for (size_t k = 0; k < sticks.size(); ++k)
            if (double d = std::abs(x(sticks[k].ppm) - e->position().x()); d < nearest) nearest = d, best = int(k);
        setCurrent(best);
    }
    void leaveEvent(QEvent*) override { setCurrent(-1); }
    void keyPressEvent(QKeyEvent* e) override {
        const int n = int(sticks.size());
        if (e->key() == Qt::Key_Right && n) setCurrent(std::min(current + 1, n - 1));
        else if (e->key() == Qt::Key_Left && n) setCurrent(std::max(current - 1, 0));
        else QWidget::keyPressEvent(e);
    }
};

static QString uiStyle(const Theme& t) {
    const Chrome c = chrome(t);
    return QString(R"(
        QMainWindow, QDialog { background: %1; color: %4; }
        QMenuBar, QStatusBar { background: %1; color: %4; border: none; }
        QMenuBar::item { padding: 5px 9px; border-radius: 6px; }
        QMenuBar::item:selected, QMenu::item:selected { background: %7; color: %6; }
        QMenu { background: %2; color: %4; border: 1px solid %3; border-radius: 10px; padding: 5px; }
        QMenu::item { padding: 5px 20px; border-radius: 6px; }
        QToolBar#tools { background: transparent; border: none; padding: 8px 4px; }
        QFrame#toolCard { background: %2; border: 1px solid %3; border-radius: 14px; }
        QFrame#toolDivider { background: %3; border: none; }
        QToolBar#tools QToolButton { color: %4; background: transparent; border: none; border-radius: 8px; padding: 5px; }
        QToolBar#tools QToolButton:hover, QToolBar#tools QToolButton:checked { background: %7; }
        QToolBar#tools QToolButton:focus { border: 2px solid %6; }
        QToolBar#tools QToolButton::menu-button { background: transparent; border: none; width: 10px; }
        QCheckBox::indicator { width: 14px; height: 14px; background: %2; border: 1px solid %5; border-radius: 4px; }
        QCheckBox::indicator:checked { background: %6; border-color: %6; }
        QToolBar#tools QToolButton#railButton { color: %5; font-size: 11px; padding: 6px 4px 4px; min-width: 50px; }
        QToolBar#tools QToolButton#railButton:checked { color: %6; background: %7; font-weight: 600; }
        QFrame#toolFlyout { background: %2; border: 1px solid %3; border-radius: 12px; }
        QFrame#toolFlyout QToolButton { color: %4; background: transparent; border: none; border-radius: 8px; padding: 5px; }
        QFrame#toolFlyout QToolButton:hover, QFrame#toolFlyout QToolButton:checked { background: %7; }
        QFrame#toolFlyout QToolButton:focus { border: 2px solid %6; }
        QFrame#toolFlyout QToolButton::menu-button { background: transparent; border: none; width: 10px; }
        QFrame#toolFlyout QToolButton#close { color: %5; font-size: 11px; padding: 2px 6px; border: none; }
        QLabel#flyoutTitle { color: %5; font-size: 11px; font-weight: 700; letter-spacing: 1px; }
        QDockWidget#properties, QDockWidget#templates { background: %1; color: %4; border: none; }
        QDockWidget::title { background: %2; color: %4; border: 1px solid %3; border-radius: 10px; padding: 8px; }
        QFrame#panelCard { background: %2; border: 1px solid %3; border-radius: 12px; }
        QFrame#welcome { background: %2; border: 1px solid %3; border-radius: 16px; }
        QLabel#welcomeTitle { font-size: 18px; font-weight: 600; }
        QLabel#welcomeHint { color: %5; }
        QToolButton#example { background: %1; color: %4; border: 1px solid %3; border-radius: 10px; padding: 8px; }
        QToolButton#example:hover { background: %7; border-color: %6; }
        QTabBar#pageTabs { background: transparent; border: none; }
        QTabBar#pageTabs::tab { color: %5; background: transparent; border: none; border-top: 2px solid transparent;
                                padding: 5px 14px; margin-right: 2px; }
        QTabBar#pageTabs::tab:hover { color: %4; background: %7; }
        QTabBar#pageTabs::tab:selected { color: %6; background: %2; border-top-color: %6; font-weight: 600; }
        QToolButton#addPage { color: %5; background: transparent; border: none; border-radius: 6px; padding: 2px 9px; font-size: 15px; }
        QToolButton#addPage:hover { background: %7; color: %6; }
        QTreeWidget, QTextEdit, QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox {
            background: %2; color: %4; border: 1px solid %3; border-radius: 7px; padding: 3px;
            selection-background-color: %7; selection-color: %6;
        }
        QPushButton { background: %2; color: %4; border: 1px solid %3; border-radius: 8px; padding: 6px 12px; }
        QPushButton:hover { border-color: %6; }
        QPushButton:default { background: %6; color: %1; border-color: %6; }
        QScrollBar:vertical { background: %1; width: 12px; margin: 0; }
        QScrollBar:horizontal { background: %1; height: 12px; margin: 0; }
        QScrollBar::handle { background: %3; border-radius: 5px; min-width: 24px; min-height: 24px; }
        QScrollBar::handle:hover { background: %5; }
        QScrollBar::add-line, QScrollBar::sub-line { width: 0; height: 0; }
        QScrollBar::add-page, QScrollBar::sub-page { background: transparent; }
    )").arg(t.window.name(), t.surface.name(), c.border.name(), t.text.name(), c.secondary.name(), t.accent.name(),
              c.accentBg.name());
}

// A tool flyout: the tools flow into as many columns as fit its width. Drag it by any part that
// isn't a tool (its title, margins and gaps) to move it off the others (#466), or by a corner
// to resize it; the tools reflow and the box follows them.
struct FlyoutFrame : QObject {
    struct Tool {
        QWidget* w;
        int span;
        bool newRow;
    };
    QWidget* fly;
    QGridLayout* grid;
    std::vector<Tool> tools;
    bool nextNewRow = false;
    static constexpr int kDefaultCols = 4;
    int cols = kDefaultCols;
    int corner = 0;  // the corner being dragged: Qt::Edges, or 0 when moving
    QPoint grab;  // where the drag began, in the window
    QRect start;  // and the box then

    FlyoutFrame(QWidget* flyout, QGridLayout* g) : QObject(flyout), fly(flyout), grid(g) {
        flyout->setCursor(Qt::SizeAllCursor);
        flyout->setMouseTracking(true);
        flyout->installEventFilter(this);
    }
    void newRow() { nextNewRow = true; }  // the next tool starts a row, at the default width (resized, the tools just flow)
    void add(QWidget* w, int span = 1) {
        tools.push_back({w, span, std::exchange(nextNewRow, false)});
        place(cols);
    }
    void place(int c) {
        cols = c;
        int slot = 0;
        for (auto& t : tools) {
            grid->removeWidget(t.w);
            const int span = std::min(t.span, c);
            if (t.newRow && c == kDefaultCols) slot = (slot + c - 1) / c * c;
            if (slot % c + span > c) slot = (slot / c + 1) * c;
            grid->addWidget(t.w, slot / c, slot % c, 1, span);
            slot += span;
        }
    }
    QSize sizeFor(int c) {
        place(c);
        fly->layout()->invalidate();
        fly->layout()->activate();
        return fly->layout()->minimumSize();
    }
    // The column count for a box dragged to `want`: as many as fit the width, or, when the drag
    // was mostly vertical, the fewest that fit the height. Never wider than the window.
    void fit(QSize want, int maxWidth, bool byHeight) {
        int most = 0;
        for (auto& t : tools) most += t.span;
        int c = 1;
        if (byHeight) {
            while (c < most && sizeFor(c).height() > want.height() && sizeFor(c + 1).width() <= maxWidth) ++c;
        } else {
            for (int n = 1; n <= most && sizeFor(n).width() <= std::min(want.width(), maxWidth); ++n) c = n;
        }
        place(c);
    }
    int cornerAt(QPoint p) const {
        constexpr int zone = 12;
        int edges = 0;
        if (p.x() < zone) edges |= Qt::LeftEdge;
        if (p.x() >= fly->width() - zone) edges |= Qt::RightEdge;
        if (p.y() < zone) edges |= Qt::TopEdge;
        if (p.y() >= fly->height() - zone) edges |= Qt::BottomEdge;
        const bool horizontal = edges & (Qt::LeftEdge | Qt::RightEdge), vertical = edges & (Qt::TopEdge | Qt::BottomEdge);
        return horizontal && vertical ? edges : 0;
    }
    static Qt::CursorShape cursorFor(int edges) {
        if (!edges) return Qt::SizeAllCursor;
        const bool topLeftOrBottomRight = bool(edges & Qt::LeftEdge) == bool(edges & Qt::TopEdge);
        return topLeftOrBottomRight ? Qt::SizeFDiagCursor : Qt::SizeBDiagCursor;
    }
    void resize(QPoint mouse) {
        QWidget* window = fly->parentWidget();
        const QRect old = start;
        const int left = corner & Qt::LeftEdge ? mouse.x() : old.left(), right = corner & Qt::RightEdge ? mouse.x() : old.right();
        const int top = corner & Qt::TopEdge ? mouse.y() : old.top(), bottom = corner & Qt::BottomEdge ? mouse.y() : old.bottom();
        const QPoint moved = mouse - grab;
        fit({right - left + 1, bottom - top + 1}, window->width(), std::abs(moved.y()) > std::abs(moved.x()));
        const QSize size = sizeFor(cols);
        // The corner opposite the one dragged stays put.
        QPoint at(corner & Qt::LeftEdge ? old.right() + 1 - size.width() : old.left(),
                  corner & Qt::TopEdge ? old.bottom() + 1 - size.height() : old.top());
        at.setX(std::clamp(at.x(), 0, std::max(0, window->width() - size.width())));
        at.setY(std::clamp(at.y(), 0, std::max(0, window->height() - size.height())));
        fly->setGeometry({at, size});
    }
    bool eventFilter(QObject*, QEvent* e) override {
        if (e->type() == QEvent::Show) {  // the tools keep the ordinary pointer
            for (auto* b : fly->findChildren<QAbstractButton*>()) b->setCursor(Qt::ArrowCursor);
            return false;
        }
        const auto type = e->type();
        if (type != QEvent::MouseButtonPress && type != QEvent::MouseMove && type != QEvent::MouseButtonRelease) return false;
        const auto* me = static_cast<QMouseEvent*>(e);
        QWidget* window = fly->parentWidget();
        const QPoint mouse = window->mapFromGlobal(me->globalPosition().toPoint());
        if (type == QEvent::MouseButtonPress) {
            corner = cornerAt(me->position().toPoint());
            grab = corner ? mouse : mouse - fly->pos();
            start = fly->geometry();
            fly->raise();  // over another one it overlaps
            return true;
        }
        if (type == QEvent::MouseButtonRelease) {
            corner = 0;
            return false;
        }
        if (!(me->buttons() & Qt::LeftButton)) {
            fly->setCursor(cursorFor(cornerAt(me->position().toPoint())));
            return false;
        }
        if (corner) {
            resize(mouse);
            return true;
        }
        const QPoint to = mouse - grab;
        fly->move(std::clamp(to.x(), 0, std::max(0, window->width() - fly->width())),
                  std::clamp(to.y(), 0, std::max(0, window->height() - fly->height())));
        return true;
    }
};

MainWindow::MainWindow() : undo_(new QUndoStack(this)), canvas_(new Canvas(undo_, this)) {
#ifdef Q_OS_MACOS
    // Registers itself with Qt, once. Qt deletes its converters when QApplication goes, so this
    // one must be on the heap: a static was freed there, aborting on quit.
    static auto* chemDraw = new ChemDrawPasteboard;
    Q_UNUSED(chemDraw);
#endif
#ifdef Q_OS_WIN
    // Registers itself with Qt's Windows plugin (only there: it asserts on any other), once. Never
    // deleted: its destructor would look for the plugin after QApplication has gone.
    static const bool emf = QGuiApplication::platformName() == "windows" && (new EmfClipboard) && (new EmbedClipboard);
    Q_UNUSED(emf);
#endif
    undoGroup_ = new QUndoGroup(this);
    undoGroup_->addStack(undo_);
    undoGroup_->setActiveStack(undo_);
    pages_ = {{tr("Page 1"), Document{}, undo_}};
    // Pages as spreadsheet-style tabs under the canvas: + adds one, drag to reorder,
    // double-click to rename, right-click to rename or delete.
    auto* central = new QWidget;
    auto* column = new QVBoxLayout(central);
    column->setContentsMargins(0, 0, 0, 0);
    column->setSpacing(0);
    column->addWidget(canvas_, 1);
    auto* strip = new QHBoxLayout;
    strip->setContentsMargins(6, 0, 6, 0);
    strip->setSpacing(2);
    pageTabs_ = new QTabBar;
    pageTabs_->setObjectName("pageTabs");
    pageTabs_->setShape(QTabBar::RoundedSouth);
    pageTabs_->setDocumentMode(true);
    pageTabs_->setExpanding(false);
    pageTabs_->setDrawBase(false);
    pageTabs_->setMovable(true);
    pageTabs_->setContextMenuPolicy(Qt::CustomContextMenu);
    pageTabs_->addTab(pages_[0].name);
    auto* addPageButton = new QToolButton;
    addPageButton->setObjectName("addPage");
    addPageButton->setText("+");
    addPageButton->setAutoRaise(true);
    addPageButton->setToolTip(tr("New Page"));
    addPageButton->setAccessibleName(tr("New Page"));
    strip->addWidget(pageTabs_);
    strip->addWidget(addPageButton);
    strip->addStretch();
    column->addLayout(strip);
    setCentralWidget(central);
    connect(addPageButton, &QToolButton::clicked, this, &MainWindow::addPage);
    connect(pageTabs_, &QTabBar::currentChanged, this, [this](int i) {
        if (i >= 0 && i != page_) showPage(i);
    });
    connect(pageTabs_, &QTabBar::tabMoved, this, [this](int from, int to) {
        PageState moved = pages_[from];
        pages_.erase(pages_.begin() + from);
        pages_.insert(pages_.begin() + to, moved);
        page_ = pageTabs_->currentIndex();
        pagesEdited_ = true;
        renumberPages();
        updateTitle();
    });
    connect(canvas_, &Canvas::documentChanged, this, &MainWindow::renumberPages);  // an edit, an undo, another page
    connect(pageTabs_, &QTabBar::tabBarDoubleClicked, this, &MainWindow::renamePage);
    connect(pageTabs_, &QWidget::customContextMenuRequested, this, [this](QPoint at) {
        const int i = pageTabs_->tabAt(at);
        if (i < 0) return;
        QMenu menu;
        menu.addAction(tr("&Rename…"), this, [=, this] { renamePage(i); });
        menu.addAction(tr("&Delete"), this, [=, this] { deletePage(i); })->setEnabled(pages_.size() > 1);
        menu.exec(pageTabs_->mapToGlobal(at));
    });
    setAcceptDrops(true);
    canvas_->viewport()->setAcceptDrops(false);  // so a dropped file reaches the window
    setWindowTitle("Penzene " PENZENE_BUILD);
    resize(1100, 750);
    // Template library (built before the menus too); filled the first time it's shown.
    templateDock_ = new QDockWidget(tr("Templates"), this);
    templateDock_->setObjectName("templates");
    templates_ = new QTreeWidget;
    templates_->setHeaderHidden(true);
    templates_->setAccessibleName(tr("Templates"));
    templates_->setIconSize({56, 40});
    templates_->setContextMenuPolicy(Qt::CustomContextMenu);
    auto* templateCard = new QFrame;
    templateCard->setObjectName("panelCard");
    auto* templateLayout = new QVBoxLayout(templateCard);
    templateLayout->setContentsMargins(8, 8, 8, 8);
    templateLayout->addWidget(templates_);
    templateDock_->setWidget(templateCard);
    addDockWidget(Qt::RightDockWidgetArea, templateDock_);
    templateDock_->hide();
    connect(templateDock_, &QDockWidget::visibilityChanged, this, [this](bool shown) {
        if (shown && !templates_->topLevelItemCount()) fillTemplates();
    });
    connect(templates_, &QTreeWidget::itemActivated, this, &MainWindow::insertTemplate);
    connect(templates_, &QTreeWidget::itemClicked, this, &MainWindow::insertTemplate);
    connect(templates_, &QTreeWidget::customContextMenuRequested, this, [this](QPoint at) {
        QTreeWidgetItem* item = templates_->itemAt(at);
        if (!item || !item->data(0, Qt::UserRole + 1).toBool()) return;  // only the user's own
        QMenu menu;
        menu.addAction(tr("Delete Template"), this, [this, item] {
            removeUserTemplate(item->text(0));
            fillTemplates();
        });
        menu.exec(templates_->viewport()->mapToGlobal(at));
    });
    // Properties panel (built before the menus, which offer its toggle): descriptors for the selection or everything.
    profileDock_ = new QDockWidget(tr("Properties"), this);
    profileDock_->setObjectName("properties");
    auto* panel = new QFrame;
    panel->setObjectName("panelCard");
    auto* panelLayout = new QVBoxLayout(panel);
    panelLayout->setContentsMargins(12, 12, 12, 12);
    profile_ = new QLabel;
    profile_->setTextInteractionFlags(Qt::TextSelectableByMouse);
    profile_->setAlignment(Qt::AlignTop | Qt::AlignLeft);
    auto* copyProfile = new QPushButton(tr("Copy as Text"));
    connect(copyProfile, &QPushButton::clicked, this, [this] { QApplication::clipboard()->setText(profileText_); });
    panelLayout->addWidget(profile_);
    panelLayout->addWidget(copyProfile);
    panelLayout->addStretch();
    panel->setMinimumWidth(300);  // room for the values beside their names
    profileDock_->setWidget(panel);
    addDockWidget(Qt::RightDockWidgetArea, profileDock_);
    profileDock_->hide();
    connect(profileDock_, &QDockWidget::visibilityChanged, this, &MainWindow::updateProfile);
    connect(canvas_, &Canvas::documentChanged, this, &MainWindow::updateProfile);
    connect(canvas_, &Canvas::selectionChanged, this, &MainWindow::updateProfile);
    // Mass Spec panel, beside Properties: the isotope pattern of the selection or everything.
    massDock_ = new QDockWidget(tr("Mass Spec"), this);
    massDock_->setObjectName("massSpec");
    auto* massCard = new QFrame;
    massCard->setObjectName("panelCard");
    auto* massLayout = new QVBoxLayout(massCard);
    massLayout->setContentsMargins(12, 12, 12, 12);
    ion_ = new QComboBox;
    ion_->addItems({"[M]⁺•", "[M+H]⁺", "[M+Na]⁺", "[M−H]⁻"});
    ion_->setCurrentIndex(1);
    ion_->setAccessibleName(tr("Ion"));  // macOS reads the current item; Windows needs a name
    spectrum_ = new SpectrumView;
    spectrum_->setAccessibleName(tr("Isotope pattern"));
    massLayout->addWidget(ion_);
    massLayout->addWidget(spectrum_, 1);
    eiIons_ = new QLabel;
    eiIons_->setWordWrap(true);
    eiIons_->setTextInteractionFlags(Qt::TextSelectableByMouse);
    massLayout->addWidget(eiIons_);
    auto* massButtons = new QHBoxLayout;
    auto* copyHrms = new QPushButton(tr("Copy HRMS Line"));
    copyHrms->setToolTip(tr("The ion's calculated mass, for the supporting information"));
    connect(copyHrms, &QPushButton::clicked, this, [this] {
        const QString line = chem::hrmsLine(canvas_->selectedSubset(), chem::Ion(ion_->currentIndex()));
        if (line.isEmpty()) return statusBar()->showMessage(tr("No valid structure for an HRMS line"), 4000);  // the clipboard kept
        QGuiApplication::clipboard()->setText(line);
    });
    auto* exportPattern = new QPushButton(tr("Export CSV…"));
    exportPattern->setToolTip(tr("The isotope pattern as m/z and intensity"));
    connect(exportPattern, &QPushButton::clicked, this, [this] {
        if (spectrum_->peaks.empty()) return;
        const QString path = QFileDialog::getSaveFileName(this, tr("Export Isotope Pattern"), "isotope-pattern.csv", tr("CSV (*.csv)"));
        QFile f(path);
        if (path.isEmpty() || !f.open(QIODevice::WriteOnly | QIODevice::Text)) return;
        QTextStream out(&f);
        out << "mz,intensity\n";
        for (const auto& k : spectrum_->peaks) out << QString::number(k.mz, 'f', 4) << ',' << QString::number(k.intensity, 'f', 2) << '\n';
    });
    massButtons->addWidget(copyHrms);
    massButtons->addWidget(exportPattern);
    massLayout->addLayout(massButtons);
    massCard->setMinimumWidth(300);
    massDock_->setWidget(massCard);
    addDockWidget(Qt::RightDockWidgetArea, massDock_);
    massDock_->hide();
    connect(massDock_, &QDockWidget::visibilityChanged, this, &MainWindow::updateMassSpec);
    connect(ion_, &QComboBox::currentIndexChanged, this, &MainWindow::updateMassSpec);
    connect(canvas_, &Canvas::documentChanged, this, &MainWindow::updateMassSpec);
    connect(canvas_, &Canvas::selectionChanged, this, &MainWindow::updateMassSpec);
    // NMR panel (#444): the predicted spectrum of the selection or everything.
    nmrDock_ = new QDockWidget(tr("NMR"), this);
    nmrDock_->setObjectName("nmr");
    auto* nmrCard = new QFrame;
    nmrCard->setObjectName("panelCard");
    auto* nmrLayout = new QVBoxLayout(nmrCard);
    nmrLayout->setContentsMargins(12, 12, 12, 12);
    nucleus_ = new QComboBox;
    nucleus_->addItems({"¹³C", "¹H"});
    nucleus_->setAccessibleName(tr("Nucleus"));
    nmr_ = new NmrView;
    nmr_->setAccessibleName(tr("Predicted spectrum"));
    nmr_->setToolTip(tr("Point at a stick, or press Left and Right, to light its atoms"));
    nmr_->light = [this](const std::vector<int>& atoms) { canvas_->setHighlight(QSet<int>(atoms.begin(), atoms.end())); };
    auto* notice = new QLabel(tr("Predicted; ~ in a stick's label marks a weaker match.") + " " + nmrshiftdbNotice().join(" "));
    notice->setWordWrap(true);
    notice->setTextInteractionFlags(Qt::TextSelectableByMouse);
    notice->setStyleSheet("font-size: 10px;");
    nmrLayout->addWidget(nucleus_);
    nmrLayout->addWidget(nmr_, 1);
    auto* copyNmr = new QPushButton(tr("Copy SI Line"));
    copyNmr->setToolTip(tr("The predicted shifts as a supporting-information line, to replace with measured ones"));
    connect(copyNmr, &QPushButton::clicked, this, [this] {
        const QString line = chem::nmrLine(nmr_->doc, nmr_->proton, nmr_->atoms);
        if (line.isEmpty()) return statusBar()->showMessage(tr("No predicted shifts to copy"), 4000);  // the clipboard kept
        QGuiApplication::clipboard()->setText(line);
    });
    nmrLayout->addWidget(copyNmr);
    nmrLayout->addWidget(notice);
    nmrCard->setMinimumWidth(300);
    nmrDock_->setWidget(nmrCard);
    addDockWidget(Qt::RightDockWidgetArea, nmrDock_);
    nmrDock_->hide();
    connect(nmrDock_, &QDockWidget::visibilityChanged, this, &MainWindow::updateNmr);
    connect(nucleus_, &QComboBox::currentIndexChanged, this, &MainWindow::updateNmr);
    connect(canvas_, &Canvas::documentChanged, this, &MainWindow::updateNmr);
    connect(canvas_, &Canvas::selectionChanged, this, &MainWindow::updateNmr);
    connect(canvas_, &Canvas::hotspotAtomChanged, this, [this](int atom) {  // an atom under the pointer: its stick
        const auto& s = nmr_->sticks;
        const auto k = std::find_if(s.begin(), s.end(), [&](const chem::NmrStick& k) { return std::count(k.atoms.begin(), k.atoms.end(), atom); });
        if (nmrDock_->isVisible()) nmr_->setCurrent(k == s.end() ? -1 : int(k - s.begin()));
    });
    buildTools();
    buildMenus();
    connect(undoGroup_, &QUndoGroup::cleanChanged, this, &MainWindow::updateTitle);
    updateTitle();
    info_ = new QLabel;
    info_->setTextInteractionFlags(Qt::TextSelectableByMouse);
    statusBar()->addPermanentWidget(info_);
    connect(canvas_, &Canvas::documentChanged, this, &MainWindow::updateInfo);
    connect(canvas_, &Canvas::selectionChanged, this, &MainWindow::updateInfo);

    auto* autosaver = new QTimer(this);
    connect(autosaver, &QTimer::timeout, this, &MainWindow::autosave);
    autosaver->start(60 * 1000);
    buildWelcome();
    for (size_t i = 0; i < flyouts_.size(); ++i) {  // the tool layout left at the last quit, after the theme sizes it
        const QVariantList v = QSettings().value(QString("toolLayout/%1").arg(i)).toList();
        if (v.size() != 3) continue;
        FlyoutFrame* f = flyouts_[i];
        f->fly->ensurePolished();
        f->fly->setGeometry({v[1].toPoint(), f->sizeFor(std::max(1, v[2].toInt()))});  // as a corner drag sizes it
        f->fly->setVisible(v[0].toBool());
    }
}

// A card over the empty canvas: examples to open, and where to learn the keys. It goes away
// once there's a drawing, or at the first click on the canvas.
void MainWindow::buildWelcome() {
    welcome_ = new QFrame;
    welcome_->setObjectName("welcome");
    auto* layout = new QVBoxLayout(welcome_);
    layout->setContentsMargins(22, 20, 22, 18);
    layout->setSpacing(10);
    auto* title = new QLabel(tr("Start drawing"));
    title->setObjectName("welcomeTitle");
    auto* hint = new QLabel(tr("Click anywhere on the page, or open an example:"));
    hint->setObjectName("welcomeHint");
    layout->addWidget(title);
    layout->addWidget(hint);
    auto* row = new QHBoxLayout;
    row->setSpacing(10);
    for (const auto& [name, doc] : exampleDocuments()) {
        auto* b = new QToolButton;
        b->setObjectName("example");
        b->setText(name);
        b->setToolButtonStyle(Qt::ToolButtonTextUnderIcon);
        b->setIconSize(kExampleIcon);
        welcomeExamples_.push_back({b, doc});
        connect(b, &QToolButton::clicked, this, [this, doc = doc] {
            if (!maybeSave()) return;
            setPages({{tr("Page 1"), doc}});
            canvas_->fitToDocument();
            path_.clear();
            updateTitle();
        });
        row->addWidget(b);
    }
    layout->addLayout(row);
    auto* links = new QLabel(QString("<a href='keys'>%1</a> &nbsp;·&nbsp; <a href='%2'>%3</a>")
                                 .arg(tr("Keyboard shortcuts"), kDocsUrl, tr("Documentation")));
    links->setObjectName("welcomeLinks");
    connect(links, &QLabel::linkActivated, this, [this](const QString& link) {
        if (link == "keys") shortcutsAction_->trigger();
        else QDesktopServices::openUrl(QUrl(link));
    });
    layout->addWidget(links);

    auto* centre = new QGridLayout(canvas_->viewport());
    centre->addWidget(welcome_, 0, 0, Qt::AlignCenter);
    canvas_->viewport()->installEventFilter(this);
    qApp->installEventFilter(this);  // macOS hands over files opened from Finder as QFileOpenEvents (#332)
    connect(canvas_, &Canvas::documentChanged, welcome_, [this] {
        if (!canvas_->document().empty()) welcome_->hide();
    });
    welcome_->setVisible(canvas_->document().empty());
    paintExamples();
}

// The examples' previews in the theme's ink (again after a theme change).
// A drawing as an icon in the interface's ink, so it reads on light and dark panels alike.
static QIcon drawingIcon(const Document& doc, QSize size, QColor ink, qreal ratio) {
    QImage img(size * ratio, QImage::Format_ARGB32_Premultiplied);
    img.setDevicePixelRatio(ratio);
    img.fill(Qt::transparent);
    QPainter p(&img);
    p.setRenderHint(QPainter::Antialiasing);
    const QRectF r = documentBounds(doc);
    const double k = std::min(size.width() / r.width(), size.height() / r.height()) * 0.92;
    p.translate(size.width() / 2.0, size.height() / 2.0);
    p.scale(k, k);
    p.translate(-r.center());
    paintDocument(p, doc, {ink, ink});
    p.end();
    return QIcon(QPixmap::fromImage(img));
}

void MainWindow::paintExamples() {
    const QColor ink = palette().color(QPalette::WindowText);
    for (auto& [button, doc] : welcomeExamples_)
        button->setIcon(drawingIcon(doc, kExampleIcon, ink, devicePixelRatioF()));
}

// A smaller window pulls the open flyouts back inside it, so none is stranded off-screen.
void MainWindow::resizeEvent(QResizeEvent* e) {
    QMainWindow::resizeEvent(e);
    for (auto* f : flyouts_)
        if (!f->fly->isHidden())  // restored ones too, before the window is first shown
            f->fly->move(std::clamp(f->fly->x(), 0, std::max(0, width() - f->fly->width())),
                         std::clamp(f->fly->y(), 0, std::max(0, height() - f->fly->height())));
}

bool MainWindow::eventFilter(QObject* watched, QEvent* e) {
    if (watched == qApp && e->type() == QEvent::FileOpen) {
        if (maybeSave()) openFile(static_cast<QFileOpenEvent*>(e)->file());
        return true;
    }
    if (watched == canvas_->viewport() && e->type() == QEvent::MouseButtonPress) {
        if (welcome_->isVisible() && !welcome_->geometry().contains(static_cast<QMouseEvent*>(e)->position().toPoint()))
            welcome_->hide();  // and the click goes on to draw
    }
    return QMainWindow::eventFilter(watched, e);
}

QString formulaHtml(const std::string& formula) {
    static const QRegularExpression charge("([+-])(\\d*)$");
    QString f = QString::fromStdString(formula);
    const auto m = charge.match(f);
    QString sup;
    if (m.hasMatch()) sup = m.captured(2) + (m.captured(1) == "+" ? "+" : "−"), f.chop(m.capturedLength());
    // Counts subscripted; an isotope written apart ("[13C]") as its mass number up front, ¹³C.
    static const QRegularExpression token("\\[(\\d+)([A-Z][a-z]?)\\]|(\\d+)");
    f = f.toHtmlEscaped();
    QString html;
    qsizetype last = 0;
    for (auto it = token.globalMatch(f); it.hasNext();) {
        const auto m = it.next();
        html += f.mid(last, m.capturedStart() - last);
        html += m.hasCaptured(1) ? "<sup>" + m.captured(1) + "</sup>" + m.captured(2) : "<sub>" + m.captured(3) + "</sub>";
        last = m.capturedEnd();
    }
    html += f.mid(last);
    return sup.isEmpty() ? html : html + "<sup>" + sup + "</sup>";
}

void MainWindow::updateProfile() {
    if (!profileDock_->isVisible()) return;
    const auto p = chem::profile(canvas_->selectedSubset());
    if (!p) {
        profile_->setText(tr("Draw or select a valid structure."));
        profileText_.clear();
        return;
    }
    const QString formula = formulaHtml(p->basic.formula);
    QStringList analysis, plain;
    for (const auto& [el, pct] : p->elemental) {
        analysis << QString("%1 %2").arg(QString::fromStdString(el)).arg(pct, 0, 'f', 2);
    }
    const QList<std::pair<QString, QString>> rows{
        {tr("Formula"), formula},
        {tr("MW"), QString::number(p->basic.mw, 'f', 2)},
        {tr("Exact mass"), QString::number(p->basic.exactMass, 'f', 4)},
        {tr("Elemental (%)"), analysis.join(", ")},
        {tr("cLogP"), QString::number(p->logP, 'f', 2)},
        {tr("TPSA (Å²)"), QString::number(p->tpsa, 'f', 1)},
        {tr("H-bond donors"), QString::number(p->hbd)},
        {tr("H-bond acceptors"), QString::number(p->hba)},
        {tr("Rotatable bonds"), QString::number(p->rotatable)},
        {tr("Heavy atoms"), QString::number(p->heavyAtoms)},
        {tr("Lipinski (Ro5)"), p->lipinskiViolations ? tr("%n violation(s)", "", p->lipinskiViolations) : tr("passes")},
        {tr("Veber"), p->veber ? tr("passes") : tr("fails")},
    };
    QString html = "<table cellspacing='4'>";
    for (const auto& [k, v] : rows) {
        html += QString("<tr><td><b>%1</b></td><td>%2</td></tr>").arg(k, v);
        // The formula as RDKit writes it: stripped of tags, "O4S2−" would read as two sulfurs.
        plain << k + "\t" + (k == tr("Formula") ? QString::fromStdString(p->basic.formula) : QString(v).remove(QRegularExpression("<[^>]*>")));
    }
    profile_->setText(html + "</table>");
    profileText_ = plain.join("\n");
}

void MainWindow::updateMassSpec() {
    if (!massDock_->isVisible()) return;
    spectrum_->peaks = chem::isotopePattern(canvas_->selectedSubset(), chem::Ion(ion_->currentIndex()));
    QStringList sticks;  // for screen readers
    for (const auto& k : spectrum_->peaks)
        sticks << QString("%1 (%2%)").arg(k.mz, 0, 'f', 4).arg(k.intensity, 0, 'f', 1);
    spectrum_->setAccessibleDescription(sticks.join(", "));
    spectrum_->update();
    const auto ions = ion_->currentIndex() == int(chem::Ion::M) ? chem::eiIons(canvas_->selectedSubset()) : std::vector<chem::EiIon>{};
    QString rows;
    for (const auto& i : ions)
        rows += QString("<tr><td>%1</td><td>%2</td><td>%3</td></tr>").arg(i.mz, 0, 'f', 4).arg(i.formula, i.from.toHtmlEscaped());
    eiIons_->setText(tr("<b>EI ions to look for</b> (from the groups drawn; no intensities)") +
                     "<table cellspacing=\"4\">" + rows + "</table>");
    eiIons_->setVisible(!ions.empty());
}

void MainWindow::updateNmr() {
    nmr_->current = -1;
    canvas_->setHighlight({});
    if (!nmrDock_->isVisible()) return;
    // One molecule's spectrum: the one holding most of the selection, else the largest.
    const Document& d = canvas_->document();
    const QSet<int>& selected = canvas_->selection();
    std::vector<bool> seen(d.atoms.size());
    std::pair<int, size_t> best{-1, 0};
    nmr_->atoms.clear();
    for (int i = 0; i < int(d.atoms.size()); ++i) {
        if (seen[i]) continue;
        auto mol = edit::moleculeOf(d, i);
        int picked = 0;
        for (int a : mol) seen[a] = true, picked += selected.contains(a);
        if (std::pair{picked, mol.size()} > best) best = {picked, mol.size()}, nmr_->atoms = std::move(mol);
    }
    std::sort(nmr_->atoms.begin(), nmr_->atoms.end());
    nmr_->proton = nucleus_->currentIndex() == 1;
    nmr_->doc = d;
    nmr_->molecule = d;
    std::vector<int> others;
    for (int i = 0; i < int(d.atoms.size()); ++i)
        if (!std::binary_search(nmr_->atoms.begin(), nmr_->atoms.end(), i)) others.push_back(i);
    nmr_->molecule.removeAtoms(others);
    nmr_->molecule.arrows.clear(), nmr_->molecule.texts.clear();
    nmr_->sticks = nmr_->atoms.empty() ? std::vector<chem::NmrStick>{} : chem::nmrSticks(d, nmr_->proton, nmr_->atoms);
    nmr_->setAccessibleDescription(nmr_->reading());
    nmr_->update();
}

// Formula and masses of the selection, or of everything.
void MainWindow::updateInfo() {
    auto p = chem::properties(canvas_->selectedSubset());
    if (!p) return info_->clear();
    const QString f = formulaHtml(p->formula);
    info_->setText(tr("%1 &nbsp;·&nbsp; MW %2 &nbsp;·&nbsp; exact mass %3")
                       .arg(f)
                       .arg(p->mw, 0, 'f', 2)
                       .arg(p->exactMass, 0, 'f', 4));
}

// Children are deleted after this destructor has run, and some signal on the way
// out (the undo stack's cleanChanged, a dock's visibilityChanged as it hides):
// none of that may reach a half-destroyed window.
MainWindow::~MainWindow() {
    for (QObject* child : findChildren<QObject*>()) child->disconnect(this);
}

void MainWindow::updateTitle() {
    QString name = path_.isEmpty() ? tr("Untitled") : QFileInfo(path_).fileName();
#ifdef Q_OS_WIN
    if (embeddedSave_) name = embeddedIn_.isEmpty() ? tr("Embedded drawing") : tr("Drawing in %1").arg(embeddedIn_);
#endif
    setWindowTitle(name + "[*] — Penzene " PENZENE_BUILD);
    setWindowModified(!isClean());
}

bool MainWindow::isClean() const {
    return !pagesEdited_ && std::all_of(pages_.begin(), pages_.end(), [](const PageState& p) { return p.undo->isClean(); });
}

std::vector<Sheet> MainWindow::sheets() const {
    std::vector<Sheet> out;
    for (int i = 0; i < int(pages_.size()); ++i)
        out.push_back({pages_[i].name, i == page_ ? canvas_->document() : pages_[i].doc});
    return out;
}

void MainWindow::setPages(const std::vector<Sheet>& sheets) {
    QUndoStack* first = pages_[0].undo;
    // Off the old pages before their stacks go: deleting the group's active stack signals
    // cleanChanged, and updateTitle would read the stacks already deleted.
    undoGroup_->setActiveStack(first);
    canvas_->setUndoStack(first);
    const auto old = std::exchange(pages_, {});
    for (size_t i = 1; i < old.size(); ++i) delete old[i].undo;
    first->clear();
    for (const Sheet& s : sheets) {
        QUndoStack* undo = pages_.empty() ? first : new QUndoStack(this);
        undoGroup_->addStack(undo);
        pages_.push_back({s.name, s.doc, undo});
    }
    page_ = 0;
    pagesEdited_ = false;
    penz1_ = false;
    undo_ = first;
    canvas_->setUndoStack(undo_);
    undoGroup_->setActiveStack(undo_);
    {
        QSignalBlocker quiet(pageTabs_);
        while (pageTabs_->count()) pageTabs_->removeTab(0);
        for (const auto& p : pages_) pageTabs_->addTab(p.name);
        pageTabs_->setCurrentIndex(0);
    }
    canvas_->setSelection({});  // the old drawing's indices mean nothing in this one
    canvas_->setDocumentSilently(pages_[0].doc);
    updateTitle();
}

// Each page's numbers start where the page before's stop. Only the page on the canvas is edited,
// so the others are renumbered here, outside their undo history: numbers follow the pages.
void MainWindow::renumberPages() {
    edit::CompoundCount count;
    for (int i = 0; i < int(pages_.size()); ++i) {
        if (i != page_) {
            count = edit::renumberCompounds(pages_[i].doc, count);
            continue;
        }
        canvas_->setCompoundStart(count);
        Document doc = canvas_->document();
        count = edit::renumberCompounds(doc, count);
        if (!(doc == canvas_->document())) canvas_->setDocumentSilently(doc);  // stale after an undo or an earlier page's edit
    }
}

void MainWindow::showPage(int i) {
    if (i < 0 || i >= int(pages_.size())) return;
    pages_[page_].doc = canvas_->document();
    page_ = i;
    undo_ = pages_[i].undo;
    canvas_->setUndoStack(undo_);
    undoGroup_->setActiveStack(undo_);
    canvas_->setSelection({});
    canvas_->setDocumentSilently(pages_[i].doc);
    QSignalBlocker quiet(pageTabs_);
    pageTabs_->setCurrentIndex(i);
}

void MainWindow::addPage() {
    int n = 1;
    auto taken = [this](const QString& name) {
        return std::any_of(pages_.begin(), pages_.end(), [&](const PageState& p) { return p.name == name; });
    };
    while (taken(tr("Page %1").arg(n))) ++n;
    Document blank;
    blank.style = QSettings().value("defaultStyle").toString();  // Edit > Preferences
    auto* undo = new QUndoStack(this);
    undoGroup_->addStack(undo);
    pages_.push_back({tr("Page %1").arg(n), blank, undo});
    {
        QSignalBlocker quiet(pageTabs_);
        pageTabs_->addTab(pages_.back().name);
    }
    pagesEdited_ = true;
    showPage(int(pages_.size()) - 1);
    updateTitle();
}

void MainWindow::renamePage(int i) {
    if (i < 0 || i >= int(pages_.size())) return;
    bool ok = false;
    const QString name = QInputDialog::getText(this, tr("Rename Page"), tr("Name:"), QLineEdit::Normal,
                                               pages_[i].name, &ok).trimmed();
    if (!ok || name.isEmpty() || name == pages_[i].name) return;
    pages_[i].name = name;
    pageTabs_->setTabText(i, name);
    pagesEdited_ = true;
    updateTitle();
}

void MainWindow::deletePage(int i) {
    if (pages_.size() < 2 || i < 0 || i >= int(pages_.size())) return;
    const Document& doc = i == page_ ? canvas_->document() : pages_[i].doc;
    if (!doc.empty() && QMessageBox::question(this, tr("Delete Page"),
                                              tr("Delete %1 and its drawing? This can't be undone.").arg(pages_[i].name)) !=
                            QMessageBox::Yes)
        return;
    if (i == page_) showPage(i == 0 ? 1 : i - 1);  // step off it first
    delete pages_[i].undo;
    pages_.erase(pages_.begin() + i);
    if (page_ > i) --page_;
    {
        QSignalBlocker quiet(pageTabs_);
        pageTabs_->removeTab(i);
        pageTabs_->setCurrentIndex(page_);
    }
    pagesEdited_ = true;
    renumberPages();
    updateTitle();
}

// Cut from this page, pasted on that one: an undo step on each.
void MainWindow::moveSelectionToPage(int i) {
    if (i == page_ || (canvas_->selection().isEmpty() && canvas_->selectedArrows().isEmpty() &&
                       canvas_->selectedTexts().isEmpty()))
        return;
    const Document part = canvas_->selectedSubset();
    canvas_->deleteSelection();
    showPage(i);
    canvas_->insert(part, tr("Move to %1").arg(pages_[i].name));
}

bool MainWindow::openFile(const QString& path) {
    const QString ext = QFileInfo(path).suffix().toLower();
    std::vector<Sheet> sheets;
    if (ext == "penz") {
        QFile f(path);
        if (f.open(QIODevice::ReadOnly)) sheets = sheetsFromJson(f.readAll());
    } else if (auto doc = chem::readFile(path)) {
        sheets = {{tr("Page 1"), *doc}};
    }
    if (sheets.empty()) {
        const QFileInfo info(path);
        if (!info.exists() && info.dir().exists()) {  // gone from its folder, not on a drive that isn't mounted (#488)
            QStringList files = recentFiles();
            files.removeAll(info.absoluteFilePath());
            QSettings().setValue("recentFiles", files);
        }
        QMessageBox::warning(this, tr("Open"),
                             !info.exists()       ? tr("%1 doesn't exist. It may have been moved or deleted.").arg(path)
                             : !info.isReadable() ? tr("%1 can't be opened for reading.").arg(path)
                                                  : tr("%1 is not a structure file I can read.").arg(path));
        return false;
    }
    setPages(sheets);
    canvas_->fitToDocument();
    // Save goes back only to a file that holds everything saved: never over a
    // ChemDraw file, or a library (SDF, SMILES) with one MOL block (#315).
    path_ = ext == "penz" || ext == "mol" ? path : QString();
    remember(path);
    updateTitle();
    return true;
}

QStringList MainWindow::recentFiles() const { return QSettings().value("recentFiles").toStringList(); }

void MainWindow::remember(const QString& path) {
    QStringList files = recentFiles();
    files.removeAll(QFileInfo(path).absoluteFilePath());
    files.prepend(QFileInfo(path).absoluteFilePath());
    QSettings().setValue("recentFiles", files.mid(0, 10));
}

// One autosave per process, locked while the process runs: a second Penzene
// (Windows opens each double-clicked file in its own) leaves it alone (#318).
QString MainWindow::autosavePath() {
    // With the start time too: a later process can get a crashed one's pid (#361).
    static const QString name =
        QString("/autosave-%1-%2.penz").arg(QCoreApplication::applicationPid()).arg(QDateTime::currentMSecsSinceEpoch());
    return QStandardPaths::writableLocation(QStandardPaths::AppDataLocation) + name;
}

void MainWindow::autosave() {
    if (isClean()) return QFile::remove(autosavePath()), void();
    QDir().mkpath(QFileInfo(autosavePath()).path());
    static QLockFile lock(autosavePath() + ".lock");  // released when this process exits, or goes stale if it crashes
    lock.tryLock(0);
    QSaveFile f(autosavePath());  // a crash mid-write keeps the previous autosave
    if (f.open(QIODevice::WriteOnly)) f.write(sheetsToJson(sheets())), f.commit();
}

// Autosaves left by Penzenes that are no longer running, newest first; the
// first one accepted is recovered, and each one offered is then removed.
void MainWindow::offerRecovery() {
    const QDir dir(QFileInfo(autosavePath()).path());
    for (const QFileInfo& file : dir.entryInfoList({"autosave*.penz"}, QDir::Files, QDir::Time)) {
        if (file.fileName() == QFileInfo(autosavePath()).fileName()) continue;
        QLockFile lock(file.filePath() + ".lock");
        if (!lock.tryLock(0)) continue;  // its Penzene is still running
        QFile f(file.filePath());
        if (!f.open(QIODevice::ReadOnly)) continue;
        auto sheets = sheetsFromJson(f.readAll());
        f.close();
        const bool drawn = std::any_of(sheets.begin(), sheets.end(), [](const Sheet& s) { return !s.doc.empty(); });
        const bool recover =
            drawn && QMessageBox::question(this, tr("Recover"),
                                           tr("Penzene closed without saving your last drawing. Recover it?")) ==
                         QMessageBox::Yes;
        QFile::remove(file.filePath());
        if (!recover) continue;
        const Document first = sheets[0].doc;
        sheets[0].doc = Document{};
        setPages(sheets);
        canvas_->commit(first, tr("Recover"));  // unsaved, so Save asks where to put it
        pagesEdited_ = sheets.size() > 1;
        path_.clear();  // not the file opened at launch, which Save would overwrite (#491)
        updateTitle();
        return;
    }
}

bool MainWindow::saveTo(const QString& path, bool v3000) {
    const auto& doc = canvas_->document();
    QByteArray data;
    if (path.endsWith(".penz", Qt::CaseInsensitive)) {
        data = penz1_ ? sheetsToJsonV1(sheets()) : sheetsToJson(sheets());
    } else if (pages_.size() > 1) {
        QMessageBox::warning(this, tr("Save"),
                             tr("%1 holds one page. Save as a Penzene document (.penz) to keep all %2 pages, or export this page.")
                                 .arg(QFileInfo(path).suffix().toUpper()).arg(pages_.size()));
        return false;
    } else if (path.endsWith(".cdxml", Qt::CaseInsensitive)) {
        data = chem::toCdxml(doc);
    } else if (path.endsWith(".cdx", Qt::CaseInsensitive)) {
        data = chem::toCdx(doc);
        if (data.isEmpty()) {
            QMessageBox::warning(this, tr("Save"), tr("This build can't write binary CDX; save as ChemDraw XML (.cdxml)."));
            return false;
        }
    } else if (path.endsWith(".rxn", Qt::CaseInsensitive) || path.endsWith(".rdf", Qt::CaseInsensitive)) {
        const auto steps = chem::reactionsOf(doc);
        if (steps.empty()) {
            QMessageBox::warning(this, tr("Save"), tr("A reaction file needs a reaction arrow in the drawing."));
            return false;
        }
        if (path.endsWith(".rdf", Qt::CaseInsensitive)) {
            data = QByteArray::fromStdString(chem::toRdf(steps));
        } else if (steps.size() > 1) {  // #334
            QMessageBox::warning(this, tr("Save"),
                                 tr("An Rxnfile holds one step. Save as an RD file (.rdf) to keep all %1.").arg(steps.size()));
            return false;
        } else {
            data = QByteArray::fromStdString(chem::toRxn(steps[0]));
        }
    } else if (path.endsWith(".sdf", Qt::CaseInsensitive)) {  // one record per molecule (#330)
        data = QByteArray::fromStdString(chem::toSdf(doc, v3000));
        if (data.isEmpty()) {
            QMessageBox::warning(this, tr("Save"), tr("An SD file holds molecules, and this page has none."));
            return false;
        }
    } else if (const QString ext = QFileInfo(path).suffix().toLower(); ext.isEmpty() || ext == "mol") {
        data = QByteArray::fromStdString(chem::toMolBlock(doc, v3000));
    } else {  // an image or library path would get MOL text (#315)
        QMessageBox::warning(this, tr("Save"),
                             tr("Penzene can't save a drawing as .%1. Save as .penz, .mol, .sdf, .rxn, .rdf, .cdxml or .cdx, "
                                "or use File → Export… for images.")
                                 .arg(ext));
        return false;
    }
    // Written aside and swapped in, so a failed save never truncates the file already there.
    QSaveFile f(path);
    f.setDirectWriteFallback(true);  // a writable file in a read-only folder (sandbox portals)
    if (!f.open(QIODevice::WriteOnly) || f.write(data) != data.size() || !f.commit()) {
        QMessageBox::warning(this, tr("Save"), tr("Cannot write %1").arg(path));
        return false;
    }
    path_ = path;
    for (auto& p : pages_) p.undo->setClean();
    pagesEdited_ = false;
    QFile::remove(autosavePath());
    remember(path);
    updateTitle();
    return true;
}

bool MainWindow::save() {
#ifdef Q_OS_WIN
    if (embeddedSave_) {  // back into the Word or PowerPoint document (#229)
        if (!embeddedSave_()) return false;
        for (auto& p : pages_) p.undo->setClean();
        pagesEdited_ = false;
        updateTitle();
        return true;
    }
#endif
    // MOL can't hold everything .penz will (text, arrows), so only .penz saves silently.
    return path_.endsWith(".penz", Qt::CaseInsensitive) ? saveTo(path_) : saveAs();
}

bool MainWindow::saveAs() {
    const QString v3000 = tr("MDL Molfile V3000 (*.mol)");
    const QString penz1 = tr("Penzene 1.x (*.penz)");  // version 1, which Penzene 1.x opens (#404)
    QString filter;
    QString path = QFileDialog::getSaveFileName(this, tr("Save As"), path_,
                                                tr("Penzene document (*.penz);;") + penz1 + tr(";;MDL Molfile (*.mol);;") + v3000 +
                                                    tr(";;MDL SD file, one record per molecule (*.sdf);;MDL Rxnfile (*.rxn);;MDL RD file, every reaction step (*.rdf);;ChemDraw XML (*.cdxml);;"
                                                       "ChemDraw, molecules only (*.cdx)"),
                                                &filter);
    if (path.isEmpty()) return false;
    penz1_ = filter == penz1;  // Save keeps writing what Save As chose
    return saveTo(path, filter == v3000);
}

bool MainWindow::maybeSave() {
    if (isClean()) return true;
    auto r = QMessageBox::question(this, tr("Unsaved changes"), tr("Save changes to this document?"),
                                   QMessageBox::Save | QMessageBox::Discard | QMessageBox::Cancel);
    return r == QMessageBox::Discard || (r == QMessageBox::Save && save());
}

static QString droppedFile(const QMimeData* mime) {
    for (const QUrl& url : mime->urls())
        if (url.isLocalFile()) return url.toLocalFile();
    return {};
}

void MainWindow::dragEnterEvent(QDragEnterEvent* e) {
    if (!droppedFile(e->mimeData()).isEmpty()) e->acceptProposedAction();
}

void MainWindow::dropEvent(QDropEvent* e) {
    const QString path = droppedFile(e->mimeData());
    if (path.isEmpty()) return;
    e->acceptProposedAction();
    if (maybeSave()) openFile(path);
}

void MainWindow::closeEvent(QCloseEvent* e) {
    if (!maybeSave()) return e->ignore();
    QSettings s;  // the tool layout, restored at the next launch
    for (size_t i = 0; i < flyouts_.size(); ++i) {
        const FlyoutFrame* f = flyouts_[i];
        s.setValue(QString("toolLayout/%1").arg(i), QVariantList{f->fly->isVisible(), f->fly->pos(), f->cols});
    }
    QFile::remove(autosavePath());  // a deliberate quit: nothing to recover
    e->accept();
}

// Export preferences (Edit > Preferences), used by Export and Copy.
static double exportDpi() { return QSettings().value("exportDpi", 300).toDouble(); }
static QColor exportBackground() {
    return QSettings().value("exportBackground").toString() == "white" ? QColor(Qt::white) : QColor(Qt::transparent);
}
static ExportOptions exportOptions() {
    QSettings s;
    return {exportDpi(), exportBackground(), s.value("exportScale", 100).toDouble() / 100,
            s.value("exportMargin", 0).toDouble()};
}

QStringList MainWindow::languages() {
    QStringList out;
    for (const QString& f : QDir(":/i18n").entryList({"penzene_*.qm"}))
        out << f.mid(QString("penzene_").size()).chopped(QString(".qm").size());
    return out;
}

bool MainWindow::installTranslations(const QString& language) {
    static QTranslator* own = nullptr;
    static QTranslator* qt = nullptr;  // Qt's own dialogs and buttons
    for (QTranslator** t : {&own, &qt})
        if (*t) QCoreApplication::removeTranslator(*t), delete *t, *t = nullptr;
    if (language == "en") return true;
    const QLocale locale = language.isEmpty() ? QLocale() : QLocale(language);
    own = new QTranslator;
    const bool found = language.isEmpty() ? own->load(locale, "penzene", "_", ":/i18n")
                                          : own->load(":/i18n/penzene_" + language);
    if (found) QCoreApplication::installTranslator(own);
    qt = new QTranslator;
    if (qt->load(locale, "qtbase", "_", QLibraryInfo::path(QLibraryInfo::TranslationsPath)))
        QCoreApplication::installTranslator(qt);
    return found;
}

void MainWindow::showPreferences() {
    QDialog dialog(this);
    dialog.setWindowTitle(tr("Preferences"));
    auto* form = new QFormLayout(&dialog);
    auto* themeBox = new QComboBox;
    for (const auto& t : themes()) themeBox->addItem(t.name);
    themeBox->setCurrentText(theme(QSettings().value("theme", "System").toString()).name);
    auto* styleBox = new QComboBox;
    styleBox->setObjectName("defaultStyle");
    for (const auto& s : drawingStyles()) styleBox->addItem(s.name);
    styleBox->setCurrentText(drawingStyle(QSettings().value("defaultStyle").toString()).name);
    auto* dpiBox = new QSpinBox;
    dpiBox->setRange(72, 1200);
    dpiBox->setSingleStep(50);
    dpiBox->setSuffix(tr(" dpi"));
    dpiBox->setValue(int(exportDpi()));
    auto* backgroundBox = new QComboBox;
    backgroundBox->setObjectName("exportBackground");
    backgroundBox->addItems({tr("Clear"), tr("White")});
    backgroundBox->setCurrentIndex(exportBackground().alpha() ? 1 : 0);
    auto* languageBox = new QComboBox;
    languageBox->setObjectName("language");
    languageBox->addItem(tr("System default"), QString());
    languageBox->addItem("English", "en");
    for (const QString& code : languages()) {
        const QString name = QLocale(code).nativeLanguageName();
        languageBox->addItem(name.isEmpty() ? code : name, code);
    }
    languageBox->setCurrentIndex(std::max(0, languageBox->findData(QSettings().value("language").toString())));
    form->addRow(tr("Theme:"), themeBox);
    form->addRow(tr("Language:"), languageBox);
    form->addRow(QString(), new QLabel(tr("A new language takes effect when Penzene restarts.")));
    form->addRow(tr("Drawing style for new documents:"), styleBox);
    form->addRow(tr("PNG resolution:"), dpiBox);
    form->addRow(tr("Export and copy background:"), backgroundBox);
    auto* scaleBox = new QSpinBox;
    scaleBox->setObjectName("exportScale");
    scaleBox->setRange(10, 400);
    scaleBox->setSingleStep(5);
    scaleBox->setSuffix("%");
    scaleBox->setValue(int(exportOptions().scale * 100 + 0.5));
    scaleBox->setToolTip(tr("e.g. 85% for a journal's column width"));
    auto* marginBox = new QSpinBox;
    marginBox->setObjectName("exportMargin");
    marginBox->setRange(0, 72);
    marginBox->setSuffix(tr(" pt"));
    marginBox->setValue(int(exportOptions().margin));
    form->addRow(tr("Export and copy scale:"), scaleBox);
    form->addRow(tr("Margin around exports:"), marginBox);
    auto* checkBefore = new QCheckBox(tr("Check structures before export and copy"));
    checkBefore->setObjectName("checkBeforeExport");
    checkBefore->setChecked(QSettings().value("checkBeforeExport", false).toBool());
    checkBefore->setToolTip(tr("Lists what Check Structure finds (valence errors, stereocentres without a wedge…) first"));
    form->addRow(QString(), checkBefore);
    auto* updates = new QCheckBox(tr("Check for a new version once a week (asks GitHub)"));
    updates->setObjectName("autoUpdates");
    updates->setChecked(QSettings().value("updates/auto", false).toBool());
    form->addRow(tr("Updates:"), updates);
    auto* buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    buttons->button(QDialogButtonBox::Ok)->setDefault(true);
    form->addRow(buttons);
    connect(buttons, &QDialogButtonBox::accepted, &dialog, &QDialog::accept);
    connect(buttons, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);
    if (dialog.exec() != QDialog::Accepted) return;
    QSettings settings;
    settings.setValue("defaultStyle", styleBox->currentIndex() ? styleBox->currentText() : QString());
    settings.setValue("exportDpi", dpiBox->value());
    settings.setValue("exportBackground", backgroundBox->currentIndex() ? "white" : "clear");
    settings.setValue("exportScale", scaleBox->value());
    settings.setValue("exportMargin", marginBox->value());
    settings.setValue("checkBeforeExport", checkBefore->isChecked());
    settings.setValue("updates/auto", updates->isChecked());
    settings.setValue("language", languageBox->currentData().toString());
    applyTheme(themeBox->currentText());
    for (auto* a : themeGroup_->actions()) a->setChecked(a->text() == themeBox->currentText());
}

QWidget* MainWindow::checkStructure() {
    auto* dialog = new QDialog(this);
    dialog->setAttribute(Qt::WA_DeleteOnClose);
    dialog->setWindowTitle(tr("Check Structure"));
    auto* layout = new QVBoxLayout(dialog);
    auto* list = new QListWidget;
    layout->addWidget(list);
    const auto problems = chem::checkStructure(canvas_->document());
    for (const auto& p : problems) {
        auto* item = new QListWidgetItem(p.message, list);
        QVariantList atoms;
        for (int i : p.atoms) atoms << i;
        item->setData(Qt::UserRole, atoms);
    }
    if (problems.empty()) list->addItem(tr("No problems found."));
    connect(list, &QListWidget::currentItemChanged, this, [this](QListWidgetItem* item) {
        if (!item) return;
        QSet<int> atoms;
        for (const QVariant& v : item->data(Qt::UserRole).toList()) atoms.insert(v.toInt());
        canvas_->setSelection(atoms);  // show where the problem is
    });
    auto* close = new QDialogButtonBox(QDialogButtonBox::Close);
    connect(close, &QDialogButtonBox::rejected, dialog, &QDialog::close);
    layout->addWidget(close);
    dialog->resize(460, 260);
    dialog->show();
    return dialog;
}

bool printDocument(QPrinter& printer, const Document& doc) {
    if (doc.empty()) return false;
    const QRectF r = outputBounds(doc);
    QPainter p;
    if (!p.begin(&printer)) return false;
    const QRectF page = printer.pageRect(QPrinter::DevicePixel);
    const double perPoint = printer.resolution() / 72.0;
    const double s = std::min({exportScale(doc) * perPoint, page.width() / r.width(), page.height() / r.height()});
    p.translate(page.width() / 2, page.height() / 2);  // the painter's origin is the printable area's corner
    p.scale(s, s);
    p.translate(-r.center());
    paintDocument(p, doc);
    return p.end();
}

void MainWindow::fillTemplates() {
    templates_->clear();
    const QColor ink = palette().color(QPalette::WindowText);
    auto templateIcon = [&](const Document& doc) {
        return drawingIcon(doc, templates_->iconSize(), ink, devicePixelRatioF());
    };
    QHash<QString, QTreeWidgetItem*> groups;
    auto group = [&](const QString& name) {
        if (!groups.contains(name)) groups[name] = new QTreeWidgetItem(templates_, {name});
        return groups[name];
    };
    const auto mine = userTemplates();
    if (!mine.empty()) {
        for (const auto& [name, doc] : mine) {
            auto* item = new QTreeWidgetItem(group(tr("My templates")), {name});
            item->setIcon(0, templateIcon(doc));
            item->setData(0, Qt::UserRole, doc.toJson());
            item->setData(0, Qt::UserRole + 1, true);
        }
    }
    for (const auto& t : builtinTemplates()) {
        const Document doc = templateDocument(t);
        if (doc.empty()) continue;
        auto* item = new QTreeWidgetItem(group(t.category), {t.name});
        item->setIcon(0, templateIcon(doc));
        item->setData(0, Qt::UserRole, doc.toJson());
    }
    if (!mine.empty()) groups[tr("My templates")]->setExpanded(true);
}

void MainWindow::insertTemplate(QTreeWidgetItem* item) {
    if (!item || item->data(0, Qt::UserRole).isNull()) return;  // a category
    if (auto doc = Document::fromJson(item->data(0, Qt::UserRole).toByteArray())) {
        doc->style = canvas_->document().style;
        canvas_->insert(*doc, tr("Insert %1").arg(item->text(0)));
    }
}

void MainWindow::saveTemplate() {
    const Document doc = canvas_->selectedSubset();
    if (doc.empty()) return;
    bool ok = false;
    const QString name = QInputDialog::getText(this, tr("Save as Template"), tr("Template name:"), QLineEdit::Normal, {}, &ok);
    if (!ok || name.trimmed().isEmpty()) return;
    if (!saveUserTemplate(name, doc)) {
        QMessageBox::warning(this, tr("Save as Template"), tr("Could not save the template."));
        return;
    }
    fillTemplates();
    templateDock_->show();
}

void MainWindow::checkForUpdates(bool quietly) {
    auto* net = new QNetworkAccessManager(this);
    QNetworkRequest request(online::latestReleaseUrl());
    request.setHeader(QNetworkRequest::UserAgentHeader, "Penzene/" PENZENE_VERSION);  // GitHub's API wants one
    request.setTransferTimeout(15000);
    QNetworkReply* reply = net->get(request);
    connect(reply, &QNetworkReply::finished, this, [this, reply, net, quietly] {
        const auto release = online::parseRelease(reply->readAll());
        const QString error = reply->error() ? reply->errorString() : QString();
        net->deleteLater();
        if (release.tag.isEmpty()) {
            if (!quietly) QMessageBox::warning(this, tr("Check for Updates"), tr("Couldn't reach GitHub: %1").arg(error));
            return;
        }
        if (online::isNewer(release.tag, PENZENE_VERSION)) {
            if (QMessageBox::question(this, tr("Check for Updates"),
                                      tr("Penzene %1 is available (you have %2). Open the download page?")
                                          .arg(QString(release.tag).remove(QRegularExpression("^v")), PENZENE_VERSION)) ==
                QMessageBox::Yes)
                QDesktopServices::openUrl(QUrl(release.url));
        } else if (!quietly) {
            QMessageBox::information(this, tr("Check for Updates"), tr("Penzene %1 is the latest version.").arg(PENZENE_VERSION));
        }
    });
}

void MainWindow::showWhatsNew() {
    QFile f(":/CHANGELOG.md");
    showWhatsNewDialog(this, f.open(QIODevice::ReadOnly) ? QString::fromUtf8(f.readAll()) : QString(), PENZENE_VERSION);
}

void MainWindow::maybeShowWhatsNew() {
    QSettings settings;
    const QString last = settings.value("lastVersion").toString();
    settings.setValue("lastVersion", PENZENE_VERSION);
    if (last.isEmpty() || !online::isNewer(PENZENE_VERSION, last)) return;  // fresh install, or not newer
    QFile f(":/CHANGELOG.md");
    if (f.open(QIODevice::ReadOnly) && !online::releaseNotes(QString::fromUtf8(f.readAll()), PENZENE_VERSION).isEmpty())
        showWhatsNew();
}

void MainWindow::maybeCheckForUpdates() {
    QSettings settings;
    if (!settings.value("updates/auto", false).toBool()) return;
    const QDateTime last = settings.value("updates/last").toDateTime();
    if (last.isValid() && last.daysTo(QDateTime::currentDateTime()) < 7) return;
    settings.setValue("updates/last", QDateTime::currentDateTime());
    checkForUpdates(true);
}

void MainWindow::print() {
    QPrinter printer(QPrinter::HighResolution);
    QPrintDialog dialog(&printer, this);
    if (dialog.exec() == QDialog::Accepted && !printDocument(printer, canvas_->selectedSubset()))
        QMessageBox::warning(this, tr("Print"), tr("Nothing to print."));
}

bool MainWindow::confirmStructure(const Document& doc, const QString& title, const QString& proceed) {
    if (!QSettings().value("checkBeforeExport", false).toBool()) return true;
    QStringList messages;
    for (const auto& p : chem::checkStructure(doc)) messages << p.message;
    if (messages.isEmpty()) return true;
    QMessageBox box(QMessageBox::Warning, title, tr("Check Structure found problems in this drawing:"),
                    QMessageBox::Cancel, this);
    box.setInformativeText(messages.join('\n'));
    box.setDefaultButton(box.addButton(proceed, QMessageBox::AcceptRole));
    box.exec();
    return box.buttonRole(box.clickedButton()) == QMessageBox::AcceptRole;
}

void MainWindow::exportImage() {
    if (!confirmStructure(canvas_->selectedSubset(), tr("Export"), tr("Export Anyway"))) return;
    QString base = path_.isEmpty() ? QString("structure") : QFileInfo(path_).completeBaseName();
    QString path = QFileDialog::getSaveFileName(this, tr("Export"), base + ".svg",
                                                tr("SVG (*.svg);;PNG image (*.png);;PDF (*.pdf)"));
    if (path.isEmpty()) return;
    if (!exportDocument(canvas_->selectedSubset(), path, exportOptions()))
        QMessageBox::warning(this, tr("Export"), tr("Nothing to export, or cannot write %1").arg(path));
}

// One CSV row per molecule in the selection, or on the page (the columns are in docs/cli.md).
void MainWindow::exportDescriptors() {
    const QString base = path_.isEmpty() ? QString("structure") : QFileInfo(path_).completeBaseName();
    const QString path = QFileDialog::getSaveFileName(this, tr("Export Descriptors"), base + ".csv", tr("CSV (*.csv)"));
    if (path.isEmpty()) return;
    std::vector<chem::Record> records;
    for (auto& m : chem::molecules(canvas_->selectedSubset()))
        records.push_back({QString("%1-%2").arg(base).arg(records.size() + 1), std::move(m)});
    const QByteArray csv = QByteArray::fromStdString(chem::descriptorsCsv(records));
    if (records.empty() || !writeWhole(path, csv))
        QMessageBox::warning(this, tr("Export Descriptors"), tr("Nothing to export, or cannot write %1").arg(path));
}

void MainWindow::importSmiles() {
    bool ok = false;
    QString s = QInputDialog::getText(this, tr("Import SMILES"), tr("SMILES:"), QLineEdit::Normal, {}, &ok);
    if (!ok || s.trimmed().isEmpty()) return;
    if (auto doc = chem::fromSmiles(s.trimmed().toStdString())) canvas_->insert(*doc, tr("Import SMILES"));
    else QMessageBox::warning(this, tr("Import SMILES"), tr("Not a valid SMILES string."));
}

void MainWindow::importSequence() {
    bool ok = false;
    const QString s = QInputDialog::getText(this, tr("Import Peptide Sequence"),
                                            tr("One-letter (GFLS; lower case for D) or three-letter (Gly-Phe-Leu-Ser, H-Gly-D-Phe-OH):"),
                                            QLineEdit::Normal, {}, &ok);
    if (!ok || s.trimmed().isEmpty()) return;
    if (auto doc = chem::fromSequence(s)) canvas_->insert(*doc, tr("Import %1").arg(s.trimmed()));
    else QMessageBox::warning(this, tr("Import Peptide Sequence"), tr("Not a peptide sequence."));
}

void MainWindow::importName() {
    bool ok = false;
    const QString name = QInputDialog::getText(this, tr("Import Name"),
                                               tr("Compound name (looked up on PubChem, online):"),
                                               QLineEdit::Normal, {}, &ok);
    if (!ok || name.trimmed().isEmpty()) return;
    QString error;
    QApplication::setOverrideCursor(Qt::WaitCursor);
    const QString smiles = pubchem::fetch(pubchem::nameToSmilesUrl(name), "SMILES", &error);
    QApplication::restoreOverrideCursor();
    if (auto doc = chem::fromSmiles(smiles.toStdString()); doc && !smiles.isEmpty())
        canvas_->insert(*doc, tr("Import %1").arg(name.trimmed()));
    else
        QMessageBox::warning(this, tr("Import Name"), tr("Could not look up “%1”: %2").arg(name.trimmed(), error));
}

bool MainWindow::copy() {
    Document doc = canvas_->selectedSubset();
    if (doc.empty() || !confirmStructure(doc, tr("Copy"), tr("Copy Anyway"))) return false;
    auto* mime = new QMimeData;
    const QByteArray penz = doc.toJson();
#ifdef Q_OS_WIN
    // First, so Ctrl+V in Word and PowerPoint can embed an object that opens here (#229): EmbedClipboard
    // offers the drawing as "Embed Source" once the installer has registered Penzene.
    mime->setData(kPenzMime, penz);
    mime->setData(kEmfMime, renderEmf(doc, exportOptions()));  // vector for Word and PowerPoint, ahead of the bitmap
#endif
    // The PNG as exported, so the drawing in its text chunk survives: Qt re-encodes an image
    // it converts itself and drops it. The image stays for apps that only read a bitmap.
    const QByteArray png = renderPng(doc, exportOptions());
    mime->setData("image/png", png);
#ifdef Q_OS_WIN
    mime->setData(kWinPngMime, png);  // Qt offers no "PNG" for an image; elsewhere image/png is enough
#endif
    mime->setImageData(QImage::fromData(png, "PNG"));
    mime->setData("image/svg+xml", renderSvg(doc, exportOptions()));
    mime->setData("application/pdf", renderPdf(doc, exportOptions()));  // vector, for Office and Keynote
    mime->setData(kPenzMime, penz);  // already first on Windows, where this keeps its place
    // For pasting into ChemDraw (macOS maps this to its pasteboard type through ChemDrawPasteboard).
    if (const QByteArray cdx = doc.atoms.empty() ? QByteArray() : chem::toCdx(doc); !cdx.isEmpty())
        mime->setData("chemical/x-cdx", cdx);
    if (!doc.atoms.empty()) {
        // No plain text: Office's ⌘V/Ctrl+V takes text over a picture (Copy as SMILES gives it).
        const std::string mol = chem::toMolBlock(doc);
        mime->setData(kMolMime, QByteArray::fromStdString(mol));
    }
    QApplication::clipboard()->setMimeData(mime);
    return true;
}

#ifdef Q_OS_MACOS
// Qt only shows pasteboard types it has a converter for: ChemDraw's copies
// carry binary CDX under this UTI, and PDF (Penzene's vector copy) has none built in.
ChemDrawPasteboard::ChemDrawPasteboard() = default;
QString ChemDrawPasteboard::mimeForUti(const QString& uti) const {
    if (uti == "com.adobe.pdf") return "application/pdf";
    if (uti == "public.png") return "image/png";
    return uti == "com.perkinelmer.chemdraw.cdx-clipboard" ? "chemical/x-cdx" : QString();
}
QString ChemDrawPasteboard::utiForMime(const QString& mime) const {
    if (mime == "application/pdf") return "com.adobe.pdf";
    if (mime == "image/png") return "public.png";
    return mime == "chemical/x-cdx" ? "com.perkinelmer.chemdraw.cdx-clipboard" : QString();
}
QVariant ChemDrawPasteboard::convertToMime(const QString&, const QList<QByteArray>& data, const QString&) const {
    return data.value(0);
}
QList<QByteArray> ChemDrawPasteboard::convertFromMime(const QString&, const QVariant& data, const QString&) const {
    return {data.toByteArray()};
}
#endif

void MainWindow::paste() {
    const QMimeData* mime = QApplication::clipboard()->mimeData();
    // Our own first: Penzene's copy also offers CDX, which drops brackets and fills.
    if (auto doc = Document::fromJson(mime->data(kPenzMime)); doc && !doc->empty())
        return canvas_->insert(*doc, tr("Paste"));
    // ChemDraw: CDX as chemical/x-cdx (macOS, via ChemDrawPasteboard) or its
    // Windows clipboard format; CDXML where an app offers that.
    for (const QString& type : mime->formats())
        if (type == "chemical/x-cdx" || type == "chemical/x-cdxml" || type.contains("ChemDraw Interchange Format"))
            if (auto doc = chem::fromChemDraw(mime->data(type)); doc && !doc->empty())
                return canvas_->insert(*doc, tr("Paste"));
    // A figure Penzene exported, copied from another app or as a file: the drawing inside it.
    for (const char* type : {"image/svg+xml", "application/pdf", "image/png", kWinPngMime})
        if (auto doc = Document::fromEmbedded(mime->data(type)); doc && !doc->empty())
            return canvas_->insert(*doc, tr("Paste"));
    for (const QUrl& url : mime->urls())
        if (auto doc = url.isLocalFile() ? chem::readFile(url.toLocalFile()) : std::nullopt; doc && !doc->empty())
            return canvas_->insert(*doc, tr("Paste"));
    std::string text = mime->hasFormat(kMolMime) ? mime->data(kMolMime).toStdString()
                                                 : mime->text().trimmed().toStdString();
    if (text.empty()) return;
    auto doc = text.find("M  END") != std::string::npos ? chem::fromMolBlock(text)
               : text.starts_with("InChI=")         ? chem::fromInchi(text)
               : text.find('>') != std::string::npos ? chem::fromReactionSmiles(text)
                                                    : chem::fromSmiles(text);
    if (doc) canvas_->insert(*doc, tr("Paste"));
    else statusBar()->showMessage(tr("Clipboard has no structure or SMILES"), 4000);
}

// Tool icons are drawn with the same renderer as the canvas, in the palette's ink.
// Returned as makers so they can be repainted when the theme changes.
using IconMaker = std::function<QIcon()>;
static IconMaker paintedIcon(std::function<void(QPainter&, QColor)> paint) {
    return [paint] {
    QPixmap pm(48, 48);
    pm.setDevicePixelRatio(2);
    pm.fill(Qt::transparent);
    QPainter p(&pm);
    p.setRenderHint(QPainter::Antialiasing);
    paint(p, QApplication::palette().color(QPalette::WindowText));
    return QIcon(pm);
    };
}

QPixmap inkGlyph(const QPixmap& icon, const QColor& ink) {
    const QImage img = icon.toImage();
    for (int y = 0; y < img.height(); ++y)
        for (int x = 0; x < img.width(); ++x)
            if (const QColor c = img.pixelColor(x, y); c.alpha() > 200 && c.value() > 100)
                return icon;  // coloured or light: it shows on a dark theme as it is
    QPixmap pm = icon;  // same alpha and anti-aliasing, in the ink, like the tool icons
    QPainter p(&pm);
    p.setCompositionMode(QPainter::CompositionMode_SourceIn);
    p.fillRect(pm.rect(), ink);
    return pm;
}

QStyle* themedStyle() {
    // QMessageBox takes its icon from the icon theme before the style, so fix the pixmap it chose.
    struct Style : QProxyStyle {
        using QProxyStyle::polish;
        void polish(QWidget* w) override {
            QProxyStyle::polish(w);
            if (auto* box = qobject_cast<QMessageBox*>(w); box && !box->iconPixmap().isNull())
                box->setIconPixmap(inkGlyph(box->iconPixmap(), box->palette().color(QPalette::WindowText)));
        }
    };
    return new Style;
}

static IconMaker docIcon(const Document& d) {
    return paintedIcon([d](QPainter& p, QColor ink) {
        QRectF r;
        for (const auto& a : d.atoms) r |= QRectF(a.pos, QSizeF(0.01, 0.01));
        for (const auto& a : d.arrows) r |= arrowPath(a).boundingRect().adjusted(-3, -3, 3, 3);
        for (const auto& t : d.texts) r |= textPath(t).boundingRect();
        const double s = 18 / std::max(r.width(), r.height());
        p.translate(12, 12);
        p.scale(s, s);
        p.translate(-r.center());
        paintDocument(p, d, {ink, ink, 1.3 / s});  // constant stroke whatever the scale
    });
}

static Document chainDoc(std::vector<QPointF> pts, int order = 1, BondStereo stereo = BondStereo::None) {
    Document d;
    for (QPointF q : pts) d.addAtom(q * kBondLength);
    for (int i = 1; i < int(pts.size()); ++i) d.bonds.push_back({i - 1, i, i == 1 ? order : 1, i == 1 ? stereo : BondStereo::None});
    return d;
}

static Document ringDoc(int n, bool aromatic) {
    Document d;
    const double r = kBondLength / (2 * std::sin(M_PI / n));
    const double turn = n % 4 == 0 ? M_PI / n : 0;  // squares sit flat, not as diamonds
    for (int k = 0; k < n; ++k)
        d.addAtom(r * QPointF(std::sin(2 * M_PI * k / n + turn), -std::cos(2 * M_PI * k / n + turn)));
    for (int k = 0; k < n; ++k) d.bonds.push_back({k, (k + 1) % n, aromatic && k % 2 == 0 ? 2 : 1});
    return d;
}

static Document arrowDoc(ArrowKind kind, double bend = 0, bool dashed = false, bool crossed = false) {
    Document d;
    const bool area = isShape(kind) && kind != ArrowKind::Line;
    d.arrows.push_back({{0, area ? -5.0 : 0.0}, {16, area ? 5.0 : 0.0}, kind, bend, {}, dashed});
    d.arrows.back().crossed = crossed;
    return d;
}

// An orbital pointing up from the icon's centre.
static Document orbitalDoc(ArrowKind kind, OrbitalLook look) {
    Document d;
    Arrow a{{0, 0}, {0, kind == ArrowKind::SOrbital ? -6.0 : -14.0}, kind};
    a.look = look;
    d.arrows.push_back(a);
    return d;
}

static Document textDoc(const QString& s) {
    Document d;
    d.texts.push_back({{0, 0}, s});
    return d;
}

// Screen readers, and tools driving the accessibility API, toggle a checkable button where a user
// clicks it, and Qt's toggle only flips the check: the tool was never picked (#535). A pressable
// button answers both press and toggle with a real click.
class PressableButton : public QAccessibleWidget {
public:
    explicit PressableButton(QAbstractButton* b) : QAccessibleWidget(b, QAccessible::CheckBox) {}
    QAccessible::State state() const override {
        QAccessible::State s = QAccessibleWidget::state();
        s.checkable = button()->isCheckable();
        s.checked = button()->isChecked();
        return s;
    }
    QString text(QAccessible::Text t) const override {
        const QString s = QAccessibleWidget::text(t);
        return t == QAccessible::Name && s.isEmpty() ? button()->text() : s;
    }
    QStringList actionNames() const override {
        QStringList names{pressAction(), toggleAction()};
        if (auto* t = qobject_cast<QToolButton*>(button()); t && t->menu()) names << showMenuAction();
        return names + QAccessibleWidget::actionNames();
    }
    void doAction(const QString& name) override {
        if (name == pressAction() || name == toggleAction()) button()->click();
        else if (name == showMenuAction()) static_cast<QToolButton*>(button())->showMenu();
        else QAccessibleWidget::doAction(name);
    }

private:
    QAbstractButton* button() const { return static_cast<QAbstractButton*>(widget()); }
};

static void pressable(QAbstractButton* b) {
    static const bool installed = [] {
        QAccessible::installFactory([](const QString&, QObject* o) -> QAccessibleInterface* {
            auto* b = qobject_cast<QAbstractButton*>(o);
            return b && b->property("pressable").toBool() ? new PressableButton(b) : nullptr;
        });
        return true;
    }();
    Q_UNUSED(installed);
    b->setProperty("pressable", true);
}

// Periodic table: main block by group and period, lanthanides and actinides
// underneath. Organic elements are bold, since they're the ones drawn most.
// A drop-down of colour swatches plus "Custom…", for the colour and ring fill tools.
// The current colour's swatch is shown pressed; `names` label the swatches.
static QMenu* colourMenu(QWidget* parent, const QList<QColor>& presets, std::function<QColor()> current,
                         std::function<void(QColor)> picked, const QStringList& names = {}) {
    auto* menu = new QMenu(parent);
    auto* w = new QWidget;
    auto* grid = new QGridLayout(w);
    grid->setSpacing(3);
    grid->setContentsMargins(6, 6, 6, 6);
    for (int i = 0; i < presets.size(); ++i) {
        auto* b = new QToolButton;
        pressable(b);
        b->setFixedSize(24, 24);
        b->setAutoRaise(true);
        b->setToolTip(names.value(i, presets[i].name()));
        b->setAccessibleName(b->toolTip());
        b->setCheckable(true);
        QPixmap swatch(16, 16);
        swatch.fill(presets[i]);
        b->setIcon(QIcon(swatch));
        QObject::connect(b, &QToolButton::clicked, menu, [=] { picked(presets[i]), menu->close(); });
        QObject::connect(menu, &QMenu::aboutToShow, b, [=] { b->setChecked(presets[i] == current()); });
        grid->addWidget(b, i / 4, i % 4);
    }
    auto* custom = new QToolButton;
    custom->setText(QObject::tr("Custom…"));
    custom->setAutoRaise(true);
    QObject::connect(custom, &QToolButton::clicked, menu, [=] {
        menu->close();
        QColor c = QColorDialog::getColor(current(), parent);
        if (c.isValid()) picked(c);
    });
    grid->addWidget(custom, (presets.size() + 3) / 4, 0, 1, 4);
    auto* action = new QWidgetAction(menu);
    action->setDefaultWidget(w);
    menu->addAction(action);
    return menu;
}

namespace {
// Arrow keys move between the buttons of a grid layout, skipping its gaps.
struct GridArrows : QObject {
    QGridLayout* grid;
    explicit GridArrows(QGridLayout* g) : QObject(g), grid(g) {}
    bool eventFilter(QObject* o, QEvent* e) override {
        if (e->type() != QEvent::KeyPress) return false;
        const int key = static_cast<QKeyEvent*>(e)->key();
        const int dr = key == Qt::Key_Down ? 1 : key == Qt::Key_Up ? -1 : 0;
        const int dc = key == Qt::Key_Right ? 1 : key == Qt::Key_Left ? -1 : 0;
        if (!dr && !dc) return false;
        int row, col, rs, cs;
        const int index = grid->indexOf(static_cast<QWidget*>(o));
        if (index < 0) return false;
        grid->getItemPosition(index, &row, &col, &rs, &cs);
        for (int r = row + dr, c = col + dc; r >= 0 && r < grid->rowCount() && c >= 0 && c < grid->columnCount();
             r += dr, c += dc)
            if (auto* item = grid->itemAtPosition(r, c); item && item->widget()) {
                item->widget()->setFocus(Qt::TabFocusReason);
                break;
            }
        return true;
    }
};

}  // namespace

static QWidget* periodicTable(const std::function<void(int)>& picked) {
    auto* w = new QWidget;
    auto* grid = new QGridLayout(w);
    auto* arrows = new GridArrows(grid);
    grid->setSpacing(2);
    grid->setContentsMargins(6, 6, 6, 6);
    auto place = [&](int z, int row, int col) {
        const QString sym = QString::fromStdString(chem::symbol(z));
        auto* b = new QToolButton;
        b->setText(sym);
        b->setToolTip(QString("%1 (%2)").arg(sym).arg(z));
        b->setAccessibleName(QString::fromStdString(chem::elementName(z)));
        b->setFixedSize(30, 26);
        b->setAutoRaise(true);
        b->setFocusPolicy(Qt::StrongFocus);
        b->installEventFilter(arrows);
        static const QSet<int> organic{1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53};
        if (organic.contains(z)) {
            QFont f = b->font();
            f.setBold(true);
            b->setFont(f);
        }
        QObject::connect(b, &QToolButton::clicked, w, [picked, z] { picked(z); });
        grid->addWidget(b, row, col);
    };
    place(1, 0, 0), place(2, 0, 17);
    for (int p = 1, z = 3; p <= 2; ++p) {  // periods 2-3: s block, then p block
        place(z++, p, 0), place(z++, p, 1);
        for (int c = 12; c < 18; ++c) place(z++, p, c);
    }
    for (int p = 3, z = 19; p <= 4; ++p)  // periods 4-5 are full
        for (int c = 0; c < 18; ++c) place(z++, p, c);
    for (int p = 5, z = 55; p <= 6; ++p, z += 32) {  // periods 6-7: Cs/Fr, Ba/Ra, then Hf/Rf onwards
        place(z, p, 0), place(z + 1, p, 1);
        for (int c = 3; c < 18; ++c) place(z + 14 + c, p, c);  // 72 (Hf) at column 3
        for (int k = 0; k < 15; ++k) place(z + 2 + k, p + 3, 2 + k);  // La-Lu, Ac-Lr below
    }
    grid->setRowMinimumHeight(7, 8);  // gap above the f block
    return w;
}

void MainWindow::buildTools() {
    // A rail of tool groups; each group's tools open in a flyout beside it (#220).
    auto* bar = new QToolBar(tr("Tools"), this);
    bar->setObjectName("tools");
    addToolBar(Qt::LeftToolBarArea, bar);
    bar->setMovable(false);
    bar->setFloatable(false);
    auto* card = new QFrame(bar);
    card->setObjectName("toolCard");
    auto* rail = new QVBoxLayout(card);
    rail->setContentsMargins(4, 6, 4, 6);
    rail->setSpacing(3);
    rail->setAlignment(Qt::AlignTop);
    bar->addWidget(card);
    auto* railGroup = new QActionGroup(card);  // exclusive: the group of the current tool is marked
    struct Group {
        QToolButton* railButton;
        QFrame* flyout;
        QGridLayout* grid;
        QAction* last = nullptr;  // the tool a click on the rail button picks again
    };
    auto groups = std::make_shared<std::vector<Group>>();
    auto showFlyout = [this, groups](int i) {
        const Group& g = (*groups)[i];
        if (!g.flyout->isVisible()) {  // an open one stays where it was dragged
            g.flyout->ensurePolished();  // its stylesheet decides its size: measure it after that, not on the first show
            g.flyout->adjustSize();
            QPoint at = g.railButton->mapTo(this, QPoint(g.railButton->width() + 12, 0));
            at.setY(std::min(at.y(), height() - g.flyout->height() - 8));
            // Not over an open one: the first free spot to the right of one, then below one (#408).
            QRect box(at, g.flyout->size());
            std::vector<QRect> open;
            for (auto* f : flyouts_)
                if (f->fly != g.flyout && f->fly->isVisible()) open.push_back(f->fly->geometry());
            auto free = [&](const QRect& r) {
                return rect().contains(r) && std::none_of(open.begin(), open.end(), [&](const QRect& o) { return o.intersects(r); });
            };
            if (!free(box)) {
                std::vector<QPoint> spots;
                for (const QRect& o : open) spots.push_back({o.right() + 9, o.top()});
                for (const QRect& o : open) spots.push_back({at.x(), o.bottom() + 9});
                for (QPoint p : spots)
                    if (free(box.translated(p - at))) {
                        at = p;
                        break;
                    }  // none free: its usual spot
            }
            g.flyout->move(at);
        }
        g.flyout->raise();
        g.flyout->show();
    };
    FlyoutFrame* frame = nullptr;  // the current group's flyout
    QWidget* palette = nullptr;
    auto startGroup = [&](const QString& name, const IconMaker& icon) {
        auto* fly = new QFrame(this);
        fly->setObjectName("toolFlyout");
        fly->hide();
        auto* layout = new QVBoxLayout(fly);
        layout->setContentsMargins(10, 8, 10, 10);
        layout->setSpacing(6);
        auto* head = new QHBoxLayout;
        auto* title = new QLabel(name.toUpper());
        title->setObjectName("flyoutTitle");
        auto* close = new QToolButton;
        close->setObjectName("close");
        close->setFocusPolicy(Qt::StrongFocus);
        close->setText("✕");
        close->setToolTip(tr("Close"));
        close->setAccessibleName(tr("Close %1").arg(name));
        connect(close, &QToolButton::clicked, fly, &QWidget::hide);
        head->addWidget(title);
        head->addStretch();
        head->addWidget(close);
        layout->addLayout(head);
        auto* grid = new QGridLayout;
        grid->setSpacing(2);
        layout->addLayout(grid);
        frame = new FlyoutFrame(fly, grid);
        auto* escape = new QAction(fly);
        escape->setShortcut(Qt::Key_Escape);
        escape->setShortcutContext(Qt::WidgetWithChildrenShortcut);
        connect(escape, &QAction::triggered, fly, &QWidget::hide);
        fly->addAction(escape);
        flyouts_.push_back(frame);

        auto* railAction = new QAction(icon(), name, this);
        railAction->setCheckable(true);
        railGroup->addAction(railAction);
        icons_.push_back({railAction, icon});
        auto* b = new QToolButton;
        pressable(b);
        b->setObjectName("railButton");
        b->setDefaultAction(railAction);
        b->setToolButtonStyle(Qt::ToolButtonTextUnderIcon);
        b->setIconSize({24, 24});
        b->setFocusPolicy(Qt::StrongFocus);
        b->setAccessibleName(name);
        rail->addWidget(b);
        const int index = int(groups->size());
        groups->push_back({b, fly, frame->grid});
        connect(railAction, &QAction::triggered, this, [groups, index, showFlyout] {
            (*groups)[index].railButton->defaultAction()->setChecked(true);
            if (auto* last = (*groups)[index].last) last->trigger();
            showFlyout(index);
        });
        palette = fly;
    };
    auto section = [&] { frame->newRow(); };
    auto* group = new QActionGroup(this);
    auto button = [&](QAction* a) {
        auto* b = new QToolButton;
        pressable(b);
        b->setDefaultAction(a);
        b->setIconSize({26, 26});
        b->setAutoRaise(true);
        b->setFocusPolicy(Qt::StrongFocus);  // Tab reaches every tool; Space picks it
        // Screen readers: the tool's name ("Benzene"), with the whole tip as its description.
        static const QRegularExpression end(R"(\s*(:| \(| —).*$)");
        b->setAccessibleName(QString(a->toolTip()).remove(end));
        b->setAccessibleDescription(a->toolTip());
        frame->add(b);
        const int index = int(groups->size()) - 1;
        if (!(*groups)[index].last) (*groups)[index].last = a;
        connect(a, &QAction::triggered, this, [groups, index, a] {
            Group& g = (*groups)[index];
            g.last = a;
            g.railButton->defaultAction()->setChecked(true);
        });
        return b;
    };
    auto add = [&](const IconMaker& icon, const QString& tip, auto setup) {
        auto* a = new QAction(icon(), {}, this);
        icons_.push_back({a, icon});
        a->setToolTip(tip);
        a->setStatusTip(tip);  // the status bar explains the tool while it's chosen
        a->setCheckable(true);
        group->addAction(a);
        connect(a, &QAction::triggered, this, setup);
        connect(a, &QAction::triggered, this, [this, tip] { statusBar()->showMessage(tip); });
        button(a);
        return a;
    };
    using T = Canvas::Tool;
    auto tool = [this](T t) { return [this, t] { canvas_->setTool(t); }; };
    auto bond = [this](int order) {
        return [this, order] { canvas_->setTool(T::Bond), canvas_->setBondOrder(order); };
    };
    // Keys that pick a tool when no atom or bond is the hotspot (ChemDraw).
    QHash<QString, QAction*> keys;
    const IconMaker select = paintedIcon([](QPainter& p, QColor ink) {
        p.setPen(QPen(ink, 1.2, Qt::DashLine));
        p.drawRect(QRectF(4.5, 5.5, 15, 13));
    });
    startGroup(tr("Select"), select);
    keys[" "] = add(select, tr("Select (drag to move, Alt+drag to rotate, or to lasso from empty space; double-click for fragment) — Space"),
                    tool(T::Select));
    const IconMaker rotate3D = paintedIcon([](QPainter& p, QColor ink) {
        p.setPen(QPen(ink, 1.3));
        p.drawEllipse(QRectF(5, 8, 14, 8));
        p.drawArc(QRectF(3, 3, 18, 18), 40 * 16, 270 * 16);
        p.drawLine(QPointF(17, 5), QPointF(21, 6));
        p.drawLine(QPointF(17, 5), QPointF(19, 9));
    });
    add(rotate3D, tr("Rotate in 3D: drag a selected molecule out of the page; keeps stereochemistry"), tool(T::Rotate3D));
    const IconMaker eraser = paintedIcon([](QPainter& p, QColor ink) {
        p.translate(12, 12);
        p.rotate(-40);
        p.setPen(QPen(ink, 1.3));
        p.drawRoundedRect(QRectF(-8, -4, 16, 8), 1.5, 1.5);
        p.drawLine(QPointF(-2, -4), QPointF(-2, 4));
    });
    add(eraser, tr("Eraser (click an atom, bond, arrow or text)"), tool(T::Erase));
    // Colour tool: paints atoms, bonds, arrows and text. Clicking the swatch opens the
    // colours (CPK first); picking one chooses the tool.
    auto colour = std::make_shared<QColor>(canvas_->colour());
    const IconMaker colourIcon = paintedIcon([colour](QPainter& p, QColor ink) {
        p.setPen(QPen(ink, 1));
        p.setBrush(*colour);
        p.drawRoundedRect(QRectF(5, 5, 14, 14), 3, 3);
    });
    auto* colourTool = add(colourIcon, tr("Colour: click the swatch to pick a colour, then click atoms, bonds, arrows "
                                          "or text to paint them (again to clear)"),
                           tool(T::Colour));
    for (auto* b : palette->findChildren<QToolButton*>())
        if (b->defaultAction() == colourTool) {
            b->setPopupMode(QToolButton::InstantPopup);
            b->setStyleSheet("QToolButton::menu-indicator { image: none; width: 0; }");
            // CPK (Jmol) colours; sulfur darkened from #FFFF30 so it reads on white paper.
            b->setMenu(colourMenu(b,
                                  {Qt::black, QColor(0x30, 0x50, 0xF8), QColor(0xFF, 0x0D, 0x0D), QColor(0xD4, 0xB0, 0x00),
                                   QColor(0xFF, 0x80, 0x00), QColor(0x90, 0xE0, 0x50), QColor(0x1F, 0xF0, 0x1F),
                                   QColor(0xA6, 0x29, 0x29), QColor(0x94, 0x00, 0x94), QColor(0xE0, 0x66, 0x33),
                                   QColor(0x90, 0x90, 0x90), QColor(0xFF, 0xB5, 0xB5)},
                                  [this] { return canvas_->colour(); },
                                  [=, this](QColor c) {
                                      *colour = c;
                                      canvas_->setColour(c);
                                      canvas_->setTool(T::Colour);
                                      colourTool->setChecked(true);
                                      colourTool->setIcon(colourIcon());
                                  },
                                  {tr("Carbon"), tr("Nitrogen"), tr("Oxygen"), tr("Sulfur"), tr("Phosphorus"),
                                   tr("Fluorine"), tr("Chlorine"), tr("Bromine"), tr("Iodine"), tr("Iron"),
                                   tr("Carbon (grey)"), tr("Boron")}));
        }
    const QPointF bondPts[] = {{0, 0}, {0.87, -0.5}};
    auto bondIcon = [&](int order, BondStereo st = BondStereo::None) {
        return docIcon(chainDoc({std::begin(bondPts), std::end(bondPts)}, order, st));
    };
    startGroup(tr("Bonds"), bondIcon(1));
    keys["x"] = add(bondIcon(1), tr("Single bond — x: click empty space or an atom to add a bond; drag to aim it; click a bond to change it"), bond(1));
    keys["x"]->setChecked(true);
    add(bondIcon(2), tr("Double bond: click an atom to add one, or a bond to make it double"), bond(2));
    add(bondIcon(3), tr("Triple bond: click an atom to add one, or a bond to make it triple"), bond(3));
    add(bondIcon(1, BondStereo::Wedge), tr("Wedge bond: points from the atom you start at; click a wedge again to flip it"), tool(T::Wedge));
    add(bondIcon(1, BondStereo::Hash), tr("Hashed bond: points from the atom you start at; click a hash again to flip it"), tool(T::Hash));
    auto styled = [this](int order, BondStereo style) {
        return [this, order, style] { canvas_->setTool(T::Bond), canvas_->setBondOrder(order, style); };
    };
    add(bondIcon(1, BondStereo::Interaction),
        tr("Interaction bond (H-bond, contact, coordination): dotted, not a covalent bond — i on a bond. "
           "Drag between atoms, also of different molecules"),
        styled(1, BondStereo::Interaction));
    add(bondIcon(1, BondStereo::Partial),
        tr("Partial bond, forming or breaking (transition states): dashed, not counted — p on a bond; P for a partial double"),
        styled(1, BondStereo::Partial));
    keys["X"] = add(docIcon(chainDoc({{0, 0}, {0.87, -0.5}, {1.73, 0}, {2.6, -0.5}})), tr("Chain — X: drag to draw a zig-zag chain; it grows with the drag"), tool(T::Chain));

    startGroup(tr("Rings"), docIcon(ringDoc(6, false)));

    auto ring = [this](int n, bool arom) {
        return [this, n, arom] { canvas_->setTool(T::Ring), canvas_->setRing(n, arom); };
    };
    keys["j"] = add(docIcon(ringDoc(6, true)), tr("Benzene — j: click empty space for a ring, an atom to attach one, or a bond to fuse one"), ring(6, true));
    for (int n = 3; n <= 8; ++n) add(docIcon(ringDoc(n, false)), tr("%1-membered ring: click empty space, an atom (spiro/attached) or a bond (fused)").arg(n), ring(n, false));
    // Ring fill: the icon shows the current fill colour; the arrow picks it.
    auto fill = std::make_shared<QColor>(canvas_->fillColor());
    const IconMaker fillIcon = [fill] {
        Document filled = ringDoc(6, false);
        filled.fills.push_back({{0, 1, 2, 3, 4, 5}, *fill});
        return docIcon(filled)();
    };
    auto* fillTool = add(fillIcon, tr("Ring fill: click inside a ring to shade it (again to clear); pick the colour from the arrow"),
                         tool(T::Fill));
    for (auto* b : palette->findChildren<QToolButton*>())
        if (b->defaultAction() == fillTool) {
            b->setPopupMode(QToolButton::MenuButtonPopup);
            // Light tints, so bonds and labels stay readable on top.
            b->setMenu(colourMenu(b,
                                  {QColor(207, 227, 255), QColor(255, 214, 214), QColor(212, 240, 210), QColor(255, 236, 196),
                                   QColor(232, 218, 250), QColor(255, 222, 240), QColor(220, 220, 220), QColor(255, 250, 200)},
                                  [this] { return canvas_->fillColor(); },
                                  [=, this](QColor c) {
                                      *fill = c;
                                      canvas_->setFillColor(c);
                                      canvas_->setTool(T::Fill);
                                      fillTool->setChecked(true);
                                      fillTool->setIcon(fillIcon());
                                  },
                                  {tr("Blue"), tr("Rose"), tr("Green"), tr("Amber"), tr("Lavender"), tr("Pink"),
                                   tr("Grey"), tr("Yellow")}));
        }


    // Element: the button shows the current element and draws it; its arrow
    // opens the periodic table, and picking one switches to the atom tool.
    // Until an element has been picked (ever: it's remembered), the button looks
    // like a small periodic table, so what it opens is obvious.
    auto element = std::make_shared<QString>(QSettings().value("element").toString());
    auto* atom = new QAction(this);
    const IconMaker tableIcon = paintedIcon([](QPainter& p, QColor ink) {
        p.setPen(Qt::NoPen);
        p.setBrush(ink);
        auto cell = [&](int col, int row) { p.drawRect(QRectF(3 + col * 2.6, 6 + row * 2.6, 2, 2)); };
        for (int row = 0; row < 4; ++row) cell(0, row), cell(6, row);  // groups 1 and 18
        for (int row = 1; row < 4; ++row) cell(1, row), cell(4, row), cell(5, row);
        for (int row = 2; row < 4; ++row) cell(2, row), cell(3, row);  // the d block
        for (int col = 1; col < 6; ++col) cell(col, 5);                 // f block underneath
    });
    startGroup(tr("Atoms"), tableIcon);
    const IconMaker atomIcon = [element, tableIcon] {
        return element->isEmpty() ? tableIcon() : docIcon(textDoc(*element))();
    };
    if (!element->isEmpty()) canvas_->setElement(chem::atomicNumber(element->toStdString()));
    atom->setIcon(atomIcon());
    icons_.push_back({atom, atomIcon});
    atom->setToolTip(tr("Atom: click to place or relabel (element from the arrow's periodic table; "
                        "or point at an atom and type N, O, S…)"));
    atom->setCheckable(true);
    group->addAction(atom);
    connect(atom, &QAction::toggled, this, [this, atom](bool on) {
        if (on) statusBar()->showMessage(atom->toolTip());
    });
    connect(atom, &QAction::triggered, this, [this] { canvas_->setTool(T::Atom); });
    auto* atomButton = button(atom);
    frame->tools.back().span = 2;  // two cells wide
    frame->place(frame->cols);
    atomButton->setPopupMode(QToolButton::MenuButtonPopup);
    atomButton->setToolButtonStyle(Qt::ToolButtonIconOnly);
    auto* menu = new QMenu(atomButton);
    auto* table = new QWidgetAction(menu);
    table->setDefaultWidget(periodicTable([=, this](int z) {
        *element = QString::fromStdString(chem::symbol(z));
        QSettings().setValue("element", *element);
        canvas_->setElement(z);
        canvas_->setTool(T::Atom);
        atom->setChecked(true);
        atom->setIcon(atomIcon());
        menu->close();
    }));
    menu->addAction(table);
    atomButton->setMenu(menu);
    auto charge = [](bool plus) {
        return paintedIcon([plus](QPainter& p, QColor ink) {
            p.setPen(QPen(ink, 1.3));
            p.drawEllipse(QPointF(12, 12), 7, 7);
            p.drawLine(QPointF(8.5, 12), QPointF(15.5, 12));
            if (plus) p.drawLine(QPointF(12, 8.5), QPointF(12, 15.5));
        });
    };
    add(charge(true), tr("Positive charge: click an atom to add +1"), tool(T::ChargePlus));
    add(charge(false), tr("Negative charge: click an atom to add −1"), tool(T::ChargeMinus));
    startGroup(tr("Arrows"), docIcon(arrowDoc(ArrowKind::Reaction)));
    auto arrow = [this](ArrowKind k, bool curved, bool dashed = false, bool crossed = false) {
        return [this, k, curved, dashed, crossed] {
            canvas_->setTool(T::Arrow), canvas_->setArrow(k, curved, dashed, OrbitalLook::Outline, crossed);
        };
    };
    const QString drag = tr(" (drag to draw; click an arrow to restyle it)");
    keys["e"] = add(docIcon(arrowDoc(ArrowKind::Reaction)), tr("Reaction arrow — e") + drag,
                    arrow(ArrowKind::Reaction, false));
    add(docIcon(arrowDoc(ArrowKind::Equilibrium)), tr("Equilibrium arrow") + drag, arrow(ArrowKind::Equilibrium, false));
    add(docIcon(arrowDoc(ArrowKind::Resonance)), tr("Resonance arrow") + drag, arrow(ArrowKind::Resonance, false));
    add(docIcon(arrowDoc(ArrowKind::Retro)), tr("Retrosynthesis arrow") + drag, arrow(ArrowKind::Retro, false));
    add(docIcon(arrowDoc(ArrowKind::Reaction, 0, false, true)), tr("No reaction (crossed arrow)") + drag,
        arrow(ArrowKind::Reaction, false, false, true));
    add(docIcon(arrowDoc(ArrowKind::Reaction, 10)), tr("Curved arrow, electron pair (click it again to flip the curve)"),
        arrow(ArrowKind::Reaction, true));
    add(docIcon(arrowDoc(ArrowKind::Fishhook, 10)), tr("Fishhook arrow, single electron (click it again to flip)"),
        arrow(ArrowKind::Fishhook, true));

    startGroup(tr("Shapes"), docIcon(arrowDoc(ArrowKind::RoundedBox)));
    const QString shape = tr(" (drag to draw; Shift for a square or circle; click one to restyle it)");
    add(docIcon(arrowDoc(ArrowKind::Line)), tr("Line (drag to draw)"), arrow(ArrowKind::Line, false));
    add(docIcon(arrowDoc(ArrowKind::Line, 0, true)), tr("Dashed line (drag to draw)"), arrow(ArrowKind::Line, false, true));
    // Solid on the left, dashed on the right.
    add(docIcon(arrowDoc(ArrowKind::RoundedBox)), tr("Rounded box") + shape, arrow(ArrowKind::RoundedBox, false));
    add(docIcon(arrowDoc(ArrowKind::RoundedBox, 0, true)), tr("Dashed rounded box") + shape, arrow(ArrowKind::RoundedBox, false, true));
    add(docIcon(arrowDoc(ArrowKind::Ellipse)), tr("Ellipse") + shape, arrow(ArrowKind::Ellipse, false));
    add(docIcon(arrowDoc(ArrowKind::Ellipse, 0, true)), tr("Dashed ellipse") + shape, arrow(ArrowKind::Ellipse, false, true));
    add(docIcon(arrowDoc(ArrowKind::Box)), tr("Box") + shape, arrow(ArrowKind::Box, false));
    // Orbitals: s, p, lobe and hybrid across, a row for each look; the colour tool colours them.
    section();
    const QString orbital = tr(" (click an atom to centre one on it, drag to point it; click one to restyle it)");
    using OL = OrbitalLook;
    for (auto [look, lookName] : {std::pair{OL::Outline, tr("outline")}, {OL::Shaded, tr("shaded")}, {OL::Gradient, tr("gradient")}})
        for (auto [kind, name] : {std::pair{ArrowKind::SOrbital, tr("s orbital")}, {ArrowKind::POrbital, tr("p orbital")},
                                  {ArrowKind::Lobe, tr("lobe")}, {ArrowKind::HybridOrbital, tr("hybrid orbital")}})
            add(docIcon(orbitalDoc(kind, look)), name + ", " + lookName + orbital,
                [this, kind, look] { canvas_->setTool(T::Arrow), canvas_->setArrow(kind, false, false, look); });
    section();
    keys["t"] = add(docIcon(textDoc("T")), tr("Text (click to add or edit; H2O is set as H₂O) — t"), tool(T::Text));
    (*groups)[1].railButton->defaultAction()->setChecked(true);  // the single bond, chosen at start
    connect(canvas_, &Canvas::toolKey, this, [keys](const QString& k) {
        if (auto* a = keys.value(k)) a->trigger();
    });
}

// "System" follows the OS; the others force light or dark, and Catppuccin
// also recolours the UI. Canvas colours always come from the theme.
void MainWindow::applyTheme(const QString& name) {
    const Theme& chosen = theme(name);
    auto* hints = QGuiApplication::styleHints();
    hints->setColorScheme(chosen.name == "System" ? Qt::ColorScheme::Unknown
                          : chosen.dark           ? Qt::ColorScheme::Dark
                                                  : Qt::ColorScheme::Light);
    const auto scheme = hints->colorScheme();
    const bool dark = chosen.name == "System"
                          ? scheme == Qt::ColorScheme::Dark ||
                                (scheme == Qt::ColorScheme::Unknown &&
                                 QApplication::style()->standardPalette().color(QPalette::Window).lightness() < 128)
                          : chosen.dark;
    const Theme& active = chosen.name == "System" ? theme(dark ? "Dark" : "Light") : chosen;
    QPalette pal = QApplication::style()->standardPalette();
    pal.setColor(QPalette::Window, active.window);
    pal.setColor(QPalette::Button, active.surface);
    for (auto role : {QPalette::WindowText, QPalette::Text, QPalette::ButtonText, QPalette::ToolTipText})
        pal.setColor(role, active.text);
    pal.setColor(QPalette::Base, active.paper);
    pal.setColor(QPalette::AlternateBase, active.surface);
    pal.setColor(QPalette::ToolTipBase, active.surface);
    pal.setColor(QPalette::Highlight, active.accent);
    pal.setColor(QPalette::HighlightedText, active.paper);
    pal.setColor(QPalette::Mid, active.surface);
    pal.setColor(QPalette::Link, active.accent);
    QApplication::setPalette(pal);
    qApp->setStyleSheet(uiStyle(active));
    canvas_->setTheme(active);
    for (auto& [action, make] : icons_) action->setIcon(make());
    paintExamples();
    if (templates_->topLevelItemCount()) fillTemplates();  // repaint the thumbnails in the new ink
    QSettings().setValue("theme", chosen.name);
}

void MainWindow::buildMenus() {
    auto* file = menuBar()->addMenu(tr("&File"));
    file->addAction(tr("&New"), QKeySequence::New, this, [this] {
        if (!maybeSave()) return;
        Document blank;
        blank.style = QSettings().value("defaultStyle").toString();  // Edit > Preferences
        setPages({{tr("Page 1"), blank}});
        path_.clear();
        updateTitle();
        welcome_->show();
    });
    file->addAction(tr("&Open…"), QKeySequence::Open, this, [this] {
        if (!maybeSave()) return;
        QString p = QFileDialog::getOpenFileName(this, tr("Open"), {},
                                                 tr("Structures (*.penz *.mol *.sdf *.smi *.inchi *.rxn *.rdf *.cdxml *.cdx);;Penzene figures (*.svg *.png *.pdf);;All files (*)"));
        if (!p.isEmpty()) openFile(p);
    });
    auto* recent = file->addMenu(tr("Open &Recent"));
    connect(recent, &QMenu::aboutToShow, this, [this, recent] {
        recent->clear();
        for (const QString& p : recentFiles())
            recent->addAction(QFileInfo(p).fileName(), this, [this, p] {
                if (maybeSave()) openFile(p);
            })->setToolTip(p);
        if (recent->isEmpty()) recent->addAction(tr("No recent files"))->setEnabled(false);
        recent->addSeparator();
        recent->addAction(tr("Clear Menu"), this, [] { QSettings().remove("recentFiles"); });
    });
    file->addAction(tr("&Save"), QKeySequence::Save, this, &MainWindow::save);
    file->addAction(tr("Save &As…"), QKeySequence::SaveAs, this, &MainWindow::saveAs);
    file->addSeparator();
    auto* importMenu = file->addMenu(tr("&Import"));
    importMenu->addAction(tr("&SMILES…"), QKeySequence(tr("Ctrl+Shift+I")), this, &MainWindow::importSmiles);
    importMenu->addAction(tr("&Name…"), this, &MainWindow::importName)->setStatusTip(tr("Look a name up on PubChem"));
    importMenu->addAction(tr("&Peptide Sequence…"), this, &MainWindow::importSequence);
    file->addAction(tr("&Export…"), QKeySequence(tr("Ctrl+E")), this, &MainWindow::exportImage);
    file->addAction(tr("Export &Descriptors…"), this, &MainWindow::exportDescriptors);
    file->addAction(tr("&Print…"), QKeySequence::Print, this, &MainWindow::print);
    file->addSeparator();
    file->addAction(tr("&Quit"), QKeySequence::Quit, this, &QWidget::close);

    auto* edit = menuBar()->addMenu(tr("&Edit"));
    auto* u = undoGroup_->createUndoAction(this);
    u->setShortcut(QKeySequence::Undo);
    auto* r = undoGroup_->createRedoAction(this);
    r->setShortcut(QKeySequence::Redo);
    edit->addAction(u);
    edit->addAction(r);
    edit->addSeparator();
    edit->addAction(tr("Cu&t"), QKeySequence::Cut, this, [this] {
        if (copy()) canvas_->deleteSelection();  // a cancelled copy keeps the drawing
    });
    edit->addAction(tr("&Copy"), QKeySequence::Copy, this, &MainWindow::copy);
    auto* copyAs = edit->addMenu(tr("Copy A&s"));
    copyAs->addAction(tr("&SMILES"), QKeySequence(tr("Ctrl+Alt+C")), this, [this] {
        QApplication::clipboard()->setText(QString::fromStdString(chem::toSmiles(canvas_->selectedSubset())));
    });
    copyAs->addAction(tr("&InChI"), this, [this] {
        QApplication::clipboard()->setText(QString::fromStdString(chem::toInchi(canvas_->selectedSubset())));
    });
    copyAs->addAction(tr("InChI&Key"), this, [this] {
        QApplication::clipboard()->setText(QString::fromStdString(chem::toInchiKey(canvas_->selectedSubset())));
    });
    copyAs->addAction(tr("&Reaction SMILES"), this, [this] {
        if (const auto steps = chem::reactionsOf(canvas_->selectedSubset()); !steps.empty()) {
            const std::string smiles = chem::toReactionSmiles(steps);  // a line a step
            if (smiles.empty())
                statusBar()->showMessage(tr("A structure in the reaction is invalid, so nothing was copied"), 6000);
            else
                QApplication::clipboard()->setText(QString::fromStdString(smiles));
        } else
            statusBar()->showMessage(tr("No reaction arrow in the drawing"), 4000);
    });
    copyAs->addAction(tr("IUPAC &Name"), this, [this] {
        const Document doc = canvas_->selectedSubset();
        if (doc.atoms.empty()) return;
        QString error;
        QApplication::setOverrideCursor(Qt::WaitCursor);
        const QString smiles = QString::fromStdString(chem::toSmiles(doc)).section(' ', 0, 0);  // plain SMILES, no CXSMILES extension
        const QString name = pubchem::fetch(pubchem::smilesToNameUrl(), "IUPACName", &error,
                                            pubchem::smilesToNameForm(smiles));
        QApplication::restoreOverrideCursor();
        if (name.isEmpty()) {
            QMessageBox::warning(this, tr("Name from PubChem"),
                                 tr("No name for %1: %2").arg(smiles, error.isEmpty() ? tr("PubChem has no match.") : error));
            return;
        }
        QApplication::clipboard()->setText(name);
        QMessageBox::information(this, tr("Name from PubChem"), tr("%1\n\n(copied to the clipboard)").arg(name));
    });
    edit->addAction(tr("&Paste"), QKeySequence::Paste, this, &MainWindow::paste);
    edit->addAction(tr("&Delete"), canvas_, &Canvas::deleteSelection);
    edit->addSeparator();
    edit->addAction(tr("Select &All"), QKeySequence::SelectAll, canvas_, &Canvas::selectAll);
    edit->addSeparator();
    auto* prefs = edit->addAction(tr("Pr&eferences…"), QKeySequence::Preferences, this, &MainWindow::showPreferences);
    prefs->setMenuRole(QAction::PreferencesRole);  // the app menu on macOS

    auto* page = menuBar()->addMenu(tr("&Page"));
    page->addAction(tr("Ne&w Page"), this, &MainWindow::addPage);
    page->addAction(tr("Rena&me Page…"), this, [this] { renamePage(page_); });
    auto* deletePageAction = page->addAction(tr("&Delete Page"), this, [this] { deletePage(page_); });
    connect(page, &QMenu::aboutToShow, this, [=, this] { deletePageAction->setEnabled(pages_.size() > 1); });
    page->addSeparator();
    auto* moveTo = page->addMenu(tr("Mo&ve Selection To"));
    connect(moveTo, &QMenu::aboutToShow, this, [=, this] {
        moveTo->clear();
        for (int i = 0; i < int(pages_.size()); ++i)
            if (i != page_) moveTo->addAction(pages_[i].name, this, [=, this] { moveSelectionToPage(i); });
        if (moveTo->isEmpty()) moveTo->addAction(tr("(add a page first)"))->setEnabled(false);
    });
    page->addSeparator();
    page->addAction(tr("&Next Page"), QKeySequence(tr("Ctrl+PgDown")), this,
                    [this] { showPage((page_ + 1) % int(pages_.size())); });
    page->addAction(tr("&Previous Page"), QKeySequence(tr("Ctrl+PgUp")), this,
                    [this] { showPage((page_ + int(pages_.size()) - 1) % int(pages_.size())); });
    page->addSeparator();
    auto* pageMenu = page->addMenu(tr("Page &Size"));
    auto* pageGroup = new QActionGroup(pageMenu);
    QStringList pages{""};
    for (const auto& p : pageSizes()) pages << p.name;
    for (const QString& name : pages) {
        auto* a = pageMenu->addAction(name.isEmpty() ? tr("None") : name);
        a->setCheckable(true);
        a->setData(name);
        a->setChecked(name.isEmpty());
        pageGroup->addAction(a);
        connect(a, &QAction::triggered, this, [this, name] {
            Document next = canvas_->document();
            if (next.page == name) return;
            const QRectF was = pageRect(next);
            next.page = name;
            // Keep the page where it was, or centre a new one on the drawing.
            const QSizeF size = pageRect(next).size();
            next.pageOrigin = was.isEmpty() ? documentBounds(next).center() - QPointF(size.width(), size.height()) / 2
                                            : was.topLeft();
            canvas_->commit(next, name.isEmpty() ? tr("No page") : tr("Page: %1").arg(name));
        });
    }
    connect(canvas_, &Canvas::documentChanged, pageGroup, [this, pageGroup] {
        for (auto* a : pageGroup->actions()) a->setChecked(a->data().toString() == canvas_->document().page);
    });

    auto* structure = menuBar()->addMenu(tr("&Structure"));
    structure->addAction(tr("&Clean Structure"), QKeySequence(tr("Ctrl+Shift+K")), this, [this] {
        const auto& sel = canvas_->selection();  // selected molecules only, else everything
        canvas_->commit(chem::clean2D(canvas_->document(), {sel.begin(), sel.end()}), tr("Clean"));
    });
    structure->addAction(tr("Chec&k Structure…"), QKeySequence(tr("Ctrl+Alt+K")), this, [this] { checkStructure(); });
    auto* pubchem = structure->addAction(tr("Look Up on &PubChem"), this, [this] {  // in the browser; Penzene stays offline
        const std::string key = chem::toInchiKey(canvas_->selectedSubset());
        if (key.empty()) return statusBar()->showMessage(tr("No valid structure to look up"), 4000);
        QDesktopServices::openUrl(QUrl("https://pubchem.ncbi.nlm.nih.gov/#query=" + QString::fromStdString(key)));
    });
    connect(structure, &QMenu::aboutToShow, this, [=, this] { pubchem->setEnabled(!canvas_->selectedSubset().atoms.empty()); });
    structure->addSeparator();
    structure->addAction(tr("Add Explicit &Hydrogens"), this, [this] {
        canvas_->commit(chem::addHydrogens(canvas_->document()), tr("Add hydrogens"));
    });
    structure->addAction(tr("Remove Explicit Hydro&gens"), this, [this] {
        canvas_->commit(chem::removeHydrogens(canvas_->document()), tr("Remove hydrogens"));
    });
    structure->addAction(tr("&Expand Abbreviations"), QKeySequence(tr("Ctrl+Shift+E")), canvas_,
                         &Canvas::expandAbbreviations);
    structure->addAction(tr("&Invert Stereochemistry"), this, [this] { canvas_->invertStereo(); });
    structure->addSeparator();
    auto* brackets = structure->addMenu(tr("&Brackets"));
    for (bool square : {true, false})
        brackets->addAction(square ? tr("&Square Brackets Around Selection…") : tr("&Round Brackets Around Selection…"), this,
                            [this, square] {
                                bool ok = false;
                                const QString label = QInputDialog::getText(this, tr("Brackets"), tr("Subscript (e.g. n; may be empty):"),
                                                                            QLineEdit::Normal, "n", &ok);
                                if (ok) canvas_->bracketSelection(square, label.trimmed());
                            });
    brackets->addAction(tr("Remove &Brackets"), this, [this] { canvas_->removeBrackets(); });
    structure->addAction(tr("Save Selection as &Template…"), this, &MainWindow::saveTemplate);

    auto* arrangeMenu = menuBar()->addMenu(tr("&Arrange"));
    arrangeMenu->addAction(tr("Flip &Horizontal"), QKeySequence(tr("Ctrl+Shift+H")), this,
                         [this] { canvas_->flipSelection(true); });
    arrangeMenu->addAction(tr("Flip &Vertical"), QKeySequence(tr("Ctrl+Shift+V")), this,
                         [this] { canvas_->flipSelection(false); });
    arrangeMenu->addAction(tr("&Transform…"), this, [this] {
        QDialog dialog(this);
        dialog.setWindowTitle(tr("Transform"));
        auto* form = new QFormLayout(&dialog);
        auto spin = [&](const QString& label, double lo, double hi, double value, const QString& suffix) {
            auto* s = new QDoubleSpinBox;
            s->setRange(lo, hi), s->setValue(value), s->setSuffix(suffix), s->setDecimals(1);
            form->addRow(label, s);
            return s;
        };
        auto* angle = spin(tr("Rotate:"), -360, 360, 0, "°");
        auto* scale = spin(tr("Scale:"), 5, 1000, 100, "%");
        auto* sx = spin(tr("Stretch across:"), 5, 1000, 100, "%");
        auto* sy = spin(tr("Stretch up and down:"), 5, 1000, 100, "%");
        auto* buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
        form->addRow(buttons);
        connect(buttons, &QDialogButtonBox::accepted, &dialog, &QDialog::accept);
        connect(buttons, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);
        if (dialog.exec() != QDialog::Accepted) return;
        const double k = scale->value() / 100;
        canvas_->transformSelection(QTransform().rotate(angle->value()).scale(k * sx->value() / 100, k * sy->value() / 100),
                                    tr("Transform"));
    });
    arrangeMenu->addSeparator();
    auto* align = arrangeMenu->addMenu(tr("&Align and Distribute"));
    using A = Canvas::Align;
    for (auto [label, edge] : {std::pair{tr("Align &Left"), A::Left}, {tr("Align &Centres"), A::HCentre},
                               {tr("Align &Right"), A::Right}, {tr("Align &Top"), A::Top},
                               {tr("Align &Middles"), A::VCentre}, {tr("Align &Bottom"), A::Bottom}})
        align->addAction(label, this, [this, edge] { canvas_->alignSelection(edge); });
    align->addSeparator();
    align->addAction(tr("Distribute &Horizontally"), this, [this] { canvas_->distributeSelection(true); });
    align->addAction(tr("Distribute &Vertically"), this, [this] { canvas_->distributeSelection(false); });
    arrangeMenu->addAction(tr("Center on &Page"), this, [this] { canvas_->centerOnPage(); })
        ->setStatusTip(tr("Move the selection, or the whole drawing, to the middle of the page"));
    arrangeMenu->addSeparator();
    arrangeMenu->addAction(tr("Arrange &Scheme"), this, [this] { canvas_->arrangeScheme(); });
    arrangeMenu->addAction(tr("&Group"), QKeySequence(tr("Ctrl+G")), canvas_, &Canvas::groupSelection);
    arrangeMenu->addAction(tr("&Ungroup"), QKeySequence(tr("Ctrl+Shift+G")), canvas_, &Canvas::ungroupSelection);
    arrangeMenu->addAction(tr("&Number Compounds"), canvas_, &Canvas::numberCompounds)
        ->setStatusTip(tr("A bold number under each selected molecule, or every one; they renumber in scheme order as you edit"));
    arrangeMenu->addSeparator();
    auto* addArrow = arrangeMenu->addAction(tr("Add A&rrow After Selection"), canvas_, &Canvas::addArrowAfter);
    auto* addText = arrangeMenu->addAction(tr("Add Te&xt After Selection…"), canvas_, &Canvas::addTextAfter);
    for (auto* a : {addArrow, addText}) a->setStatusTip(tr("Just right of the selection, or at the hotspot"));
    auto* flipArrow = arrangeMenu->addAction(tr("&Flip Curved Arrow"), this, [this] { canvas_->bendArrow(-1); });
    flipArrow->setStatusTip(tr("Bow the selected curved arrow the other way; Alt+Up and Alt+Down bend it more or less"));
    connect(arrangeMenu, &QMenu::aboutToShow, this, [=, this] {
        addArrow->setEnabled(bool(canvas_->nextPlace()));
        addText->setEnabled(bool(canvas_->nextPlace()));
        const auto& arrows = canvas_->selectedArrows();
        flipArrow->setEnabled(arrows.size() == 1 && canvas_->document().arrows[*arrows.begin()].bend);
    });

    auto* format = menuBar()->addMenu(tr("F&ormat"));
    // Drawing style presets, like ChemDraw's document settings; stored in the .penz.
    auto* styles = format->addMenu(tr("Drawing &Style"));
    auto* styleGroup = new QActionGroup(this);
    for (const auto& st : drawingStyles()) {
        auto* act = styles->addAction(st.name);
        act->setCheckable(true);
        styleGroup->addAction(act);
        connect(act, &QAction::triggered, this, [this, name = st.name] {
            Document next = canvas_->document();
            next.style = name == drawingStyles()[0].name ? QString() : name;
            next.labelRatio = 0;  // the style's own label size, not an imported file's (#203)
            if (!(next == canvas_->document())) canvas_->commit(next, tr("Drawing style"));
        });
    }
    auto syncStyle = [this, styleGroup] {
        const QString current = drawingStyle(canvas_->document().style).name;
        for (auto* act : styleGroup->actions()) act->setChecked(act->text() == current);
    };
    connect(canvas_, &Canvas::documentChanged, this, syncStyle);
    syncStyle();
    // Display options belong to the document (saved, and in exports), so changing one is an edit.
    auto setDisplay = [this](auto change, const QString& what) {
        Document next = canvas_->document();
        change(next);
        if (!(next == canvas_->document())) canvas_->commit(next, what);
    };
    auto* carbons = format->addMenu(tr("&Carbon Labels"));
    auto* carbonGroup = new QActionGroup(carbons);
    using CL = Document::CarbonLabels;
    for (auto [text, mode] : {std::pair{tr("&None (skeletal)"), CL::None}, {tr("&Terminal CH₃"), CL::Terminal},
                              {tr("&All Carbons"), CL::All}}) {
        auto* a = carbons->addAction(text, this, [=] {
            setDisplay([mode](Document& d) { d.carbonLabels = mode; }, tr("Carbon labels"));
        });
        a->setCheckable(true);
        a->setData(int(mode));
        carbonGroup->addAction(a);
    }
    connect(canvas_, &Canvas::documentChanged, this, [this, carbonGroup] {
        for (auto* a : carbonGroup->actions()) a->setChecked(a->data().toInt() == int(canvas_->document().carbonLabels));
    });
    carbonGroup->actions().first()->setChecked(true);
    format->addSeparator();
    // Stereo labels belong to the document (saved, and in exports), so toggling is an edit.
    auto* stereo = format->addAction(tr("S&tereo Labels"));
    stereo->setCheckable(true);
    connect(stereo, &QAction::toggled, this, [this](bool on) {
        if (canvas_->document().showStereo == on) return;
        Document next = canvas_->document();
        next.showStereo = on;
        canvas_->commit(next, on ? tr("Show Stereo Labels") : tr("Hide Stereo Labels"));
    });
    connect(canvas_, &Canvas::documentChanged, stereo, [this, stereo] {
        QSignalBlocker quiet(stereo);
        stereo->setChecked(canvas_->document().showStereo);
    });
    auto* numbers = format->addAction(tr("Atom &Numbers"));
    numbers->setCheckable(true);
    numbers->setStatusTip(tr("Number every atom; ' on an atom sets its reaction map number"));
    connect(numbers, &QAction::toggled, this, [this](bool on) {
        if (canvas_->document().showAtomNumbers == on) return;
        Document next = canvas_->document();
        next.showAtomNumbers = on;
        canvas_->commit(next, on ? tr("Show atom numbers") : tr("Hide atom numbers"));
    });
    connect(canvas_, &Canvas::documentChanged, numbers, [this, numbers] {
        QSignalBlocker quiet(numbers);
        numbers->setChecked(canvas_->document().showAtomNumbers);
    });
    auto* shifts = format->addAction(tr("Predicted N&MR Shifts"));
    shifts->setCheckable(true);
    shifts->setStatusTip(tr("13C and 1H shifts beside each atom, looked up in nmrshiftdb2 by HOSE code"));
    connect(shifts, &QAction::toggled, this, [this](bool on) {
        if (canvas_->document().showShifts == on) return;
        Document next = canvas_->document();
        next.showShifts = on;
        canvas_->commit(next, on ? tr("Show predicted NMR shifts") : tr("Hide predicted NMR shifts"));
    });
    connect(canvas_, &Canvas::documentChanged, shifts, [this, shifts] {
        QSignalBlocker quiet(shifts);
        shifts->setChecked(canvas_->document().showShifts);
    });
    auto* circles = format->addAction(tr("&Aromatic Circles"));
    circles->setCheckable(true);
    connect(circles, &QAction::toggled, this, [this](bool on) {
        if (canvas_->document().aromaticCircles == on) return;
        Document next = canvas_->document();
        next.aromaticCircles = on;
        next.aromaticCircleOverrides.clear();
        canvas_->commit(next, on ? tr("Aromatic circles") : tr("Kekulé rings"));
    });
    connect(canvas_, &Canvas::documentChanged, circles, [this, circles] {
        QSignalBlocker quiet(circles);
        circles->setChecked(canvas_->document().aromaticCircles);
    });
    auto* selectedCircles = format->addAction(tr("Circles in Selected &Rings"));
    selectedCircles->setStatusTip(tr("Turn aromatic circles on or off in the rings whose atoms are all selected"));
    connect(selectedCircles, &QAction::triggered, this, [this] {
        Document next = canvas_->document();
        const auto& selected = canvas_->selection();
        for (auto ring : chem::aromaticRings(next)) {
            if (!std::all_of(ring.begin(), ring.end(), [&](int i) { return selected.contains(i); })) continue;
            std::sort(ring.begin(), ring.end());
            auto it = std::find(next.aromaticCircleOverrides.begin(), next.aromaticCircleOverrides.end(), ring);
            if (it == next.aromaticCircleOverrides.end()) next.aromaticCircleOverrides.push_back(ring);
            else next.aromaticCircleOverrides.erase(it);
        }
        if (!(next == canvas_->document())) canvas_->commit(next, tr("Toggle aromatic circles"));
    });
    format->addSeparator();
    format->addAction(tr("C&olour Selection"), this, [this] { canvas_->colourSelection(); });
    format->addAction(tr("Ring &Fill Colour…"), this, [this] {
        QColor c = QColorDialog::getColor(canvas_->fillColor(), this, tr("Ring fill colour"));
        if (c.isValid()) canvas_->setFillColor(c);
    });
    auto* heads = format->addMenu(tr("Arrow&head Size"));
    for (auto [name, size] : {std::pair{tr("&Small"), 0.6}, {tr("&Normal"), 1.0}, {tr("&Large"), 1.5}, {tr("&Extra Large"), 2.2}})
        heads->addAction(name, this, [this, size] { canvas_->setArrowHead(size); });
    connect(format, &QMenu::aboutToShow, this, [=, this] { heads->setEnabled(!canvas_->selectedArrows().isEmpty()); });

    auto* view = menuBar()->addMenu(tr("&View"));
    view->addAction(tr("Zoom &In"), QKeySequence::ZoomIn, this, [this] { canvas_->zoomBy(1.25); });
    view->addAction(tr("Zoom &Out"), QKeySequence::ZoomOut, this, [this] { canvas_->zoomBy(0.8); });
    view->addAction(tr("&Fit to Window"), QKeySequence(tr("Ctrl+0")), canvas_, &Canvas::fitToSelection)
        ->setStatusTip(tr("Zoom to the selection, or to the whole drawing"));
    view->addSeparator();
    // Guides are the user's own, not the document's: remembered, never saved or exported.
    auto* grid = view->addAction(tr("&Grid"));
    auto* rulers = view->addAction(tr("&Rulers"));
    for (auto* a : {grid, rulers}) {
        a->setCheckable(true);
        a->setChecked(QSettings().value(a == grid ? "showGrid" : "showRulers", true).toBool());
        connect(a, &QAction::toggled, this, [=, this] {
            QSettings().setValue("showGrid", grid->isChecked());
            QSettings().setValue("showRulers", rulers->isChecked());
            canvas_->setGuides(grid->isChecked(), rulers->isChecked());
        });
    }
    canvas_->setGuides(grid->isChecked(), rulers->isChecked());
    view->addAction(tr("Reset Tool &Layout"), this, [this] {
        for (auto* f : flyouts_) {
            f->fly->hide();
            f->place(FlyoutFrame::kDefaultCols);
        }
        QSettings().remove("toolLayout");
    })->setStatusTip(tr("Close the tool flyouts and put them back at their usual size and place"));
    view->addSeparator();
    auto* templatesToggle = templateDock_->toggleViewAction();
    templatesToggle->setText(tr("&Templates"));
    templatesToggle->setShortcut(QKeySequence(tr("Ctrl+Shift+T")));
    view->addAction(templatesToggle);
    auto* panelToggle = profileDock_->toggleViewAction();
    panelToggle->setText(tr("&Properties Panel"));
    panelToggle->setShortcut(QKeySequence(tr("Ctrl+I")));
    view->addAction(panelToggle);
    auto* massToggle = massDock_->toggleViewAction();
    massToggle->setText(tr("&Mass Spec Panel"));
    view->addAction(massToggle);
    auto* nmrToggle = nmrDock_->toggleViewAction();
    nmrToggle->setText(tr("&NMR Panel"));
    view->addAction(nmrToggle);
    view->addSeparator();
    auto* themeMenu = view->addMenu(tr("T&heme"));
    auto* themeGroup = themeGroup_ = new QActionGroup(this);
    const QString current = QSettings().value("theme", "System").toString();
    for (const auto& t : themes()) {
        auto* a = themeMenu->addAction(t.name, this, [this, n = t.name] { applyTheme(n); });
        a->setCheckable(true);
        a->setChecked(t.name == theme(current).name);
        themeGroup->addAction(a);
        if (t.name == "Dark") themeMenu->addSeparator();
    }
    applyTheme(current);
    // Following the OS: repaint the canvas and icons when it switches.
    connect(QGuiApplication::styleHints(), &QStyleHints::colorSchemeChanged, this, [this] {
        if (QSettings().value("theme", "System").toString() == "System") applyTheme("System");
    });

    auto* help = menuBar()->addMenu(tr("&Help"));
    help->addAction(tr("&Documentation"), this, [] { QDesktopServices::openUrl(QUrl(kDocsUrl)); });
    shortcutsAction_ = help->addAction(tr("&Keyboard Shortcuts"), QKeySequence(tr("F1")), this, [this] {
        QMessageBox box(this);
        box.setWindowTitle(tr("Keyboard Shortcuts"));
        box.setTextFormat(Qt::RichText);
        box.setText(tr(R"(<p>Point at an atom or bond to make it the <b>hotspot</b>. It stays put when the mouse
moves off, so you can keep typing.</p>
<table cellspacing="5">
<tr><th colspan="2" align="left">Moving the hotspot</th></tr>
<tr><td><b>←↑→↓</b></td><td>atom → bond → atom; with <b>Shift</b>: atom → atom, bond → bond (with a selection: nudge it, below)</td></tr>
<tr><td><b>Space</b> / <b>g</b></td><td>select the hotspot's molecule / just its atom or bond</td></tr>
<tr><td><b>G</b></td><td>add the hotspot's atom or bond to the selection; the arrow keys go on to the next (Esc when done)</td></tr>
<tr><td><b>&gt;</b> … <b>&gt;</b></td><td>a curved arrow from the first hotspot to the second, selected</td></tr>
<tr><td><b>Alt+↑</b> / <b>Alt+↓</b></td><td>bend the selected curved arrow more / less (Arrange → Flip Curved Arrow turns it over)</td></tr>
<tr><td><b>Esc</b></td><td>clear hotspot and selection</td></tr>
<tr><th colspan="2" align="left">Atom: sprout</th></tr>
<tr><td><b>1</b> / <b>0</b></td><td>single bond, linear / cyclic mode (0 is longer on 2°/3° carbons)</td></tr>
<tr><td><b>2</b></td><td>acetyl (1°), C=O (2°), CH<sub>2</sub>-acetyl (3°/aromatic)</td></tr>
<tr><td><b>3</b> or <b>a</b></td><td>phenyl</td></tr>
<tr><td><b>4</b> / <b>5</b></td><td>wedged / hashed methyl</td></tr>
<tr><td><b>6 7 u v</b></td><td>cyclohexane, cyclopentane, cyclobutane, cyclopropane (spiro on 2°)</td></tr>
<tr><td><b>8 9 z</b></td><td>methylidene, dimethyl / gem-dimethyl / isopropyl, alkyne</td></tr>
<tr><td><b>k K</b></td><td>sulfonyl, t-Bu</td></tr>
<tr><td><b>.</b> / <b>j</b> / <b>J</b></td><td>attachment point / η⁵-cyclopentadienyl / η⁶-benzene</td></tr>
<tr><th colspan="2" align="left">Atom: label and marks</th></tr>
<tr><td><b>c n/w o/q s p f l b i h d</b></td><td>C N O S P F Cl Br I H D (deuterium)</td></tr>
<tr><td><b>B S L</b></td><td>B, Si, Li</td></tr>
<tr><td><b>m e P A</b></td><td>Me, Et, Ph, Ac</td></tr>
<tr><td><b>O N F E Z</b></td><td>OMe, NO<sub>2</sub>, CF<sub>3</sub>, CO<sub>2</sub>Me, N<sub>3</sub></td></tr>
<tr><td><b>y/Y H Q M</b></td><td>Boc, Cbz, Fmoc, MgBr</td></tr>
<tr><td><b>x r</b></td><td>X, R (generic atoms)</td></tr>
<tr><td><b>+ −</b></td><td>charge</td></tr>
<tr><td><b>:</b> / <b>*</b></td><td>lone pairs (0–3) / radical dot</td></tr>
<tr><td><b>'</b></td><td>atom-map number (next free, or off)</td></tr>
<tr><td><b>Enter</b>, <b>=</b> or <b>t</b></td><td>type a label: element, abbreviation (OMe, Boc, TBS…), SMILES or any text</td></tr>
<tr><td><b>/</b> or <b>?</b></td><td>atom properties</td></tr>
<tr><td><b>Delete</b></td><td>remove label (C stays), or delete a carbon</td></tr>
<tr><th colspan="2" align="left">Bond</th></tr>
<tr><td><b>1 2 3</b></td><td>single, double, triple; <b>2</b> on a double bond swaps the side of its second line</td></tr>
<tr><td><b>w</b> / <b>h</b>, <b>W</b>, <b>H</b></td><td>wedged / hashed (press again to flip)</td></tr>
<tr><td><b>a z</b></td><td>fuse benzene / cyclopentadiene</td></tr>
<tr><td><b>v 4–8</b></td><td>fuse ring of that size (v = 3)</td></tr>
<tr><td><b>9</b> / <b>0</b></td><td>fuse chair cyclohexane (two orientations)</td></tr>
<tr><td><b>d b y</b></td><td>dashed, bold, wavy</td></tr>
<tr><td><b>D</b> / <b>B</b></td><td>dashed double / bold double</td></tr>
<tr><td><b>i</b></td><td>interaction: H-bond or contact, dotted, not a bond</td></tr>
<tr><td><b>p</b> / <b>P</b></td><td>partial bond forming or breaking / partial double (transition states)</td></tr>
<tr><td><b>l c r</b></td><td>double bond's second line left / centred / right</td></tr>
<tr><td><b>f</b></td><td>bring to front: bonds it crosses get a gap</td></tr>
<tr><th colspan="2" align="left">No hotspot (Esc)</th></tr>
<tr><td><b>x X j e t Space</b></td><td>bond, chain, benzene, arrow, text, select tool</td></tr>
<tr><th colspan="2" align="left">Selection</th></tr>
<tr><td><b>←↑→↓</b></td><td>nudge 1 pt; with <b>Shift</b> 10 pt</td></tr>
<tr><td><b>Enter</b></td><td>back to a hotspot on the selection</td></tr>
<tr><td><b>Drag onto an atom</b></td><td>merge (Select tool) &nbsp;•&nbsp; <b>Shift+drag</b> move straight; draw a bond at any angle</td></tr>
<tr><td><b>Ctrl+←↑→↓</b></td><td>duplicate across the next arrow that way (or alongside)</td></tr>
<tr><td><b>Alt+← →</b></td><td>rotate 15° &nbsp;•&nbsp; <b>Alt+drag</b> rotate freely, or lasso from empty space • <b>double-click</b> select fragment, or edit text</td></tr>
<tr><td><b>Ctrl+0</b></td><td>zoom to the selection (to everything with none)</td></tr>
<tr><td><b>Shift+Alt+←↑→↓</b></td><td>rotate 15° out of the page (3D), keeping stereo &nbsp;•&nbsp; choose the Rotate in 3D tool from Select and drag freely</td></tr>
</table>)"));
        box.exec();
    });
    help->addAction(tr("&What's New"), this, &MainWindow::showWhatsNew);
    help->addAction(tr("Check for &Updates…"), this, [this] { checkForUpdates(false); });
    help->addAction(tr("&About Penzene"), this, [this] {
        QMessageBox about(this);
        about.setWindowTitle(tr("About Penzene"));
        about.setText(tr("<h3>Penzene %1</h3><p>An open-source chemical structure editor.</p>"
                         "<p>GPL-3.0 • <a href='https://github.com/JamesOBrien2/penzene'>GitHub</a></p>"
                         "<p>Chemistry by RDKit. GUI by Qt.</p>").arg(PENZENE_BUILD));
        auto* licenses = about.addButton(tr("Third-party licenses…"), QMessageBox::ActionRole);
        about.addButton(QMessageBox::Close);
        about.exec();
        if (about.clickedButton() != licenses) return;

        QFile file(":/THIRD_PARTY_LICENSES.md");
        if (!file.open(QIODevice::ReadOnly)) return;
        QDialog dialog(this);
        dialog.setWindowTitle(tr("Third-party licenses"));
        auto* layout = new QVBoxLayout(&dialog);
        auto* text = new QTextBrowser;
        text->setMarkdown(QString::fromUtf8(file.readAll()));
        text->setOpenExternalLinks(true);
        layout->addWidget(text);
        auto* close = new QDialogButtonBox(QDialogButtonBox::Close);
        connect(close, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);
        layout->addWidget(close);
        dialog.resize(700, 550);
        dialog.exec();
    });
}

#ifdef Q_OS_WIN
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>

bool EmfClipboard::canConvertFromMime(const FORMATETC& format, const QMimeData* mime) const {
    return format.cfFormat == CF_ENHMETAFILE && (format.tymed & TYMED_ENHMF) && mime->hasFormat(kEmfMime);
}

bool EmfClipboard::convertFromMime(const FORMATETC& format, const QMimeData* mime, STGMEDIUM* medium) const {
    if (!canConvertFromMime(format, mime)) return false;
    const QByteArray emf = mime->data(kEmfMime);
    HENHMETAFILE handle = SetEnhMetaFileBits(UINT(emf.size()), reinterpret_cast<const BYTE*>(emf.constData()));
    if (!handle) return false;
    medium->tymed = TYMED_ENHMF;  // a fresh handle each time: the receiver frees it
    medium->hEnhMetaFile = handle;
    medium->pUnkForRelease = nullptr;
    return true;
}

QList<FORMATETC> EmfClipboard::formatsForMime(const QString& type, const QMimeData*) const {
    if (type != kEmfMime) return {};
    return {FORMATETC{CF_ENHMETAFILE, nullptr, DVASPECT_CONTENT, -1, TYMED_ENHMF}};
}

void MainWindow::editEmbedded(const std::vector<Sheet>& sheets, std::function<bool()> save) {
    setPages(sheets.empty() ? std::vector<Sheet>{{tr("Page 1"), Document{}}} : sheets);
    canvas_->fitToDocument();
    path_.clear();
    embeddedSave_ = std::move(save);
    updateTitle();
}

void MainWindow::setEmbeddedIn(const QString& document) {
    embeddedIn_ = document;
    updateTitle();
}

QByteArray MainWindow::embeddedPicture() const { return renderEmf(sheets()[0].doc, exportOptions()); }
#endif
