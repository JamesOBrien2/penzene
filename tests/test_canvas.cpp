#include "Canvas.h"
#include "Chem.h"
#include "Edit.h"
#include "Geometry.h"
#include "MainWindow.h"
#include "NumberBox.h"
#include "Online.h"
#include "WhatsNew.h"
#include "Templates.h"
#include "Stoichiometry.h"
#include "Render.h"
#include "Styles.h"
#include "StyleEditor.h"

#include <QApplication>
#include <QDoubleSpinBox>
#include <QTableWidget>
#include <QNativeGestureEvent>
#include <QSettings>
#include <QStatusBar>
#include <QFileDialog>
#include <QFileOpenEvent>

#include <QLockFile>
#include <QMessageBox>
#include <QPainter>
#include <QStyle>
#include <memory>
#include <QStandardPaths>
#include <QTemporaryDir>
#include <QListWidget>
#include <QToolBar>

#include <QLabel>
#include <QDockWidget>
#include <QPushButton>
#include <QFrame>
#include <QStackedWidget>
#include <QElapsedTimer>
#include <QEventLoop>
#include <QMenuBar>
#include <QAccessible>
#include <QTest>

#include <QComboBox>
#include <QDoubleSpinBox>
#include <QEventLoop>
#include <QDialog>
#include <QSpinBox>
#include <QTest>
#include <QTemporaryDir>
#include <QLineEdit>
#include <QTcpServer>
#include <QTcpSocket>
#include <QInputDialog>
#include <QTimer>
#include <QCheckBox>
#include <QTreeWidget>
#include <QFileInfo>
#include <QPrinter>
#include <QClipboard>
#include <QMimeData>
#include <qpa/qwindowsysteminterface.h>
#include <QMenu>
#include <QTabBar>
#include <QToolButton>
#include <QTextBrowser>
#include <QWidgetAction>
#include <QUndoGroup>
#include <QUndoStack>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QScrollArea>
#include <QScrollBar>
#include <QScreen>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/reporters/catch_reporter_event_listener.hpp>
#include <catch2/reporters/catch_reporter_registrars.hpp>

// Drives the real canvas with synthetic mouse events (QT_QPA_PLATFORM=offscreen).
struct App {  // base class so the QApplication exists before any widget member
    App() {
        static int argc = 1;
        static char name[] = "tests";
        static char* argv[] = {name, nullptr};
        // Match the app's style; Qt's native caches outlive local static destructors.
        static QApplication* app = [] {
            auto* app = new QApplication(argc, argv);
            app->setStyle(themedStyle());
            return app;
        }();
        QCoreApplication::sendPostedEvents(nullptr, QEvent::DeferredDelete);
        qApp->setStyleSheet({});  // New fixtures need fresh caches after native widgets have been destroyed.
        QApplication::setOrganizationName("penzene-tests");  // keep the user's settings out of it
        QStandardPaths::setTestModeEnabled(true);  // and their app data (autosave)
        if (const auto theme = qgetenv("PENZENE_PREVIEW_THEME"); !theme.isEmpty()) QSettings().setValue("theme", QString::fromUtf8(theme));
    }
};

// Catch2 runs the tests in a random order, and rendering needs the application's fonts,
// so it exists before the first test rather than when some test first asks for it.
struct AppFirst : Catch::EventListenerBase {
    using EventListenerBase::EventListenerBase;
    void testRunStarting(const Catch::TestRunInfo&) override { App(); }
};
CATCH_REGISTER_LISTENER(AppFirst)

struct Fixture : App {
    Fixture() {
        canvas.resize(800, 600);
        canvas.show();
    }
    QPoint at(QPointF scene) { return canvas.mapFromScene(scene); }
    void click(QPointF scene) { QTest::mouseClick(canvas.viewport(), Qt::LeftButton, {}, at(scene)); }
    void drag(QPointF from, QPointF to) {
        QTest::mousePress(canvas.viewport(), Qt::LeftButton, {}, at(from));
        QMouseEvent move(QEvent::MouseMove, at(to), canvas.viewport()->mapToGlobal(at(to)), Qt::NoButton,
                         Qt::LeftButton, {});
        QApplication::sendEvent(canvas.viewport(), &move);
        QTest::mouseRelease(canvas.viewport(), Qt::LeftButton, {}, at(to));
    }
    void hover(QPointF scene) {
        QMouseEvent move(QEvent::MouseMove, at(scene), canvas.viewport()->mapToGlobal(at(scene)), Qt::NoButton,
                         Qt::NoButton, {});
        QApplication::sendEvent(canvas.viewport(), &move);
    }
    void key(const QString& k) { QTest::keyClicks(canvas.viewport(), k); }
    const Document& doc() const { return canvas.document(); }
    int doubles() const {
        int n = 0;
        for (auto& b : canvas.document().bonds) n += b.order == 2;
        return n;
    }
    QUndoStack undo;
    Canvas canvas{&undo};
};

TEST_CASE("stoichiometry calculates mmol, mg and liquid volumes and copies a table (#486)") {
    App app;
    MainWindow w;
    w.resize(1200, 900);
    w.show();
    const Document doc = *chem::fromReactionSmiles("CCO.CC(=O)O>O>CCOC(C)=O");
    auto* canvas = w.findChild<Canvas*>();
    canvas->setDocumentSilently(doc);
    auto* dock = w.findChild<QDockWidget*>("stoichiometry");
    REQUIRE(dock);
    dock->show();
    QApplication::processEvents();
    auto* table = dock->findChild<QTableWidget*>("stoichTable");
    auto* scale = dock->findChild<QDoubleSpinBox*>("stoichScale");
    auto* limiting = dock->findChild<QComboBox*>("stoichLimiting");
    auto* copy = dock->findChild<QPushButton*>("stoichCopy");
    REQUIRE(table->rowCount() == 3);
    CHECK(limiting->count() == 3);
    const QString preview = qEnvironmentVariable("PENZENE_STOICHIOMETRY_PREVIEW");
    if (!preview.isEmpty()) REQUIRE(w.grab().save(preview + "/before.png"));
    scale->setValue(2.5);
    qobject_cast<QDoubleSpinBox*>(table->cellWidget(1, 3))->setValue(1.5);
    qobject_cast<QDoubleSpinBox*>(table->cellWidget(0, 6))->setValue(0.789);
    qobject_cast<QDoubleSpinBox*>(table->cellWidget(2, 3))->setValue(2);
    qobject_cast<QDoubleSpinBox*>(table->cellWidget(2, 6))->setValue(0.997);
    CHECK(table->item(0, 4)->text().toDouble() == Catch::Approx(2.5));
    CHECK(table->item(0, 5)->text().toDouble() == Catch::Approx(115.1725).margin(0.01));
    CHECK(table->item(0, 7)->text().toDouble() == Catch::Approx(145.973).margin(0.02));
    CHECK(table->item(1, 4)->text().toDouble() == Catch::Approx(3.75));
    CHECK(table->item(1, 5)->text().toDouble() == Catch::Approx(225.195).margin(0.01));
    CHECK(table->item(1, 7)->text() == "—");
    CHECK(table->item(2, 7)->text().toDouble() == Catch::Approx(90.346).margin(0.02));
    CHECK_FALSE(table->cellWidget(0, 3)->isEnabled());
    copy->click();
    CHECK(QApplication::clipboard()->text().contains("Mass (mg)"));
    CHECK(QApplication::clipboard()->text().contains("Volume (µL)"));
    CHECK(QApplication::clipboard()->text().contains("C2H6O\t46.069\t1\t2.5\t115.1725"));
    CHECK(canvas->document() == doc);
    canvas->fitToDocument();
    QApplication::processEvents();
    if (!preview.isEmpty()) REQUIRE(w.grab().save(preview + "/after.png"));
    limiting->setCurrentIndex(1);
    CHECK(table->cellWidget(0, 3)->isEnabled());
    CHECK_FALSE(table->cellWidget(1, 3)->isEnabled());
    CHECK(table->item(1, 4)->text().toDouble() == Catch::Approx(2.5));
}

TEST_CASE("stoichiometry retains page and reaction inputs but refuses unknown formulas (#486)") {
    App app;
    StoichiometryPanel panel;
    QUndoStack first, second;
    const Document a = *chem::fromReactionSmiles("CCO>O>CC=O");
    const Document b = *chem::fromReactionSmiles("CCN>N>CC=N");
    Document scheme = a;
    scheme.append(b, {0, 150});
    panel.setDocument(scheme, &first);
    auto* reaction = panel.findChild<QComboBox*>("stoichReaction");
    auto* scale = panel.findChild<QDoubleSpinBox*>("stoichScale");
    auto* table = panel.findChild<QTableWidget*>("stoichTable");
    auto* copy = panel.findChild<QPushButton*>("stoichCopy");
    REQUIRE(reaction->count() == 2);
    scale->setValue(2.5);
    reaction->setCurrentIndex(1);
    scale->setValue(4);
    reaction->setCurrentIndex(0);
    CHECK(scale->value() == Catch::Approx(2.5));
    Document swapped = b;
    swapped.append(a, {0, 150});
    panel.setDocument(swapped, &first);
    CHECK(reaction->currentIndex() == 1);
    CHECK(scale->value() == Catch::Approx(2.5));
    reaction->setCurrentIndex(0);
    CHECK(scale->value() == Catch::Approx(4));
    swapped.atoms[b.atoms.size() - 1].z = 8;  // edit B's product while keeping B's starting materials
    panel.setDocument(swapped, &first);
    CHECK(scale->value() == Catch::Approx(4));
    reaction->setCurrentIndex(1);
    Document moved = scheme;
    for (auto& a : moved.atoms) a.pos += QPointF(40, 20);
    for (auto& a : moved.arrows) a.from += QPointF(40, 20), a.to += QPointF(40, 20);
    panel.setDocument(moved, &first);
    CHECK(scale->value() == Catch::Approx(2.5));
    panel.setDocument(scheme, &second);
    CHECK(scale->value() == Catch::Approx(1));
    panel.setDocument(scheme, &first);
    CHECK(scale->value() == Catch::Approx(2.5));
    Document unknown = *chem::fromReactionSmiles("CCO>*>CC=O");
    panel.setDocument(unknown, &first);
    REQUIRE(table->rowCount() == 2);
    CHECK(table->item(1, 1)->text() == "Unknown structure");
    CHECK(table->item(1, 5)->text() == "—");
    CHECK_FALSE(table->cellWidget(1, 3)->isEnabled());
    panel.setDocument(Document{}, &first);
    CHECK(table->rowCount() == 0);
    CHECK_FALSE(scale->isEnabled());
    CHECK_FALSE(copy->isEnabled());
    panel.reset();
    panel.setDocument(scheme, &first);
    CHECK(scale->value() == Catch::Approx(1));
    auto* disposable = new QUndoStack;
    panel.setDocument(scheme, disposable);
    delete disposable;
    CHECK(table->rowCount() == 0);
    CHECK_FALSE(copy->isEnabled());
}

TEST_CASE("benzene then fused ring gives naphthalene") {
    Fixture f;
    f.canvas.setTool(Canvas::Tool::Ring);
    f.canvas.setRing(6, true);
    f.click({0, 0});
    REQUIRE(f.canvas.document().atoms.size() == 6);
    CHECK(f.canvas.document().bonds.size() == 6);
    CHECK(f.doubles() == 3);

    // Click the midpoint of a bond to fuse.
    const auto& d = f.canvas.document();
    QPointF mid = (d.atoms[d.bonds[1].a].pos + d.atoms[d.bonds[1].b].pos) / 2;
    f.click(mid);
    CHECK(f.canvas.document().atoms.size() == 10);
    CHECK(f.canvas.document().bonds.size() == 11);
    CHECK(f.doubles() == 5);

    f.undo.undo();
    CHECK(f.canvas.document().atoms.size() == 6);
    f.undo.redo();
    CHECK(f.canvas.document().atoms.size() == 10);
}

TEST_CASE("bond click, drag, cycle, chain and erase") {
    Fixture f;
    f.canvas.setTool(Canvas::Tool::Bond);
    f.click({0, 0});  // empty: new bond
    REQUIRE(f.canvas.document().bonds.size() == 1);
    f.click(f.canvas.document().atoms[1].pos);  // atom: grow zig-zag
    CHECK(f.canvas.document().bonds.size() == 2);

    auto& d = f.canvas.document();
    f.click((d.atoms[0].pos + d.atoms[1].pos) / 2);  // bond: cycle order
    CHECK(f.canvas.document().bonds[0].order == 2);

    f.canvas.setTool(Canvas::Tool::Chain);
    f.drag({0, 100}, {60, 100});
    CHECK(f.canvas.document().bonds.size() == 2 + 5);

    f.canvas.setTool(Canvas::Tool::Erase);
    f.click(f.canvas.document().atoms[0].pos);
    CHECK(f.canvas.document().bonds.size() == 6);
}

TEST_CASE("erasing an atom drops neighbours left without bonds, as erasing a bond does (#335)") {
    Fixture f;
    f.canvas.setDocumentSilently(*chem::fromSmiles("CCC"));
    f.canvas.setTool(Canvas::Tool::Erase);
    f.click(f.doc().atoms[1].pos);  // the middle carbon
    CHECK(f.doc().atoms.empty());  // not two lone methanes

    f.canvas.setDocumentSilently(*chem::fromSmiles("CCO"));
    f.click(f.doc().atoms[0].pos);  // an end: the rest stays bonded
    CHECK(f.doc().atoms.size() == 2);
    CHECK(f.doc().bonds.size() == 1);
}

TEST_CASE("drag erase sweeps skipped points, previews changes and makes one undo step (#487)") {
    Fixture f;
    Document d;
    d.atoms = {{{0, 0}}, {{30, 0}}, {{200, 0}, 8}, {{0, 40}, 7}};
    d.bonds = {{0, 1}};
    d.arrows = {{{60, -12}, {60, 12}, ArrowKind::Line}, {{40, 20}, {80, 20}}, {{180, -12}, {180, 12}, ArrowKind::Line}};
    d.arrows[1].bend = 20;
    d.texts = {{{90, 4}, "heat"}, {{180, 50}, "keep"}};
    d.page = "A4";
    f.canvas.setDocumentSilently(d);
    f.canvas.setTool(Canvas::Tool::Erase);
    const QString preview = qEnvironmentVariable("PENZENE_ERASE_PREVIEW");
    QApplication::processEvents();
    if (!preview.isEmpty()) REQUIRE(f.canvas.grab().save(preview + "/before.png"));
    QTest::mousePress(f.canvas.viewport(), Qt::LeftButton, {}, f.at({-20, 0}));
    QMouseEvent move(QEvent::MouseMove, f.at({150, 0}), f.canvas.viewport()->mapToGlobal(f.at({150, 0})),
                     Qt::NoButton, Qt::LeftButton, {});
    QApplication::sendEvent(f.canvas.viewport(), &move);
    CHECK(f.doc() == d);
    CHECK(f.undo.count() == 0);
    QApplication::processEvents();
    if (!preview.isEmpty()) REQUIRE(f.canvas.grab().save(preview + "/pending.png"));
    QTest::mouseRelease(f.canvas.viewport(), Qt::LeftButton, {}, f.at({150, 0}));
    CHECK(f.doc().atoms.size() == 2);
    CHECK(f.doc().bonds.empty());
    REQUIRE(f.doc().arrows.size() == 1);
    CHECK(f.doc().arrows[0] == d.arrows[2]);
    REQUIRE(f.doc().texts.size() == 1);
    CHECK(f.doc().texts[0].text == "keep");
    CHECK(f.doc().page == "A4");
    CHECK(f.undo.count() == 1);
    QApplication::processEvents();
    if (!preview.isEmpty()) REQUIRE(f.canvas.grab().save(preview + "/after.png"));
    const Document erased = f.doc();
    f.undo.undo();
    CHECK(f.doc() == d);
    f.undo.redo();
    CHECK(f.doc() == erased);
}

TEST_CASE("drag erase handles bond-only strokes, loopbacks and cancellation (#487)") {
    Fixture f;
    Document d;
    d.atoms = {{{-30, 0}}, {{30, 0}}, {{-30, 30}}, {{30, 30}}};
    d.bonds = {{0, 1}, {0, 2}, {1, 3}};
    f.canvas.setDocumentSilently(d);
    f.canvas.setTool(Canvas::Tool::Erase);
    f.drag({0, -15}, {0, 15});
    CHECK(f.doc().atoms.size() == 4);
    CHECK(f.doc().bonds.size() == 2);
    CHECK(f.undo.count() == 1);
    f.undo.undo();
    CHECK(f.doc() == d);
    QTest::mousePress(f.canvas.viewport(), Qt::LeftButton, {}, f.at({0, -15}));
    for (QPointF p : {QPointF{0, 15}, QPointF{0, -15}}) {
        QMouseEvent move(QEvent::MouseMove, f.at(p), f.canvas.viewport()->mapToGlobal(f.at(p)), Qt::NoButton, Qt::LeftButton, {});
        QApplication::sendEvent(f.canvas.viewport(), &move);
    }
    QTest::mouseRelease(f.canvas.viewport(), Qt::LeftButton, {}, f.at({0, -15}));
    CHECK(f.doc().bonds.size() == 2);  // returning to the press position is still a drag
    f.undo.undo();
    QTest::mousePress(f.canvas.viewport(), Qt::LeftButton, {}, f.at({0, -15}));
    QTest::keyClick(f.canvas.viewport(), Qt::Key_Escape);
    QTest::mouseRelease(f.canvas.viewport(), Qt::LeftButton, {}, f.at({0, 15}));
    CHECK(f.doc() == d);
    QTest::mousePress(f.canvas.viewport(), Qt::LeftButton, {}, f.at({0, -15}));
    Document replacement = d;
    replacement.addAtom({80, 80}, 8);
    f.canvas.setDocumentSilently(replacement);
    QTest::mouseRelease(f.canvas.viewport(), Qt::LeftButton, {}, f.at({0, 15}));
    CHECK(f.doc() == replacement);
    const int steps = f.undo.count();
    f.drag({100, -15}, {100, 15});
    CHECK(f.undo.count() == steps);
    Document box;
    box.arrows = {{{-40, -40}, {40, 40}, ArrowKind::Box}};
    f.canvas.setDocumentSilently(box);
    f.drag({-10, 0}, {10, 0});  // crossing an unfilled box's interior does not cross its outline
    CHECK(f.doc() == box);
}

TEST_CASE("select and delete, charges") {
    Fixture f;
    f.canvas.setTool(Canvas::Tool::Bond);
    f.click({0, 0});
    f.canvas.setTool(Canvas::Tool::ChargePlus);
    f.click(f.canvas.document().atoms[0].pos);
    CHECK(f.canvas.document().atoms[0].charge == 1);
    f.canvas.selectAll();
    f.canvas.deleteSelection();
    CHECK(f.canvas.document().atoms.empty());
    f.undo.undo();
    CHECK(f.canvas.document().atoms.size() == 2);
}

TEST_CASE("Alt-drag on empty space selects what a freehand loop takes in (#499)") {
    Fixture f;
    Document d;  // any rectangle round atoms 0 and 2 takes in atom 1 too
    d.atoms = {{{0, 0}}, {{30, 10}}, {{40, 40}}};
    f.canvas.setDocumentSilently(d);
    f.canvas.setTool(Canvas::Tool::Select);
    const std::vector<QPointF> loop{{-8, 0}, {0, -8}, {48, 40}, {40, 48}};  // a band along the diagonal
    QTest::mousePress(f.canvas.viewport(), Qt::LeftButton, Qt::AltModifier, f.at(loop[0]));
    for (QPointF p : loop) {
        QMouseEvent move(QEvent::MouseMove, f.at(p), f.canvas.viewport()->mapToGlobal(f.at(p)), Qt::NoButton, Qt::LeftButton,
                         Qt::AltModifier);
        QApplication::sendEvent(f.canvas.viewport(), &move);
    }
    QTest::mouseRelease(f.canvas.viewport(), Qt::LeftButton, Qt::AltModifier, f.at(loop.back()));
    CHECK(f.canvas.selection() == QSet<int>{0, 2});
}

TEST_CASE("insert centres the fragment and selects it") {
    Fixture f;
    auto frag = chem::fromSmiles("CCO");
    REQUIRE(frag);
    f.canvas.insert(*frag, "Paste");
    f.canvas.insert(*frag, "Paste");
    CHECK(f.canvas.document().atoms.size() == 6);
    CHECK(f.canvas.document().bonds[2].a == 3);
    CHECK(f.canvas.selection() == QSet<int>{3, 4, 5});
}

TEST_CASE("Find Substructure selects all matches, preserves edits and seeds from selection (#485)") {
    App app;
    MainWindow w;
    w.resize(1000, 800);
    w.show();
    auto* canvas = w.findChild<Canvas*>();
    REQUIRE(canvas);
    const Document doc = *chem::fromSmiles("CC(=O)OC.CC=O");
    canvas->setDocumentSilently(doc);
    canvas->fitToSelection();
    canvas->setSelection({1, 2});
    QAction* find = nullptr;
    for (auto* a : w.findChildren<QAction*>()) if (a->text() == "&Find Substructure…") find = a;
    REQUIRE(find);
    QString initial;
    const QString preview = qEnvironmentVariable("PENZENE_FIND_PREVIEW");
    QApplication::processEvents();
    if (!preview.isEmpty()) REQUIRE(w.grab().save(preview + "/before.png"));
    auto search = [&](const QString& pattern) {
        QTimer::singleShot(0, &w, [&, pattern] {
            auto* dialog = qobject_cast<QInputDialog*>(QApplication::activeModalWidget());
            REQUIRE(dialog);
            initial = dialog->textValue();
            dialog->setTextValue(pattern);
            dialog->accept();
        });
        find->trigger();
    };
    search("C=O");
    CHECK_FALSE(initial.isEmpty());
    CHECK(canvas->selection() == QSet<int>{1, 2, 6, 7});
    CHECK(w.statusBar()->currentMessage() == "2 matches selected");
    CHECK(canvas->document() == doc);
    QApplication::processEvents();
    if (!preview.isEmpty()) REQUIRE(w.grab().save(preview + "/after.png"));
    QTimer::singleShot(0, &w, [&] {
        QTimer::singleShot(0, &w, [] {
            auto* box = qobject_cast<QMessageBox*>(QApplication::activeModalWidget());
            REQUIRE(box);
            box->accept();
        });
    });
    search("C1");
    CHECK(canvas->selection() == QSet<int>{1, 2, 6, 7});
    search("N");
    CHECK(canvas->selection().isEmpty());
    CHECK(w.statusBar()->currentMessage() == "No matches found");
    CHECK(canvas->document() == doc);
}

// The documentation's screenshots: `pixi run screenshots` writes them to docs/_static.
TEST_CASE("main window screenshots") {
    App app;
    const QString dir = qEnvironmentVariable("PENZENE_SCREENSHOTS");
    struct Shot {
        const char *name, *example, *dock;  // a Welcome example, with a panel open
    };
    for (QString theme : {"Light", "Dark"}) {
        QSettings().setValue("theme", theme);
        QSettings().remove("element");  // the periodic-table icon, as on first run
        for (auto s : {Shot{"screenshot", "Aspirin", nullptr}, Shot{"templates", "Aspirin", "templates"},
                       Shot{"properties", "Aspirin", "properties"}, Shot{"scheme", "Reaction scheme", nullptr}}) {
            MainWindow w;
            w.resize(1000, 800);
            w.show();
            if (s.dock) {
                auto* dock = w.findChild<QDockWidget*>(s.dock);
                REQUIRE(dock);
                dock->show();
                if (auto* tree = dock->findChild<QTreeWidget*>()) tree->topLevelItem(0)->setExpanded(true);
            }
            QApplication::processEvents();  // lay out the panel before the example fits the view
            QToolButton* example = nullptr;
            for (auto* b : w.findChildren<QToolButton*>("example"))
                if (b->text() == s.example) example = b;
            REQUIRE(example);
            example->click();
            QApplication::processEvents();
            if (!dir.isEmpty()) w.grab().save(QString("%1/%2-%3.png").arg(dir, s.name, theme.toLower()));
        }
    }
    QSettings().remove("theme");
}

// #42: a ring on a terminal atom must continue straight on, so the substituent
// bond bisects the ring's outside angle.
TEST_CASE("ring on a terminal atom bisects the outside angle") {
    for (int n : {3, 6}) {
        Fixture f;
        f.canvas.setTool(Canvas::Tool::Bond);
        f.click({0, 0});
        const QPointF stem = f.doc().atoms[0].pos;
        f.canvas.setTool(Canvas::Tool::Ring);
        f.canvas.setRing(n, n == 6);
        f.click(f.doc().atoms[1].pos);
        REQUIRE(f.doc().atoms.size() == size_t(n + 1));
        QPointF p = f.doc().atoms[1].pos, sum;
        for (int nb : f.doc().neighbors(1))
            if (f.doc().atoms[nb].pos != stem) {
                QPointF v = f.doc().atoms[nb].pos - p;
                sum += v / std::hypot(v.x(), v.y());
            }
        QPointF s = stem - p;
        double cosang = (sum.x() * s.x() + sum.y() * s.y()) / (std::hypot(sum.x(), sum.y()) * std::hypot(s.x(), s.y()));
        CHECK(cosang < -0.999);
    }
}

TEST_CASE("hotkeys: chain, labels, groups, bonds") {
    Fixture f;
    f.canvas.setTool(Canvas::Tool::Bond);
    f.click({0, 0});
    f.hover(f.doc().atoms[1].pos);
    f.key("111");  // hotspot follows each new atom
    CHECK(f.doc().atoms.size() == 5);
    CHECK(f.doc().bonds.size() == 4);

    f.key("O");  // Shift+o: OMe, as an abbreviation
    REQUIRE(f.doc().atoms.size() == 5);
    CHECK(f.doc().atoms[4].z == 8);
    CHECK(f.doc().atoms[4].label == "OMe");
    f.key("F");  // Shift+f: CF3 replaces it
    CHECK(f.doc().atoms[4].z == 6);
    CHECK(f.doc().atoms[4].label == "CF3");

    // Shift+arrow jumps atom to atom back along the chain.
    f.hover(f.doc().atoms[1].pos);
    QPointF dir = f.doc().atoms[0].pos - f.doc().atoms[1].pos;
    QTest::keyClick(f.canvas.viewport(), dir.x() < 0 ? Qt::Key_Left : Qt::Key_Right, Qt::ShiftModifier);
    f.key("n");
    CHECK(f.doc().atoms[0].z == 7);

    // Bond hotkeys: 2 on a bond makes it double (on an atom it sprouts a carbonyl).
    const auto& d = f.doc();
    f.hover((d.atoms[d.bonds[1].a].pos + d.atoms[d.bonds[1].b].pos) / 2);
    f.key("2");
    CHECK(f.doc().bonds[1].order == 2);
    f.key("w");
    CHECK(f.doc().bonds[1].stereo == BondStereo::Wedge);
    int a = f.doc().bonds[1].a;
    f.key("w");
    CHECK(f.doc().bonds[1].b == a);  // flipped
}

static std::string noStereo(std::string smi) {
    std::erase(smi, '@');
    return chem::toSmiles(*chem::fromSmiles(smi));
}

// The worked example from ChemDraw's cheat sheet: from H2N-CH3, "42n152o" builds Ala-Ala.
TEST_CASE("hotkeys: ChemDraw dipeptide example") {
    Fixture f;
    f.canvas.setTool(Canvas::Tool::Bond);
    f.click({0, 0});
    f.hover(f.doc().atoms[0].pos);
    f.key("n");
    f.hover(f.doc().atoms[1].pos);
    f.key("42n152o");
    CHECK(noStereo(chem::toSmiles(f.doc())) == noStereo("CC(N)C(=O)NC(C)C(=O)O"));
    int wedges = 0, hashes = 0;
    for (auto& b : f.doc().bonds) wedges += b.stereo == BondStereo::Wedge, hashes += b.stereo == BondStereo::Hash;
    CHECK(wedges == 1);
    CHECK(hashes == 1);
}

TEST_CASE("hotkeys: context-dependent sprouts") {
    Fixture f;
    f.canvas.setTool(Canvas::Tool::Chain);
    f.drag({0, 0}, {40, 0});  // propane-ish chain
    REQUIRE(f.doc().atoms.size() >= 3);
    f.hover(f.doc().atoms[1].pos);  // secondary carbon
    f.key("2");
    CHECK(f.doc().atoms.back().z == 8);  // ketone on the hotspot
    CHECK(f.doc().bonds.back().order == 2);

    Fixture g;
    g.canvas.setTool(Canvas::Tool::Ring);
    g.canvas.setRing(6, false);
    g.click({0, 0});
    g.hover(g.doc().atoms[0].pos);  // secondary ring carbon
    g.key("9");                     // gem-dimethyl
    CHECK(g.doc().atoms.size() == 8);
    CHECK(g.doc().neighbors(0).size() == 4);

    Fixture t;  // tertiary ring carbon: "6" adds a C-C bond, then a cyclohexane
    t.canvas.setTool(Canvas::Tool::Ring);
    t.canvas.setRing(6, false);
    t.click({0, 0});
    t.hover(t.doc().atoms[0].pos);
    t.key("1");
    t.hover(t.doc().atoms[0].pos);
    t.key("6");
    CHECK(t.doc().atoms.size() == 6 + 1 + 1 + 5);
    CHECK(t.doc().neighbors(0).size() == 4);

    Fixture h;  // phenyl on an aromatic carbon goes via a C-C bond (biphenyl)
    h.canvas.setTool(Canvas::Tool::Ring);
    h.canvas.setRing(6, true);
    h.click({0, 0});
    h.hover(h.doc().atoms[0].pos);
    h.key("a");
    CHECK(noStereo(chem::toSmiles(h.doc())) == noStereo("c1ccc(-c2ccccc2)cc1"));
}

TEST_CASE("hotspot is sticky and arrows walk atom -> bond -> atom") {
    Fixture f;
    f.canvas.setTool(Canvas::Tool::Bond);
    f.click({0, 0});
    f.hover(f.doc().atoms[1].pos);
    f.hover({300, 300});  // drift off into empty space
    CHECK(f.canvas.hotspotAtom() == 1);
    f.key("1");
    CHECK(f.doc().atoms.size() == 3);
    CHECK(f.canvas.hotspotAtom() == 2);

    QPointF back = f.doc().atoms[1].pos - f.doc().atoms[2].pos;
    auto toward = [&](QPointF v) {
        return std::abs(v.x()) > std::abs(v.y()) ? (v.x() < 0 ? Qt::Key_Left : Qt::Key_Right)
                                                 : (v.y() < 0 ? Qt::Key_Up : Qt::Key_Down);
    };
    QTest::keyClick(f.canvas.viewport(), toward(back));
    CHECK(f.canvas.hotspotBond() == 1);
    QTest::keyClick(f.canvas.viewport(), toward(back));
    CHECK(f.canvas.hotspotAtom() == 1);
    QTest::keyClick(f.canvas.viewport(), Qt::Key_Escape);
    CHECK(f.canvas.hotspotAtom() == -1);
}

TEST_CASE("a drawing beyond ±5000 pt can be brought into view (#326)") {
    Fixture f;
    Document d;
    d.addAtom({8000, 0}, 7);
    f.canvas.setDocumentSilently(d);
    f.canvas.fitToDocument();
    CHECK(f.canvas.viewport()->rect().contains(f.at({8000, 0})));
}

TEST_CASE("hotkeys: fused ring on a bond") {
    Fixture f;
    f.canvas.setTool(Canvas::Tool::Bond);
    f.click({0, 0});
    f.hover((f.doc().atoms[0].pos + f.doc().atoms[1].pos) / 2);
    f.key("a");
    CHECK(f.doc().atoms.size() == 6);  // benzene shares the bond's 2 atoms
    CHECK(f.doubles() == 3);
}

TEST_CASE("applyLabel understands elements, groups and SMILES") {
    Document d;
    d.atoms = {{{0, 0}}, {{kBondLength, 0}}};
    d.bonds = {{0, 1}};
    CHECK(edit::applyLabel(d, 1, "Br"));
    CHECK(d.atoms[1].z == 35);
    CHECK(edit::applyLabel(d, 1, "NO2"));  // abbreviation: one labelled atom
    CHECK(d.atoms[1].label == "NO2");
    CHECK(d.atoms[1].charge == 1);
    CHECK(d.atoms.size() == 2);
    CHECK(chem::toSmiles(d) == "C[N+](=O)[O-]");  // chemistry sees the full group
    CHECK(edit::applyLabel(d, 1, "OH"));
    CHECK(d.atoms[1].z == 8);
    CHECK(d.atoms[1].label.isEmpty());
    CHECK(edit::applyLabel(d, 1, "C(=O)Cl"));  // SMILES: drawn out
    CHECK(d.atoms.size() == 4);
    CHECK_FALSE(edit::applyLabel(d, 0, "notachem!!"));
}

TEST_CASE("abbreviations: valence, clean, expand") {
    Fixture f;
    Document d;
    d.atoms = {{{0, 0}}, {{kBondLength, 0}}, {{2 * kBondLength, 5}}};
    d.bonds = {{0, 1}, {1, 2}};
    REQUIRE(edit::applyLabel(d, 2, "Boc"));
    CHECK_FALSE(chem::atomInfo(d)[2].valenceError);
    CHECK(chem::toSmiles(d) == "CCC(=O)OC(C)(C)C");
    Document clean = chem::clean2D(d);
    CHECK(clean.atoms.size() == 3);  // still abbreviated
    CHECK(clean.bonds.size() == 2);
    CHECK(clean.atoms[2].label == "Boc");
    f.canvas.setDocumentSilently(d);
    f.canvas.expandAbbreviations();
    CHECK(f.doc().atoms.size() == 9);
    CHECK(f.doc().atoms[2].label.isEmpty());
    CHECK(chem::toSmiles(f.doc()) == "CCC(=O)OC(C)(C)C");
    auto back = Document::fromJson(d.toJson());
    REQUIRE(back);
    CHECK(back->atoms[2].label == "Boc");
}

TEST_CASE("expanding a charged group keeps its charge, from the menu and the atom's own action (#589)") {
    Fixture f;
    Document d;
    d.addAtom({0, 0}, 7);
    REQUIRE(edit::applyLabel(d, 0, "N3-", true));
    REQUIRE(chem::properties(d));
    const std::string azide = chem::properties(d)->formula;
    CHECK(azide.find("H") == std::string::npos);
    auto charge = [](const Document& doc) {
        int sum = 0;
        for (const Atom& a : doc.atoms) sum += a.charge;
        return sum;
    };
    f.canvas.setDocumentSilently(d);
    f.canvas.expandAbbreviations();
    CHECK(f.doc().atoms.size() == 3);
    CHECK(charge(f.doc()) == -1);
    CHECK(chem::properties(f.doc())->formula == azide);
    Document one = d;
    REQUIRE(chem::expandLabel(one, 0));  // what the atom's Expand Abbreviation action calls
    CHECK(charge(one) == -1);
}

static double angleAt(const Document& d, int centre, int x, int y) {
    QPointF u = d.atoms[x].pos - d.atoms[centre].pos, v = d.atoms[y].pos - d.atoms[centre].pos;
    return std::acos((u.x() * v.x() + u.y() * v.y()) / (std::hypot(u.x(), u.y()) * std::hypot(v.x(), v.y()))) * 180 / M_PI;
}

// #43: sp centres (allenes, alkynes) are linear.
TEST_CASE("allene and alkyne centres are linear") {
    Fixture f;  // drawing double bonds onto a double bond
    f.canvas.setTool(Canvas::Tool::Bond);
    f.canvas.setBondOrder(2);
    f.click({0, 0});
    f.click(f.doc().atoms[1].pos);
    REQUIRE(f.doc().atoms.size() == 3);
    CHECK(angleAt(f.doc(), 1, 0, 2) > 179);

    Fixture g;  // zig-zag chain, then make both bonds double: terminal atom swings into line
    g.canvas.setTool(Canvas::Tool::Chain);
    g.drag({0, 0}, {40, 0});
    REQUIRE(g.doc().atoms.size() == 4);
    for (int bi : {1, 2}) {
        const auto& d = g.doc();
        g.hover((d.atoms[d.bonds[bi].a].pos + d.atoms[d.bonds[bi].b].pos) / 2);
        g.key("2");
    }
    CHECK(angleAt(g.doc(), 2, 1, 3) > 179);

    Fixture h;  // growing from an alkyne carbon continues straight
    h.canvas.setTool(Canvas::Tool::Bond);
    h.canvas.setBondOrder(3);
    h.click({0, 0});
    h.hover(h.doc().atoms[1].pos);
    h.key("1");
    CHECK(angleAt(h.doc(), 1, 0, 2) > 179);
}

TEST_CASE("a curved arrow stays on the atoms and bonds it's drawn between (#498)") {
    Fixture f;
    Document d;  // H3C–C=O, a little off the grid so Clean moves it
    const QPointF c(kBondLength, 2), o = c + QPointF(0.5, -0.866) * kBondLength;
    d.atoms = {{{0, 0}}, {c}, {o, 8}};
    d.bonds = {{0, 1}, {1, 2, 2}};
    f.canvas.setDocumentSilently(d);
    f.canvas.setTool(Canvas::Tool::Arrow);
    f.canvas.setArrow(ArrowKind::Reaction, true);
    f.drag(c + (o - c) * 0.4, o + QPointF(0.6, -0.2) * kBondLength);  // from along the C=O bond to beside the O
    f.canvas.setArrow(ArrowKind::Reaction, false);
    f.drag({0, 60}, {40, 60});
    REQUIRE(f.doc().arrows.size() == 2);
    CHECK(f.doc().arrows[0].fromAt == std::array{1, 2});
    CHECK(f.doc().arrows[0].toAt == std::array{2, -1});
    CHECK(len(f.doc().arrows[0].from - (c + o) / 2) < 1e-6);  // snapped to the bond's middle
    CHECK(std::abs(len(f.doc().arrows[0].to - o) - 0.55 * kBondLength) < 1e-6);  // and to just off the O's label
    CHECK(f.doc().arrows[1].fromAt == std::array{-1, -1});  // a reaction arrow rests on nothing
    f.canvas.setArrow(ArrowKind::Reaction, true);
    f.click(arrowPath(f.doc().arrows[0]).pointAtPercent(0.5));  // flipped, it still rests on them
    CHECK(f.doc().arrows[0].toAt == std::array{2, -1});
    CHECK(f.doc().arrows[0].fromAt == std::array{1, 2});

    // Drag the O: the end on it moves as far, the one on the bond half as far.
    f.canvas.setTool(Canvas::Tool::Select);
    Document was = f.doc();
    f.drag(d.atoms[2].pos, d.atoms[2].pos + QPointF(10, 0));
    QPointF moved = f.doc().atoms[2].pos - was.atoms[2].pos;
    REQUIRE(len(moved) > 5);
    CHECK(len(f.doc().arrows[0].to - (was.arrows[0].to + moved)) < 1e-6);
    CHECK(len(f.doc().arrows[0].from - (was.arrows[0].from + moved / 2)) < 1e-6);
    CHECK(f.doc().arrows[1] == was.arrows[1]);
    was = f.doc();
    f.drag(c, c + QPointF(0, 8));  // now the C: the O's end stays
    moved = f.doc().atoms[1].pos - was.atoms[1].pos;
    REQUIRE(len(moved) > 5);
    CHECK(len(f.doc().arrows[0].from - (was.arrows[0].from + moved / 2)) < 1e-6);
    CHECK(f.doc().arrows[0].to == was.arrows[0].to);

    // Clean moves the atoms too.
    was = f.doc();
    f.canvas.commit(chem::clean2D(f.doc()), "Clean");
    moved = f.doc().atoms[2].pos - was.atoms[2].pos;
    REQUIRE(len(moved) > 0.5);
    CHECK(len(f.doc().arrows[0].to - (was.arrows[0].to + moved)) < 1e-6);

    // Saved and opened again, it still knows.
    CHECK(Document::fromJson(f.doc().toJson())->arrows == f.doc().arrows);

    // Delete the O: the arrow stays, with both ends free (the C=O bond went too).
    f.canvas.setSelection({2});
    f.canvas.deleteSelection();
    REQUIRE(f.doc().arrows.size() == 2);
    CHECK(f.doc().arrows[0].fromAt == std::array{-1, -1});
    CHECK(f.doc().arrows[0].toAt == std::array{-1, -1});
}

TEST_CASE("a lone selected curved arrow reshapes by its handles; heads come in sizes (#534)") {
    Fixture f;
    Document d;
    d.atoms = {{{0, 0}}, {{60, 0}, 8}, {{0, 60}, 7}};
    Arrow a{{10, 0}, {50, 0}};
    a.bend = 10, a.fromAt = {0, -1}, a.toAt = {1, -1};
    d.arrows = {a};
    f.canvas.setDocumentSilently(d);
    f.canvas.setTool(Canvas::Tool::Select);
    f.canvas.setSelection({}, {0});
    auto top = [&] { return arrowPath(f.doc().arrows[0]).boundingRect().top(); };
    REQUIRE(std::abs(top() + 10) < 0.5);

    f.drag({30, -10}, {32, -30});  // the curve's handle: its top goes where it's dropped
    CHECK(std::abs(top() + 30) < 0.5);
    f.drag({30, -30}, {30, 12});  // and through the chord, to the other side
    CHECK(std::abs(arrowPath(f.doc().arrows[0]).boundingRect().bottom() - 12) < 0.5);

    f.drag({50, 0}, {4, 64});  // an end dropped by the N rests on it, just off it
    CHECK(f.doc().arrows[0].toAt == std::array{2, -1});
    CHECK(std::abs(len(f.doc().arrows[0].to - QPointF(0, 60)) - 0.55 * kBondLength) < 1e-6);
    f.drag(f.doc().arrows[0].to, {-5, 62});  // moved round the N, it stays on it
    CHECK(f.doc().arrows[0].toAt == std::array{2, -1});
    const QPointF end = f.doc().arrows[0].to;
    f.drag({0, 60}, {20, 60});  // and follows it
    CHECK(len(f.doc().arrows[0].to - (end + QPointF(20, 0))) < 1e-6);

    f.canvas.setSelection({}, {0});
    f.canvas.setArrowHead(1.5);
    CHECK(f.doc().arrows[0].head == 1.5);
    CHECK(Document::fromJson(f.doc().toJson())->arrows == f.doc().arrows);
    QJsonObject file = QJsonDocument::fromJson(f.doc().toJson()).object();  // a file asking for a huge head
    QJsonArray arrows = file["arrows"].toArray();
    QJsonObject huge = arrows[0].toObject();
    huge["head"] = 1e9;
    arrows[0] = huge, file["arrows"] = arrows;
    CHECK(Document::fromJson(QJsonDocument(file).toJson())->arrows[0].head == 4);
}

TEST_CASE("arrows: draw, restyle, select, move, delete; text subscripts") {
    Fixture f;
    f.canvas.setTool(Canvas::Tool::Arrow);
    f.canvas.setArrow(ArrowKind::Reaction, false);
    f.drag({0, 0}, {50, 3});  // snaps to horizontal
    REQUIRE(f.doc().arrows.size() == 1);
    CHECK(std::abs(f.doc().arrows[0].to.y()) < 0.5);

    f.canvas.setArrow(ArrowKind::Equilibrium, false);
    f.click({25, 0});  // restyles instead of adding
    REQUIRE(f.doc().arrows.size() == 1);
    CHECK(f.doc().arrows[0].kind == ArrowKind::Equilibrium);

    f.canvas.setArrow(ArrowKind::Reaction, true);
    f.drag({0, 40}, {40, 40});
    REQUIRE(f.doc().arrows.size() == 2);
    const double bend = f.doc().arrows[1].bend;
    CHECK(bend != 0);
    QPointF mid = arrowPath(f.doc().arrows[1]).pointAtPercent(0.5);
    f.click(mid);  // same tool again: flips the curve
    CHECK(f.doc().arrows[1].bend == -bend);

    Document withText = f.doc();
    withText.texts.push_back({{0, -30}, "CH2Cl2"});
    f.canvas.setDocumentSilently(withText);
    f.canvas.setTool(Canvas::Tool::Select);
    f.drag({-10, -50}, {60, 10});  // rubber band: straight arrow and text, not the curve
    CHECK(f.canvas.selectedArrows() == QSet<int>{0});
    CHECK(f.canvas.selectedTexts() == QSet<int>{0});
    f.drag({25, 0}, {25, 20});  // drag the arrow: text moves with it
    CHECK(std::abs(f.doc().arrows[0].from.y() - 20) < 1);
    CHECK(std::abs(f.doc().texts[0].pos.y() + 10) < 1);
    f.canvas.deleteSelection();
    CHECK(f.doc().arrows.size() == 1);
    CHECK(f.doc().texts.empty());

    // Formula subscripts sit below the baseline; digits after a space don't.
    auto bottom = [](const QString& s) { return textPath({{0, 0}, s}).boundingRect().bottom(); };
    CHECK(bottom("H2") > bottom("H") + 1);
    CHECK(bottom("80 C") <= bottom("H") + 0.5);
    CHECK(bottom("(2 equiv)") <= bottom("(") + 0.5);

    // Tabs jump to stops, so columns line up; leading spaces indent.
    auto left = [](const QString& s) { return textPath({{0, 0}, s}).boundingRect().left(); };
    auto right = [](const QString& s) { return textPath({{0, 0}, s}).boundingRect().right(); };
    CHECK(std::abs(right("\tA") - right("ab\tA")) < 0.01);
    CHECK(left("\tA") > left("A") + 5);
    CHECK(left("  A") > left("A") + 2);
    CHECK(bottom("CH2Cl2") > bottom("CHCl") + 1);  // run layout keeps subscripts
}

TEST_CASE("hotkeys: bond styles, positions, chair; duplicate across an arrow; tool keys") {
    Fixture f;
    f.canvas.setTool(Canvas::Tool::Chain);
    f.drag({0, 0}, {40, 0});
    REQUIRE(f.doc().bonds.size() >= 2);
    auto hoverBond = [&](int i) {
        const auto& d = f.doc();
        f.hover((d.atoms[d.bonds[i].a].pos + d.atoms[d.bonds[i].b].pos) / 2);
    };
    hoverBond(0);
    f.key("y");
    CHECK(f.doc().bonds[0].stereo == BondStereo::Wavy);
    f.key("e");
    CHECK(f.doc().bonds[0].stereo == BondStereo::DoubleDashed);
    CHECK(f.doc().bonds[0].order == 1);
    f.key("2");
    CHECK(f.doc().bonds[0].order == 2);
    CHECK(f.doc().bonds[0].stereo == BondStereo::None);
    f.key("e");
    CHECK(f.doc().bonds[0].stereo == BondStereo::DoubleDashed);
    CHECK(f.doc().bonds[0].order == 2);  // keeps a double bond double
    f.key("Q");
    CHECK(f.doc().bonds[0].order == 4);
    CHECK(f.doc().bonds[0].stereo == BondStereo::None);
    f.key("B");
    CHECK(f.doc().bonds[0].order == 2);
    CHECK(f.doc().bonds[0].stereo == BondStereo::Bold);
    f.key("r");
    CHECK(f.doc().bonds[0].position == BondPosition::Right);
    auto back = Document::fromJson(f.doc().toJson());
    REQUIRE(back);
    CHECK(*back == f.doc());

    const size_t before = f.doc().atoms.size();
    hoverBond(1);
    f.key("9");  // chair: four new atoms, bonds all about one bond long
    CHECK(f.doc().atoms.size() == before + 4);
    for (const auto& b : f.doc().bonds) {
        QPointF v = f.doc().atoms[b.a].pos - f.doc().atoms[b.b].pos;
        CHECK(std::abs(std::hypot(v.x(), v.y()) - kBondLength) < 0.15 * kBondLength);
    }

    // Duplicate across an arrow on the right.
    Document d;
    d.atoms = {{{0, 0}}, {{kBondLength, 0}}};
    d.bonds = {{0, 1}};
    d.arrows = {{{40, 0}, {80, 0}}};
    f.canvas.setDocumentSilently(d);
    f.canvas.setSelection({0, 1});
    f.canvas.duplicateSelection({1, 0});
    REQUIRE(f.doc().atoms.size() == 4);
    CHECK(f.doc().atoms[2].pos.x() > 80);  // beyond the arrow head
    CHECK(f.canvas.selection() == QSet<int>{2, 3});

    QString picked;
    QObject::connect(&f.canvas, &Canvas::toolKey, [&](const QString& k) { picked = k; });
    QTest::keyClick(f.canvas.viewport(), Qt::Key_Escape);
    f.key("e");
    CHECK(picked == "e");
}

TEST_CASE("drawing style presets: JDP from its ChemDraw stationery") {
    App app;
    const auto& jdp = drawingStyle("JDP");
    CHECK(jdp.name == "JDP");
    CHECK(std::abs(jdp.lineWidth * jdp.bondLength / kBondLength - 0.879) < 1e-9);  // native pt in exports
    CHECK(drawingStyle("").name == "ACS 1996");
    CHECK(drawingStyle("no such style").name == "ACS 1996");

    // The style is part of the document: it changes the rendering.
    Document d = *chem::fromSmiles("CC(=O)O");
    Document j = d;
    j.style = "JDP";
    CHECK(renderSvg(d) != renderSvg(j));
}

TEST_CASE("themes: Catppuccin palettes; exports stay black") {
    CHECK(theme("Catppuccin Mocha").paper == QColor("#1e1e2e"));
    CHECK(theme("Catppuccin Latte").ink == QColor("#4c4f69"));
    CHECK(theme("no such theme").name == "System");
    App app;
    Fixture f;
    f.canvas.setTheme(theme("Catppuccin Mocha"));
    auto doc = chem::fromSmiles("CO");
    REQUIRE(doc);
    QImage img = renderImage(*doc, {72});
    bool dark = false;
    for (int y = 0; y < img.height(); ++y)
        for (int x = 0; x < img.width(); ++x)
            if (QColor c = img.pixelColor(x, y); c.alpha() > 200 && c.lightness() < 60) dark = true;
    CHECK(dark);  // black ink, not the theme's pale text
}

TEST_CASE("imported themes persist, appear in the menu and keep exports black (#475)") {
    App app;
    Theme custom = theme("Light");
    custom.name = "Shared & ../theme";
    custom.paper = QColor("#eff1f5");
    const QString dir = QStandardPaths::writableLocation(QStandardPaths::AppDataLocation) + "/themes";
    const QStringList before = QDir(dir).entryList({"*.json"}, QDir::Files);
    QTemporaryDir files;
    REQUIRE(files.isValid());
    const QString source = files.filePath("shared.json"), exported = files.filePath("exported");
    REQUIRE(writeWhole(source, custom.toJson()));
    MainWindow w;
    w.resize(1000, 800);
    w.show();
    auto* canvas = w.findChild<Canvas*>();
    REQUIRE(canvas);
    canvas->setDocumentSilently(*chem::fromSmiles("CCO"));
    auto action = [&](const QString& name) {
        for (auto* a : w.findChildren<QAction*>()) if (a->text() == name || a->data().toString() == name) return a;
        return static_cast<QAction*>(nullptr);
    };
    const QString preview = qEnvironmentVariable("PENZENE_THEME_PREVIEW");
    if (!preview.isEmpty()) {
        action("Light")->trigger();
        QApplication::processEvents();
        REQUIRE(w.grab().save(QString(preview).replace(".png", "-before.png")));
    }
    QApplication::setAttribute(Qt::AA_DontUseNativeDialogs);
    auto chooseFile = [&](const QString& path, bool cancelReplacement = false) {
        QTimer::singleShot(0, &w, [&, path, cancelReplacement] {
            auto* dialog = qobject_cast<QFileDialog*>(QApplication::activeModalWidget());
            REQUIRE(dialog);
            if (cancelReplacement) QTimer::singleShot(0, &w, [] {
                auto* box = qobject_cast<QMessageBox*>(QApplication::activeModalWidget());
                REQUIRE(box);
                box->button(QMessageBox::No)->click();
            });
            dialog->selectFile(path);
            static_cast<QDialog*>(dialog)->accept();
        });
    };
    REQUIRE(action("Import Theme…"));
    chooseFile(source);
    action("Import Theme…")->trigger();
    QApplication::processEvents();
    CHECK(theme(custom.name).toJson() == custom.toJson());
    CHECK(QSettings().value("theme").toString() == custom.name);
    REQUIRE(action(custom.name));
    CHECK(action(custom.name)->isChecked());
    CHECK(action(custom.name)->text() == "Shared && ../theme");
    CHECK(QApplication::palette().color(QPalette::Base) == custom.paper);
    CHECK(renderSvg(canvas->document()).contains("#000000"));
    if (!preview.isEmpty()) REQUIRE(w.grab().save(preview));
    chooseFile(exported);
    action("Export Theme…")->trigger();
    QFile output(exported + ".json");
    REQUIRE(output.open(QIODevice::ReadOnly));
    CHECK(output.readAll() == custom.toJson());
    Theme replacement = custom;
    replacement.paper = QColor("#fafafa");
    REQUIRE(writeWhole(source, replacement.toJson()));
    chooseFile(source, true);
    action("Import Theme…")->trigger();
    CHECK(theme(custom.name).toJson() == custom.toJson());
    CHECK(QApplication::palette().color(QPalette::Base) == custom.paper);
    QApplication::setAttribute(Qt::AA_DontUseNativeDialogs, false);
    MainWindow reopened;
    bool found = false;
    for (auto* a : reopened.findChildren<QAction*>()) if (a->data().toString() == custom.name) found = a->isChecked();
    CHECK(found);
    Theme preset = theme("Light");
    preset.ink = Qt::red;
    CHECK_FALSE(saveUserTheme(preset));
    CHECK(theme("Light").ink != QColor(Qt::red));
    for (const QString& file : QDir(dir).entryList({"*.json"}, QDir::Files))
        if (!before.contains(file)) REQUIRE(QFile::remove(dir + "/" + file));
    QSettings().remove("theme");
}

TEST_CASE("message box icons: a black glyph takes the theme's text colour (#418)") {
    App app;
    const QColor text = theme("Catppuccin Mocha").text;
    QPixmap glyph(32, 32);  // like macOS's question mark: black on clear
    glyph.fill(Qt::transparent);
    QPainter(&glyph).fillRect(12, 4, 8, 24, Qt::black);
    const QImage inked = inkGlyph(glyph, text).toImage();
    CHECK(inked.pixelColor(16, 16) == text);
    CHECK(inked.pixelColor(2, 2).alpha() == 0);
    QPixmap coloured = glyph;  // Fusion's blue circle keeps its own colours
    QPainter(&coloured).fillRect(0, 0, 8, 8, QColor(40, 120, 255));
    CHECK(inkGlyph(coloured, text).toImage() == coloured.toImage());
    QPixmap light = glyph;  // macOS's grey-and-white information bubble
    light.fill(Qt::lightGray);
    CHECK(inkGlyph(light, text).toImage() == light.toImage());
    std::unique_ptr<QStyle> style(themedStyle());
    QMessageBox box;
    QPalette pal = box.palette();
    pal.setColor(QPalette::WindowText, text);
    box.setPalette(pal);
    box.setIconPixmap(glyph);
    box.setStyle(style.get());
    box.ensurePolished();  // as showing it does
    CHECK(box.iconPixmap().toImage().pixelColor(16, 16) == text);
}

TEST_CASE("ring fill: click inside toggles; survives delete, copy and save") {
    Fixture f;
    auto nap = chem::fromSmiles("c1ccc2ccccc2c1");
    REQUIRE(nap);
    f.canvas.setDocumentSilently(*nap);
    // Centre of one ring: mean of the atoms in its fill.
    auto rings = chem::rings(f.doc());
    REQUIRE(rings.size() == 2);
    QPointF c;
    for (int i : rings[0]) c += f.doc().atoms[i].pos / 6;
    f.canvas.setTool(Canvas::Tool::Fill);
    f.click(c);
    REQUIRE(f.doc().fills.size() == 1);
    CHECK(f.doc().fills[0].atoms.size() == 6);
    auto back = Document::fromJson(f.doc().toJson());
    REQUIRE(back);
    CHECK(*back == f.doc());
    f.click(c);  // same colour again: cleared
    CHECK(f.doc().fills.empty());
    f.click(c);
    REQUIRE(f.doc().fills.size() == 1);
    Document d = f.doc();
    d.removeAtoms({d.fills[0].atoms[0]});
    CHECK(d.fills.empty());  // a ring missing an atom loses its fill
}

TEST_CASE("2 on a double bond swaps the side of its second line") {
    Fixture f;
    auto d = chem::fromSmiles("CC=CC");  // trans-2-butene: an offset (not centred) double bond
    REQUIRE(d);
    f.canvas.setDocumentSilently(*d);
    int db = -1;
    for (int i = 0; i < int(f.doc().bonds.size()); ++i)
        if (f.doc().bonds[i].order == 2) db = i;
    REQUIRE(db >= 0);
    const auto& doc = f.doc();
    f.hover((doc.atoms[doc.bonds[db].a].pos + doc.atoms[doc.bonds[db].b].pos) / 2);
    f.key("2");
    CHECK(f.doc().bonds[db].order == 2);
    const BondPosition first = f.doc().bonds[db].position;
    CHECK(first != BondPosition::Auto);
    f.key("2");
    CHECK(f.doc().bonds[db].position != first);  // and back again
    CHECK(f.doc().bonds[db].position != BondPosition::Auto);
}

TEST_CASE("picking from the periodic table switches to the atom tool") {
    App app;
    QSettings().remove("element");  // first run: nothing picked yet
    MainWindow w;
    w.show();
    auto* canvas = w.findChild<Canvas*>();
    REQUIRE(canvas);
    QToolButton* nitrogen = nullptr;  // in the periodic table popup
    for (auto* wa : w.findChildren<QWidgetAction*>())
        if (wa->defaultWidget())
            for (auto* b : wa->defaultWidget()->findChildren<QToolButton*>())
                if (b->text() == "N") nitrogen = b;
    REQUIRE(nitrogen);
    nitrogen->click();  // no need to pick the atom tool first
    QTest::mouseClick(canvas->viewport(), Qt::LeftButton, {}, canvas->viewport()->rect().center());
    REQUIRE(canvas->document().atoms.size() == 1);
    CHECK(canvas->document().atoms[0].z == 7);
    CHECK(QSettings().value("element").toString() == "N");  // remembered for next time (#128)
    MainWindow again;  // a new window starts with the picked element
    auto* c2 = again.findChild<Canvas*>();
    c2->setTool(Canvas::Tool::Atom);
    again.show();
    QTest::mouseClick(c2->viewport(), Qt::LeftButton, {}, c2->viewport()->rect().center());
    REQUIRE(c2->document().atoms.size() == 1);
    CHECK(c2->document().atoms[0].z == 7);
}

TEST_CASE("Ac, Pr and Ts are groups, not actinium, praseodymium and tennessine (#125)") {
    for (auto [label, formula] : {std::pair{"Ac", "C3H6O"}, {"Pr", "C4H10"}, {"Ts", "C8H10O2S"}}) {
        Document d;
        d.atoms = {{{0, 0}}, {{kBondLength, 0}}};
        d.bonds = {{0, 1}};
        REQUIRE(edit::applyLabel(d, 1, label));
        CHECK(d.atoms[1].label == label);
        CHECK(chem::properties(d)->formula == formula);
    }
    Fixture f;  // and the Shift+A hotkey
    Document d;
    d.atoms = {{{0, 0}}, {{kBondLength, 0}}};
    d.bonds = {{0, 1}};
    f.canvas.setDocumentSilently(d);
    f.hover(f.doc().atoms[1].pos);
    f.key("A");
    CHECK(f.doc().atoms[1].label == "Ac");
}

TEST_CASE("curved arrows are circular arcs, exact past 180 degrees (#86)") {
    // Chord 40, bend 30: more than a semicircle, radius (20² + 30²) / 60 = 21.7.
    const QRectF r = arrowPath({{0, 0}, {40, 0}, ArrowKind::Reaction, 30}).boundingRect();
    CHECK(std::abs(r.top() + 30) < 0.2);          // the arc's midpoint sits 30 above the chord
    CHECK(r.left() < -1.4);                          // and it bulges past both ends
    CHECK(r.right() > 41.4);
    // A small bend still matches the old midpoint.
    QPointF mid = arrowPath({{0, 0}, {40, 0}, ArrowKind::Reaction, 8}).pointAtPercent(0.5);
    CHECK(std::abs(mid.y() + 8) < 0.2);
}

TEST_CASE("a no-reaction arrow is crossed, draws its cross, and keeps it through .penz and CDXML") {
    Fixture f;
    f.canvas.setDocumentSilently({});
    f.canvas.setTool(Canvas::Tool::Arrow);
    f.canvas.setArrow(ArrowKind::Reaction, false, false, OrbitalLook::Outline, true);
    f.drag({0, 0}, {60, 0});
    REQUIRE(f.doc().arrows.size() == 1);
    CHECK(f.doc().arrows[0].crossed);

    Document d = *chem::fromSmiles("CCO");
    d.arrows.push_back({{40, 0}, {100, 0}});
    Document crossed = d;
    crossed.arrows[0].crossed = true;
    CHECK(renderImage(crossed, {300}) != renderImage(d, {300}));  // the ✕ is drawn
    CHECK(Document::fromJson(crossed.toJson())->arrows[0].crossed);
    CHECK_FALSE(Document::fromJson(d.toJson())->arrows[0].crossed);
    const auto back = chem::fromChemDraw(chem::toCdxml(crossed));
    REQUIRE(back);
    REQUIRE(back->arrows.size() == 1);
    CHECK(back->arrows[0].crossed);
    CHECK(chem::toCdxml(crossed).contains("NoGo=\"Cross\""));
}

TEST_CASE("arrowhead size goes out to ChemDraw and comes back (#538)") {
    Document d = *chem::fromSmiles("CCO");
    d.arrows.push_back({{40, 0}, {100, 0}});
    d.arrows[0].head = 1.5;
    const QByteArray xml = chem::toCdxml(d);
    CHECK(xml.contains("LineWidth=\"0.6\""));
    CHECK(xml.contains("HeadSize=\"1500\""));  // 9 pt at the style's 0.6 pt line: 1.5 × our 6
    CHECK(xml.contains("ArrowheadCenterSize=\"1313\""));
    CHECK(xml.contains("ArrowheadWidth=\"375\""));
    CHECK(std::abs(chem::fromChemDraw(xml)->arrows[0].head - 1.5) < 1e-9);
    CHECK(std::abs(chem::fromChemDraw(chem::toCdx(d))->arrows[0].head - 1.5) < 1e-9);

    // ChemDraw draws a head HeadSize % of the line width long: 1000 on a 0.6 pt line is our 6
    auto head = [](const char* root, const char* arrow) {
        const QString xml = QString(R"(<CDXML BondLength="14.4" %1><page><arrow Tail3D="0 0 0" Head3D="100 0 0" )"
                                    R"(ArrowheadHead="Full" ArrowheadType="Solid" %2/></page></CDXML>)").arg(root, arrow);
        return chem::fromChemDraw(xml.toUtf8())->arrows.at(0).head;
    };
    CHECK(std::abs(head(R"(LineWidth="0.6")", R"(HeadSize="1000")") - 1) < 1e-9);
    CHECK(std::abs(head(R"(LineWidth="0.6")", "") - 1) < 1e-9);
    CHECK(std::abs(head(R"(LineWidth="0.6")", R"(HeadSize="1000" LineWidth="1.2")") - 2) < 1e-9);
    CHECK(std::abs(head(R"(LineWidth="0.6")", R"(HeadSize="2250")") - 2.25) < 1e-9);
}

TEST_CASE("a huge arrow bend is drawn with a bounded number of points (#314)") {
    Arrow a{{0, 0}, {10, 0}};
    a.bend = 1e7;  // from a file: this took 2.9 GB to draw
    CHECK(arrowPath(a).elementCount() <= 1000);
}

TEST_CASE("flip mirrors (enantiomer with wedges kept), align and distribute (#87)") {
    Fixture f;
    auto ala = chem::fromSmiles("C[C@H](N)C(=O)O");
    REQUIRE(ala);
    Document d = *ala;
    d.arrows.push_back({{60, 0}, {100, 0}, ArrowKind::Reaction, 10});
    int dbl = -1;
    for (int i = 0; i < int(d.bonds.size()); ++i)
        if (d.bonds[i].order == 2) dbl = i;
    d.bonds[dbl].position = BondPosition::Left;
    f.canvas.setDocumentSilently(d);
    const std::string before = chem::toSmiles(f.doc());
    f.canvas.flipSelection(true);
    CHECK(chem::toSmiles(f.doc()) == chem::toSmiles(*chem::fromSmiles("C[C@@H](N)C(=O)O")));  // mirror image
    CHECK(f.doc().bonds[dbl].position == BondPosition::Right);
    CHECK(f.doc().arrows[0].bend == -10);
    f.canvas.flipSelection(true);  // flipping back restores everything
    CHECK(chem::toSmiles(f.doc()) == before);
    for (size_t i = 0; i < d.atoms.size(); ++i)
        CHECK(std::hypot(f.doc().atoms[i].pos.x() - d.atoms[i].pos.x(), f.doc().atoms[i].pos.y() - d.atoms[i].pos.y()) < 1e-6);

    // Three methanols at uneven spacing and heights.
    Document three;
    for (double x : {0.0, 30.0, 100.0}) {
        int c = three.addAtom({x, x / 5});
        int o = three.addAtom({x + kBondLength, x / 5}, 8);
        three.bonds.push_back({c, o});
    }
    f.canvas.setDocumentSilently(three);
    f.canvas.alignSelection(Canvas::Align::Top);
    for (int i = 0; i < 6; ++i) CHECK(std::abs(f.doc().atoms[i].pos.y() - f.doc().atoms[0].pos.y()) < 1e-6);
    f.canvas.distributeSelection(true);
    const double gap1 = f.doc().atoms[2].pos.x() - f.doc().atoms[1].pos.x();
    const double gap2 = f.doc().atoms[4].pos.x() - f.doc().atoms[3].pos.x();
    CHECK(std::abs(gap1 - gap2) < 1e-6);
    CHECK(std::abs(f.doc().atoms[0].pos.x()) < 1e-9);  // outermost objects stay put
    CHECK(std::abs(f.doc().atoms[4].pos.x() - 100) < 1e-9);
}

TEST_CASE("invert stereochemistry swaps wedges and hashes in place (#374)") {
    Fixture f;
    auto d = chem::fromSmiles("C[C@@H](N)C(=O)O.C/C=C/C");  // (R)-alanine and (E)-but-2-ene
    REQUIRE(d);
    f.canvas.setDocumentSilently(*d);
    auto labels = [&] {
        QStringList out;
        for (const auto& l : chem::stereoLabels(f.doc())) out << l.text;
        out.sort();
        return out;
    };
    const QStringList before = labels();
    REQUIRE(before.join(' ').contains('R'));
    REQUIRE(before.join(' ').contains('E'));
    f.canvas.invertStereo();
    QStringList after = labels();
    CHECK(after.join(' ').contains('S'));
    CHECK_FALSE(after.join(' ').contains('R'));
    CHECK(after.join(' ').contains('E'));  // the double bond is untouched
    for (size_t i = 0; i < d->atoms.size(); ++i) CHECK(f.doc().atoms[i].pos == d->atoms[i].pos);
    CHECK(f.undo.count() == 1);
    f.undo.undo();
    CHECK(labels() == before);
}

TEST_CASE("exports come out at the drawing style's own bond length (#92)") {
    App app;
    auto d = chem::fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
    REQUIRE(d);
    const QImage acs = renderImage(*d, {72});
    d->style = "RSC";
    const QImage rsc = renderImage(*d, {72});
    CHECK(std::abs(exportScale(*d) - 12.2 / 14.4) < 1e-12);
    // Same drawing, same model bounds (RSC's labels are smaller, so compare against its own bounds).
    const QRectF model = outputBounds(*d);
    CHECK(std::abs(rsc.width() - model.width() * 12.2 / 14.4) <= 1);
    CHECK(rsc.width() < acs.width());
}

TEST_CASE("canvas drawing style uses physical bond length without changing model coordinates (#624)") {
    Fixture f;
    const Document acs = *chem::fromSmiles("CCO");
    f.canvas.setDocumentSilently(acs);
    f.canvas.centerOn(20, 30);
    const QPointF centre = f.canvas.mapToScene(f.canvas.viewport()->rect().center());
    const double zoom = f.canvas.transform().m11();
    Document rsc = acs;
    rsc.style = "RSC";
    f.canvas.commit(rsc, "Drawing style");
    CHECK(std::abs(f.canvas.transform().m11() / zoom - 12.2 / 14.4) < 1e-10);
    CHECK(len(f.canvas.mapToScene(f.canvas.viewport()->rect().center()) - centre) < 1);
    for (size_t i = 0; i < acs.atoms.size(); ++i) CHECK(f.doc().atoms[i].pos == acs.atoms[i].pos);
    f.undo.undo();
    CHECK(std::abs(f.canvas.transform().m11() - zoom) < 1e-10);
    f.undo.redo();
    CHECK(std::abs(f.canvas.transform().m11() / exportScale(rsc) - zoom) < 1e-10);
    f.canvas.setTool(Canvas::Tool::Select);
    f.click(rsc.atoms[2].pos);
    CHECK(f.canvas.selection().contains(2));
    f.canvas.zoomBy(0.21 / zoom);
    CHECK(std::abs(f.canvas.transform().m11() / exportScale(rsc) - 0.21) < 1e-10);
    const double low = f.canvas.transform().m11();
    f.canvas.zoomBy(0.5);
    CHECK(f.canvas.transform().m11() == low);
    f.canvas.zoomBy(39 / 0.21);
    const double high = f.canvas.transform().m11();
    f.canvas.zoomBy(2);
    CHECK(f.canvas.transform().m11() == high);
    f.canvas.fitToSelection();
    CHECK(f.canvas.transform().m11() / exportScale(rsc) <= 10.001);
}

TEST_CASE("style switching, pages and opening update the canvas scale (#624)") {
    App app;
    const QVariant oldTheme = QSettings().value("theme");
    QSettings().setValue("theme", "Light");
    MainWindow w;
    w.resize(1000, 700);
    w.show();
    auto* canvas = w.findChild<Canvas*>();
    const Document acs = *chem::fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
    canvas->setDocumentSilently(acs);
    canvas->setGuides(false, false);
    canvas->zoomBy(2);
    canvas->centerOn(documentBounds(acs).center());
    QApplication::processEvents();
    const double zoom = canvas->transform().m11();
    const QString preview = qEnvironmentVariable("PENZENE_STYLE_PREVIEW_DIR");
    if (!preview.isEmpty()) REQUIRE(w.grab().save(preview + "/acs.png"));
    QAction* rsc = nullptr;
    QAction* undo = nullptr;
    QAction* redo = nullptr;
    for (auto* a : w.findChildren<QAction*>()) {
        if (a->text() == "RSC") rsc = a;
        if (a->shortcut() == QKeySequence::Undo) undo = a;
        if (a->shortcut() == QKeySequence::Redo) redo = a;
    }
    REQUIRE(rsc);
    REQUIRE(undo);
    REQUIRE(redo);
    rsc->trigger();
    QApplication::processEvents();
    CHECK(std::abs(canvas->transform().m11() / zoom - 12.2 / 14.4) < 1e-10);
    if (!preview.isEmpty()) {
        REQUIRE(w.grab().save(preview + "/rsc.png"));
        REQUIRE(writeWhole(preview + "/aspirin-rsc.penz", canvas->document().toJson()));
    }
    undo->trigger();
    CHECK(std::abs(canvas->transform().m11() - zoom) < 1e-10);
    redo->trigger();
    CHECK(std::abs(canvas->transform().m11() / zoom - 12.2 / 14.4) < 1e-10);
    auto* tabs = w.findChild<QTabBar*>("pageTabs");
    w.findChild<QToolButton*>("addPage")->click();
    canvas->setDocumentSilently(acs);
    const double second = canvas->transform().m11();
    tabs->setCurrentIndex(0);
    CHECK(std::abs(canvas->transform().m11() / second - 12.2 / 14.4) < 1e-10);
    tabs->setCurrentIndex(1);
    CHECK(std::abs(canvas->transform().m11() - second) < 1e-10);
    QTemporaryDir dir;
    Document saved = acs;
    saved.style = "RSC";
    REQUIRE(writeWhole(dir.filePath("rsc.penz"), saved.toJson()));
    REQUIRE(w.openFile(dir.filePath("rsc.penz")));
    CHECK(canvas->document().style == "RSC");
    const double opened = canvas->transform().m11() / exportScale(saved);
    canvas->setDocumentSilently(acs);
    CHECK(std::abs(canvas->transform().m11() - opened) < 1e-10);
    if (oldTheme.isValid()) QSettings().setValue("theme", oldTheme);
    else QSettings().remove("theme");
}

// Finds a menu entry by its text, looking inside submenus.
static QAction* findAction(QMenu* menu, const QString& text) {
    for (QAction* a : menu->actions()) {
        if (a->text() == text) return a;
        if (a->menu())
            if (QAction* sub = findAction(a->menu(), text)) return sub;
    }
    return nullptr;
}

TEST_CASE("terminal bond deletion drops the end atom through each UI path") {
    for (int path = 0; path < 3; ++path) {
        Fixture f;
        Document d;
        d.atoms = {{{0, 0}}, {{kBondLength, 0}}, {{2 * kBondLength, 0}}};
        d.bonds = {{0, 1}, {1, 2}};
        f.canvas.setDocumentSilently(d);
        const QPointF mid(1.5 * kBondLength, 0);
        if (path == 0) {
            f.hover(mid);
            QTest::keyClick(f.canvas.viewport(), Qt::Key_Delete);
        } else if (path == 1) {
            f.canvas.setTool(Canvas::Tool::Erase);
            f.click(mid);
        } else {
            QMenu* menu = f.canvas.contextMenuAt(mid);
            QAction* action = findAction(menu, "Delete Bond");
            REQUIRE(action);
            action->trigger();
        }
        REQUIRE(f.doc().atoms.size() == 2);
        CHECK(f.doc().bonds.size() == 1);
        CHECK(f.doc().atoms[1].pos == QPointF(kBondLength, 0));
        if (path == 0)
            if (auto out = qgetenv("PENZENE_DELETE_SHOT"); !out.isEmpty()) f.canvas.grab().save(out);
        f.undo.undo();
        CHECK(f.doc() == d);
    }
}

TEST_CASE("Delete returns the cursor to the previous point (#675)") {
    Document d;
    d.atoms = {{{0, 0}}, {{kBondLength, 0}}, {{2 * kBondLength, 0}}};
    d.bonds = {{0, 1}, {1, 2}};
    SECTION("after deleting the atom under the cursor") {
        Fixture f;
        f.canvas.setDocumentSilently(d);
        f.hover({2 * kBondLength, 0});
        REQUIRE(f.canvas.hotspotAtom() == 2);
        QTest::keyClick(f.canvas.viewport(), Qt::Key_Delete);
        CHECK(f.doc().atoms.size() == 2);
        CHECK(f.canvas.hotspotAtom() == 1);
        f.key("1");  // drawing carries on from there
        CHECK(f.doc().atoms.size() == 3);
        CHECK(f.doc().bonds.size() == 2);
    }
    SECTION("after deleting the terminal bond under the cursor") {
        Fixture f;
        f.canvas.setDocumentSilently(d);
        f.hover({1.5 * kBondLength, 0});
        QTest::keyClick(f.canvas.viewport(), Qt::Key_Delete);
        CHECK(f.doc().atoms.size() == 2);
        CHECK(f.canvas.hotspotAtom() == 1);
        f.key("1");
        CHECK(f.doc().atoms.size() == 3);
        CHECK(f.doc().bonds.size() == 2);
    }
}

TEST_CASE("no two items in a menu share an access key") {
    App app;
    MainWindow w;
    for (auto* menu : w.findChildren<QMenu*>()) {
        QMap<QChar, QString> taken;
        for (auto* a : menu->actions()) {
            const QString t = a->text();
            const int i = t.indexOf('&');
            if (i < 0 || i + 1 >= t.size() || t[i + 1] == '&') continue;
            INFO(menu->title().toStdString() << ": " << t.toStdString());
            CHECK_FALSE(taken.contains(t[i + 1].toLower()));
            taken[t[i + 1].toLower()] = t;
        }
    }
}

TEST_CASE("right-click menus for atoms, bonds, selection and canvas (#89)") {
    Fixture f;
    Document d;
    d.atoms = {{{0, 0}}, {{kBondLength, 0}}};
    d.bonds = {{0, 1}};
    f.canvas.setDocumentSilently(d);

    QMenu* atomMenu = f.canvas.contextMenuAt(f.doc().atoms[1].pos);
    REQUIRE(findAction(atomMenu, "N"));
    findAction(atomMenu, "N")->trigger();
    CHECK(f.doc().atoms[1].z == 7);
    findAction(f.canvas.contextMenuAt(f.doc().atoms[1].pos), "Boc")->trigger();
    CHECK(f.doc().atoms[1].label == "Boc");

    QMenu* bondMenu = f.canvas.contextMenuAt({kBondLength / 2, 0});
    REQUIRE(findAction(bondMenu, "Double"));
    findAction(bondMenu, "Double")->trigger();
    CHECK(f.doc().bonds[0].order == 2);
    CHECK(findAction(f.canvas.contextMenuAt({kBondLength / 2, 0}), "Right"));  // double: position submenu

    f.canvas.selectAll();
    CHECK(findAction(f.canvas.contextMenuAt(f.doc().atoms[0].pos), "Flip Horizontal"));
    f.canvas.setSelection({});
    CHECK(findAction(f.canvas.contextMenuAt({200, 200}), "Select All"));
}

TEST_CASE("a stereocentre's menu tags it And 1, then offers And 2 (#389)") {
    Fixture f;
    auto d = chem::fromSmiles("C[C@H](N)C(=O)O");  // alanine
    REQUIRE(d);
    f.canvas.setDocumentSilently(*d);
    CHECK_FALSE(findAction(f.canvas.contextMenuAt(f.doc().atoms[0].pos), "And 1"));  // not a stereocentre
    findAction(f.canvas.contextMenuAt(f.doc().atoms[1].pos), "And 1")->trigger();
    CHECK(stereoGroupTag(f.doc().atoms[1]) == "&1");
    CHECK(findAction(f.canvas.contextMenuAt(f.doc().atoms[1].pos), "And 2"));
    findAction(f.canvas.contextMenuAt(f.doc().atoms[1].pos), "None")->trigger();
    CHECK(f.doc().atoms[1].stereoGroup == StereoGroup::None);
}

TEST_CASE("colour atoms, bonds, arrows and text; exports keep the colour (#82)") {
    Fixture f;
    Document d;
    d.atoms = {{{0, 0}}, {{kBondLength, 0}, 8}};
    d.bonds = {{0, 1}};
    d.arrows = {{{40, 0}, {80, 0}}};
    d.texts = {{{40, -10}, "heat"}};
    f.canvas.setDocumentSilently(d);
    const QColor red(214, 39, 40);
    f.canvas.setColour(red);
    f.canvas.setTool(Canvas::Tool::Colour);
    f.click(f.doc().atoms[1].pos);
    CHECK(f.doc().atoms[1].color == red);
    f.click(f.doc().atoms[1].pos);  // same colour again: cleared
    CHECK_FALSE(f.doc().atoms[1].color.isValid());

    f.canvas.selectAll();
    f.canvas.colourSelection();
    CHECK(f.doc().bonds[0].color == red);
    CHECK(f.doc().arrows[0].color == red);
    CHECK(f.doc().texts[0].color == red);
    auto back = Document::fromJson(f.doc().toJson());
    REQUIRE(back);
    CHECK(*back == f.doc());

    Document bond;  // a red bond exports red (colours are the user's, unlike theme ink)
    bond.atoms = {{{0, 0}}, {{kBondLength * 3, 0}}};
    bond.bonds = {{0, 1}};
    bond.bonds[0].color = red;
    const QImage img = renderImage(bond, {150});
    bool sawRed = false;
    for (int y = 0; y < img.height(); ++y)
        for (int x = 0; x < img.width(); ++x)
            if (QColor c = img.pixelColor(x, y); c.alpha() > 200 && c.red() > 150 && c.green() < 90) sawRed = true;
    CHECK(sawRed);
}

TEST_CASE("copy and paste within Penzene keeps what ChemDraw files don't, such as brackets") {
    App app;
    MainWindow w;
    auto* canvas = w.findChild<Canvas*>();
    Document d = *chem::fromSmiles("CCO");
    d.brackets.push_back({{0, 1}, true, "n"});
    canvas->setDocumentSilently(d);
    canvas->selectAll();
    for (const QString text : {"&Copy", "&Paste"})
        for (auto* a : w.findChildren<QAction*>())
            if (a->text() == text) a->trigger();
    REQUIRE(canvas->document().atoms.size() == 6);
    CHECK(canvas->document().brackets.size() == 2);
}

TEST_CASE("choosing a tool explains it in the status bar (#93)") {
    App app;
    MainWindow w;
    w.show();
    QAction* chain = nullptr;
    for (auto* a : w.findChildren<QAction*>())
        if (a->toolTip().startsWith("Chain")) chain = a;
    REQUIRE(chain);
    chain->trigger();
    CHECK(w.statusBar()->currentMessage().startsWith("Chain"));
    CHECK(w.statusBar()->currentMessage().contains("drag"));
}

TEST_CASE("Look Up on PubChem needs a structure (#500)") {
    App app;
    MainWindow w;
    QAction* lookUp = nullptr;
    for (auto* a : w.findChildren<QAction*>())
        if (a->text() == "Look Up on &PubChem") lookUp = a;
    REQUIRE(lookUp);
    QMenu* structure = nullptr;
    for (auto* m : w.findChildren<QMenu*>())
        if (m->title() == "&Structure") structure = m;
    REQUIRE(structure);
    emit structure->aboutToShow();
    CHECK_FALSE(lookUp->isEnabled());  // an empty drawing
    w.findChild<Canvas*>()->setDocumentSilently(*chem::fromSmiles("c1ccccc1"));
    emit structure->aboutToShow();
    CHECK(lookUp->isEnabled());
}

TEST_CASE("drop an atom on another to merge; Shift for free angles and straight moves (#88)") {
    Fixture f;
    auto drag = [&](QPointF from, QPointF to, Qt::KeyboardModifiers mods) {
        QTest::mousePress(f.canvas.viewport(), Qt::LeftButton, {}, f.at(from));
        QMouseEvent move(QEvent::MouseMove, f.at(to), f.canvas.viewport()->mapToGlobal(f.at(to)), Qt::NoButton,
                         Qt::LeftButton, mods);
        QApplication::sendEvent(f.canvas.viewport(), &move);
        QTest::mouseRelease(f.canvas.viewport(), Qt::LeftButton, mods, f.at(to));
    };
    // Two separate bonds; drag the second's end atom onto the first's.
    Document d;
    d.atoms = {{{0, 0}}, {{kBondLength, 0}}, {{60, 30}}, {{60 + kBondLength, 30}}};
    d.bonds = {{0, 1}, {2, 3}};
    f.canvas.setDocumentSilently(d);
    f.canvas.setTool(Canvas::Tool::Select);
    drag({60, 30}, {kBondLength, 0}, {});
    CHECK(f.doc().atoms.size() == 3);
    CHECK(f.doc().bonds.size() == 2);
    CHECK(f.doc().neighbors(1).size() == 2);  // one connected chain now
    f.undo.undo();
    CHECK(f.doc().atoms.size() == 4);

    // Shift-drag: moves along one axis only.
    f.canvas.setSelection({2, 3});
    drag({60, 30}, {90, 36}, Qt::ShiftModifier);
    CHECK(std::abs(f.doc().atoms[2].pos.y() - 30) < 1e-6);
    CHECK(std::abs(f.doc().atoms[2].pos.x() - 90) < 0.5);

    // Shift while drawing a bond: any angle, not snapped to 30 degrees.
    f.canvas.setDocumentSilently({});
    f.canvas.setTool(Canvas::Tool::Bond);
    drag({0, 0}, {20, 20}, Qt::ShiftModifier);
    REQUIRE(f.doc().bonds.size() == 1);
    QPointF v = f.doc().atoms[1].pos - f.doc().atoms[0].pos;
    CHECK(std::abs(std::atan2(v.y(), v.x()) * 180 / M_PI - 45) < 2);
}

TEST_CASE("drag the ring tool to size a ring (#127)") {
    Fixture f;
    f.canvas.setTool(Canvas::Tool::Ring);
    f.canvas.setRing(6, false);
    // One atom per half bond length of drag, from 3.
    Document none;
    f.canvas.setDocumentSilently(none);
    f.drag({0, 0}, {5.2 * kBondLength / 2, 0});  // 5 steps: an 8-membered ring
    CHECK(f.doc().atoms.size() == 8);
    CHECK(f.doc().bonds.size() == 8);
    f.canvas.setDocumentSilently(none);
    f.drag({0, 0}, {0.5, 0});  // barely moved: a click, so the chosen 6-ring
    CHECK(f.doc().atoms.size() == 6);
    // Dragging from a bond fuses the sized ring onto it.
    Document bond;
    bond.atoms = {{{0, 0}}, {{kBondLength, 0}}};
    bond.bonds = {{0, 1}};
    f.canvas.setDocumentSilently(bond);
    f.drag({kBondLength / 2, 0}, {kBondLength / 2, 2.2 * kBondLength / 2});  // 2 steps: a 5-ring
    CHECK(f.doc().atoms.size() == 5);  // 2 shared + 3 new
    CHECK(f.doc().bonds.size() == 5);
}

TEST_CASE("ring fill colour is picked from the fill button (#126)") {
    App app;
    MainWindow w;
    w.show();
    auto* canvas = w.findChild<Canvas*>();
    QToolButton* fillButton = nullptr;
    for (auto* b : w.findChildren<QToolButton*>())
        if (b->defaultAction() && b->defaultAction()->toolTip().startsWith("Ring fill")) fillButton = b;
    REQUIRE(fillButton);
    REQUIRE(fillButton->menu());
    QToolButton* pink = nullptr;  // a swatch in the drop-down
    for (auto* wa : fillButton->menu()->findChildren<QWidgetAction*>())
        for (auto* b : wa->defaultWidget()->findChildren<QToolButton*>())
            if (b->toolTip() == "Rose") pink = b;  // QColor(255, 214, 214)
    REQUIRE(pink);
    pink->click();
    CHECK(canvas->fillColor() == QColor(255, 214, 214));
    CHECK(fillButton->defaultAction()->isChecked());  // and the fill tool is chosen
}

TEST_CASE("recent files, autosave and crash recovery (#91)") {
    App app;
    QSettings().remove("recentFiles");
    QFile::remove(MainWindow::autosavePath());
    MainWindow w;
    auto* canvas = w.findChild<Canvas*>();
    REQUIRE(w.openFile(QString(PENZENE_TEST_DATA) + "/aspirin.mol"));
    CHECK(w.recentFiles().value(0).endsWith("aspirin.mol"));

    canvas->commit(*chem::fromSmiles("CCO"), "edit");  // unsaved changes
    w.autosave();
    REQUIRE(QFile::exists(MainWindow::autosavePath()));
    // "Crash": the autosave is now another, no longer running, Penzene's.
    const QString crashed = QFileInfo(MainWindow::autosavePath()).dir().filePath("autosave-crashed.penz");
    QFile::remove(crashed);
    REQUIRE(QFile::rename(MainWindow::autosavePath(), crashed));

    // A fresh window after the crash, opened on a file, offers the autosave back; answer Yes.
    MainWindow after;
    REQUIRE(after.openFile(QString(PENZENE_TEST_DATA) + "/aspirin.mol"));
    QTimer::singleShot(0, &after, [] {
        if (auto* box = qobject_cast<QMessageBox*>(QApplication::activeModalWidget()))
            box->button(QMessageBox::Yes)->click();
    });
    after.offerRecovery();
    auto* c2 = after.findChild<Canvas*>();
    CHECK(chem::toSmiles(c2->document()) == "CCO");
    CHECK(after.windowTitle().startsWith("Untitled"));  // so Save doesn't overwrite aspirin.mol (#491)
    CHECK_FALSE(QFile::exists(crashed));  // offered once only

    after.autosave();  // once changes are saved (the stack is clean), autosave removes its copy
    REQUIRE(QFile::exists(MainWindow::autosavePath()));
    after.findChild<QUndoStack*>()->setClean();
    after.autosave();
    CHECK_FALSE(QFile::exists(MainWindow::autosavePath()));
}

TEST_CASE("a file handed over by the system (Finder's Open) opens in the window (#332)") {
    App app;
    MainWindow w;
    QFileOpenEvent open(QString(PENZENE_TEST_DATA) + "/aspirin.mol");
    QApplication::sendEvent(qApp, &open);
    CHECK(chem::toSmiles(w.findChild<Canvas*>()->document()) == chem::toSmiles(*chem::fromSmiles("CC(=O)Oc1ccccc1C(=O)O")));
}

TEST_CASE("preferences: default style for new documents, export resolution and background (#90)") {
    App app;
    QSettings().remove("defaultStyle");
    QSettings().remove("exportBackground");
    MainWindow w;
    auto* canvas = w.findChild<Canvas*>();
    QTimer::singleShot(0, [] {
        auto* dialog = qobject_cast<QDialog*>(QApplication::activeModalWidget());
        REQUIRE(dialog);
        auto* style = dialog->findChild<QComboBox*>("defaultStyle");
        auto* background = dialog->findChild<QComboBox*>("exportBackground");
        REQUIRE(style);
        REQUIRE(background);
        style->setCurrentText("RSC");
        background->setCurrentIndex(1);  // white
        dialog->findChild<QSpinBox*>()->setValue(150);
        dialog->accept();
    });
    w.showPreferences();
    CHECK(QSettings().value("defaultStyle").toString() == "RSC");
    CHECK(QSettings().value("exportDpi").toInt() == 150);
    CHECK(QSettings().value("exportBackground").toString() == "white");

    for (auto* a : w.findChildren<QAction*>())
        if (a->text() == "&New") a->trigger();
    CHECK(canvas->document().style == "RSC");

    // White background: the corner pixel of an export is opaque white, not clear.
    QTemporaryDir dir;
    auto doc = chem::fromSmiles("CCO");
    REQUIRE(exportDocument(*doc, dir.filePath("x.png"), {150, Qt::white}));
    QImage img(dir.filePath("x.png"));
    CHECK(img.pixelColor(0, 0) == QColor(Qt::white));
    REQUIRE(exportDocument(*doc, dir.filePath("y.png")));
    CHECK(QImage(dir.filePath("y.png")).pixelColor(0, 0).alpha() == 0);
    QSettings().remove("defaultStyle");
    QSettings().remove("exportBackground");
    QSettings().remove("exportDpi");
}

TEST_CASE("valence error colouring is a persisted view preference (#620)") {
    App app;
    QSettings().remove("markValenceErrors");
    const QVariant oldTheme = QSettings().value("theme");
    QSettings().setValue("theme", "Light");
    auto hasColour = [](const QImage& img, bool red) {
        for (int y = 0; y < img.height(); ++y)
            for (int x = 0; x < img.width(); ++x) {
                const QColor c = img.pixelColor(x, y);
                if (c.alpha() > 200 && (red ? c.red() > 150 && c.green() < 90 && c.blue() < 100
                                          : c.blue() > 150 && c.red() < 100 && c.green() < 150)) return true;
            }
        return false;
    };
    Document d;
    d.addAtom({0, 0});
    d.atoms[0].charge = 1;
    d.atoms[0].color = QColor(32, 90, 214);
    for (QPointF v : {QPointF(1, 0), QPointF(0, 1), QPointF(-1, 0), QPointF(0, -1)}) {
        d.addAtom(v * kBondLength);
        d.bonds.push_back({0, int(d.atoms.size()) - 1});
    }
    REQUIRE(chem::atomInfo(d)[0].valenceError);
    const auto problems = chem::checkStructure(d);
    REQUIRE(!problems.empty());
    ExportOptions off;
    off.markValenceErrors = false;
    for (auto labels : {Document::CarbonLabels::None, Document::CarbonLabels::All}) {
        d.carbonLabels = labels;
        CHECK(hasColour(renderImage(d), true));
        CHECK_FALSE(hasColour(renderImage(d, off), true));
        CHECK(hasColour(renderImage(d, off), false));
        CHECK(renderSvg(d).contains("#dc2828"));
        CHECK_FALSE(renderSvg(d, off).contains("#dc2828"));
        CHECK(hasColour(QImage::fromData(renderPng(d, off)), false));
        CHECK_FALSE(hasColour(QImage::fromData(renderPng(d, off)), true));
    }
    MainWindow w;
    w.resize(1000, 700);
    w.show();
    auto* canvas = w.findChild<Canvas*>();
    REQUIRE(canvas);
    canvas->setDocumentSilently(d);
    canvas->setGuides(false, false);
    canvas->zoomBy(2);
    canvas->centerOn(0, 0);
    QAction* toggle = nullptr;
    QAction* copy = nullptr;
    for (auto* a : w.findChildren<QAction*>()) {
        if (a->text() == "Mark &Valence Errors") toggle = a;
        if (a->shortcut() == QKeySequence::Copy) copy = a;
    }
    REQUIRE(toggle);
    REQUIRE(copy);
    CHECK(toggle->isChecked());
    const auto* undo = w.findChild<QUndoGroup*>()->activeStack();
    const int history = undo->count();
    const QByteArray original = canvas->document().toJson();
    QApplication::processEvents();
    CHECK(hasColour(canvas->viewport()->grab().toImage(), true));
    const QString preview = qEnvironmentVariable("PENZENE_VALENCE_PREVIEW_DIR");
    if (!preview.isEmpty()) {
        REQUIRE(writeWhole(preview + "/charged-carbon.penz", d.toJson()));
        REQUIRE(w.grab().save(preview + "/on.png"));
    }
    toggle->trigger();
    QApplication::processEvents();
    CHECK_FALSE(toggle->isChecked());
    CHECK_FALSE(QSettings().value("markValenceErrors").toBool());
    CHECK(canvas->document().toJson() == original);
    CHECK(undo->count() == history);
    CHECK(undo->isClean());
    CHECK_FALSE(hasColour(canvas->viewport()->grab().toImage(), true));
    CHECK(hasColour(canvas->viewport()->grab().toImage(), false));
    CHECK(hasColour(renderImage(d), true));  // standalone rendering retains its default
    QWidget* dialog = w.checkStructure();
    auto* list = dialog->findChild<QListWidget*>();
    REQUIRE(list);
    CHECK(list->count() == int(problems.size()));
    dialog->close();
    if (!preview.isEmpty()) REQUIRE(w.grab().save(preview + "/off.png"));
    canvas->selectAll();
    copy->trigger();
    CHECK_FALSE(hasColour(QImage::fromData(QApplication::clipboard()->mimeData()->data("image/png")), true));
    CHECK_FALSE(QApplication::clipboard()->mimeData()->data("image/svg+xml").contains("#dc2828"));
    MainWindow reopened;
    for (auto* a : reopened.findChildren<QAction*>())
        if (a->text() == "Mark &Valence Errors") CHECK_FALSE(a->isChecked());
    QSettings().remove("markValenceErrors");
    if (oldTheme.isValid()) QSettings().setValue("theme", oldTheme);
    else QSettings().remove("theme");
}

TEST_CASE("Check Structure dialog selects the problem's atoms (#95)") {
    App app;
    MainWindow w;
    auto* canvas = w.findChild<Canvas*>();
    canvas->setDocumentSilently(*chem::fromSmiles("CCC(C)O"));
    QWidget* dialog = w.checkStructure();
    auto* list = dialog->findChild<QListWidget*>();
    REQUIRE(list);
    REQUIRE(list->count() == 1);
    list->setCurrentRow(0);
    CHECK(canvas->selection() == QSet<int>{2});  // the unassigned stereocentre
    if (auto out = qgetenv("PENZENE_CHECK_SHOT"); !out.isEmpty()) dialog->grab().save(out);
    dialog->close();
}

TEST_CASE("export checks the structure first, when turned on in Preferences (#390)") {
    App app;
    MainWindow w;
    auto* canvas = w.findChild<Canvas*>();
    QAction* exportAction = nullptr;
    for (auto* a : w.findChildren<QAction*>())
        if (a->text() == "&Export…") exportAction = a;
    REQUIRE(exportAction);
    QApplication::setAttribute(Qt::AA_DontUseNativeDialogs);  // a native dialog can't be answered from here
    // Exports and says which modal came up first: the warning (cancelled) or the file dialog (dismissed).
    auto firstDialog = [&] {
        QString seen;
        QTimer::singleShot(0, &w, [&] {
            QWidget* modal = QApplication::activeModalWidget();
            if (auto* box = qobject_cast<QMessageBox*>(modal)) {
                seen = box->informativeText();
                if (auto out = qgetenv("PENZENE_EXPORT_CHECK_SHOT"); !out.isEmpty()) box->grab().save(out);
                box->button(QMessageBox::Cancel)->click();
            } else if (auto* files = qobject_cast<QFileDialog*>(modal)) {
                seen = "file dialog";
                files->reject();
            }
        });
        exportAction->trigger();
        return seen;
    };

    // The preference, off by default, is a checkbox in Preferences.
    QSettings().remove("checkBeforeExport");
    QTimer::singleShot(0, &w, [] {
        auto* dialog = qobject_cast<QDialog*>(QApplication::activeModalWidget());
        REQUIRE(dialog);
        auto* check = dialog->findChild<QCheckBox*>("checkBeforeExport");
        REQUIRE(check);
        CHECK_FALSE(check->isChecked());
        check->setChecked(true);
        if (auto out = qgetenv("PENZENE_PREFS_SHOT"); !out.isEmpty()) dialog->grab().save(out);
        dialog->accept();
    });
    w.showPreferences();
    CHECK(QSettings().value("checkBeforeExport").toBool());

    canvas->setDocumentSilently(*chem::fromSmiles("CC(O)CC"));  // an unassigned stereocentre
    CHECK(firstDialog().contains("Stereocentre"));
    canvas->setDocumentSilently(*chem::fromSmiles("CC(=O)Oc1ccccc1C(=O)O"));  // aspirin: nothing to report
    CHECK(firstDialog() == "file dialog");
    // Cut is a copy first: cancelling the warning leaves the drawing where it was (#428).
    QAction* cut = nullptr;
    for (auto* a : w.findChildren<QAction*>())
        if (a->text() == "Cu&t") cut = a;
    REQUIRE(cut);
    const Document butanol = *chem::fromSmiles("CC(O)CC");
    canvas->setDocumentSilently(butanol);
    canvas->selectAll();
    QTimer::singleShot(0, &w, [] {
        if (auto* box = qobject_cast<QMessageBox*>(QApplication::activeModalWidget())) box->button(QMessageBox::Cancel)->click();
    });
    cut->trigger();
    CHECK(canvas->document() == butanol);
    QSettings().setValue("checkBeforeExport", false);
    canvas->setDocumentSilently(*chem::fromSmiles("CC(O)CC"));
    CHECK(firstDialog() == "file dialog");
    QApplication::setAttribute(Qt::AA_DontUseNativeDialogs, false);
    QSettings().remove("checkBeforeExport");
}

TEST_CASE("selected aromatic rings can use circles independently (#98)") {
    App app;
    MainWindow w;
    auto* canvas = w.findChild<Canvas*>();
    auto doc = chem::fromSmiles("c1ccc2ccccc2c1");
    REQUIRE(doc);
    const auto before = chem::toSmiles(*doc);
    auto rings = chem::aromaticRings(*doc);
    REQUIRE(rings.size() == 2);
    canvas->setDocumentSilently(*doc);
    QSet<int> selected(rings[0].begin(), rings[0].end());
    canvas->setSelection(selected);
    QAction* oneRing = nullptr;
    QAction* allRings = nullptr;
    for (auto* action : w.findChildren<QAction*>()) {
        if (action->text() == "Circles in Selected &Rings") oneRing = action;
        if (action->text() == "&Aromatic Circles") allRings = action;
    }
    REQUIRE(oneRing);
    REQUIRE(allRings);
    oneRing->trigger();
    CHECK_FALSE(canvas->document().aromaticCircles);
    CHECK(canvas->document().aromaticCircleOverrides.size() == 1);
    CHECK(chem::toSmiles(canvas->document()) == before);
    auto back = Document::fromJson(canvas->document().toJson());
    REQUIRE(back);
    CHECK(back->aromaticCircleOverrides == canvas->document().aromaticCircleOverrides);
    if (auto out = qgetenv("PENZENE_ONE_CIRCLE_SHOT"); !out.isEmpty())
        CHECK(exportDocument(canvas->document(), QString::fromUtf8(out), {150, Qt::white}));
    canvas->setSelection(selected);
    oneRing->trigger();
    CHECK(canvas->document().aromaticCircleOverrides.empty());
    allRings->trigger();
    CHECK(canvas->document().aromaticCircles);
    CHECK(canvas->document().aromaticCircleOverrides.empty());
}

TEST_CASE("properties panel shows descriptors for the selection (#96)") {
    App app;
    MainWindow w;
    w.resize(1100, 700);
    w.show();
    auto* canvas = w.findChild<Canvas*>();
    canvas->setDocumentSilently(*chem::fromSmiles("CC(=O)Oc1ccccc1C(=O)O"));
    QAction* toggle = nullptr;
    for (auto* a : w.findChildren<QAction*>())
        if (a->text() == "&Properties Panel") toggle = a;
    REQUIRE(toggle);
    toggle->trigger();
    QApplication::processEvents();
    QLabel* panel = nullptr;
    for (auto* l : w.findChildren<QLabel*>())
        if (l->text().contains("cLogP")) panel = l;
    REQUIRE(panel);
    CHECK(panel->text().contains("C<sub>9</sub>H<sub>8</sub>O<sub>4</sub>"));
    CHECK(panel->text().contains("63.6"));
    if (auto out = qgetenv("PENZENE_PANEL_SHOT"); !out.isEmpty()) w.grab().save(out);
}

TEST_CASE("mass spec panel shows the isotope pattern of the selection (#397)") {
    App app;
    for (QString t : {"Light", "Dark"}) {
        QSettings().setValue("theme", t);
        MainWindow w;
        w.resize(1100, 700);
        w.show();
        auto* canvas = w.findChild<Canvas*>();
        canvas->setDocumentSilently(*chem::fromSmiles("CC(=O)Oc1ccccc1C(=O)O"));
        QAction* toggle = nullptr;
        for (auto* a : w.findChildren<QAction*>())
            if (a->text() == "&Mass Spec Panel") toggle = a;
        REQUIRE(toggle);
        toggle->trigger();
        QApplication::processEvents();
        auto* dock = w.findChild<QDockWidget*>("massSpec");
        REQUIRE(dock);
        QWidget* spectrum = nullptr;
        for (auto* c : dock->findChildren<QWidget*>())
            if (c->accessibleName() == "Isotope pattern") spectrum = c;
        REQUIRE(spectrum);
        CHECK(spectrum->accessibleDescription().startsWith("181.0495 (100.0%), 182.05"));
        auto* ion = dock->findChild<QComboBox*>();
        ion->setCurrentIndex(0);  // [M]: EI, with the ions to look for
        QLabel* ei = nullptr;
        for (auto* l : dock->findChildren<QLabel*>())
            if (l->text().contains("EI ions")) ei = l;
        REQUIRE(ei);
        CHECK(ei->isVisible());
        CHECK(ei->text().contains("43.0178"));  // aspirin's acetyl
        ion->setCurrentIndex(3);  // [M−H]⁻
        CHECK_FALSE(ei->isVisible());
        CHECK(spectrum->accessibleDescription().startsWith("179.0350 (100.0%)"));
        canvas->setDocumentSilently(*chem::fromSmiles("Clc1ccccc1"));
        emit canvas->documentChanged();
        CHECK(spectrum->accessibleDescription().contains("(32.0%)"));
        ion->setCurrentIndex(1);
        canvas->setDocumentSilently(*chem::fromSmiles("CC(=O)Oc1ccccc1C(=O)O"));
        emit canvas->documentChanged();
        QApplication::processEvents();
        for (auto* b : dock->findChildren<QPushButton*>())
            if (b->text() == "Copy HRMS Line") b->click();
        CHECK(QGuiApplication::clipboard()->text() == "HRMS (ESI) m/z: [M+H]+ calcd for C9H9O4 181.0495");
        canvas->setDocumentSilently(*chem::fromSmiles("*CC"));  // an R group has no mass: the clipboard keeps what it had
        for (auto* b : dock->findChildren<QPushButton*>())
            if (b->text() == "Copy HRMS Line") b->click();
        CHECK(QGuiApplication::clipboard()->text() == "HRMS (ESI) m/z: [M+H]+ calcd for C9H9O4 181.0495");
        canvas->setDocumentSilently(*chem::fromSmiles("CC(=O)Oc1ccccc1C(=O)O"));
        if (auto out = qEnvironmentVariable("PENZENE_MASS_SHOT"); !out.isEmpty())
            w.grab().save(QString(out).replace(".png", "-" + t.toLower() + ".png"));
    }
    QSettings().remove("theme");
}

TEST_CASE("Arrange → Number Compounds numbers each molecule once, in scheme order as they move (#504)") {
    App app;
    MainWindow w;
    auto* canvas = w.findChild<Canvas*>();
    canvas->setDocumentSilently(*chem::fromSmiles("CCO.c1ccccc1"));
    QAction* number = nullptr;
    for (auto* a : w.findChildren<QAction*>())
        if (a->text() == "&Number Compounds") number = a;
    REQUIRE(number);
    number->trigger();
    const auto& d = canvas->document();
    REQUIRE(d.texts.size() == 2);
    CHECK(d.texts[0].text == "1");
    CHECK(d.texts[1].text == "2");
    CHECK(d.texts[0].compound);
    CHECK(d.atoms[d.texts[1].anchor].z == 6);
    double lo = 1e9, hi = -1e9;
    for (int i = 3; i < 9; ++i) lo = std::min(lo, d.atoms[i].pos.x()), hi = std::max(hi, d.atoms[i].pos.x());
    CHECK(std::abs(textPath(d.texts[1], documentStyle(d)).boundingRect().center().x() - (lo + hi) / 2) < 0.1 * kBondLength);  // centred under it
    number->trigger();
    CHECK(d.texts.size() == 2);  // numbered already
    canvas->setSelection({3, 4, 5, 6, 7, 8});
    canvas->transformSelection(QTransform::fromTranslate(-30 * kBondLength, 0), "Move");  // benzene first now
    CHECK(d.texts[1].text == "1");
    CHECK(d.texts[0].text == "2");
    CHECK(textPath(d.texts[1], documentStyle(d)).boundingRect().center().x() < d.atoms[0].pos.x());  // it came along
    canvas->setSelection({});
    canvas->arrangeScheme();  // numbers stay under their molecules, not lined up between them
    for (const Text& t : d.texts) {
        const auto mol = edit::moleculeOf(d, t.anchor);
        double lo = 1e9, hi = -1e9, bottom = -1e9;
        for (int i : mol) lo = std::min(lo, d.atoms[i].pos.x()), hi = std::max(hi, d.atoms[i].pos.x()), bottom = std::max(bottom, d.atoms[i].pos.y());
        const QRectF box = textPath(t, documentStyle(d)).boundingRect();
        CHECK(box.center().x() > lo - 0.1 * kBondLength);
        CHECK(box.center().x() < hi + 0.1 * kBondLength);
        CHECK(box.top() > bottom);
    }
    if (auto out = qEnvironmentVariable("PENZENE_NUMBERS_SHOT"); !out.isEmpty()) {
        w.resize(1000, 600), w.show(), canvas->setSelection({}), canvas->fitToDocument();
        w.grab().save(out);
    }
}

TEST_CASE("NMR panel: a stick per set of equivalent atoms, lighting them on the canvas (#444)") {
    App app;
    MainWindow w;
    w.resize(1100, 700);
    w.show();
    auto* canvas = w.findChild<Canvas*>();
    canvas->setDocumentSilently(*chem::fromSmiles("CCO.Cc1ccccc1"));  // ethanol 0-2, toluene 3-9
    for (auto* a : w.findChildren<QAction*>())
        if (a->text() == "&NMR Panel") a->trigger();
    QApplication::processEvents();
    auto* dock = w.findChild<QDockWidget*>("nmr");
    REQUIRE(dock);
    REQUIRE(dock->isVisible());
    QWidget* view = nullptr;
    for (auto* c : dock->findChildren<QWidget*>())
        if (c->accessibleName() == "Predicted spectrum") view = c;
    REQUIRE(view);
    CHECK(view->accessibleDescription().split("; ").size() == 5);  // one molecule: the largest, toluene
    canvas->setSelection({1});
    CHECK(view->accessibleDescription().split("; ").size() == 2);  // the selected one, all of it: ethanol
    canvas->setSelection({0, 1, 4});
    CHECK(view->accessibleDescription().split("; ").size() == 2);  // the one with most of the selection
    canvas->setSelection({3, 4, 5, 6, 7, 8, 9});
    CHECK(view->accessibleDescription().split("; ").size() == 5);
    CHECK(view->accessibleDescription().contains("2 C: C7, C9"));  // the two meta carbons, one stick
    view->setFocus();
    QTest::keyClick(view, Qt::Key_Right);
    CHECK(view->accessibleDescription().startsWith("δ 13"));  // the ipso carbon first, as a spectrum reads
    CHECK(canvas->highlight() == QSet<int>{4});
    if (auto out = qEnvironmentVariable("PENZENE_NMR_SHOT"); !out.isEmpty()) w.grab().save(out);
    QTest::keyClick(view, Qt::Key_Right);
    CHECK(canvas->highlight().size() == 2);
    QEvent leave(QEvent::Leave);
    QApplication::sendEvent(view, &leave);  // the pointer gone: every stick is read again
    CHECK(view->accessibleDescription().split("; ").size() == 5);
    QTest::keyClick(view, Qt::Key_Right);
    QTest::keyClick(view, Qt::Key_Left);
    CHECK(canvas->highlight() == QSet<int>{4});
    dock->findChild<QPushButton*>("copyNmrSi")->click();
    CHECK(QGuiApplication::clipboard()->text().startsWith("13C NMR (predicted) δ 13"));  // toluene's, not ethanol's
    CHECK(QGuiApplication::clipboard()->text().count(", ") == 4);
    canvas->setHotspot(5);  // and the other way: an ortho carbon under the pointer marks its stick
    CHECK(canvas->highlight() == QSet<int>{5, 9});
    CHECK(view->accessibleDescription().startsWith("δ 12"));
    canvas->setHotspot(-1);
    CHECK(canvas->highlight().isEmpty());
    dock->findChild<QComboBox*>()->setCurrentIndex(1);  // 1H
    CHECK(canvas->highlight().isEmpty());
    CHECK(view->accessibleDescription().split("; ").last().contains("3 H, s: C4"));  // the methyl, furthest upfield
    QTest::keyClick(view, Qt::Key_Right);
    CHECK_FALSE(canvas->highlight().isEmpty());
    canvas->commit(*chem::fromSmiles("CCO"), "Replace");  // an edit (or undo) that drops atoms: nothing stale stays lit
    CHECK(canvas->highlight().isEmpty());
    CHECK(view->accessibleDescription().contains("3 H, t: C1"));  // first order (#552)
    if (auto out = qEnvironmentVariable("PENZENE_NMR_SHOT"); !out.isEmpty()) w.grab().save(QString(out).replace(".png", "-1h.png"));
    QTest::keyClick(view, Qt::Key_Right);
    REQUIRE_FALSE(canvas->highlight().isEmpty());
    for (int i : canvas->highlight()) CHECK(i < 3);
    canvas->commit(*chem::fromSmiles("[Na+].[Cl-]"), "Replace");
    dock->findChild<QPushButton*>("copyNmrSi")->click();  // nothing predicted: the clipboard keeps what it had
    CHECK(QGuiApplication::clipboard()->text().startsWith("13C NMR"));
    dock->hide();
    CHECK(canvas->highlight().isEmpty());
}

TEST_CASE("Experimental JCAMP overlay: import, nucleus, removal and failed replacement (#553)") {
    App app;
    QTemporaryDir dir;
    MainWindow w;
    w.resize(1100, 700);
    w.show();
    auto* canvas = w.findChild<Canvas*>();
    canvas->setDocumentSilently(*chem::fromSmiles("CCO"));
    for (auto* action : w.findChildren<QAction*>())
        if (action->text() == "&NMR Panel") action->trigger();
    QApplication::processEvents();
    auto* dock = w.findChild<QDockWidget*>("nmr");
    auto* view = dock->findChild<QWidget*>("nmrSpectrum");
    auto* label = dock->findChild<QLabel*>("experimentalSpectrum");
    auto* remove = dock->findChild<QPushButton*>("removeExperimental");
    auto* nucleus = dock->findChild<QComboBox*>();
    REQUIRE(view);
    REQUIRE(label);
    REQUIRE(remove);
    CHECK_FALSE(remove->isEnabled());
    const QString preview = qEnvironmentVariable("PENZENE_JCAMP_PREVIEW");
    if (!preview.isEmpty()) {
        QDir().mkpath(preview);
        canvas->fitToDocument();
        REQUIRE(w.grab().save(preview + "/before.png"));
    }
    QByteArray measured = "##TITLE=synthetic ethanol\n##DATA TYPE=NMR SPECTRUM\n##.OBSERVE NUCLEUS=1H\n"
                          "##XUNITS=PPM\n##NPOINTS=1201\n##XYPOINTS=(XY..XY)\n";
    for (int i = 0; i <= 1200; ++i) {
        const double ppm = 12 - i / 100.0;
        const double intensity = std::exp(-std::pow((ppm - 3.55) / 0.10, 2)) +
                                 1.5 * std::exp(-std::pow((ppm - 1.12) / 0.08, 2));
        measured += QByteArray::number(ppm) + "," + QByteArray::number(intensity) + "\n";
    }
    measured += "##END=\n";
    auto write = [&](const QString& name, const QByteArray& data) {
        const QString path = dir.filePath(name);
        QFile f(path);
        REQUIRE(f.open(QIODevice::WriteOnly));
        REQUIRE(f.write(data) == data.size());
        return path;
    };
    QApplication::setAttribute(Qt::AA_DontUseNativeDialogs);
    auto import = [&](const QString& path, bool invalid = false, int choice = -1) {
        QTimer::singleShot(0, &w, [&, path, invalid, choice] {
            auto* files = qobject_cast<QFileDialog*>(QApplication::activeModalWidget());
            REQUIRE(files);
            if (path.isEmpty()) return files->reject();
            if (invalid || choice >= 0) QTimer::singleShot(0, &w, [&, invalid, choice] {
                if (invalid) {
                    auto* warning = qobject_cast<QMessageBox*>(QApplication::activeModalWidget());
                    REQUIRE(warning);
                    CHECK_FALSE(warning->text().isEmpty());
                    warning->accept();
                } else {
                    auto* dialog = qobject_cast<QInputDialog*>(QApplication::activeModalWidget());
                    REQUIRE(dialog);
                    dialog->findChild<QComboBox*>()->setCurrentIndex(choice);
                    dialog->accept();
                }
            });
            files->selectFile(path);
            static_cast<QDialog*>(files)->accept();
        });
        dock->findChild<QPushButton*>("importExperimental")->click();
        QApplication::processEvents();
    };
    import(write("ethanol.jdx", measured));
    CHECK(nucleus->currentIndex() == 1);
    CHECK(label->text().contains("ethanol.jdx"));
    CHECK(view->accessibleDescription().contains("Experimental synthetic ethanol: 1201 samples"));
    CHECK(remove->isEnabled());
    const auto imported = view->grab().toImage();
    if (!preview.isEmpty()) REQUIRE(w.grab().save(preview + "/after.png"));
    dock->findChild<QPushButton*>("copyNmrSi")->click();
    CHECK(QGuiApplication::clipboard()->text().startsWith("1H NMR (predicted)"));
    const auto identity = label->text();
    import(write("broken.dx", "##TITLE=broken\n"), true);
    CHECK(label->text() == identity);
    CHECK(view->grab().toImage() == imported);
    import({});
    CHECK(label->text() == identity);
    nucleus->setCurrentIndex(0);
    CHECK(label->text().contains("matching nucleus"));
    CHECK_FALSE(view->accessibleDescription().contains("Experimental"));
    nucleus->setCurrentIndex(1);
    canvas->commit({}, "Clear drawing");
    CHECK(view->accessibleDescription().contains("Experimental synthetic ethanol"));
    CHECK_FALSE(view->grab().isNull());
    remove->click();
    CHECK(label->text().isEmpty());
    CHECK_FALSE(remove->isEnabled());
    const QByteArray carbon = "##TITLE=carbon\n##DATA TYPE=NMR PEAK TABLE\n##.OBSERVE NUCLEUS=13C\n"
                              "##XUNITS=PPM\n##PEAK TABLE=(XY..XY)\n77,1;60,-0.5;20,?;15,2\n##END=\n";
    import(write("linked.jcamp", "##TITLE=two\n##DATA TYPE=LINK\n" + measured + carbon + "##END=\n"), false, 1);
    CHECK(nucleus->currentIndex() == 0);
    CHECK(view->accessibleDescription().contains("Experimental carbon: 4 samples"));
    CHECK_FALSE(view->grab().isNull());
    const QString drawing = write("fresh.penz", Document().toJson());
    REQUIRE(w.openFile(drawing));
    CHECK(label->text().isEmpty());
    CHECK_FALSE(remove->isEnabled());
    QApplication::setAttribute(Qt::AA_DontUseNativeDialogs, false);
}

TEST_CASE("PubChem name lookup: URL and response parsing (#28)") {
    CHECK(pubchem::nameToSmilesUrl(" acetylsalicylic acid ").toString(QUrl::FullyEncoded) ==
          "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/acetylsalicylic%20acid/property/SMILES/JSON");
    const QByteArray found = R"({"PropertyTable": {"Properties": [{"CID": 2244, "SMILES": "CC(=O)OC1=CC=CC=C1C(=O)O"}]}})";
    CHECK(pubchem::property(found, "SMILES") == "CC(=O)OC1=CC=CC=C1C(=O)O");
    CHECK(pubchem::property(found, "IUPACName").isEmpty());
    CHECK(pubchem::property(R"({"Fault": {"Code": "PUGREST.NotFound"}})", "SMILES").isEmpty());
    CHECK(pubchem::property("not json", "SMILES").isEmpty());
}

TEST_CASE("PubChem structure lookup: POSTed SMILES and the IUPAC name (#27)") {
    CHECK(pubchem::smilesToNameUrl().toString() ==
          "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/property/IUPACName/JSON");
    CHECK(pubchem::smilesToNameForm("C#N/C=C/O") == "smiles=C%23N%2FC%3DC%2FO");
    CHECK(pubchem::property(R"({"PropertyTable": {"Properties": [{"CID": 2244, "IUPACName": "2-acetyloxybenzoic acid"}]}})",
                            "IUPACName") == "2-acetyloxybenzoic acid");
}

TEST_CASE("x and r label an atom X and R; free-text labels keep unspecified chemistry (#156)") {
    Document d = *chem::fromSmiles("CCO");
    CHECK(edit::hotkey(d, {0, -1}, "x").atom == 0);
    CHECK(edit::hotkey(d, {2, -1}, "r").atom == 2);
    CHECK(d.atoms[0].label == "X");
    CHECK(d.atoms[2].label == "R");
    CHECK(d.atoms[0].z == 0);
    CHECK(d.atoms.size() == 3);  // nothing invented
    CHECK(chem::toSmiles(d) == "*C*");

    CHECK_FALSE(edit::applyLabel(d, 1, "MgEt"));  // strict unless asked (the Python API)
    REQUIRE(edit::applyLabel(d, 1, "MgEt", true));
    CHECK(d.atoms[1].label == "MgEt");
    auto back = Document::fromJson(d.toJson());
    REQUIRE(back);
    CHECK(back->atoms[1].label == "MgEt");
    CHECK(edit::applyLabel(d, 1, "OMe", true));  // real groups still win
    CHECK(d.atoms[1].z == 8);

    // With no hotspot, x is still the bond tool; on a bond, r still places the double bond.
    Document e = *chem::fromSmiles("C=C");
    CHECK(edit::hotkey(e, {-1, 0}, "r").valid());
    CHECK(e.atoms[0].label.isEmpty());
}

TEST_CASE("the window title names the build (#165)") {
    App app;
    MainWindow w;
    CHECK(w.windowTitle().endsWith("Penzene " PENZENE_BUILD));
    CHECK(QString(PENZENE_BUILD).startsWith(PENZENE_VERSION));
}

TEST_CASE("a PubChem lookup holds user input until it returns (#171)") {
    App app;
    QTcpServer server;  // accepts, stays silent, then hangs up
    REQUIRE(server.listen(QHostAddress::LocalHost));
    QObject::connect(&server, &QTcpServer::newConnection, [&] {
        QTcpSocket* s = server.nextPendingConnection();
        QTimer::singleShot(300, s, [s] { s->close(); });
    });
    QLineEdit typing;
    typing.show();
    QTimer::singleShot(50, [&] {
        // As the window system delivers it (a posted QKeyEvent would bypass the filter).
        QWindowSystemInterface::handleKeyEvent(typing.windowHandle(), QEvent::KeyPress, Qt::Key_A, Qt::NoModifier, "a");
    });
    QString error;
    const QUrl url(QString("http://127.0.0.1:%1/").arg(server.serverPort()));
    CHECK(pubchem::fetch(url, "SMILES", &error).isEmpty());
    CHECK(typing.text().isEmpty());  // not handled mid-request
    QApplication::processEvents();
    CHECK(typing.text() == "a");  // delivered afterwards
}

TEST_CASE("exported SVG and PNG reopen as the editable drawing (#100)") {
    App app;
    Document doc = *chem::fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
    doc.atoms[0].color = Qt::red;
    doc.atoms[3].map = 4;
    doc.texts.push_back({{0, 60}, "aspirin"});
    QTemporaryDir dir;
    for (const char* ext : {"svg", "png"}) {
        const QString path = dir.filePath(QString("aspirin.") + ext);
        REQUIRE(exportDocument(doc, path));
        auto back = chem::readFile(path);
        REQUIRE(back);
        CHECK(*back == doc);
    }
    CHECK_FALSE(Document::fromEmbedded(QByteArray("<svg xmlns='http://www.w3.org/2000/svg'/>")));
    CHECK_FALSE(Document::fromEmbedded(QByteArray()));

    // A copied figure pastes back as structure, not as a picture.
    MainWindow w;
    auto* mime = new QMimeData;
    mime->setData("image/svg+xml", renderSvg(doc));
    QApplication::clipboard()->setMimeData(mime);
    w.findChild<Canvas*>()->setDocumentSilently(Document{});
    for (auto* a : w.findChildren<QAction*>())
        if (a->shortcut() == QKeySequence::Paste) a->trigger();
    CHECK(w.findChild<Canvas*>()->document().atoms.size() == doc.atoms.size());
}

TEST_CASE("copied PNG keeps the embedded drawing and pastes back (#399)") {
    App app;
    MainWindow w;
    auto* canvas = w.findChild<Canvas*>();
    Document doc = *chem::fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
    doc.texts.push_back({{0, 60}, "aspirin"});
    canvas->setDocumentSilently(doc);
    canvas->selectAll();
    const Document copied = canvas->selectedSubset();
    for (auto* a : w.findChildren<QAction*>())
        if (a->shortcut() == QKeySequence::Copy) a->trigger();
    const QByteArray png = QApplication::clipboard()->mimeData()->data("image/png");
    REQUIRE(png.startsWith("\x89PNG"));
    auto embedded = Document::fromEmbedded(png);
    REQUIRE(embedded);
    CHECK(*embedded == copied);
    CHECK(QApplication::clipboard()->mimeData()->hasImage());  // still a bitmap for apps that want one
    CHECK(!QApplication::clipboard()->mimeData()->hasText());  // else Office pastes the SMILES, not the picture
#ifdef Q_OS_WIN
    CHECK(QApplication::clipboard()->mimeData()->data("application/x-qt-windows-mime;value=\"PNG\"") == png);
#endif
#ifdef Q_OS_MACOS
    ChemDrawPasteboard uti;
    CHECK(uti.utiForMime("image/png") == "public.png");
    CHECK(uti.mimeForUti("public.png") == "image/png");
#endif

    // The PNG alone, as another app (Office) hands it back, pastes as the drawing.
    for (const char* type : {"image/png", "application/x-qt-windows-mime;value=\"PNG\""}) {
        INFO(type);
        auto* mime = new QMimeData;
        mime->setData(type, png);
        QApplication::clipboard()->setMimeData(mime);
        canvas->setDocumentSilently(Document{});
        for (auto* a : w.findChildren<QAction*>())
            if (a->shortcut() == QKeySequence::Paste) a->trigger();
        CHECK(canvas->document().atoms.size() == doc.atoms.size());
        CHECK(canvas->document().texts.size() == 1);
    }
}

TEST_CASE("copied pictures measure what is painted, not the label room (#681)") {
    App app;
    QSettings settings;
    const QVariant oldMargin = settings.value("exportMargin");
    const QVariant oldBackground = settings.value("exportBackground");
    const QVariant oldScale = settings.value("exportScale");
    const QVariant oldCheck = settings.value("checkBeforeExport");
    settings.setValue("exportMargin", 0);
    settings.setValue("exportBackground", "clear");
    settings.setValue("exportScale", 100);
    settings.setValue("checkBeforeExport", false);
    MainWindow w;
    auto* canvas = w.findChild<Canvas*>();
    canvas->setDocumentSilently(*chem::fromSmiles("CC(=O)Oc1ccccc1C(=O)O"));
    canvas->selectAll();
    auto copied = [&] {
        for (auto* a : w.findChildren<QAction*>())
            if (a->shortcut() == QKeySequence::Copy) a->trigger();
        return QImage::fromData(QApplication::clipboard()->mimeData()->data("image/png")).convertToFormat(QImage::Format_ARGB32);
    };
    const QImage tight = copied();
    REQUIRE(!tight.isNull());
    const double dpi = tight.dotsPerMeterX() * 0.0254;
    int l = tight.width(), t = tight.height(), r = -1, b = -1;
    for (int y = 0; y < tight.height(); ++y)
        for (int x = 0; x < tight.width(); ++x)
            if (qAlpha(tight.pixel(x, y))) l = std::min(l, x), r = std::max(r, x), t = std::min(t, y), b = std::max(b, y);
    REQUIRE(r >= 0);
    const int slack = int(std::ceil(2 * dpi / 72));  // the painted box is whole model units, padded by one
    CHECK(l <= slack);
    CHECK(t <= slack);
    CHECK(tight.width() - 1 - r <= slack);
    CHECK(tight.height() - 1 - b <= slack);
    settings.setValue("exportMargin", 10);
    const QImage padded = copied();
    const int add = int(20 * dpi / 72);
    CHECK(std::abs(padded.width() - tight.width() - add) <= 1);  // the setting still adds its margin
    CHECK(std::abs(padded.height() - tight.height() - add) <= 1);
    auto restore = [&](const char* key, const QVariant& value) {
        if (value.isValid()) settings.setValue(key, value);
        else settings.remove(key);
    };
    restore("exportMargin", oldMargin);
    restore("exportBackground", oldBackground);
    restore("exportScale", oldScale);
    restore("checkBeforeExport", oldCheck);
}

TEST_CASE("paste from ChemDraw: CDX/CDXML clipboard formats (#104)") {
    App app;
    const QString path = QString(PENZENE_TEST_DATA) + "/scheme.cdxml";
    QFile f(path);
    REQUIRE(f.open(QIODevice::ReadOnly));
    const QByteArray cdxml = f.readAll();
    const auto expected = chem::readFile(path);
    REQUIRE(expected);
    // macOS (through ChemDrawPasteboard) and Windows name the format differently.
    for (const char* type : {"chemical/x-cdx", "application/x-qt-windows-mime;value=\"ChemDraw Interchange Format\""}) {
        INFO(type);
        MainWindow w;
        auto* mime = new QMimeData;
        mime->setData(type, cdxml);
        mime->setText("not a structure");  // ChemDraw also offers text; the CDX wins
        QApplication::clipboard()->setMimeData(mime);
        for (auto* a : w.findChildren<QAction*>())
            if (a->shortcut() == QKeySequence::Paste) a->trigger();
        CHECK(w.findChild<Canvas*>()->document().atoms.size() == expected->atoms.size());
    }
#ifdef Q_OS_MACOS
    ChemDrawPasteboard uti;
    CHECK(uti.mimeForUti("com.perkinelmer.chemdraw.cdx-clipboard") == "chemical/x-cdx");
    CHECK(uti.mimeForUti("public.utf8-plain-text").isEmpty());
    CHECK(uti.convertToMime("chemical/x-cdx", {cdxml}, {}).toByteArray() == cdxml);
#endif
}

TEST_CASE("export scale and margin (#103)") {
    App app;
    const Document doc = *chem::fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
    const QImage full = renderImage(doc, {72});
    const QImage half = renderImage(doc, {72, Qt::transparent, 0.5});
    const QImage padded = renderImage(doc, {72, Qt::white, 1, 10});
    CHECK(std::abs(half.width() - full.width() / 2) <= 1);
    CHECK(std::abs(half.height() - full.height() / 2) <= 1);
    CHECK(std::abs(padded.width() - (full.width() + 20)) <= 1);  // 10 pt a side at 72 dpi
    CHECK(padded.pixelColor(2, 2) == QColor(Qt::white));        // the margin takes the background
    const QByteArray svg = renderSvg(doc, {300, Qt::transparent, 0.85});
    CHECK(svg.contains("<svg"));
    CHECK(svg.size() > 100);
}

TEST_CASE("exports hold atom numbers, lone pairs, δ and stereo labels whole (#325)") {
    App app;
    Document d = *chem::fromSmiles("C[C@H](N)C(=O)O");
    d.showAtomNumbers = d.showStereo = true;
    for (auto& a : d.atoms)
        if (a.z == 8) a.lonePairs = 2, a.partial = -1;
    QTemporaryDir dir;
    const QString path = qEnvironmentVariable("PENZENE_BOUNDS_SHOT", dir.filePath("x.png"));
    REQUIRE(exportDocument(d, path));
    const QImage img = QImage(path).convertToFormat(QImage::Format_ARGB32);
    REQUIRE(!img.isNull());
    int edge = 0;  // anything painted on the image's border was cut off
    for (int x = 0; x < img.width(); ++x) edge += qAlpha(img.pixel(x, 0)) + qAlpha(img.pixel(x, img.height() - 1));
    for (int y = 0; y < img.height(); ++y) edge += qAlpha(img.pixel(0, y)) + qAlpha(img.pixel(img.width() - 1, y));
    CHECK(edge == 0);
}

TEST_CASE("printing: to PDF, at export size, centred (#33)") {
    App app;
    QTemporaryDir dir;
    QPrinter printer(QPrinter::HighResolution);
    printer.setOutputFormat(QPrinter::PdfFormat);
    printer.setOutputFileName(dir.filePath("aspirin.pdf"));
    CHECK(printDocument(printer, *chem::fromSmiles("CC(=O)Oc1ccccc1C(=O)O")));
    CHECK(QFileInfo(dir.filePath("aspirin.pdf")).size() > 1000);
    if (auto out = qgetenv("PENZENE_PRINT_SHOT"); !out.isEmpty()) QFile::copy(dir.filePath("aspirin.pdf"), out);
    CHECK_FALSE(printDocument(printer, Document{}));
}

TEST_CASE("page mode: shown, saved, exported at page size (#105)") {
    App app;
    MainWindow w;
    auto* canvas = w.findChild<Canvas*>();
    canvas->setDocumentSilently(*chem::fromSmiles("CC(=O)Oc1ccccc1C(=O)O"));
    QAction* column = nullptr;
    for (auto* a : w.findChildren<QAction*>())
        if (a->text() == "ACS single column") column = a;
    REQUIRE(column);
    column->trigger();
    const Document& doc = canvas->document();
    CHECK(doc.page == "ACS single column");
    const QRectF page = pageRect(doc);
    CHECK(page.contains(documentBounds(doc).center()));
    auto back = Document::fromJson(doc.toJson());
    REQUIRE(back);
    CHECK(*back == doc);
    const QImage img = renderImage(doc, {72});  // the page: 240 × 684 pt
    CHECK(std::abs(img.width() - 240) <= 1);
    CHECK(std::abs(img.height() - 684) <= 1);
    canvas->selectAll();
    CHECK(canvas->selectedSubset().page.isEmpty());  // a selection exports just the drawing
    canvas->setSelection({0});
    canvas->deleteSelection();
    CHECK(canvas->document().page == "ACS single column");  // but deleting from it keeps the page (#492)
    if (auto out = qgetenv("PENZENE_PAGE_SHOT"); !out.isEmpty()) {
        w.resize(900, 900);
        w.show();
        canvas->setSelection({});
        QApplication::processEvents();
        w.grab().save(out);
    }
}

TEST_CASE("stretch, squash and transform a selection (#174)") {
    Fixture f;
    auto size = [](const Document& d) {
        QPolygonF pts;
        for (const auto& a : d.atoms) pts << a.pos;
        return pts.boundingRect().size();
    };
    const Document benzene = *chem::fromSmiles("c1ccccc1");
    f.canvas.setDocumentSilently(benzene);
    const QSizeF before = size(f.doc());
    f.canvas.transformSelection(QTransform::fromScale(1, 0.5), "Squash");  // nothing selected: everything
    CHECK(std::abs(size(f.doc()).height() - before.height() / 2) < 0.01);
    CHECK(std::abs(size(f.doc()).width() - before.width()) < 0.01);
    CHECK(chem::toSmiles(f.doc()) == chem::toSmiles(benzene));
    f.undo.undo();
    CHECK(f.doc() == benzene);

    // Handles: drag the bottom edge down to stretch, then a corner out to scale.
    f.canvas.setTool(Canvas::Tool::Select);
    f.canvas.selectAll();
    QPolygonF pts;
    for (const auto& a : f.doc().atoms) pts << a.pos;
    const QRectF box = pts.boundingRect().adjusted(-6, -6, 6, 6);
    const QPointF bottom(box.center().x(), box.bottom());
    f.drag(bottom, bottom + QPointF(0, box.height()));
    CHECK(std::abs(size(f.doc()).width() - before.width()) < 0.5);
    CHECK(size(f.doc()).height() > 1.8 * before.height());
    f.undo.undo();
    f.canvas.selectAll();
    f.drag(box.bottomRight(), box.bottomRight() + (box.bottomRight() - box.topLeft()));  // twice the diagonal
    CHECK(std::abs(size(f.doc()).width() / before.width() - size(f.doc()).height() / before.height()) < 0.05);
    CHECK(size(f.doc()).width() > 1.8 * before.width());
}

TEST_CASE("lone pairs, radicals and brackets from the keyboard and menu (#106)") {
    Fixture f;
    f.canvas.setDocumentSilently(*chem::fromSmiles("CCO"));
    Document d = f.doc();
    REQUIRE(edit::hotkey(d, {2, -1}, ":").valid());
    edit::hotkey(d, {2, -1}, ":");
    CHECK(d.atoms[2].lonePairs == 2);
    edit::hotkey(d, {0, -1}, "*");
    CHECK(d.atoms[0].radicals == 1);
    CHECK(chem::toSmiles(d) == "[CH2]CO");

    f.canvas.setSelection({0, 1});
    f.canvas.bracketSelection(true, "n");
    REQUIRE(f.doc().brackets.size() == 1);
    CHECK(f.doc().brackets[0].atoms == std::vector<int>{0, 1});
    Document moved = f.doc();
    moved.removeAtoms({0});  // a bracket keeps around what's left
    CHECK(moved.brackets[0].atoms == std::vector<int>{0});
    f.canvas.removeBrackets();
    CHECK(f.doc().brackets.empty());
}

TEST_CASE("template library: built-ins insert, selections save as templates (#32)") {
    App app;
    REQUIRE(builtinTemplates().size() >= 50);
    for (const auto& t : builtinTemplates()) {
        INFO(t.name.toStdString());
        CHECK(chem::fromSmiles(t.smiles.toStdString()));
    }
    MainWindow w;
    w.resize(1000, 700);
    w.show();
    auto* canvas = w.findChild<Canvas*>();
    canvas->setDocumentSilently(Document{});
    for (auto* a : w.findChildren<QAction*>())
        if (a->text() == "&Templates") a->trigger();
    QApplication::processEvents();
    auto* tree = w.findChild<QTreeWidget*>();
    REQUIRE(tree);
    auto found = tree->findItems("L-Alanine", Qt::MatchExactly | Qt::MatchRecursive);
    REQUIRE(found.size() == 1);
    emit tree->itemClicked(found[0], 0);
    CHECK(chem::toSmiles(canvas->document()) == chem::toSmiles(*chem::fromSmiles("C[C@@H](C(=O)O)N")));

    REQUIRE(saveUserTemplate("My alanine", canvas->document()));
    const auto mine = userTemplates();
    CHECK(std::any_of(mine.begin(), mine.end(), [](const auto& t) { return t.first == "My alanine"; }));
    CHECK(removeUserTemplate("My alanine"));
    if (auto out = qgetenv("PENZENE_TEMPLATES_SHOT"); !out.isEmpty()) {
        tree->findItems("Sugars", Qt::MatchExactly)[0]->setExpanded(true);
        QApplication::processEvents();
        w.grab().save(out);
    }
}

TEST_CASE("templates are reachable from the welcome card and the tool rail (#679)") {
    App app;
    MainWindow w;
    w.resize(1000, 700);
    w.show();
    auto* dock = w.findChild<QDockWidget*>("templates");
    REQUIRE(dock);
    CHECK_FALSE(dock->isVisible());
    QToolButton* library = nullptr;
    for (auto* b : w.findChild<QFrame*>("toolCard")->findChildren<QToolButton*>("railButton"))
        if (b->text() == "Templates") library = b;
    REQUIRE(library);
    CHECK(library->isCheckable());
    library->click();
    QApplication::processEvents();
    CHECK(dock->isVisible());
    CHECK(library->isChecked());
    auto* tree = dock->findChild<QTreeWidget*>();
    REQUIRE(tree);
    CHECK(tree->topLevelItemCount() > 0);  // filled the first time it opens
    library->click();
    CHECK_FALSE(dock->isVisible());
    CHECK_FALSE(library->isChecked());

    auto* welcome = w.findChild<QLabel*>("welcomeLinks");
    REQUIRE(welcome);
    CHECK(welcome->text().contains("Templates"));
    emit welcome->linkActivated("templates");
    QApplication::processEvents();
    CHECK(dock->isVisible());
    CHECK(library->isChecked());  // the View toggle and the rail button agree

    QAction* viewToggle = nullptr;
    for (auto* a : w.findChildren<QAction*>())
        if (a->text() == "&Templates") viewToggle = a;
    REQUIRE(viewToggle);
    CHECK(viewToggle->shortcut() == QKeySequence("Ctrl+Shift+T"));
    viewToggle->trigger();
    CHECK_FALSE(dock->isVisible());
    CHECK_FALSE(library->isChecked());
}

TEST_CASE("Added templates retain their chemical identities and bond geometry (#474)") {
    App app;
    QFile file(QString(PENZENE_TEST_DATA) + "/template-identities.tsv");
    REQUIRE(file.open(QIODevice::ReadOnly | QIODevice::Text));
    Document gallery;
    int checked = 0;
    for (const auto& line : QString::fromUtf8(file.readAll()).split('\n', Qt::SkipEmptyParts)) {
        if (line.startsWith('#')) continue;
        const auto columns = line.split('\t');
        REQUIRE(columns.size() == 5);
        INFO(columns[1].toStdString());
        const auto& templates = builtinTemplates();
        const auto found = std::find_if(templates.begin(), templates.end(), [&](const Template& t) { return t.name == columns[1]; });
        REQUIRE(found != templates.end());
        CHECK(found->category == columns[0]);
        CHECK(std::count_if(templates.begin(), templates.end(), [&](const Template& t) { return t.name == columns[1]; }) == 1);
        const Document doc = templateDocument(*found);
        const auto properties = chem::properties(doc);
        REQUIRE(properties);
        CHECK(properties->formula == columns[2].toStdString());
        CHECK(chem::toInchiKey(doc) == columns[3].toStdString());
        for (const auto& bond : doc.bonds) {
            const QLineF edge(doc.atoms[bond.a].pos, doc.atoms[bond.b].pos);
            CHECK(std::abs(edge.length() - kBondLength) < 0.1 * kBondLength);
            for (const auto& other : doc.bonds)
                if (bond.a != other.a && bond.a != other.b && bond.b != other.a && bond.b != other.b)
                    CHECK(edge.intersects(QLineF(doc.atoms[other.a].pos, doc.atoms[other.b].pos)) != QLineF::BoundedIntersection);
        }
        const QPointF at((checked % 4) * 340, (checked / 4) * 240);
        gallery.append(doc, at - documentBounds(doc).center());
        gallery.texts.push_back({at + QPointF(-150, 100), found->name, 0.8});
        ++checked;
    }
    CHECK(checked == 24);
    const QString preview = qEnvironmentVariable("PENZENE_MORE_TEMPLATES_PREVIEW");
    if (!preview.isEmpty()) {
        QDir().mkpath(preview);
        REQUIRE(renderImage(gallery, {110, Qt::white}).save(preview + "/gallery.png"));
    }
    MainWindow w;
    w.resize(1100, 750);
    w.show();
    auto* canvas = w.findChild<Canvas*>();
    canvas->setDocumentSilently({});
    for (auto* action : w.findChildren<QAction*>())
        if (action->text() == "&Templates") action->trigger();
    QApplication::processEvents();
    auto* tree = w.findChild<QDockWidget*>("templates")->findChild<QTreeWidget*>();
    REQUIRE(tree);
    tree->collapseAll();
    const auto macrocycles = tree->findItems("Macrocycles", Qt::MatchExactly);
    REQUIRE(macrocycles.size() == 1);
    macrocycles[0]->setExpanded(true);
    tree->scrollToItem(macrocycles[0]);
    if (!preview.isEmpty()) REQUIRE(w.grab().save(preview + "/before.png"));
    const auto found = tree->findItems("Cyclen", Qt::MatchExactly | Qt::MatchRecursive);
    REQUIRE(found.size() == 1);
    emit tree->itemClicked(found[0], 0);
    CHECK(chem::properties(canvas->document())->formula == "C8H20N4");
    canvas->fitToDocument();
    if (!preview.isEmpty()) REQUIRE(w.grab().save(preview + "/after.png"));
}

TEST_CASE("crown ethers sit in the ring's standard drawing (#680)") {
    App app;
    for (const auto& [n, name] : {std::pair{4, "12-Crown-4"}, std::pair{5, "15-Crown-5"}, std::pair{6, "18-Crown-6"}}) {
        INFO(name);
        const auto& templates = builtinTemplates();
        const auto t = std::find_if(templates.begin(), templates.end(), [&](const Template& t) { return t.name == name; });
        REQUIRE(t != templates.end());
        const Document d = templateDocument(*t);
        if (auto dir = qgetenv("PENZENE_CROWN_PREVIEW"); !dir.isEmpty())
            REQUIRE(renderImage(d, {600, Qt::white}).save(QString::fromUtf8(dir) + "/" + name + ".png"));
        CHECK(chem::toInchiKey(d) == chem::toInchiKey(*chem::fromSmiles(t->smiles.toStdString())));
        std::vector<int> o, c;
        for (int i = 0; i < int(d.atoms.size()); ++i) (d.atoms[i].z == 8 ? o : c).push_back(i);
        REQUIRE(o.size() == size_t(n));
        REQUIRE(c.size() == size_t(2 * n));
        QPointF centre;
        for (int i : o) centre += d.atoms[i].pos;
        centre /= n;
        auto angle = [&](int i) { return std::atan2(d.atoms[i].pos.y() - centre.y(), d.atoms[i].pos.x() - centre.x()); };
        auto radius = [&](int i) { return QLineF(centre, d.atoms[i].pos).length(); };
        std::sort(o.begin(), o.end(), [&](int a, int b) { return angle(a) < angle(b); });
        // The oxygens are evenly spaced on one circle.
        for (int i : o) CHECK(std::abs(radius(i) - radius(o.front())) < 1e-6);
        for (size_t k = 0; k < o.size(); ++k) {
            double gap = angle(o[(k + 1) % o.size()]) - angle(o[k]);
            if (gap < 0) gap += 2 * M_PI;
            CHECK(std::abs(gap - 2 * M_PI / n) < 1e-6);
        }
        // Each CH2CH2 bridge is two carbons outside the oxygens' circle, one O and one C
        // neighbour each, and its C-O bonds turn inwards to meet the oxygen.
        for (int i : c) {
            std::vector<int> co, cc;
            for (const Bond& b : d.bonds) {
                const int j = b.a == i ? b.b : b.b == i ? b.a : -1;
                if (j < 0) continue;
                (d.atoms[j].z == 8 ? co : cc).push_back(j);
            }
            CHECK(co.size() == 1);
            REQUIRE(cc.size() == 1);
            CHECK(radius(i) > radius(o.front()));
            CHECK(QLineF(centre, (d.atoms[i].pos + d.atoms[cc[0]].pos) / 2).length() > radius(o.front()));
            const QPointF inward = centre - d.atoms[i].pos, toO = d.atoms[co[0]].pos - d.atoms[i].pos;
            CHECK(QPointF::dotProduct(inward, toO) > 0);
        }
    }
}

TEST_CASE("Haworth and Fischer projections read with PubChem's stereo (#108)") {
    App app;
    int checked = 0;
    for (const auto& t : builtinTemplates()) {
        if (t.category != "Projections" || t.smiles.isEmpty()) continue;
        INFO(t.name.toStdString());
        const Document doc = templateDocument(t);
        REQUIRE_FALSE(doc.empty());
        // Through SMILES and back, which drops the Fischer drawings' explicit H.
        CHECK(chem::toSmiles(*chem::fromSmiles(chem::toSmiles(doc))) == chem::toSmiles(*chem::fromSmiles(t.smiles.toStdString())));
        ++checked;
    }
    CHECK(checked == 5);
    // A Fischer centre needs its four bonds exactly on the axes; tilt one and it's unspecified again.
    Document glyceraldehyde = templateDocument(*std::find_if(builtinTemplates().begin(), builtinTemplates().end(),
                                                             [](const Template& t) { return t.name.contains("glyceraldehyde"); }));
    for (auto& a : glyceraldehyde.atoms)
        if (a.z == 8 && a.pos.x() > 1) a.pos += QPointF(0, 5);
    CHECK(chem::toSmiles(glyceraldehyde).find('@') == std::string::npos);
    // An ordinary ring drawn with a bold bond isn't a Haworth ring.
    Document cyclohexanol = *chem::fromSmiles("OC1CCCCC1");
    cyclohexanol.bonds[2].stereo = BondStereo::Bold;
    CHECK(chem::toSmiles(cyclohexanol) == "OC1CCCCC1");
}

TEST_CASE("The anti-butane Newman template puts its methyls 180 degrees apart (#604)") {
    App app;
    const auto& templates = builtinTemplates();
    const auto t = std::find_if(templates.begin(), templates.end(), [](const Template& t) { return t.name == "Newman: butane, anti (staggered)"; });
    REQUIRE(t != templates.end());
    const auto doc = templateDocument(*t);
    std::vector<QPointF> methyls;
    for (const auto& text : doc.texts)
        if (text.text == "CH3") methyls.push_back(text.pos + QPointF(2.9 * text.text.size(), -3.5));
    REQUIRE(methyls.size() == 2);
    CHECK(QLineF(methyls[0], -methyls[1]).length() < 1e-6);
    CHECK(doc.atoms.empty());
    CHECK(chem::toSmiles(doc).empty());
}

TEST_CASE("Arrange Scheme lines a reaction up (#109)") {
    Fixture f;
    // Salicylic acid + acetic anhydride → aspirin, drawn untidily, pyridine over a tilted arrow.
    Document d = *chem::fromSmiles("OC(=O)c1ccccc1O");
    d.append(*chem::fromSmiles("CC(=O)OC(C)=O"), {70, 25});
    d.texts.push_back({{35, 5}, "+"});
    d.arrows.push_back({{110, 10}, {150, 0}});
    d.append(*chem::fromSmiles("c1ccncc1"), {128, -35});  // the agent, over the arrow
    d.append(*chem::fromSmiles("CC(=O)Oc1ccccc1C(=O)O"), {200, 30});
    f.canvas.setDocumentSilently(d);
    f.canvas.arrangeScheme();
    const Document& out = f.doc();
    const Arrow& arrow = out.arrows[0];
    CHECK(std::abs(arrow.from.y() - arrow.to.y()) < 1e-9);  // straightened
    auto centre = [&](const QString& smiles) {  // a molecule's centre, found by its SMILES
        QPolygonF pts;
        const auto want = chem::toSmiles(*chem::fromSmiles(smiles.toStdString()));
        const auto parts = *chem::reactionOf(out);
        for (const auto& r : {parts.reactants, parts.agents, parts.products})
            for (const auto& m : r)
                if (chem::toSmiles(m) == want) {
                    for (const auto& a : m.atoms) pts << a.pos;
                    return pts.boundingRect().center();
                }
        return QPointF(1e9, 1e9);
    };
    const auto r = chem::reactionOf(out);
    REQUIRE(r);
    CHECK(r->reactants.size() == 2);
    CHECK(r->agents.size() == 1);
    CHECK(r->products.size() == 1);
    CHECK(std::abs(centre("c1ccncc1").x() - (arrow.from.x() + arrow.to.x()) / 2) < 1);  // agent centred on the arrow
    CHECK(centre("c1ccncc1").y() < arrow.from.y());                                     // and above it
    CHECK(std::abs(centre("CC(=O)Oc1ccccc1C(=O)O").y() - arrow.from.y()) < 1);           // product on the baseline
}

TEST_CASE("Arrange Scheme moves an orbital with the atom it is drawn on (#409)") {
    Fixture f;
    Document d = *chem::fromSmiles("CO");
    Arrow orbital{d.atoms[1].pos, d.atoms[1].pos + QPointF(0, -kBondLength)};
    orbital.kind = ArrowKind::POrbital;
    d.arrows.push_back(orbital);
    d.arrows.push_back({{60, 0}, {100, 0}});
    d.append(*chem::fromSmiles("CC"), {150, 30});
    f.canvas.setDocumentSilently(d);
    f.canvas.arrangeScheme();
    const Document& out = f.doc();
    REQUIRE(out.arrows.size() == 2);
    CHECK(QLineF(out.arrows[0].from, out.atoms[1].pos).length() < 1e-6);  // still on the oxygen
    CHECK(QLineF(out.arrows[0].to - out.arrows[0].from, QPointF(0, -kBondLength)).length() < 1e-6);
}

TEST_CASE("3D rotation from the mouse and keyboard keeps stereo (#173)") {
    Fixture f;
    const Document menthol = *chem::fromSmiles("CC(C)[C@@H]1CC[C@@H](C)C[C@H]1O");
    const std::string smiles = chem::toSmiles(menthol);
    f.canvas.setDocumentSilently(menthol);
    f.canvas.setTool(Canvas::Tool::Select);
    f.canvas.selectAll();
    const QPointF on = f.doc().atoms[3].pos;
    QTest::mousePress(f.canvas.viewport(), Qt::LeftButton, Qt::ShiftModifier | Qt::AltModifier, f.at(on));
    QMouseEvent move(QEvent::MouseMove, f.at(on + QPointF(30, 12)), f.canvas.viewport()->mapToGlobal(f.at(on + QPointF(30, 12))),
                     Qt::NoButton, Qt::LeftButton, Qt::ShiftModifier | Qt::AltModifier);
    QApplication::sendEvent(f.canvas.viewport(), &move);
    QTest::mouseRelease(f.canvas.viewport(), Qt::LeftButton, Qt::ShiftModifier | Qt::AltModifier, f.at(on + QPointF(30, 12)));
    CHECK_FALSE(f.doc() == menthol);  // it turned
    CHECK(chem::toSmiles(f.doc()) == smiles);
    QTest::keyClick(f.canvas.viewport(), Qt::Key_Up, Qt::ShiftModifier | Qt::AltModifier);
    CHECK(chem::toSmiles(f.doc()) == smiles);
    if (auto out = qgetenv("PENZENE_3D_SHOT"); !out.isEmpty()) {
        exportDocument(menthol, QString::fromUtf8(out) + ".before.png", {150, Qt::white});
        exportDocument(f.doc(), QString::fromUtf8(out), {150, Qt::white});
    }
}

TEST_CASE("the Select flyout offers 3D rotation without modifier keys") {
    App app;
    MainWindow w;
    w.resize(1000, 700);
    w.show();
    auto* canvas = w.findChild<Canvas*>();
    REQUIRE(canvas);
    QAction* rotate = nullptr;
    for (auto* action : w.findChildren<QAction*>())
        if (action->toolTip().startsWith("Rotate in 3D:")) rotate = action;
    REQUIRE(rotate);
    rotate->trigger();

    const Document menthol = *chem::fromSmiles("CC(C)[C@@H]1CC[C@@H](C)C[C@H]1O");
    canvas->setDocumentSilently(menthol);
    canvas->selectAll();
    const QPoint on = canvas->mapFromScene(menthol.atoms[3].pos);
    const QPoint to = on + QPoint(75, 30);
    QTest::mousePress(canvas->viewport(), Qt::LeftButton, {}, on);
    QMouseEvent move(QEvent::MouseMove, to, canvas->viewport()->mapToGlobal(to), Qt::NoButton, Qt::LeftButton, {});
    QApplication::sendEvent(canvas->viewport(), &move);
    QTest::mouseRelease(canvas->viewport(), Qt::LeftButton, {}, to);
    CHECK_FALSE(canvas->document() == menthol);
    CHECK(chem::toSmiles(canvas->document()) == chem::toSmiles(menthol));
}

TEST_CASE("a triple bond is gapped where a bond in front crosses it (#429)") {
    Document cross;
    cross.atoms = {{{-20, 0}}, {{20, 0}}, {{0, -20}}, {{0, 20}}};
    cross.bonds = {{0, 1}, {2, 3}};
    cross.bonds[0].order = 3;
    auto inked = [&](QPointF at) {
        const QImage img = renderImage(cross, {300, Qt::white});
        const QPoint px = ((at - outputBounds(cross).topLeft()) * exportScale(cross) * 300 / 72).toPoint();
        return img.pixelColor(px).lightness() < 160;
    };
    CHECK_FALSE(inked({1.6, 0}));  // beside the crossing, on the middle line
    CHECK(inked({12, 0}));         // but drawn away from it
}

TEST_CASE("attachment points, π-ligands, bring to front and atom properties (#60)") {
    Fixture f;
    Document d = *chem::fromSmiles("CC");
    edit::hotkey(d, {1, -1}, ".");
    CHECK(chem::toSmiles(d) == "*CC");
    Document fe = *chem::fromSmiles("[Fe]");
    edit::hotkey(fe, {0, -1}, "j");
    edit::hotkey(fe, {0, -1}, "J");
    INFO(chem::properties(fe)->formula << " " << chem::toSmiles(fe));
    CHECK(chem::properties(fe)->formula.find("C11H11") != std::string::npos);  // Cp⁻ (C5H5) + benzene (C6H6)

    // Two crossing bonds: the earlier one gets a gap where the later one crosses; f swaps which.
    Document cross;
    cross.atoms = {{{-20, 0}}, {{20, 0}}, {{0, -20}}, {{0, 20}}};
    cross.bonds = {{0, 1}, {2, 3}};
    auto inked = [](const Document& doc, QPointF at) {
        const QImage img = renderImage(doc, {300, Qt::white});
        const QRectF r = outputBounds(doc);
        const QPoint px = ((at - r.topLeft()) * exportScale(doc) * 300 / 72).toPoint();
        return img.pixelColor(px).lightness() < 160;
    };
    const QPointF onHorizontal(1.6, 0);  // just beside the crossing, on the horizontal bond
    CHECK_FALSE(inked(cross, onHorizontal));  // the horizontal bond (earlier) has the gap
    REQUIRE(edit::hotkey(cross, {-1, 0}, "f").bond == 1);
    CHECK(inked(cross, onHorizontal));  // now in front, drawn through

    f.canvas.setDocumentSilently(*chem::fromSmiles("CO"));
    QTimer::singleShot(0, [] {
        for (auto* w : QApplication::topLevelWidgets())
            if (auto* dialog = qobject_cast<QDialog*>(w); dialog && dialog->isVisible()) {
                dialog->findChildren<QSpinBox*>()[0]->setValue(-1);  // charge
                dialog->findChildren<QSpinBox*>()[2]->setValue(2);   // lone pairs
                dialog->accept();
            }
    });
    f.canvas.editAtomProperties(1);
    CHECK(f.doc().atoms[1].charge == -1);
    CHECK(f.doc().atoms[1].lonePairs == 2);
    CHECK(chem::toSmiles(f.doc()) == "C[O-]");
}

TEST_CASE("ChemDraw shortcut parity: y, W, g, ?, Space, Enter, nudging (#157)") {
    Fixture f;
    Document d = *chem::fromSmiles("CCC");
    REQUIRE(edit::hotkey(d, {2, -1}, "y").valid());
    CHECK(d.atoms[2].label == "Boc");
    REQUIRE(edit::hotkey(d, {-1, 0}, "W").valid());
    CHECK(d.bonds[0].stereo == BondStereo::Hash);

    f.canvas.setDocumentSilently(*chem::fromSmiles("CCO.CC"));
    f.canvas.setTool(Canvas::Tool::Select);
    QTest::mouseMove(f.canvas.viewport(), f.at(f.doc().atoms[1].pos));
    QTest::keyClick(f.canvas.viewport(), Qt::Key_G);
    CHECK(f.canvas.selection() == QSet<int>{1});  // g grabs the hotspot atom
    QTest::keyClick(f.canvas.viewport(), Qt::Key_Return);  // back to a hotspot
    CHECK(f.canvas.selection().isEmpty());
    QTest::keyClick(f.canvas.viewport(), Qt::Key_Space);  // hotspot → its molecule
    CHECK(f.canvas.selection() == (QSet<int>{0, 1, 2}));
    const QPointF before = f.doc().atoms[0].pos;
    QTest::keyClick(f.canvas.viewport(), Qt::Key_Right, Qt::ShiftModifier);  // nudge 10
    QTest::keyClick(f.canvas.viewport(), Qt::Key_Down);                      // nudge 1
    CHECK(QLineF(f.doc().atoms[0].pos, before + QPointF(10, 1)).length() < 1e-9);
    CHECK(QLineF(f.doc().atoms[3].pos, chem::fromSmiles("CCO.CC")->atoms[3].pos).length() < 1e-9);  // the other molecule stays
}

TEST_CASE("Structure → Join Fragments bonds two selected atoms from different fragments (#683)") {
    App app;
    MainWindow w;
    w.resize(900, 600);
    w.show();
    auto* canvas = w.findChild<Canvas*>();
    auto* stack = w.findChild<QUndoGroup*>()->activeStack();
    QAction* join = nullptr;
    for (auto* a : w.findChildren<QAction*>())
        if (a->shortcut() == QKeySequence(Qt::CTRL | Qt::Key_J)) join = a;
    REQUIRE(join);
    CHECK(join->shortcut() == QKeySequence(Qt::CTRL | Qt::Key_J));

    canvas->setDocumentSilently(*chem::fromSmiles("C.C"));
    REQUIRE(canvas->document().bonds.empty());
    canvas->setSelection({0, 1});
    join->trigger();
    REQUIRE(canvas->document().bonds.size() == 1);
    CHECK(canvas->document().bonds[0].order == 1);
    CHECK(canvas->document().bondBetween(0, 1) >= 0);
    CHECK(stack->undoText() == "Join");

    const int steps = stack->index();
    canvas->setSelection({0, 1});
    join->trigger();
    CHECK(canvas->document().bonds.size() == 1);  // already joined
    CHECK(stack->index() == steps);

    canvas->setDocumentSilently(*chem::fromSmiles("CCCC"));
    canvas->setSelection({0, 3});
    join->trigger();
    CHECK(canvas->document().bonds.size() == 3);  // one fragment: nothing to join
}

TEST_CASE("keyboard: G picks atoms one by one, > draws a curved arrow, no mouse needed (#541)") {
    Fixture f;
    f.canvas.setDocumentSilently(*chem::fromSmiles("CC=O"));
    f.canvas.setFocus();
    auto key = [&](int k, Qt::KeyboardModifiers m = {}) { QTest::keyClick(f.canvas.viewport(), Qt::Key(k), m); };
    auto type = [&](char c) { QTest::keyClick(f.canvas.viewport(), c); };
    f.canvas.setHotspot(0);
    type('G');
    CHECK(f.canvas.selection() == QSet<int>{0});
    key(Qt::Key_Right, Qt::ShiftModifier);  // the hotspot moves on to the next atom instead of nudging atom 0
    CHECK(f.canvas.hotspotAtom() == 1);
    type('G');
    CHECK(f.canvas.selection() == (QSet<int>{0, 1}));
    const QPointF before = f.doc().atoms[0].pos;
    key(Qt::Key_Escape);  // done picking: the selection stays and the arrows nudge it
    CHECK(f.canvas.selection() == (QSet<int>{0, 1}));
    key(Qt::Key_Down);
    CHECK(QLineF(f.doc().atoms[0].pos, before + QPointF(0, 1)).length() < 1e-9);

    // From the C=O bond to the oxygen, as an electron-pushing arrow.
    key(Qt::Key_Escape);
    f.canvas.setHotspot(-1, 1);
    type('>');
    CHECK(f.canvas.accessibleDescription().contains("curved arrow from"));
    const QPointF toO = f.doc().atoms[2].pos - f.doc().atoms[1].pos;
    key(std::abs(toO.x()) > std::abs(toO.y()) ? (toO.x() > 0 ? Qt::Key_Right : Qt::Key_Left) : (toO.y() > 0 ? Qt::Key_Down : Qt::Key_Up));
    CHECK(f.canvas.hotspotAtom() == 2);
    type('>');
    REQUIRE(f.doc().arrows.size() == 1);
    const Arrow& a = f.doc().arrows[0];
    CHECK(a.kind == ArrowKind::Reaction);
    CHECK(a.bend != 0);
    CHECK(a.fromAt == (std::array<int, 2>{1, 2}));
    CHECK(a.toAt == (std::array<int, 2>{2, -1}));

    // A hotspot the pointer set, with a selection, still nudges: picking is for the keyboard only.
    f.canvas.setSelection({0});
    QTest::mouseMove(f.canvas.viewport(), f.at(f.doc().atoms[2].pos));
    const QPointF c0 = f.doc().atoms[0].pos;
    key(Qt::Key_Right);
    CHECK(QLineF(f.doc().atoms[0].pos, c0 + QPointF(1, 0)).length() < 1e-9);
}

TEST_CASE("keyboard: bend and flip a curved arrow, and place an arrow or text after the selection (#549)") {
    Fixture f;
    f.canvas.setDocumentSilently(*chem::fromSmiles("CC=O"));
    f.canvas.setTool(Canvas::Tool::Select);  // as the window does when the canvas asks for it
    auto key = [&](int k, Qt::KeyboardModifiers m = {}) { QTest::keyClick(f.canvas.viewport(), Qt::Key(k), m); };
    auto toward = [](QPointF v) {
        return std::abs(v.x()) > std::abs(v.y()) ? (v.x() > 0 ? Qt::Key_Right : Qt::Key_Left) : (v.y() > 0 ? Qt::Key_Down : Qt::Key_Up);
    };
    f.canvas.setHotspot(-1, 1);
    QTest::keyClick(f.canvas.viewport(), '>');
    key(toward(f.doc().atoms[2].pos - f.doc().atoms[1].pos));
    QTest::keyClick(f.canvas.viewport(), '>');
    REQUIRE(f.doc().arrows.size() == 1);
    CHECK(f.canvas.selectedArrows() == QSet<int>{0});  // the new arrow, ready to reshape
    const Arrow drawn = f.doc().arrows[0];
    key(toward(f.doc().atoms[1].pos - f.doc().atoms[2].pos));  // the hotspot goes on, the arrow stays put
    CHECK(f.canvas.hotspotBond() == 1);
    CHECK(f.doc().arrows[0] == drawn);

    key(Qt::Key_Up, Qt::AltModifier);
    CHECK(std::abs(f.doc().arrows[0].bend / drawn.bend - 1.25) < 1e-9);
    key(Qt::Key_Down, Qt::AltModifier);
    CHECK(std::abs(f.doc().arrows[0].bend - drawn.bend) < 1e-9);
    f.canvas.setTool(Canvas::Tool::Bond);  // as selected by a screen reader's Press, with another tool on (#561)
    f.canvas.bendArrow(-1);  // Arrange → Flip Curved Arrow
    CHECK(std::abs(f.doc().arrows[0].bend + drawn.bend) < 1e-9);
    f.canvas.setTool(Canvas::Tool::Select);

    // After the molecule: a reaction arrow, then text at the new arrow's right.
    f.canvas.setSelection({0, 1, 2});
    double right = -1e9;
    for (const Atom& a : f.doc().atoms) right = std::max(right, a.pos.x());
    f.canvas.addArrowAfter();
    REQUIRE(f.doc().arrows.size() == 2);
    const Arrow& after = f.doc().arrows[1];
    CHECK(after.from.x() > right);
    CHECK(after.bend == 0);
    CHECK(std::abs(QLineF(after.from, after.to).length() - 3 * kBondLength) < 1e-9);
    CHECK(f.canvas.selectedArrows() == QSet<int>{1});
    QTimer::singleShot(0, [] {
        if (auto* d = qobject_cast<QInputDialog*>(QApplication::activeModalWidget())) d->setTextValue("heat"), d->accept();
    });
    f.canvas.addTextAfter();
    REQUIRE(f.doc().texts.size() == 1);
    CHECK(f.doc().texts[0].text == "heat");
    CHECK(f.doc().texts[0].pos.x() > f.doc().arrows[1].to.x());

    key(Qt::Key_Escape);
    f.canvas.setHotspot(2);
    REQUIRE(f.canvas.nextPlace());
    CHECK(len(*f.canvas.nextPlace() - f.doc().atoms[2].pos) >= kBondLength);  // beside the hotspot, not on top of it
    f.canvas.setHotspot(-1);
    CHECK_FALSE(f.canvas.nextPlace());  // nothing to place it by
}

TEST_CASE("update check: parsing GitHub's answer, comparing versions, off by default (#115)") {
    App app;
    const auto r = online::parseRelease(R"({"tag_name": "v0.9.0", "html_url": "https://github.com/JamesOBrien2/penzene/releases/tag/v0.9.0"})");
    CHECK(r.tag == "v0.9.0");
    CHECK(r.url.endsWith("/v0.9.0"));
    CHECK(online::isNewer("v0.9.0", "0.8.0"));
    CHECK(online::isNewer("v1.0.0", "0.10.2"));  // numerically, not as text
    CHECK_FALSE(online::isNewer("v0.8.0", "0.8.0"));
    CHECK_FALSE(online::isNewer("v0.7.1", "0.8.0"));
    CHECK_FALSE(online::isNewer("", "0.8.0"));
    CHECK(online::parseRelease("not json").tag.isEmpty());
    QSettings().remove("updates");
    MainWindow w;
    QTimer::singleShot(0, [] {
        for (auto* top : QApplication::topLevelWidgets())
            if (auto* dialog = qobject_cast<QDialog*>(top); dialog && dialog->isVisible()) {
                CHECK_FALSE(dialog->findChild<QCheckBox*>("autoUpdates")->isChecked());  // off unless chosen
                dialog->reject();
            }
    });
    w.showPreferences();
}

TEST_CASE("What's New: a release's notes, once after an update, never on a fresh install (#208)") {
    App app;
    const QString log = "# Changelog\n\n## Unreleased\n\n- soon\n\n## 0.9.0 (2026-10-01)\n\n- New thing\n- Other\n\n## 0.8.0 (2026-09-25)\n\n- Old\n";
    CHECK(online::releaseNotes(log, "0.9.0") == "- New thing\n- Other");
    CHECK(online::releaseNotes(log, "0.8.0") == "- Old");
    CHECK(online::releaseNotes(log, "0.9").isEmpty());  // not a prefix of 0.9.0
    CHECK(online::releaseNotes(log, "1.0.0").isEmpty());
    QFile bundled(":/CHANGELOG.md");
    REQUIRE(bundled.open(QIODevice::ReadOnly));  // shipped in the app
    CHECK(bundled.readAll().startsWith("# Changelog"));

    MainWindow w;
    int shown = 0;
    auto count = [&] {
        QTimer::singleShot(0, [&] {
            for (auto* top : QApplication::topLevelWidgets())
                if (auto* d = qobject_cast<QDialog*>(top); d && d->isVisible() && d->windowTitle().startsWith("What's New"))
                    ++shown, d->accept();
        });
    };
    QSettings().remove("lastVersion");
    count();
    w.maybeShowWhatsNew();  // fresh install: nothing
    QApplication::processEvents();
    CHECK(shown == 0);
    CHECK(QSettings().value("lastVersion").toString() == PENZENE_VERSION);
    w.maybeShowWhatsNew();  // same version again: nothing
    QApplication::processEvents();
    CHECK(shown == 0);
}

TEST_CASE("What's New separates highlighted additions from smaller changes (#215)") {
    const auto notes = online::parseReleaseNotes(
        "- **Drawing tools**: New shapes. <!-- icon: shapes -->\n"
        "- **Templates**: Ready-made structures.\n"
        "- Faster opening with **large** files. <!-- internal note -->");
    REQUIRE(notes.highlights.size() == 2);
    CHECK(notes.highlights[0].icon == "shapes");
    CHECK(notes.highlights[0].title == "Drawing tools");
    CHECK(notes.highlights[0].detail == "New shapes.");
    CHECK(notes.highlights[1].icon == "sparkles");
    CHECK(notes.highlights[1].title == "Templates");
    CHECK(notes.highlights[1].detail == "Ready-made structures.");
    CHECK(notes.others == QStringList{"Faster opening with large files."});
}

TEST_CASE("2.0 release notes are short and pick 4-10 highlights (#691)") {
    QFile f(":/CHANGELOG.md");
    REQUIRE(f.open(QIODevice::ReadOnly));
    const auto notes = online::parseReleaseNotes(online::releaseNotes(QString::fromUtf8(f.readAll()), "2.0.0"));
    CHECK(notes.highlights.size() >= 4);
    CHECK(notes.highlights.size() <= 10);
    auto words = [](const QString& s) { return int(s.simplified().split(' ', Qt::SkipEmptyParts).size()); };
    for (const auto& h : notes.highlights) {
        const int count = words(h.title) + words(h.detail);
        CHECK(count >= 2);
        CHECK(count <= 12);
    }
    for (const auto& other : notes.others) {
        CHECK(words(other) >= 2);
        CHECK(words(other) <= 12);
    }
}

TEST_CASE("What's New shows one card per highlighted addition (#215)") {
    App app;
    const QString changelog = "## 0.9.0 (2026-10-01)\n\n"
                              "- **Drawing tools**: Lines, boxes, rounded boxes, ellipses, arcs, brackets and more, all snapping to the grid. <!-- icon: shapes -->\n"
                              "- **Templates**: Ready-made structures.\n"
                              "- **Colours**: A new palette. <!-- icon: palette -->\n"
                              "- **Rotate in 3D**: Shift+Alt+drag turns a structure out of the page and keeps its stereo.\n"
                              "- **Projections**: Haworth, Fischer and Newman drawings with the right stereo.\n"
                              "- Faster opening.\n"
                              "- Lone pairs, radicals, and brackets with a subscript.\n"
                              "- Stretch and squash with handles; Structure → Transform.\n"
                              "- R-groups and generic atoms; attachment points and η-bonded rings.\n"
                              "- Arrow keys nudge a selection.\n";
    bool inspected = false;
    QTimer::singleShot(0, [&] {
        for (auto* top : QApplication::topLevelWidgets()) {
            auto* dialog = qobject_cast<QDialog*>(top);
            if (!dialog || dialog->objectName() != "whatsNew") continue;
            inspected = true;
            CHECK(dialog->findChildren<QFrame*>("highlight").size() == 5);
            for (auto* detail : dialog->findChildren<QLabel*>("detail"))  // changelog text isn't markup (#246)
                CHECK(detail->textFormat() == Qt::PlainText);
            for (auto* detail : dialog->findChildren<QLabel*>("detail"))  // wrapped text isn't clipped
                CHECK(detail->height() >= detail->heightForWidth(detail->width()));
            auto* others = dialog->findChild<QLabel*>("others");
            CHECK(others);
            if (others) {
                CHECK(others->textFormat() == Qt::PlainText);
                CHECK(others->text().startsWith("• Faster opening."));
                CHECK(dialog->findChildren<QLabel*>("others").size() == 5);
                auto* scroll = dialog->findChild<QScrollArea*>("moreScroll");
                CHECK(scroll);
                if (scroll) CHECK((scroll->widget() && scroll->widget()->isAncestorOf(others)));
            }
            auto* many = dialog->findChild<QLabel*>("manyMore");
            CHECK(many);
            if (many) CHECK(many->text().contains("many more"));
            dialog->accept();
        }
    });
    showWhatsNewDialog(nullptr, changelog, "0.9.0");
    CHECK(inspected);
}

TEST_CASE("What's New stays on screen and scrolls its smaller changes (#691)") {
    App app;
    QString changelog = "## 0.9.0 (2026-10-01)\n\n";
    for (int i = 0; i < 6; ++i)
        changelog += QString("- **Feature %1**: Short detail. <!-- icon: sparkles -->\n").arg(i);
    for (int i = 0; i < 40; ++i) changelog += QString("- Small change %1.\n").arg(i);
    bool inspected = false;
    QTimer::singleShot(0, [&] {
        for (auto* top : QApplication::topLevelWidgets()) {
            auto* dialog = qobject_cast<QDialog*>(top);
            if (!dialog || dialog->objectName() != "whatsNew") continue;
            inspected = true;
            if (const QScreen* screen = dialog->screen()) {
                const QRect available = screen->availableGeometry();
                CHECK(available.contains(dialog->geometry()));
            }
            auto* ok = dialog->findChild<QPushButton*>("continue");
            REQUIRE(ok);
            CHECK(ok->isVisible());
            CHECK(dialog->rect().contains(ok->geometry().center()));
            auto* scroll = dialog->findChild<QScrollArea*>("moreScroll");
            REQUIRE(scroll);
            CHECK(scroll->verticalScrollBar()->maximum() > 0);
            dialog->accept();
        }
    });
    showWhatsNewDialog(nullptr, changelog, "0.9.0");
    CHECK(inspected);
}

TEST_CASE("the colour swatch opens CPK colours; the pick paints atoms and bonds (#225)") {
    App app;
    MainWindow w;
    w.resize(1000, 700);
    w.show();
    auto* canvas = w.findChild<Canvas*>();
    QToolButton* colourButton = nullptr;
    for (auto* b : w.findChildren<QToolButton*>())
        if (b->defaultAction() && b->defaultAction()->toolTip().startsWith("Colour")) colourButton = b;
    REQUIRE(colourButton);
    CHECK(colourButton->popupMode() == QToolButton::InstantPopup);  // the swatch itself opens the colours
    REQUIRE(colourButton->menu());
    QToolButton* nitrogen = nullptr;
    for (auto* wa : colourButton->menu()->findChildren<QWidgetAction*>())
        for (auto* b : wa->defaultWidget()->findChildren<QToolButton*>())
            if (b->toolTip() == "Nitrogen") nitrogen = b;
    REQUIRE(nitrogen);
    nitrogen->click();
    const QColor blue(0x30, 0x50, 0xF8);
    CHECK(canvas->colour() == blue);
    CHECK(colourButton->defaultAction()->isChecked());  // picking chooses the tool

    Document d;
    d.atoms = {{{0, 0}}, {{kBondLength, 0}, 8}};
    d.bonds = {{0, 1}};
    canvas->setDocumentSilently(d);
    auto click = [&](QPointF scene) {
        QTest::mouseClick(canvas->viewport(), Qt::LeftButton, {}, canvas->mapFromScene(scene));
    };
    click(canvas->document().atoms[1].pos);
    click((canvas->document().atoms[0].pos + canvas->document().atoms[1].pos) / 2);
    CHECK(canvas->document().atoms[1].color == blue);
    CHECK(canvas->document().bonds[0].color == blue);
}

TEST_CASE("custom style editor saves a portable style and canvas undo restores its scale (#625)") {
    App app;
    QApplication::processEvents();
    MainWindow w;
    w.resize(1000, 700);
    w.show();
    auto* canvas = w.findChild<Canvas*>();
    canvas->setDocumentSilently(*chem::fromSmiles("[NH3]->[Pt+2](<-[NH3])(Cl)Cl"));
    canvas->fitToDocument();
    const QString name = "Penzene editor regression";
    deleteDrawingStyle(name);
    struct Cleanup { QString name; ~Cleanup() { deleteDrawingStyle(name); } } cleanup{name};
    QTimer::singleShot(0, [&] {
        auto* dialog = QApplication::activeModalWidget();
        REQUIRE(dialog);
        auto* length = dialog->findChild<QDoubleSpinBox*>("bondLength");
        auto* save = dialog->findChild<QPushButton*>("saveStyle");
        REQUIRE(length); REQUIRE(save);
        CHECK_FALSE(length->isEnabled());
        CHECK_FALSE(save->isEnabled());
        dialog->findChild<QPushButton*>("duplicateStyle")->click();
        CHECK(length->isEnabled());
        dialog->findChild<QLineEdit*>("styleName")->setText(name);
        length->setValue(20);
        dialog->findChild<QDoubleSpinBox*>("labelMargin")->setValue(2);
        save->click();
        CHECK(dialog->findChild<QPushButton*>("deleteStyle")->isEnabled());
        const QString preview = qEnvironmentVariable("PENZENE_STYLES_PREVIEW");
        if (!preview.isEmpty()) {
            REQUIRE(dialog->grab().save(preview + "/style-editor.png"));
            if (qEnvironmentVariableIsSet("PENZENE_FEEDBACK_PREVIEW_LIVE")) {
                QEventLoop loop;
                QObject::connect(qobject_cast<QDialog*>(dialog), &QDialog::finished, &loop, &QEventLoop::quit);
                loop.exec();
                return;
            }
        }
        qobject_cast<QDialog*>(dialog)->accept();
    });
    editDrawingStyles(&w, canvas->document());
    if (qEnvironmentVariableIsSet("PENZENE_FEEDBACK_PREVIEW_LIVE")) return;
    QTimer::singleShot(0, [&] {
        auto* preferences = qobject_cast<QDialog*>(QApplication::activeModalWidget());
        REQUIRE(preferences);
        for (auto* box : preferences->findChildren<QComboBox*>())
            if (box->findText("Dark") >= 0) box->setCurrentText("Dark");
        preferences->accept();
    });
    w.showPreferences();
    const Document before = canvas->document();
    const double scale = canvas->transform().m11();
    Document after = before;
    REQUIRE(applyDrawingStyle(after, name));
    canvas->commit(after, "Drawing style");
    CHECK(canvas->transform().m11() / scale == Catch::Approx(20 / 14.4));
    CHECK(canvas->selectedSubset().styleSettings == after.styleSettings);
    w.findChild<QUndoGroup*>()->undo();
    CHECK(canvas->document() == before);
    CHECK(canvas->transform().m11() == Catch::Approx(scale));
    w.findChild<QUndoGroup*>()->redo();
    CHECK(canvas->document().styleSettings == after.styleSettings);
}

TEST_CASE("Circled Charges is an undoable document setting (#623)") {
    App app;
    QApplication::processEvents();
    MainWindow w;
    w.resize(1100, 750);
    w.show();
    auto* canvas = w.findChild<Canvas*>();
    Document d = *chem::fromSmiles("[NH3]->[Pt+2](<-[NH3])(Cl)Cl");
    d.append(*chem::fromSmiles("C[CH2+]"), {65, 0});
    d.append(*chem::fromSmiles("[Cl-]"), {100, 0});
    canvas->setDocumentSilently(d);
    canvas->fitToDocument();
    QAction* charges = nullptr;
    for (auto* a : w.findChildren<QAction*>()) if (a->text() == "Circled Charges") charges = a;
    REQUIRE(charges);
    CHECK_FALSE(charges->isChecked());
    const QString preview = qEnvironmentVariable("PENZENE_CHARGES_PREVIEW");
    QApplication::processEvents();
    if (!preview.isEmpty()) REQUIRE(w.grab().save(preview + "/charges-off.png"));
    charges->trigger();
    CHECK(canvas->document().circledCharges);
    QApplication::processEvents();
    if (!preview.isEmpty()) {
        REQUIRE(w.grab().save(preview + "/charges-on.png"));
        REQUIRE(writeWhole(preview + "/charges.penz", canvas->document().toJson()));
    }
    w.findChild<QUndoGroup*>()->undo();
    CHECK_FALSE(charges->isChecked());
    CHECK(canvas->document() == d);
    w.findChild<QUndoGroup*>()->redo();
    CHECK(charges->isChecked());
    CHECK(canvas->selectedSubset().circledCharges);
}

TEST_CASE("changing charge magnitude removes and restores its circle (#623)") {
    Fixture f;
    for (const auto& smiles : {"[Na+]", "[Cl-]", "C[CH2+]"}) {
        Document d = *chem::fromSmiles(smiles);
        d.circledCharges = true;
        const int at = int(d.atoms.size()) - 1;
        f.canvas.setDocumentSilently(d);
        f.canvas.setHotspot(at);
        const int single = d.atoms[at].charge;
        const QImage before = renderImage(d);
        f.key(single > 0 ? "+" : "-");
        CHECK(f.doc().atoms[at].charge == single * 2);
        CHECK(f.doc().circledCharges);
        CHECK_FALSE(renderSvg(f.doc()).contains("<circle"));
        Document plain = f.doc();
        plain.circledCharges = false;
        CHECK(renderImage(f.doc()) == renderImage(plain));
        f.undo.undo();
        CHECK(renderImage(f.doc()) == before);
        CHECK(renderSvg(f.doc()).contains("<circle"));
        f.undo.redo();
        CHECK_FALSE(renderSvg(f.doc()).contains("<circle"));
    }
}

TEST_CASE("journal presets appear in the menu and change physical canvas size (#473)") {
    App app;
    QApplication::processEvents();
    MainWindow w;
    w.resize(1100, 750);
    w.show();
    auto* canvas = w.findChild<Canvas*>();
    Document d = *chem::fromSmiles("Oc1ccccc1");
    d.append(*chem::fromSmiles("COc1ccccc1"), {95, 0});
    d.arrows.push_back({{40, 0}, {60, 0}});
    d.texts.push_back({{37, -6}, "MeI, base", 0.7});
    d.append(*chem::fromSmiles("[NH3]->[Pt+2](<-[NH3])(Cl)Cl"), {40, 55});
    canvas->setDocumentSilently(d);
    canvas->fitToDocument();
    canvas->zoomBy(0.8);  // Leave room for the larger Angewandte preset.
    const double zoom = canvas->transform().m11();
    const QString preview = qEnvironmentVariable("PENZENE_JOURNALS_PREVIEW");
    for (const QString& name : {QString("Nature"), QString("Angewandte"), QString("SYNLETT/SYNTHESIS")}) {
        QAction* action = nullptr;
        for (auto* a : w.findChildren<QAction*>()) if (a->text() == name) action = a;
        REQUIRE(action);
        action->trigger();
        CHECK(canvas->document().style == name);
        CHECK(canvas->transform().m11() / zoom == Catch::Approx(drawingStyle(name).bondLength / kBondLength));
        CHECK(action->isChecked());
        QApplication::processEvents();
        if (!preview.isEmpty()) {
            const QString file = name == "SYNLETT/SYNTHESIS" ? "synlett-synthesis" : name.toLower();
            REQUIRE(w.grab().save(preview + "/" + file + ".png"));
            REQUIRE(writeWhole(preview + "/" + file + ".penz", canvas->document().toJson()));
        }
    }
}

TEST_CASE("Text has its own rail button and keeps its hotkey (#622)") {
    App app;
    // Cocoa treats Catch's test filter as a file to open; consume launch events before the window listens.
    QApplication::processEvents();
    MainWindow w;
    w.resize(1000, 700);
    w.show();
    auto* canvas = w.findChild<Canvas*>();
    QHash<QString, QToolButton*> rail;
    for (auto* b : w.findChildren<QToolButton*>("railButton")) rail[b->text()] = b;
    REQUIRE(rail.contains("Text"));
    auto* text = rail["Text"];
    REQUIRE(text->isVisible());
    CHECK(text->accessibleName() == "Text");
    CHECK(text->toolTip().endsWith("— t"));
    CHECK(text->focusPolicy() == Qt::StrongFocus);
    CHECK(text->menu() == nullptr);
    for (auto* f : w.findChildren<QFrame*>("toolFlyout"))
        for (auto* b : f->findChildren<QToolButton*>()) CHECK(b->defaultAction() != text->defaultAction());
    const QString preview = qEnvironmentVariable("PENZENE_TEXT_PREVIEW_DIR");
    canvas->setDocumentSilently(*chem::fromSmiles("CCO"));
    canvas->fitToDocument();
    canvas->zoomBy(0.5);
    QApplication::processEvents();
    if (!preview.isEmpty()) REQUIRE(w.grab().save(preview + "/resting.png"));
    auto drawText = [&](QPointF pos, const QString& value) {
        QTimer::singleShot(0, [value] {
            if (auto* d = qobject_cast<QInputDialog*>(QApplication::activeModalWidget()))
                d->setTextValue(value), d->accept();
        });
        QTest::mouseClick(canvas->viewport(), Qt::LeftButton, {}, canvas->mapFromScene(pos));
    };
    text->click();
    CHECK(text->isChecked());
    for (auto* b : rail) if (b != text) CHECK_FALSE(b->isChecked());
    for (auto* f : w.findChildren<QFrame*>("toolFlyout")) CHECK_FALSE(f->isVisible());
    drawText({50, 30}, "H2O");
    REQUIRE(canvas->document().texts.size() == 1);
    CHECK(canvas->document().texts[0].text == "H2O");
    QApplication::processEvents();
    if (!preview.isEmpty()) {
        REQUIRE(w.grab().save(preview + "/text-added.png"));
        REQUIRE(writeWhole(preview + "/sample.penz", canvas->document().toJson()));
        if (qEnvironmentVariableIsSet("PENZENE_TEXT_PREVIEW_LIVE")) {
            w.findChild<QUndoGroup*>()->activeStack()->setClean();
            QEventLoop loop;
            QObject::connect(qApp, &QApplication::lastWindowClosed, &loop, &QEventLoop::quit);
            loop.exec();
        }
    }
    rail["Shapes"]->click();
    CHECK(rail["Shapes"]->isChecked());
    CHECK_FALSE(text->isChecked());
    for (auto* f : w.findChildren<QFrame*>("toolFlyout")) f->hide();
    canvas->setHotspot(-1);
    QTest::keyClick(canvas->viewport(), Qt::Key_T);
    CHECK(text->isChecked());
    CHECK_FALSE(rail["Shapes"]->isChecked());
    drawText({50, 55}, "heat");
    REQUIRE(canvas->document().texts.size() == 2);
    CHECK(canvas->document().texts[1].text == "heat");
    rail["Bonds"]->click();
    CHECK(rail["Bonds"]->isChecked());
    CHECK_FALSE(text->isChecked());
}

TEST_CASE("White and teal theme; tools on a rail whose groups open beside it (#273, #220)") {
    CHECK(theme("Light").paper == QColor("#FFFFFF"));
    CHECK(theme("Light").window == QColor("#F3F5F6"));
    CHECK(theme("Light").accent == QColor("#0F6E56"));
    CHECK(theme("Dark").paper == QColor("#15171A"));
    CHECK(theme("Dark").accent == QColor("#4CC9A0"));
    App app;
    MainWindow w;
    w.resize(1000, 700);
    w.show();
    CHECK(w.findChild<QFrame*>("toolCard"));
    CHECK_FALSE(w.findChild<QToolBar*>("modeBar"));  // no Draw / Chemistry / Figure switch
    QHash<QString, QToolButton*> rail;
    for (auto* b : w.findChildren<QToolButton*>("railButton")) rail[b->text()] = b;
    CHECK(rail.keys().size() == 8);
    REQUIRE(rail.contains("Text"));
    REQUIRE(rail.contains("Templates"));
    for (const char* g : {"Select", "Bonds", "Rings", "Atoms", "Arrows", "Shapes"}) REQUIRE(rail.contains(g));
    CHECK(rail["Bonds"]->isChecked());  // the single bond, chosen at start
    auto flyout = [&](const QString& g) {
        for (auto* f : w.findChildren<QFrame*>("toolFlyout"))
            if (f->findChild<QLabel*>("flyoutTitle")->text() == g.toUpper()) return f;
        return static_cast<QFrame*>(nullptr);
    };
    auto* rings = flyout("Rings");
    REQUIRE(rings);
    CHECK_FALSE(rings->isVisible());
    rail["Rings"]->click();
    CHECK(rings->isVisible());
    CHECK(rail["Rings"]->isChecked());
    CHECK(rings->findChildren<QToolButton*>().size() >= 8);
    // Clicking a rail button again leaves its flyout where it opened, whichever group it is.
    for (const char* g : {"Select", "Bonds", "Rings", "Atoms", "Arrows", "Shapes"}) {
        rail[g]->click();
        const QPoint opened = flyout(g)->pos();
        rail[g]->click();
        CHECK(flyout(g)->pos() == opened);
        CHECK(flyout(g)->geometry().bottom() <= w.height());
    }
    // Flyouts stay open until closed: picking a tool, opening another group or drawing doesn't hide them.
    auto pick = [](QFrame* f, int i) {
        auto tools = f->findChildren<QToolButton*>();
        tools.removeIf([](QToolButton* b) { return b->objectName() == "close"; });
        tools[i]->click();
    };
    rail["Arrows"]->click();
    pick(rings, 2);
    pick(flyout("Arrows"), 1);
    CHECK(rings->isVisible());
    CHECK(flyout("Arrows")->isVisible());
    // The X closes just its own flyout.
    flyout("Arrows")->findChild<QToolButton*>("close")->click();
    CHECK_FALSE(flyout("Arrows")->isVisible());
    CHECK(rings->isVisible());
    // A flyout can be dragged by its title, and stays there (#466).
    auto* title = rings->findChild<QLabel*>("flyoutTitle");
    const QPoint before = rings->pos(), grab = title->rect().center();
    QTest::mousePress(title, Qt::LeftButton, {}, grab);
    QTest::mouseMove(title, grab + QPoint(60, 40));
    QTest::mouseRelease(title, Qt::LeftButton, {}, grab + QPoint(60, 40));
    CHECK(rings->pos() == before + QPoint(60, 40));
    rail["Rings"]->click();
    CHECK(rings->pos() == before + QPoint(60, 40));
    // ... and by any bare part of it, but not from a tool.
    const QPoint bare(rings->width() / 2, 2);  // the top margin; the corners resize
    QTest::mousePress(rings, Qt::LeftButton, {}, bare);
    QTest::mouseMove(rings, bare + QPoint(-20, 10));
    QTest::mouseRelease(rings, Qt::LeftButton, {}, bare + QPoint(-20, 10));
    CHECK(rings->pos() == before + QPoint(40, 50));
    auto* tool = rings->findChildren<QToolButton*>(Qt::FindDirectChildrenOnly).back();
    QTest::mousePress(tool, Qt::LeftButton);
    QTest::mouseMove(tool, tool->rect().center() + QPoint(30, 30));
    QTest::mouseRelease(tool, Qt::LeftButton, {}, tool->rect().center() + QPoint(500, 500));  // off it: not a click
    CHECK(rings->pos() == before + QPoint(40, 50));
    // A corner resizes it: the tools reflow into as many columns as fit, and the box follows them.
    auto columns = [](QFrame* f) {
        QSet<int> xs;
        for (auto* b : f->findChildren<QToolButton*>(Qt::FindDirectChildrenOnly))  // not the menus' swatches
            if (b->objectName() != "close") xs.insert(b->x());
        return int(xs.size());
    };
    auto holdsTools = [](QFrame* f) {
        for (auto* b : f->findChildren<QToolButton*>(Qt::FindDirectChildrenOnly))
            if (!f->rect().contains(b->geometry())) return false;
        return true;
    };
    auto dragCorner = [&](QPoint from, QPoint by, QFrame* f = nullptr) {
        QFrame* target = f ? f : rings;
        QTest::mousePress(target, Qt::LeftButton, {}, from);
        QTest::mouseMove(target, from + by);
        QTest::mouseRelease(target, Qt::LeftButton, {}, from + by);
    };
    const int usual = columns(rings), tall = rings->height();
    dragCorner({rings->width() - 3, rings->height() - 3}, {200, 0});
    CHECK(columns(rings) > usual);
    CHECK(rings->height() < tall);
    CHECK(holdsTools(rings));
    const QPoint topLeft = rings->pos();
    dragCorner({rings->width() - 3, rings->height() - 3}, {40 - rings->width(), 0});
    CHECK(columns(rings) < usual);
    CHECK(rings->height() > tall);
    CHECK(rings->pos() == topLeft);  // the far corner is the one that moves
    CHECK(holdsTools(rings));
    const int right = rings->geometry().right();
    dragCorner({3, 3}, {-60, 0});  // from the top left, the right edge stays
    CHECK(rings->geometry().right() == right);
    CHECK(holdsTools(rings));
    const int fewer = columns(rings);
    dragCorner({rings->width() - 3, rings->height() - 3}, {0, -rings->height() / 2});  // a shorter box needs more columns
    CHECK(columns(rings) > fewer);
    CHECK(holdsTools(rings));
    // Shapes lays its orbitals in rows at its usual width, but resized it flows like the rest.
    rail["Shapes"]->click();
    auto* shapes = flyout("Shapes");
    dragCorner({shapes->width() - 3, shapes->height() - 3}, {500, 0}, shapes);
    QSet<int> ys;
    for (auto* b : shapes->findChildren<QToolButton*>(Qt::FindDirectChildrenOnly))
        if (b->objectName() != "close") ys.insert(b->y());
    CHECK(ys.size() <= 2);  // not held to a row for each section
    CHECK(holdsTools(shapes));
    // Shrinking the window keeps an open flyout inside it.
    rail["Select"]->click();
    auto* select = flyout("Select");
    select->move(w.width() - select->width() - 10, w.height() - select->height() - 10);
    w.resize(600, 400);
    CHECK(w.rect().contains(select->geometry()));
    w.resize(1000, 700);
    rings->findChild<QToolButton*>("close")->click();
    CHECK_FALSE(rings->isVisible());
    // Properties and Templates are in the View menu, not the palette.
    CHECK(w.findChild<QDockWidget*>("properties"));
    CHECK(w.findChild<QDockWidget*>("templates"));
    CHECK(w.findChildren<QFrame*>("panelCard").size() == 5);  // Properties, Templates, Mass Spec, NMR
    for (auto* f : w.findChildren<QFrame*>("toolFlyout"))
        for (auto* b : f->findChildren<QToolButton*>()) CHECK(b->toolButtonStyle() != Qt::ToolButtonTextOnly);
}

TEST_CASE("PDF: copied and exported as vectors, with the drawing attached (#228)") {
    App app;
    Document doc = *chem::fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
    doc.atoms[0].color = Qt::red;
    doc.texts.push_back({{0, 60}, "aspirin"});
    const QByteArray pdf = renderPdf(doc);
    CHECK(pdf.startsWith("%PDF"));
    CHECK(pdf.contains("/EmbeddedFiles"));
    auto back = Document::fromEmbedded(pdf);
    REQUIRE(back);
    CHECK(*back == doc);
    QTemporaryDir dir;
    REQUIRE(exportDocument(doc, dir.filePath("aspirin.pdf")));
    back = chem::readFile(dir.filePath("aspirin.pdf"));
    REQUIRE(back);
    CHECK(*back == doc);
    CHECK_FALSE(Document::fromEmbedded(QByteArray("%PDF-1.4\n1 0 obj\n<<>>\nstream\nnot a drawing\nendstream\n")));
    // Only attached files are read (#246): a drawing in a plain stream isn't, so a big PDF from
    // elsewhere isn't inflated stream by stream.
    const QByteArray json = Document(*chem::fromSmiles("CCO")).toJson();
    const QByteArray plain = "%PDF-1.4\n1 0 obj\n<<>>\nstream\n" + json + "\nendstream\nendobj\n";
    CHECK_FALSE(Document::fromEmbedded(plain));
    CHECK(Document::fromEmbedded(plain + "2 0 obj\n<< /Type/Filespec /EF <</F 1 0 R>> >>\nendobj\n"));

    // Copy offers the PDF; a PDF alone pastes back as the drawing.
    MainWindow w;
    auto* canvas = w.findChild<Canvas*>();
    canvas->setDocumentSilently(doc);
    canvas->selectAll();
    for (auto* a : w.findChildren<QAction*>())
        if (a->shortcut() == QKeySequence::Copy) a->trigger();
    const QByteArray copied = QApplication::clipboard()->mimeData()->data("application/pdf");
    CHECK(copied.startsWith("%PDF"));
    auto* mime = new QMimeData;
    mime->setData("application/pdf", copied);
    QApplication::clipboard()->setMimeData(mime);
    canvas->setDocumentSilently(Document{});
    for (auto* a : w.findChildren<QAction*>())
        if (a->shortcut() == QKeySequence::Paste) a->trigger();
    CHECK(canvas->document().atoms.size() == doc.atoms.size());
    CHECK(canvas->document().texts.size() == 1);
}

TEST_CASE("welcome card: examples on an empty page, gone once drawing starts (#116)") {
    App app;
    MainWindow w;
    w.resize(1100, 750);
    w.show();
    auto* canvas = w.findChild<Canvas*>();
    auto* welcome = w.findChild<QFrame*>("welcome");
    REQUIRE(welcome);
    CHECK(welcome->isVisible());
    QList<QToolButton*> examples;
    for (auto* b : welcome->findChildren<QToolButton*>("example")) examples << b;
    REQUIRE(examples.size() == 3);  // aspirin, a reaction scheme, a mechanism
    examples[2]->click();
    CHECK_FALSE(welcome->isVisible());
    CHECK_FALSE(canvas->document().arrows.empty());  // the mechanism's curved arrows
    for (auto* a : w.findChildren<QAction*>())
        if (a->shortcut() == QKeySequence::New) a->trigger();
    CHECK(canvas->document().empty());
    CHECK(welcome->isVisible());
    // A click on the page, away from the card, dismisses it and still draws.
    canvas->setTool(Canvas::Tool::Bond);
    QTest::mouseClick(canvas->viewport(), Qt::LeftButton, {}, QPoint(40, 40));
    CHECK_FALSE(welcome->isVisible());
    CHECK_FALSE(canvas->document().empty());
}

TEST_CASE("bold bonds meet without a notch at a skeletal atom (#217)") {
    App app;
    // Two bold bonds in a V from the carbon at the origin; bold is 2 pt wide, so each half is 1 pt.
    const double h = kBondLength * std::sqrt(3.0) / 2, v = kBondLength / 2;
    for (bool boldBoth : {true, false}) {
        Document doc;
        doc.addAtom({0, 0});
        doc.addAtom({-h, v});
        doc.addAtom({h, v});
        doc.bonds.push_back({0, 1});
        doc.bonds.push_back({0, 2});
        doc.bonds[0].stereo = BondStereo::Bold;
        if (boldBoth) doc.bonds[1].stereo = BondStereo::Bold;
        QImage img(400, 400, QImage::Format_ARGB32_Premultiplied);
        img.fill(Qt::transparent);
        QPainter p(&img);
        p.translate(200, 200);
        p.scale(40, 40);
        paintDocument(p, doc);
        p.end();
        auto inked = [&](QPointF at) { return qAlpha(img.pixel((at * 40 + QPointF(200, 200)).toPoint())) > 128; };
        INFO((boldBoth ? "two bold bonds" : "bold and plain"));
        // Just past the bold bond's flat end, on the outside of the corner.
        CHECK(inked({-0.25, -0.55}));
        if (boldBoth) CHECK(inked({0, -0.7}));
        CHECK_FALSE(inked({0, -1.3}));  // and nothing sticks out beyond the bonds' width
    }
}

TEST_CASE("the logo is the teal ring on a white tile (#273)") {
    App app;
    QFile f(":/logo.svg");
    REQUIRE(f.open(QIODevice::ReadOnly));
    const QString svg = QString::fromUtf8(f.readAll()).toLower();
    CHECK(svg.contains(theme("Light").accent.name()));  // teal ink
    CHECK(svg.contains(theme("Light").paper.name()));   // on a white tile
}

TEST_CASE("the window and What's New share one palette (#246, #273)") {
    const Chrome light = chrome(theme("Light")), dark = chrome(theme("Dark"));
    CHECK(light.border == QColor("#DDE1E4"));
    CHECK(light.secondary == QColor("#5E6770"));
    CHECK(light.accentBg == QColor("#E1F5EE"));
    CHECK(dark.border == QColor("#353A40"));
    CHECK(dark.accentBg == QColor("#0B3B30"));
    CHECK(chrome(theme("Catppuccin Mocha")).secondary == theme("Catppuccin Mocha").text);  // keeps its own
}

TEST_CASE("template thumbnails are drawn in the theme's ink (#261)") {
    App app;
    auto lightest = [](const QIcon& icon) {
        const QImage img = icon.pixmap(56, 40).toImage();
        int most = 0;
        for (int y = 0; y < img.height(); ++y)
            for (int x = 0; x < img.width(); ++x)
                if (qAlpha(img.pixel(x, y)) > 128) most = std::max(most, QColor(img.pixel(x, y)).lightness());
        return most;
    };
    for (QString t : {"Dark", "Light"}) {
        QSettings().setValue("theme", t);
        MainWindow w;
        w.show();
        auto* dock = w.findChild<QDockWidget*>("templates");
        REQUIRE(dock);
        dock->show();
        auto* tree = dock->findChild<QTreeWidget*>();
        REQUIRE(tree->topLevelItemCount() > 0);
        const int ink = lightest(tree->topLevelItem(0)->child(0)->icon(0));
        if (t == "Dark") CHECK(ink > 128);  // light strokes on the dark panel (black before)
        else CHECK(ink < 128);
    }
    QSettings().remove("theme");
}

TEST_CASE("drawing and chemistry time grow linearly with the drawing (#114)") {
    App app;
    Document one = *chem::fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
    auto page = [&](int copies) {
        Document doc;
        for (int i = 0; i < copies; ++i) doc.append(one, QPointF((i % 20) * 150, (i / 20) * 150));
        return doc;
    };
    // The best of a few runs, so a busy machine doesn't decide it.
    auto best = [](const Document& doc) {
        qint64 ns = std::numeric_limits<qint64>::max();
        for (int k = 0; k < 3; ++k) {
            QElapsedTimer t;
            t.start();
            QPicture pic;
            QPainter p(&pic);
            paintDocument(p, doc);
            p.end();
            chem::properties(doc);
            ns = std::min(ns, t.nsecsElapsed());
        }
        return double(ns);
    };
    best(page(10));  // warm up fonts and RDKit
    const double small = best(page(160)), large = best(page(640));  // 2080 and 8320 atoms
    INFO("4x the atoms took " << large / small << "x as long");
    CHECK(large / small < 8);  // linear is ~4, quadratic ~16
}

TEST_CASE("translations: a built-in language is offered, and applies to the menus (#113)") {
    App app;
    CHECK(MainWindow::languages().contains("xx"));  // the test build's made-up language
    REQUIRE(MainWindow::installTranslations("xx"));
    {
        MainWindow w;
        CHECK(w.menuBar()->actions().value(0)->text() == "&Fichier-xx");
    }
    REQUIRE(MainWindow::installTranslations("en"));  // back to the source strings
    MainWindow w;
    CHECK(w.menuBar()->actions().value(0)->text() == "&File");
    // Preferences lists it under Language.
    bool offered = false;
    QTimer::singleShot(0, [&] {
        if (auto* dialog = qobject_cast<QDialog*>(QApplication::activeModalWidget())) {
            if (auto* box = dialog->findChild<QComboBox*>("language")) offered = box->findData("xx") >= 0;
            dialog->reject();
        }
    });
    w.showPreferences();
    CHECK(offered);
}

static double contrast(QColor a, QColor b) {  // WCAG 2 contrast ratio
    auto lum = [](QColor c) {
        auto ch = [](double v) { return v <= 0.03928 ? v / 12.92 : std::pow((v + 0.055) / 1.055, 2.4); };
        return 0.2126 * ch(c.redF()) + 0.7152 * ch(c.greenF()) + 0.0722 * ch(c.blueF());
    };
    const double x = lum(a), y = lum(b);
    return (std::max(x, y) + 0.05) / (std::min(x, y) + 0.05);
}

TEST_CASE("every tool button does its job pressed through accessibility, as screen readers press it (#535)") {
    App app;
    MainWindow w;
    w.resize(1100, 750);
    w.show();
    auto press = [](QToolButton* b) {  // the button's own action, as VoiceOver offers it
        QAccessibleActionInterface* act = QAccessible::queryAccessibleInterface(b)->actionInterface();
        REQUIRE(act);
        const QStringList names = act->actionNames();
        REQUIRE((names.contains(QAccessibleActionInterface::toggleAction()) || names.contains(QAccessibleActionInterface::pressAction())));
        act->doAction(names.contains(QAccessibleActionInterface::toggleAction()) ? QAccessibleActionInterface::toggleAction()
                                                                                 : QAccessibleActionInterface::pressAction());
        QCoreApplication::processEvents();
    };
    for (auto* rail : w.findChildren<QToolButton*>("railButton")) {
        INFO(rail->text().toStdString());
        QFrame* fly = nullptr;
        for (auto* f : w.findChildren<QFrame*>("toolFlyout"))
            if (f->findChild<QLabel*>("flyoutTitle")->text() == rail->text().toUpper()) fly = f;
        if (rail->text() == "Text") {
            CHECK_FALSE(fly);
            press(rail);
            CHECK(rail->isChecked());
            CHECK(w.statusBar()->currentMessage() == rail->toolTip());
            continue;
        }
        if (rail->text() == "Templates") {
            CHECK_FALSE(fly);
            auto* dock = w.findChild<QDockWidget*>("templates");
            REQUIRE(dock);
            press(rail);
            CHECK(dock->isVisible());
            CHECK(rail->isChecked());
            press(rail);
            CHECK_FALSE(dock->isVisible());
            CHECK_FALSE(rail->isChecked());
            continue;
        }
        REQUIRE(fly);
        fly->hide();
        press(rail);
        CHECK(fly->isVisible());  // its tools open
        for (auto* b : fly->findChildren<QToolButton*>()) {
            if (!b->defaultAction() || !b->defaultAction()->isCheckable()) continue;  // the close button
            if (b->popupMode() == QToolButton::InstantPopup) continue;  // a press opens its menu: the swatches, below
            INFO(b->toolTip().toStdString());
            const bool current = b->defaultAction()->isChecked();  // its group's tool, picked by the rail
            w.statusBar()->clearMessage();
            press(b);
            CHECK(b->defaultAction()->isChecked());
            if (!current) CHECK(w.statusBar()->currentMessage() == b->toolTip());  // the tool was picked, not just the button lit
        }
    }
    // A colour swatch picks its colour.
    QToolButton* oxygen = nullptr;
    for (auto* wa : w.findChildren<QWidgetAction*>())
        for (auto* b : wa->defaultWidget()->findChildren<QToolButton*>())
            if (b->isCheckable() && b->accessibleName() == "Oxygen") oxygen = b;
    REQUIRE(oxygen);
    press(oxygen);
    CHECK(w.findChild<Canvas*>()->colour() == QColor(0xFF, 0x0D, 0x0D));
}

TEST_CASE("accessibility: named, focusable tools; arrow keys in the periodic table; contrast (#112)") {
    // Every theme: text 4.5:1 or better, marks drawn on the page (hotspot, selection, errors) 3:1.
    for (const auto& t : themes()) {
        INFO(t.name.toStdString());
        CHECK(contrast(t.ink, t.paper) >= 4.5);
        if (t.text.isValid()) CHECK(contrast(t.text, t.window) >= 4.5);
        CHECK(contrast(t.accent, t.paper) >= 3);
        CHECK(contrast(t.hotspot, t.paper) >= 3);
        CHECK(contrast(t.error, t.paper) >= 3);
    }

    App app;
    MainWindow w;
    w.resize(1100, 750);
    w.show();
    QToolButton* benzene = nullptr;
    QList<QToolButton*> tools = w.findChild<QFrame*>("toolCard")->findChildren<QToolButton*>();
    for (auto* f : w.findChildren<QFrame*>("toolFlyout")) tools += f->findChildren<QToolButton*>();
    for (auto* b : tools) {
        INFO(b->toolTip().toStdString());
        const QString name = QAccessible::queryAccessibleInterface(b)->text(QAccessible::Name);
        CHECK_FALSE(name.isEmpty());  // what a screen reader says
        if (!qobject_cast<QMenu*>(b->window())) CHECK(b->focusPolicy() == Qt::StrongFocus);  // Tab reaches it
        if (name == "Benzene") benzene = b;
    }
    CHECK(benzene);

    // The periodic table: arrow keys move between elements, over the table's gaps.
    QHash<QString, QToolButton*> elements;
    for (auto* wa : w.findChildren<QWidgetAction*>())
        for (auto* b : wa->defaultWidget()->findChildren<QToolButton*>()) elements[b->text()] = b;
    REQUIRE(elements.contains("C"));
    CHECK(elements["C"]->accessibleName() == "Carbon");
    QWidget* table = elements["C"]->parentWidget();
    table->show();
    auto press = [&](QToolButton* from, Qt::Key key) {
        from->setFocus();
        QKeyEvent e(QEvent::KeyPress, key, Qt::NoModifier);
        QApplication::sendEvent(from, &e);
        return table->window()->focusWidget();
    };
    CHECK(press(elements["C"], Qt::Key_Right) == elements["N"]);
    CHECK(press(elements["C"], Qt::Key_Down) == elements["Si"]);
    CHECK(press(elements["Be"], Qt::Key_Right) == elements["B"]);  // across the d-block gap
    CHECK(press(elements["He"], Qt::Key_Up) == elements["He"]);    // nothing above: stays

}

TEST_CASE("a file opened before the window shows is fitted to the window once it does (#545)") {
    App app;
    MainWindow w;
    w.resize(1100, 750);
    auto* canvas = w.findChild<Canvas*>();
    REQUIRE(w.openFile(QString(PENZENE_TEST_DATA) + "/aspirin.mol"));  // as `penzene aspirin.mol` does
    w.show();
    QApplication::processEvents();
    const QRect drawn = canvas->mapFromScene(documentBounds(canvas->document())).boundingRect();
    const QSize view = canvas->viewport()->size();
    CHECK((drawn.width() > view.width() / 2 || drawn.height() > view.height() / 2));  // it fills the view, not a corner of it
    CHECK(canvas->viewport()->rect().contains(drawn));
}

TEST_CASE("accessibility: the hotspot is announced to screen readers (#112)") {
    Fixture f;
    f.canvas.setDocumentSilently(*chem::fromSmiles("CO"));
    f.hover(f.doc().atoms[1].pos);
    f.canvas.viewport()->repaint();
    CHECK(f.canvas.accessibleDescription() == "Hotspot: atom O2, 1 bond");
    f.hover((f.doc().atoms[0].pos + f.doc().atoms[1].pos) / 2);
    f.canvas.viewport()->repaint();
    CHECK(f.canvas.accessibleDescription() == "Hotspot: single bond, C1 to O2");
}

TEST_CASE("accessibility: the drawing's atoms, bonds, arrows and text can be read, found and selected (#536)") {
    Fixture f;
    Document d = *chem::fromSmiles("CCO");
    Arrow curved{d.atoms[2].pos + QPointF(0, -8), (d.atoms[0].pos + d.atoms[1].pos) / 2};
    curved.bend = 6, curved.fromAt = {2, -1}, curved.toAt = {0, 1};
    d.arrows = {curved, Arrow{{0, 40}, {60, 40}}};
    d.texts.push_back({{0, 70}, "heat"});
    f.canvas.setDocumentSilently(d);
    QAccessibleInterface* canvas = QAccessible::queryAccessibleInterface(&f.canvas);
    REQUIRE(canvas);
    auto find = [&](const QString& name) -> QAccessibleInterface* {
        for (int i = 0; i < canvas->childCount(); ++i)
            if (auto* c = canvas->child(i); c && c->text(QAccessible::Name) == name) return c;
        return nullptr;
    };
    auto press = [](QAccessibleInterface* c) { c->actionInterface()->doAction(QAccessibleActionInterface::pressAction()); };

    QAccessibleInterface* o = find("atom O3, 1 bond");
    REQUIRE(o);
    CHECK(o->parent() == canvas);
    CHECK(canvas->child(canvas->indexOfChild(o)) == o);
    const QPoint onScreen = f.canvas.viewport()->mapToGlobal(f.at(f.doc().atoms[2].pos));
    CHECK(o->rect().contains(onScreen));
    CHECK(canvas->childAt(onScreen.x(), onScreen.y()) == o);
    press(o);
    CHECK(f.canvas.selection() == QSet<int>{2});
    CHECK(o->state().selected);
    CHECK(o->state().focused);  // the hotspot: a hotkey now builds from the O
    f.key("n");
    CHECK(f.doc().atoms[2].z == 7);

    REQUIRE(find("single bond, C2 to N3"));
    press(find("single bond, C2 to N3"));
    CHECK(f.canvas.selection() == QSet<int>{1, 2});
    REQUIRE(find("curved arrow, from N3 to the C1–C2 bond"));
    press(find("reaction arrow"));
    CHECK(f.canvas.selectedArrows() == QSet<int>{1});
    press(find("text: heat"));
    CHECK(f.canvas.selectedTexts() == QSet<int>{0});

    // A part stays as the pointer moves, stops being valid when the drawing changes, and is then dropped.
    QAccessibleInterface* n = find("atom N3, 1 bond");
    REQUIRE(n);
    f.hover(f.doc().atoms[0].pos);
    f.canvas.setTheme(Theme{});
    CHECK(n->isValid());
    const QAccessible::Id id = QAccessible::uniqueId(n);
    const int before = canvas->childCount();
    f.canvas.setSelection({2});
    f.canvas.deleteSelection();
    CHECK_FALSE(n->isValid());
    CHECK(canvas->childCount() == before - 2);  // the N and its bond
    CHECK(QAccessible::accessibleInterface(id) == nullptr);
    CHECK_FALSE(find("atom N3, 1 bond"));
}

TEST_CASE("accessibility: every control is named, in the window, its closed flyouts and its dialogs (#536)") {
    App app;
    MainWindow w;
    w.resize(1100, 750);
    w.show();
    auto* canvas = w.findChild<Canvas*>();
    canvas->setDocumentSilently(*chem::fromSmiles("CCO"));
    canvas->selectAll();
    QStringList unnamed;
    auto walk = [&](QWidget* top, const QString& where) {
        for (QWidget* c : top->findChildren<QWidget*>()) {
            if (c->focusPolicy() == Qt::NoFocus && !qobject_cast<QAbstractButton*>(c) || qobject_cast<QLabel*>(c)) continue;  // not a control
            if (c->objectName().startsWith("qt_") || qobject_cast<QAbstractScrollArea*>(c->parentWidget())) continue;  // Qt's own parts
            bool inside = false;  // a combo box's list or a spin box's editor: read as the box
            for (QWidget* p = c->parentWidget(); p; p = p->parentWidget()) inside |= qobject_cast<QComboBox*>(p) || qobject_cast<QAbstractSpinBox*>(p);
            if (inside) continue;
            if (auto* ai = QAccessible::queryAccessibleInterface(c); ai && ai->text(QAccessible::Name).trimmed().isEmpty())
                unnamed << where + ": " + c->metaObject()->className() + " " + c->objectName() + " " + c->toolTip();
        }
    };
    walk(&w, "window");
    auto dialog = [&](const QString& what, const std::function<void()>& open) {
        QTimer::singleShot(0, &w, [&, what] {
            if (QWidget* modal = QApplication::activeModalWidget()) {
                walk(modal, what);
                if (auto* d = qobject_cast<QDialog*>(modal)) d->reject();
                else modal->close();
            }
        });
        open();
        QApplication::processEvents();
    };
    // Every dialog a menu opens ("…"), except the system's own file, print and colour pickers.
    for (QAction* a : w.findChildren<QAction*>()) {
        const QString text = a->text();
        if (text.endsWith(u'…') && a->isEnabled() && !text.contains("Open") && !text.contains("Save") && !text.contains("Print") &&
            !text.contains("Export") && !text.contains("Colour"))
            dialog(text, [a] { a->trigger(); });
    }
    dialog("atom properties", [&] { canvas->editAtomProperties(0); });
    dialog("text", [&] { canvas->editText(-1, {}); });
    INFO(unnamed.join("\n").toStdString());
    CHECK(unnamed.isEmpty());
}

TEST_CASE("accessibility: an edit at the hotspot is announced too (#244)") {
    Fixture f;
    f.canvas.setDocumentSilently(*chem::fromSmiles("CC"));
    f.hover(f.doc().atoms[1].pos);
    f.canvas.viewport()->repaint();
    CHECK(f.canvas.accessibleDescription() == "Hotspot: atom C2, 1 bond");
    f.key("n");
    f.canvas.viewport()->repaint();
    CHECK(f.canvas.accessibleDescription() == "Hotspot: atom N2, 1 bond");
}

TEST_CASE("the Figure tools include a dashed ellipse (#222)") {
    App app;
    MainWindow w;
    w.resize(1000, 700);
    w.show();
    auto* canvas = w.findChild<Canvas*>();
    QAction* dashedEllipse = nullptr;
    for (auto* b : w.findChildren<QToolButton*>())
        if (b->defaultAction() && b->defaultAction()->toolTip().startsWith("Dashed ellipse")) dashedEllipse = b->defaultAction();
    REQUIRE(dashedEllipse);
    dashedEllipse->trigger();
    QTest::mousePress(canvas->viewport(), Qt::LeftButton, {}, canvas->mapFromScene(QPointF(0, 0)));
    const QPoint to = canvas->mapFromScene(QPointF(60, 40));
    QMouseEvent move(QEvent::MouseMove, to, canvas->viewport()->mapToGlobal(to), Qt::NoButton, Qt::LeftButton, {});
    QApplication::sendEvent(canvas->viewport(), &move);
    QTest::mouseRelease(canvas->viewport(), Qt::LeftButton, {}, to);
    REQUIRE(canvas->document().arrows.size() == 1);
    CHECK(canvas->document().arrows[0].kind == ArrowKind::Ellipse);
    CHECK(canvas->document().arrows[0].dashed);
}

TEST_CASE("the selection's rotate handle turns it; Shift snaps to 15° steps, Ctrl to 45° from the page (#218)") {
    Fixture f;
    f.canvas.setTool(Canvas::Tool::Select);
    Document d;
    d.atoms = {{{0, 0}}, {{kBondLength, 0}}};
    d.bonds = {{0, 1}};
    // Drags the handle to `degrees` clockwise from straight up, about the pair's middle.
    auto turn = [&](double degrees, Qt::KeyboardModifiers mods) {
        f.canvas.setDocumentSilently(d);
        f.canvas.setSelection({0, 1});
        const auto knob = f.canvas.rotateHandle();
        REQUIRE(knob);
        const QPointF c(kBondLength / 2, 0);
        const double r = QLineF(*knob, c).length(), a = qDegreesToRadians(degrees);
        const QPoint to = f.at(c + r * QPointF(std::sin(a), -std::cos(a)));
        QTest::mousePress(f.canvas.viewport(), Qt::LeftButton, mods, f.at(*knob));
        QMouseEvent move(QEvent::MouseMove, to, f.canvas.viewport()->mapToGlobal(to), Qt::NoButton, Qt::LeftButton, mods);
        QApplication::sendEvent(f.canvas.viewport(), &move);
        QTest::mouseRelease(f.canvas.viewport(), Qt::LeftButton, mods, to);
        const QPointF v = f.doc().atoms[1].pos - f.doc().atoms[0].pos;  // was (1, 0)
        return qRadiansToDegrees(std::atan2(v.y(), v.x()));
    };
    CHECK(std::abs(turn(90, {}) - 90) < 1.5);
    CHECK(std::abs(turn(70, Qt::ShiftModifier) - 75) < 0.01);
    CHECK(std::abs(turn(70, Qt::ControlModifier) - 90) < 0.01);
    // Drawn 1° off: Shift keeps the 1° (15° steps from the start), Ctrl squares it up.
    d.atoms[1].pos = kBondLength * QPointF(std::cos(qDegreesToRadians(1.0)), std::sin(qDegreesToRadians(1.0)));
    CHECK(std::abs(turn(40, Qt::ShiftModifier) - 46) < 0.01);
    CHECK(std::abs(turn(40, Qt::ControlModifier) - 45) < 0.01);
    CHECK(f.canvas.document().bonds.size() == 1);  // turned, not redrawn
}

TEST_CASE("Fit to Window zooms to the selection") {
    Fixture f;
    Document d = *chem::fromSmiles("c1ccccc1");
    Document far = *chem::fromSmiles("CCO");
    d.append(far, {600, 400});
    f.canvas.setDocumentSilently(d);
    f.canvas.fitToDocument();
    const double all = f.canvas.transform().m11();
    f.canvas.setSelection({0, 1, 2, 3, 4, 5});  // the benzene
    f.canvas.fitToSelection();
    CHECK(f.canvas.transform().m11() > 2 * all);
    const QPoint centre = f.canvas.mapFromScene(documentBounds(f.canvas.selectedSubset()).center());
    CHECK((centre - f.canvas.viewport()->rect().center()).manhattanLength() < 4);
    CHECK_FALSE(f.canvas.viewport()->rect().contains(f.canvas.mapFromScene(d.atoms.back().pos)));  // the ethanol is off screen
    f.canvas.setSelection({0});
    f.canvas.fitToSelection();
    CHECK(f.canvas.transform().m11() <= 10.001);  // one atom: close up, not window-filling
}

TEST_CASE("a selection never points at other atoms after an erase or undo renumbers them") {
    Fixture f;
    Document d;
    for (int i = 0; i < 5; ++i) d.addAtom({i * 40.0, 0});
    f.canvas.setDocumentSilently(d);
    f.canvas.setTool(Canvas::Tool::Erase);
    f.click({0, 0});  // atom 0 goes: the rest move down one
    REQUIRE(f.doc().atoms.size() == 4);
    f.canvas.setSelection({2, 3});
    f.undo.undo();  // atom 0 is back: {2, 3} would now be the wrong atoms
    REQUIRE(f.doc().atoms.size() == 5);
    CHECK(f.canvas.selection().isEmpty());
}

TEST_CASE("opening a file drops the previous drawing's selection") {
    App app;
    QTemporaryDir dir;
    Document d;
    for (int i = 0; i < 3; ++i) d.addAtom({i * 40.0, 0});
    const QString path = dir.filePath("three.penz");
    QFile out(path);
    REQUIRE(out.open(QIODevice::WriteOnly));
    out.write(d.toJson());
    out.close();
    MainWindow w;
    auto* canvas = w.findChild<Canvas*>();
    canvas->setDocumentSilently(d);
    canvas->setSelection({0, 1});
    REQUIRE(w.openFile(path));  // same atom count: nothing else would clear it
    CHECK(canvas->selection().isEmpty());
}

TEST_CASE("the atom tool relabels an abbreviation") {
    Fixture f;
    Document d;
    d.addAtom({0, 0});
    d.atoms[0].label = "Ph";
    f.canvas.setDocumentSilently(d);
    f.canvas.setTool(Canvas::Tool::Atom);
    f.canvas.setElement(7);
    f.click({0, 0});
    CHECK(f.doc().atoms[0].z == 7);
    CHECK(f.doc().atoms[0].label.isEmpty());
}

TEST_CASE("a click that wobbles under the drag threshold moves nothing") {
    Fixture f;
    Document d;
    d.addAtom({0, 0}), d.addAtom({kBondLength, 0});
    d.bonds = {{0, 1}};
    f.canvas.setDocumentSilently(d);
    f.canvas.setTool(Canvas::Tool::Select);
    f.canvas.selectAll();
    f.drag({0, 0}, {0.8, 0});  // 2 screen pixels at the default zoom: a click
    CHECK(f.doc() == d);
}

TEST_CASE("an edit that changes nothing is not an undo step") {
    Fixture f;
    Document d;
    d.addAtom({0, 0}), d.addAtom({kBondLength, 0});
    d.bonds = {{0, 1}};
    f.canvas.setDocumentSilently(d);
    f.canvas.commit(chem::removeHydrogens(d), "Remove hydrogens");  // there are none
    CHECK(f.undo.count() == 0);
    CHECK(f.undo.isClean());
}

TEST_CASE("Save replaces the file whole, and still saves where only the file is writable") {
    App app;
    QTemporaryDir dir;
    const QString path = dir.filePath("mol.penz");
    Document d;
    d.addAtom({0, 0});
    QFile out(path);
    REQUIRE(out.open(QIODevice::WriteOnly));
    out.write(d.toJson());
    out.close();
    MainWindow w;
    REQUIRE(w.openFile(path));
    auto* canvas = w.findChild<Canvas*>();
    Document two = d;
    two.addAtom({kBondLength, 0});
    canvas->commit(two, "Add");
    // A folder that can't take the temporary file (as behind a sandbox's document portal).
    QFile::setPermissions(dir.path(), QFileDevice::ReadOwner | QFileDevice::ExeOwner);
    QTimer::singleShot(0, &w, [] {  // a failed save warns: dismiss it rather than hang
        if (auto* box = qobject_cast<QMessageBox*>(QApplication::activeModalWidget())) box->reject();
    });
    for (auto* a : w.findChildren<QAction*>())
        if (a->text() == "&Save") a->trigger();
    QFile::setPermissions(dir.path(), QFileDevice::ReadOwner | QFileDevice::WriteOwner | QFileDevice::ExeOwner);
    QFile in(path);
    REQUIRE(in.open(QIODevice::ReadOnly));
    const auto saved = Document::fromJson(in.readAll());
    REQUIRE(saved);
    CHECK(saved->atoms.size() == 2);
}

TEST_CASE("a file dropped on the window opens") {
    App app;
    QTemporaryDir dir;
    const QString path = dir.filePath("ethane.penz");
    QFile out(path);
    REQUIRE(out.open(QIODevice::WriteOnly));
    out.write(chem::fromSmiles("CC")->toJson());
    out.close();
    MainWindow w;
    w.show();
    auto* canvas = w.findChild<Canvas*>();
    CHECK_FALSE(canvas->viewport()->acceptDrops());  // the canvas lets it through to the window
    QMimeData mime;
    mime.setUrls({QUrl::fromLocalFile(path)});
    const QPointF at = w.rect().center();
    QDragEnterEvent enter(at.toPoint(), Qt::CopyAction, &mime, Qt::LeftButton, {});
    QApplication::sendEvent(&w, &enter);
    REQUIRE(enter.isAccepted());
    QDropEvent drop(at, Qt::CopyAction, &mime, Qt::LeftButton, {});
    QApplication::sendEvent(&w, &drop);
    CHECK(canvas->document().atoms.size() == 2);
    CHECK(w.windowTitle().startsWith("ethane.penz"));
}

TEST_CASE("formula HTML: counts subscripted, charge superscripted") {
    CHECK(formulaHtml("O4S-2") == "O<sub>4</sub>S<sup>2−</sup>");
    CHECK(formulaHtml("C2H3O2-") == "C<sub>2</sub>H<sub>3</sub>O<sub>2</sub><sup>−</sup>");
    CHECK(formulaHtml("C6H6") == "C<sub>6</sub>H<sub>6</sub>");
    CHECK(formulaHtml("C[13C]H6O") == "C<sup>13</sup>CH<sub>6</sub>O");  // an isotope, apart (#285)
}

TEST_CASE("the formula counts isotopes apart: 13C isn't C, D isn't H (#285)") {
    CHECK(chem::properties(*chem::fromSmiles("[13CH3]CO"))->formula == "C[13C]H6O");
    CHECK(chem::properties(*chem::fromSmiles("[2H]C(Cl)(Cl)Cl"))->formula == "CDCl3");
}

TEST_CASE("Properties → Copy as Text keeps the formula's charge unambiguous") {
    App app;
    MainWindow w;
    w.show();
    auto* canvas = w.findChild<Canvas*>();
    canvas->setDocumentSilently(*chem::fromSmiles("[O-]S(=O)(=O)[O-]"));
    w.findChild<QDockWidget*>("properties")->show();
    for (auto* b : w.findChildren<QPushButton*>())
        if (b->text() == "Copy as Text") b->click();
    CHECK(QApplication::clipboard()->text().startsWith("Formula\tO4S-2\n"));
}

TEST_CASE("Delete removes an attachment point outright, not turning it into a carbon") {
    Fixture f;
    Document d;
    d.addAtom({-kBondLength, 0}), d.addAtom({0, 0}), d.addAtom({kBondLength, 0}, 0);
    d.bonds = {{0, 1}, {1, 2, 1, BondStereo::Wavy}};
    f.canvas.setDocumentSilently(d);
    f.hover({kBondLength, 0});
    QTest::keyClick(f.canvas.viewport(), Qt::Key_Delete);
    CHECK(f.doc().atoms.size() == 2);  // the ethyl stays
}

TEST_CASE("a dashed arrow's head is solid") {
    App app;
    auto head = [](ArrowKind kind, bool dashed) {
        Document d;
        d.arrows.push_back({{0, 0}, {50, 0}, kind, 0, {}, dashed});
        const QImage img = renderImage(d, {288});  // 4 px per point; the frame starts 4 pt left of the tail
        return img.copy(QRect(4 * (4 + 46.5), 0, 4 * 4, img.height()));  // the head beyond where the shafts stop
    };
    for (auto kind : {ArrowKind::Reaction, ArrowKind::Retro}) CHECK(head(kind, true) == head(kind, false));
}

TEST_CASE("text: a formula's charge is a superscript, its counts stay subscripts") {
    App app;
    const DrawingStyle& st = drawingStyle("");
    auto box = [&](const char* s) { return textPath({{0, 0}, s}, st).boundingRect(); };
    const double cap = QFontMetricsF(labelFont(st)).capHeight();
    auto raised = [&](const char* s, const char* plain) { return box(s).top() < box(plain).top() - 0.2 * cap; };
    auto lowered = [&](const char* s, const char* plain) { return box(s).bottom() > box(plain).bottom() + 0.2 * cap; };
    CHECK(raised("Cu2+", "Cu"));  // one element: the 2 is the charge
    CHECK_FALSE(lowered("Cu2+", "Cu"));
    CHECK(raised("NH4+", "NH"));  // the 4 counts hydrogens
    CHECK(lowered("NH4+", "NH"));
    CHECK(raised("SO4^2-", "SO"));
    CHECK(box("SO4^2-").width() < box("SO4^^2-").width());  // the ^ isn't drawn
    CHECK(raised("[Fe(CN)6]3-", "[Fe(CN)]"));
    CHECK(raised("Na+, Cl-", "Na, Cl"));
    CHECK(lowered("H2O", "HO"));
    CHECK_FALSE(raised("H2O", "HO"));
    CHECK_FALSE(raised("cis- and trans-", "cis and trans"));  // words, not formulas
    CHECK_FALSE(raised("THF, -78 °C", "THF, 78 °C"));
    CHECK_FALSE(raised("CH2Cl2-MeOH", "CH2Cl2 MeOH"));
}

TEST_CASE("atom numbers stay clear of lone pairs and δ on the same atom (#349)") {
    App app;
    const Document base = *chem::fromSmiles("C[C@H](N)C(=O)O");
    // Ink painted, at 8x; a mark's ink doesn't depend on where it goes, so ink lost
    // when numbers and marks are both shown is ink they share.
    auto ink = [&](bool numbers, bool marks) {
        Document d = base;
        d.showAtomNumbers = numbers;
        for (auto& a : d.atoms)
            if (a.z == 8 && marks) a.lonePairs = 2, a.partial = -1;
        QImage img(1200, 1200, QImage::Format_ARGB32_Premultiplied);
        img.fill(Qt::transparent);
        QPainter p(&img);
        p.translate(600, 600);
        p.scale(8, 8);
        p.translate(-documentBounds(base).center());
        paintDocument(p, d);
        p.end();
        if (auto out = qgetenv("PENZENE_NUMBERS_SHOT"); !out.isEmpty() && numbers && marks) img.save(out);
        int n = 0;
        for (int y = 0; y < img.height(); ++y)
            for (int x = 0; x < img.width(); ++x) n += qAlpha(img.pixel(x, y)) > 128;
        return n;
    };
    const int plain = ink(false, false), numbers = ink(true, false) - plain, marks = ink(false, true) - plain;
    const int shared = plain + numbers + marks - ink(true, true);
    INFO("numbers " << numbers << " marks " << marks << " shared " << shared);
    // Edge pixels, not a digit on a dot: 0.04% here (macOS) and 1% on Windows' larger glyphs; 4.6% before.
    CHECK(shared < marks / 50);
}

#ifndef _WIN32
#include <csignal>
#include <sys/resource.h>

TEST_CASE("an export that can't be written in full leaves the old file whole (#304)") {
    App app;
    QTemporaryDir dir;
    auto doc = chem::fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
    REQUIRE(doc);
    for (const char* name : {"old.png", "old.svg", "old.pdf"}) {
        const QString path = dir.filePath(name);
        QFile old(path);
        REQUIRE(old.open(QIODevice::WriteOnly));
        old.write("old");
        old.close();
        // Files may grow to 64 bytes: as a full disk, every export stops part way.
        auto previous = std::signal(SIGXFSZ, SIG_IGN);
        rlimit saved{};
        getrlimit(RLIMIT_FSIZE, &saved);
        rlimit small = saved;
        small.rlim_cur = 64;
        setrlimit(RLIMIT_FSIZE, &small);
        const bool ok = exportDocument(*doc, path);
        setrlimit(RLIMIT_FSIZE, &saved);
        std::signal(SIGXFSZ, previous);
        CHECK_FALSE(ok);
        REQUIRE(old.open(QIODevice::ReadOnly));
        CHECK(old.readAll() == "old");
        old.close();
    }
}
#endif

TEST_CASE("A trackpad pinch zooms the canvas: spread in, pinch out (#452)") {
    Fixture f;
    auto pinch = [&](double value) {
        QNativeGestureEvent e(Qt::ZoomNativeGesture, QPointingDevice::primaryPointingDevice(), 2, {10, 10}, {10, 10},
                              {10, 10}, value, {}, {});
        QApplication::sendEvent(f.canvas.viewport(), &e);
    };
    const double before = f.canvas.transform().m11();
    pinch(0.2);
    CHECK(std::abs(f.canvas.transform().m11() / before - 1.2) < 1e-9);
    pinch(-0.2);
    CHECK(f.canvas.transform().m11() < before * 1.2);
}

TEST_CASE("Center on Page moves the whole drawing, or just the selection, keeping spacing (#454)") {
    Fixture f;
    Document doc = *chem::fromSmiles("CC.O");
    doc.page = "ACS single column";
    f.canvas.setDocumentSilently(doc);
    const QPointF middle = pageRect(doc).center();
    const auto gap = [&] { return f.canvas.document().atoms[0].pos - f.canvas.document().atoms[2].pos; };
    const QPointF spacing = gap();
    f.canvas.centerOnPage();
    CHECK(QLineF(documentBounds(f.canvas.document()).center(), middle).length() < 1e-6);
    CHECK(gap() == spacing);
    f.undo.undo();
    CHECK(f.canvas.document() == doc);

    f.canvas.setSelection({2}, {}, {});  // the lone oxygen
    f.canvas.centerOnPage();
    CHECK(QLineF(f.canvas.document().atoms[2].pos, middle).length() < 1e-6);
    CHECK(f.canvas.document().atoms[0].pos == doc.atoms[0].pos);
}

TEST_CASE("View shows a light grid and rulers, measured at the final size (#219)") {
    Fixture f;
    f.canvas.setTheme(theme("Light"));
    f.canvas.setDocumentSilently(Document{});
    const QPoint origin = f.at({0, 0});  // with no page, guides count from the origin
    auto pixel = [&](QPoint at) { return f.canvas.viewport()->grab().toImage().pixelColor(at); };
    const QPoint onLine(origin.x(), origin.y() + 3), band(4, 4);  // band: the rulers' corner
    CHECK(pixel(onLine) == theme("Light").paper);
    CHECK(pixel(band) == theme("Light").paper);
    f.canvas.setGuides(true, false);
    CHECK(pixel(onLine) != theme("Light").paper);  // a grid line through the origin
    CHECK(pixel(band) == theme("Light").paper);
    f.canvas.setGuides(false, true);
    CHECK(pixel(onLine) == theme("Light").paper);
    CHECK(pixel(band) == theme("Light").window);  // the rulers
    if (auto out = qgetenv("PENZENE_GUIDES_SHOT"); !out.isEmpty()) {
        Document aspirin = *chem::fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
        aspirin.page = "A4";
        f.canvas.setDocumentSilently(aspirin);
        f.canvas.setGuides(true, true);
        f.canvas.fitToDocument();
        f.canvas.grab().save(out);
    }
}

TEST_CASE("pages: tabs along the bottom, each with its own drawing and undo history (#219)") {
    App app;
    QFile::remove(MainWindow::autosavePath());
    MainWindow w;
    auto* canvas = w.findChild<Canvas*>();
    auto* tabs = w.findChild<QTabBar*>("pageTabs");
    REQUIRE(tabs);
    REQUIRE(tabs->count() == 1);
    CHECK(tabs->tabText(0) == "Page 1");
    QAction* undo = nullptr;
    for (auto* a : w.findChildren<QAction*>())
        if (a->shortcut() == QKeySequence::Undo) undo = a;
    REQUIRE(undo);

    canvas->commit(*chem::fromSmiles("CCO"), "Draw");  // page 1
    w.findChild<QToolButton*>("addPage")->click();
    REQUIRE(tabs->count() == 2);
    CHECK(tabs->currentIndex() == 1);
    CHECK(tabs->tabText(1) == "Page 2");
    CHECK(canvas->document().empty());
    CHECK_FALSE(undo->isEnabled());  // page 2 has no history of its own yet
    canvas->commit(*chem::fromSmiles("c1ccccc1"), "Draw");  // page 2

    tabs->setCurrentIndex(0);
    CHECK(chem::toSmiles(canvas->document()) == "CCO");
    undo->trigger();  // undoes page 1's edit, not page 2's
    CHECK(canvas->document().empty());
    tabs->setCurrentIndex(1);
    CHECK(chem::toSmiles(canvas->document()) == "c1ccccc1");
    CHECK(w.isWindowModified());

    // Move a selection over to page 1: cut here, pasted there.
    canvas->selectAll();
    QMenu* moveTo = nullptr;
    for (auto* m : w.findChildren<QMenu*>())
        if (m->title() == "Mo&ve Selection To") moveTo = m;
    REQUIRE(moveTo);
    emit moveTo->aboutToShow();
    REQUIRE(moveTo->actions().size() == 1);
    moveTo->actions()[0]->trigger();
    CHECK(tabs->currentIndex() == 0);
    CHECK(chem::toSmiles(canvas->document()) == "c1ccccc1");
    tabs->setCurrentIndex(1);
    CHECK(canvas->document().empty());

    // Every page is written, and a file with pages opens as tabs.
    w.autosave();
    QFile f(MainWindow::autosavePath());
    REQUIRE(f.open(QIODevice::ReadOnly));
    const auto sheets = sheetsFromJson(f.readAll());
    f.close();
    REQUIRE(sheets.size() == 2);
    CHECK(chem::toSmiles(sheets[0].doc) == "c1ccccc1");
    CHECK(sheets[1].doc.empty());
    QTemporaryDir dir;
    const QString path = dir.filePath("pages.penz");
    QFile out(path);
    REQUIRE(out.open(QIODevice::WriteOnly));
    out.write(sheetsToJson({{"Scheme", *chem::fromSmiles("CCO")}, {"Notes", *chem::fromSmiles("O")}}));
    out.close();
    MainWindow other;
    REQUIRE(other.openFile(path));
    auto* otherTabs = other.findChild<QTabBar*>("pageTabs");
    REQUIRE(otherTabs->count() == 2);
    CHECK(otherTabs->tabText(1) == "Notes");
    CHECK_FALSE(other.isWindowModified());
    QFile::remove(MainWindow::autosavePath());
    if (auto shot = qgetenv("PENZENE_PAGES_SHOT"); !shot.isEmpty()) {
        other.resize(1000, 650);
        other.show();
        other.openFile(QString(PENZENE_TEST_DATA) + "/aspirin.mol");
        other.findChild<QToolButton*>("addPage")->click();
        other.findChild<QToolButton*>("addPage")->click();
        other.findChild<QTabBar*>("pageTabs")->setCurrentIndex(0);
        QApplication::processEvents();
        other.grab().save(QString::fromUtf8(shot));
    }
}

TEST_CASE("compound numbers run on across the pages, in page order (#569)") {
    App app;
    MainWindow w;
    auto* canvas = w.findChild<Canvas*>();
    auto* tabs = w.findChild<QTabBar*>("pageTabs");
    QAction *number = nullptr, *undo = nullptr;
    for (auto* a : w.findChildren<QAction*>()) {
        if (a->text() == "&Number Compounds") number = a;
        if (a->shortcut() == QKeySequence::Undo) undo = a;
    }
    REQUIRE(number);
    REQUIRE(undo);
    auto numbers = [&](int page) {
        const auto sheets = w.sheets();
        QStringList out;
        for (const Text& t : sheets[page].doc.texts) out << t.text;
        return out.join(' ').toStdString();
    };
    canvas->setDocumentSilently(*chem::fromSmiles("CCO"));
    number->trigger();
    w.findChild<QToolButton*>("addPage")->click();
    canvas->setDocumentSilently(*chem::fromSmiles("c1ccccc1"));
    number->trigger();
    CHECK(numbers(0) == "1");
    CHECK(numbers(1) == "2");  // on from page 1

    tabs->setCurrentIndex(0);
    Document more = canvas->document();
    more.append(*chem::fromSmiles("N"), {20 * kBondLength, 0});
    canvas->commit(more, "Draw");
    number->trigger();
    CHECK(numbers(0) == "1 2");
    CHECK(numbers(1) == "3");  // an earlier page's new number moves the later ones on
    tabs->setCurrentIndex(1);
    undo->trigger();  // page 2's own history: back to before its number
    CHECK(numbers(1).empty());
    for (auto* a : w.findChildren<QAction*>())
        if (a->shortcut() == QKeySequence::Redo) a->trigger();
    CHECK(numbers(1) == "3");  // redone as it is now, not as it was saved

    QUndoStack* history = w.findChild<QUndoGroup*>()->activeStack();
    const int steps = history->count();
    canvas->commit(canvas->document(), "Nothing");
    CHECK(history->count() == steps);  // numbered on from page 1 already: no change, no step

    tabs->moveTab(1, 0);  // benzene's page first
    CHECK(numbers(0) == "1");
    CHECK(numbers(1) == "2 3");

    Document shared = w.sheets()[0].doc;  // 1a on one page, 1b on the next: still one compound
    shared.texts[0].text = "1a";
    tabs->setCurrentIndex(0);
    canvas->commit(shared, "Edit");
    tabs->setCurrentIndex(1);
    Document next = canvas->document();
    next.texts[0].text = "1b";
    canvas->commit(next, "Edit");
    CHECK(numbers(0) == "1a");
    CHECK(numbers(1) == "1b 2");
}

TEST_CASE("a series shared across pages keeps one number through an insert, undo and redo (#583)") {
    App app;
    MainWindow w;
    auto* canvas = w.findChild<Canvas*>();
    auto* tabs = w.findChild<QTabBar*>("pageTabs");
    QAction *number = nullptr, *undo = nullptr, *redo = nullptr;
    for (auto* a : w.findChildren<QAction*>()) {
        if (a->text() == "&Number Compounds") number = a;
        if (a->shortcut() == QKeySequence::Undo) undo = a;
        if (a->shortcut() == QKeySequence::Redo) redo = a;
    }
    REQUIRE(number);
    REQUIRE(undo);
    REQUIRE(redo);
    auto numbers = [&](int page) {
        QStringList out;  // the page on the canvas is the canvas's
        const Document d = page == tabs->currentIndex() ? canvas->document() : w.sheets()[page].doc;
        for (const Text& t : d.texts) out << t.text;
        out.sort();
        return out.join(' ').toStdString();
    };
    Document first = *chem::fromSmiles("CCO");
    first.append(*chem::fromSmiles("CN"), {20 * kBondLength, 0});
    canvas->setDocumentSilently(first);
    number->trigger();
    Document lettered = canvas->document();
    for (Text& t : lettered.texts)
        if (t.text == "2") t.text = "2a";
    canvas->commit(lettered, "Edit");
    w.findChild<QToolButton*>("addPage")->click();
    canvas->setDocumentSilently(*chem::fromSmiles("c1ccccc1"));
    number->trigger();
    Document shared = canvas->document();
    shared.texts[0].text = "2b";
    canvas->commit(shared, "Edit");
    REQUIRE(numbers(0) == "1 2a");
    REQUIRE(numbers(1) == "2b");

    tabs->setCurrentIndex(0);  // a step inserted between 1 and 2a
    Document inserted = canvas->document();
    inserted.append(*chem::fromSmiles("CCl"), {10 * kBondLength, 0});
    canvas->commit(inserted, "Draw");
    number->trigger();
    CHECK(numbers(0) == "1 2 3a");
    CHECK(numbers(1) == "3b");
    undo->trigger();
    undo->trigger();
    CHECK(numbers(0) == "1 2a");
    CHECK(numbers(1) == "2b");
    redo->trigger();
    redo->trigger();
    CHECK(numbers(0) == "1 2 3a");
    CHECK(numbers(1) == "3b");
    tabs->setCurrentIndex(1);
    tabs->setCurrentIndex(0);
    CHECK(numbers(1) == "3b");
}

TEST_CASE("opening a file while on a later page doesn't touch freed undo stacks (#359)") {
    App app;
    QTemporaryDir dir;
    const QString three = dir.filePath("three.penz");
    {
        QFile f(three);
        REQUIRE(f.open(QIODevice::WriteOnly));
        f.write(sheetsToJson({{"A", *chem::fromSmiles("C")}, {"B", *chem::fromSmiles("N")}, {"C", *chem::fromSmiles("O")}}));
    }
    MainWindow w;
    REQUIRE(w.openFile(three));
    auto* tabs = w.findChild<QTabBar*>("pageTabs");
    tabs->setCurrentIndex(2);  // the last page's stack is the group's active one
    REQUIRE(w.openFile(QString(PENZENE_TEST_DATA) + "/aspirin.mol"));
    CHECK(tabs->count() == 1);
}

TEST_CASE("Save after opening an SDF doesn't overwrite it, and never writes MOL into an image (#315)") {
    App app;
    QTemporaryDir dir;
    const QString sdf = dir.filePath("library.sdf");
    QFile f(sdf);
    REQUIRE(f.open(QIODevice::WriteOnly));
    const QByteArray library = QByteArray::fromStdString(chem::toMolBlock(*chem::fromSmiles("CCO")) + "$$$$\n" +
                                                         chem::toMolBlock(*chem::fromSmiles("N")) + "$$$$\n");
    f.write(library);
    f.close();
    MainWindow w;
    REQUIRE(w.openFile(sdf));
    w.findChild<Canvas*>()->commit(*chem::fromSmiles("CCC"), "edit");

    // Save asks where to save instead; answer with an image path. The refusal is dismissed too.
    QApplication::setAttribute(Qt::AA_DontUseNativeDialogs);  // a native dialog can't be answered from here
    QString offered, image;
    QTimer::singleShot(0, &w, [&] {
        auto* dialog = qobject_cast<QFileDialog*>(QApplication::activeModalWidget());
        REQUIRE(dialog);
        offered = dialog->selectedFiles().value(0);
        QTimer::singleShot(0, &w, [] {
            if (auto* box = qobject_cast<QMessageBox*>(QApplication::activeModalWidget())) box->reject();
        });
        image = dir.filePath("drawing.png");
        dialog->selectFile(image);
        static_cast<QDialog*>(dialog)->accept();
    });
    for (auto* a : w.findChildren<QAction*>())
        if (a->text() == "&Save") a->trigger();
    QApplication::setAttribute(Qt::AA_DontUseNativeDialogs, false);
    CHECK_FALSE(offered.endsWith("library.sdf"));
    QFile back(sdf);
    REQUIRE(back.open(QIODevice::ReadOnly));
    CHECK(back.readAll() == library);
    REQUIRE_FALSE(image.isEmpty());
    CHECK_FALSE(QFile::exists(image));
}

TEST_CASE("About opens bundled third-party license notices (#451)") {
    App app;
    MainWindow w;
    bool opened = false;
    QTimer::singleShot(0, &w, [&] {
        auto* about = qobject_cast<QMessageBox*>(QApplication::activeModalWidget());
        REQUIRE(about);
        QPushButton* licenses = nullptr;
        for (auto* button : about->findChildren<QPushButton*>())
            if (button->text().contains("Third-party licenses")) licenses = button;
        REQUIRE(licenses);
        QTimer::singleShot(0, &w, [&] {
            auto* dialog = qobject_cast<QDialog*>(QApplication::activeModalWidget());
            REQUIRE(dialog);
            auto* text = dialog->findChild<QTextBrowser*>();
            REQUIRE(text);
            opened = text->toPlainText().contains("Revvity ChemDraw file library") &&
                     text->toPlainText().contains("GNU Lesser General Public License") &&
                     text->toPlainText().contains("nmrshiftdb2 Database License");
            dialog->accept();
        });
        licenses->click();
    });
    for (auto* action : w.findChildren<QAction*>())
        if (action->text() == "&About Penzene") action->trigger();
    CHECK(opened);
}

TEST_CASE("another running Penzene's autosave is neither offered nor removed (#318)") {
    App app;
    // One autosave per process: a second Penzene never writes or removes this one's.
    CHECK(QFileInfo(MainWindow::autosavePath()).fileName().contains(QString::number(QCoreApplication::applicationPid())));
    const QDir dir(QFileInfo(MainWindow::autosavePath()).path());
    dir.mkpath(".");
    const QString theirs = dir.filePath("autosave-running.penz");
    QLockFile running(theirs + ".lock");  // held, as by a live process
    REQUIRE(running.tryLock(0));
    Document d = *chem::fromSmiles("CCO");
    {
        QFile f(theirs);
        REQUIRE(f.open(QIODevice::WriteOnly));
        f.write(d.toJson());
    }
    MainWindow w;
    w.autosave();  // nothing unsaved here: removes only its own copy
    int offered = 0;
    QTimer::singleShot(0, &w, [&] {
        if (auto* box = qobject_cast<QMessageBox*>(QApplication::activeModalWidget())) ++offered, box->reject();
    });
    w.offerRecovery();
    CHECK(offered == 0);
    CHECK(QFile::exists(theirs));
    running.unlock();
    QFile::remove(theirs);
}

TEST_CASE("a crashed Penzene's autosave is offered even when this process has its pid (#361)") {
    App app;
    const QDir dir(QFileInfo(MainWindow::autosavePath()).path());
    dir.mkpath(".");
    const QString crashed = dir.filePath(QString("autosave-%1.penz").arg(QCoreApplication::applicationPid()));
    {
        QFile f(crashed);
        REQUIRE(f.open(QIODevice::WriteOnly));
        f.write(chem::fromSmiles("CCO")->toJson());
    }
    MainWindow w;
    int offered = 0;
    QTimer::singleShot(0, &w, [&] {
        if (auto* box = qobject_cast<QMessageBox*>(QApplication::activeModalWidget())) ++offered, box->reject();
    });
    w.offerRecovery();
    CHECK(offered == 1);
    QFile::remove(crashed);
}

TEST_CASE("orbitals stack over the molecule, or behind it once sent to the back (#204)") {
    App app;
    Document doc;
    doc.addAtom({0, 0});
    doc.addAtom({28.8, 0});
    doc.bonds = {{0, 1}};
    Arrow lobe{{5, 0}, {25, 0}, ArrowKind::Lobe};  // along the bond
    lobe.look = OrbitalLook::Shaded;
    lobe.color = QColor(220, 20, 20);
    doc.arrows = {lobe};
    auto colourOnBond = [](const Document& d) {
        QImage img(400, 200, QImage::Format_ARGB32);
        img.fill(Qt::white);
        QPainter painter(&img);
        painter.translate(100, 100);
        painter.scale(4, 4);
        paintDocument(painter, d);
        painter.end();
        return img.pixelColor(100 + 4 * 15, 100);
    };
    CHECK(colourOnBond(doc).red() > 180);  // new orbitals go on top, as in ChemDraw
    CHECK(edit::restack(doc, {0}, edit::Restack::Back) == std::vector<int>{0});
    CHECK(doc.arrows[0].behind);
    CHECK(colourOnBond(doc).red() < 80);  // the bond shows through
}

TEST_CASE("the Shapes flyout has every orbital in every look; one clicked on an atom is centred there (#204)") {
    App app;
    MainWindow w;
    w.resize(1000, 700);
    w.show();
    QFrame* shapes = nullptr;
    for (auto* f : w.findChildren<QFrame*>("toolFlyout"))
        if (f->findChild<QLabel*>("flyoutTitle")->text() == "SHAPES") shapes = f;
    REQUIRE(shapes);
    int orbitals = 0;
    QToolButton* shadedP = nullptr;
    for (auto* b : shapes->findChildren<QToolButton*>()) {
        const QString tip = b->defaultAction() ? b->defaultAction()->toolTip() : QString();
        orbitals += tip.contains("orbital,") || tip.startsWith("lobe,");
        if (tip.startsWith("p orbital, shaded")) shadedP = b;
    }
    CHECK(orbitals == 12);  // s, p, lobe and hybrid, each outline, shaded and gradient
    REQUIRE(shadedP);
    shadedP->click();
    auto* canvas = w.findChild<Canvas*>();
    Document d;
    d.addAtom({0, 0}, 7);
    canvas->setDocumentSilently(d);
    QTest::mouseClick(canvas->viewport(), Qt::LeftButton, {}, canvas->mapFromScene(QPointF(0, 0)));
    REQUIRE(canvas->document().arrows.size() == 1);
    const Arrow& p = canvas->document().arrows[0];
    CHECK(p.kind == ArrowKind::POrbital);
    CHECK(p.look == OrbitalLook::Shaded);
    CHECK(QLineF(p.from, QPointF(0, 0)).length() < 0.5);  // centred on the nitrogen
    // A click a little off the atom's centre is still a click: the same size as one dead centre (#431).
    const double full = QLineF(p.from, p.to).length();
    canvas->setDocumentSilently(d);
    QTest::mouseClick(canvas->viewport(), Qt::LeftButton, {}, canvas->mapFromScene(QPointF(2, 2)));
    REQUIRE(canvas->document().arrows.size() == 1);
    CHECK(std::abs(QLineF(canvas->document().arrows[0].from, canvas->document().arrows[0].to).length() - full) < 1e-6);
    if (auto shot = qgetenv("PENZENE_ORBITAL_TOOLS_SHOT"); !shot.isEmpty()) {
        for (auto* b : w.findChildren<QToolButton*>("railButton"))
            if (b->text() == "Shapes") b->click();
        QApplication::processEvents();
        w.grab().save(QString::fromUtf8(shot));
    }
}



TEST_CASE("isotopes: typed, drawn, and kept through .penz, SMILES, MOL and CDXML") {
    App app;
    Document d = *chem::fromSmiles("CO");
    REQUIRE(edit::applyLabel(d, 0, "13C"));
    CHECK(d.atoms[0].z == 6);
    CHECK(d.atoms[0].isotope == 13);
    CHECK(edit::atomText(d.atoms[0]) == "13C");
    REQUIRE(edit::applyLabel(d, 1, "18OH"));
    CHECK(d.atoms[1].z == 8);
    CHECK(d.atoms[1].isotope == 18);
    CHECK_FALSE(edit::applyLabel(d, 1, "1C"));  // fewer nucleons than protons
    Document plain = d;
    REQUIRE(edit::applyLabel(plain, 0, "C"));  // typing the plain element clears it
    CHECK(plain.atoms[0].isotope == 0);

    CHECK(chem::toSmiles(d) == "[13CH3][18OH]");
    for (const Document& back : {*Document::fromJson(d.toJson()), *chem::fromSmiles(chem::toSmiles(d)),
                                 *chem::fromMolBlock(chem::toMolBlock(d)), *chem::fromChemDraw(chem::toCdxml(d))}) {
        REQUIRE(back.atoms.size() == 2);
        CHECK(back.atoms[0].isotope + back.atoms[1].isotope == 31);
    }
    CHECK(chem::properties(d)->mw > chem::properties(*chem::fromSmiles("CO"))->mw + 2.9);  // ¹³C and ¹⁸O: about 3 heavier

    Document cd4 = *chem::fromSmiles("[2H]C([2H])([2H])[2H]");
    CHECK(edit::atomText(cd4.atoms[0]) == "D");
    CHECK(std::abs(chem::properties(cd4)->mw - 20.07) < 0.02);
    CHECK(chem::removeHydrogens(cd4).atoms.size() == 5);  // deuterium isn't an explicit H to tidy away
    Document odd = *chem::fromSmiles("C");
    REQUIRE(edit::applyLabel(odd, 0, "999C"));  // an isotope RDKit has no mass for
    CHECK(chem::properties(odd));

    Document ethane = *chem::fromSmiles("CC"), labelled = ethane;  // a chain carbon isn't labelled, ¹³C is
    labelled.atoms[0].isotope = 13;
    CHECK(documentBounds(labelled).width() > documentBounds(ethane).width());  // room for H₃¹³C
    CHECK(renderImage(labelled, {300}) != renderImage(ethane, {300}));
}

TEST_CASE("d makes deuterium; Delete takes an isotope off before the atom") {
    Fixture f;
    Document d;
    d.addAtom({0, 0}), d.addAtom({kBondLength, 0});
    d.bonds = {{0, 1}};
    d.atoms[1].isotope = 13;
    f.canvas.setDocumentSilently(d);
    f.hover({0, 0});
    f.key("d");
    CHECK(f.doc().atoms[0].z == 1);
    CHECK(f.doc().atoms[0].isotope == 2);
    f.hover({kBondLength, 0});
    QTest::keyClick(f.canvas.viewport(), Qt::Key_Delete);
    REQUIRE(f.doc().atoms.size() == 2);
    CHECK(f.doc().atoms[1].isotope == 0);
}

TEST_CASE("orbital tools are named alike: s orbital, p orbital, lobe, hybrid orbital (#376)") {
    App app;
    MainWindow w;
    QStringList names;
    for (auto* a : w.findChildren<QAction*>())
        if (a->toolTip().contains(", outline (click an atom")) names << a->toolTip().section(',', 0, 0);
    names.sort();
    CHECK(names == QStringList{"hybrid orbital", "lobe", "p orbital", "s orbital"});
}

#ifdef Q_OS_WIN
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>

TEST_CASE("Windows copy offers an Enhanced Metafile for Office (#394)") {
    App app;
    MainWindow w;
    auto* canvas = w.findChild<Canvas*>();
    canvas->setDocumentSilently(*chem::fromSmiles("CC(=O)Oc1ccccc1C(=O)O"));
    for (auto* a : w.findChildren<QAction*>())
        if (a->shortcut() == QKeySequence::Copy) a->trigger();
    const QMimeData* mime = QApplication::clipboard()->mimeData();
    const QByteArray emf = mime->data("image/x-emf");
    REQUIRE(emf.size() > qsizetype(sizeof(ENHMETAHEADER)));
    HENHMETAFILE handle = SetEnhMetaFileBits(UINT(emf.size()), reinterpret_cast<const BYTE*>(emf.constData()));
    REQUIRE(handle);
    ENHMETAHEADER header{};
    REQUIRE(GetEnhMetaFileHeader(handle, sizeof header, &header));
    CHECK(header.iType == EMR_HEADER);
    CHECK(header.dSignature == ENHMETA_SIGNATURE);
    CHECK(header.nBytes == DWORD(emf.size()));
    CHECK(header.rclBounds.right > header.rclBounds.left);  // something was drawn
    CHECK(header.rclBounds.bottom > header.rclBounds.top);
    CHECK(header.rclFrame.right > 0);
    DeleteEnhMetaFile(handle);

    // What Office asks the clipboard for: CF_ENHMETAFILE, handed over as a metafile handle.
    // The converter needs Qt's Windows plugin, so this part is skipped offscreen (as on CI).
    if (QGuiApplication::platformName() != "windows") return;
    EmfClipboard converter;
    const QList<FORMATETC> formats = converter.formatsForMime("image/x-emf", mime);
    REQUIRE(formats.size() == 1);
    CHECK(formats[0].cfFormat == CF_ENHMETAFILE);
    CHECK(converter.canConvertFromMime(formats[0], mime));
    STGMEDIUM medium{};
    REQUIRE(converter.convertFromMime(formats[0], mime, &medium));
    CHECK(medium.tymed == TYMED_ENHMF);
    CHECK(GetEnhMetaFileBits(medium.hEnhMetaFile, 0, nullptr) == UINT(emf.size()));
    DeleteEnhMetaFile(medium.hEnhMetaFile);
}

#include "OleServer.h"

// COM without the registry or Office: the object as Word or PowerPoint drives it, in process.
TEST_CASE("A drawing embedded in Word or PowerPoint opens, draws and saves (#229)") {
    App app;
    REQUIRE(SUCCEEDED(OleInitialize(nullptr)));  // offscreen, Qt hasn't
    const std::vector<Sheet> sheets{{"Page 1", *chem::fromSmiles("c1ccccc1O")}, {"Two", *chem::fromSmiles("CCO")}};
    IStorage* storage = ole::embedSource(sheets);  // what Copy offers as "Embed Source"
    REQUIRE(storage);
    CLSID clsid{};
    REQUIRE(SUCCEEDED(ReadClassStg(storage, &clsid)));
    CHECK(IsEqualCLSID(clsid, ole::kClsid));

    MainWindow w;
    IUnknown* object = ole::newObject(w);
    IPersistStorage* persist = nullptr;
    IDataObject* data = nullptr;
    IOleObject* oleObject = nullptr;
    REQUIRE(SUCCEEDED(object->QueryInterface(IID_IPersistStorage, reinterpret_cast<void**>(&persist))));
    REQUIRE(SUCCEEDED(object->QueryInterface(IID_IDataObject, reinterpret_cast<void**>(&data))));
    REQUIRE(SUCCEEDED(object->QueryInterface(IID_IOleObject, reinterpret_cast<void**>(&oleObject))));
    object->Release();
    REQUIRE(persist->Load(storage) == S_OK);
    REQUIRE(w.sheets().size() == 2);
    CHECK(w.sheets()[1].name == "Two");
    CHECK(w.sheets()[0].doc.atoms.size() == 7);
    CHECK(persist->IsDirty() == S_FALSE);
    CHECK(oleObject->SetHostNames(L"PowerPoint", L"Talk.pptx") == S_OK);
    CHECK(w.windowTitle().contains("Talk.pptx"));

    // The picture the container caches: an EMF, and an old-style metafile of the same size.
    FORMATETC emfFormat{CF_ENHMETAFILE, nullptr, DVASPECT_CONTENT, -1, TYMED_ENHMF};
    STGMEDIUM medium{};
    REQUIRE(data->GetData(&emfFormat, &medium) == S_OK);
    REQUIRE(medium.tymed == TYMED_ENHMF);
    ENHMETAHEADER header{};
    REQUIRE(GetEnhMetaFileHeader(medium.hEnhMetaFile, sizeof header, &header));
    const LONG width = header.rclFrame.right - header.rclFrame.left, height = header.rclFrame.bottom - header.rclFrame.top;
    CHECK(width > 0);
    CHECK(height > 0);
    ReleaseStgMedium(&medium);
    FORMATETC wmfFormat{CF_METAFILEPICT, nullptr, DVASPECT_CONTENT, -1, TYMED_MFPICT};
    REQUIRE(data->GetData(&wmfFormat, &medium) == S_OK);
    REQUIRE(medium.tymed == TYMED_MFPICT);
    const auto* pict = static_cast<const METAFILEPICT*>(GlobalLock(medium.hMetaFilePict));
    CHECK(pict->hMF);
    CHECK(pict->xExt == width);
    GlobalUnlock(medium.hMetaFilePict);
    ReleaseStgMedium(&medium);
    SIZEL size{};
    REQUIRE(oleObject->GetExtent(DVASPECT_CONTENT, &size) == S_OK);
    CHECK(size.cx == width);
    CHECK(size.cy == height);
    FORMATETC png{CLIPFORMAT(RegisterClipboardFormatW(L"PNG")), nullptr, DVASPECT_CONTENT, -1, TYMED_HGLOBAL};
    CHECK(data->QueryGetData(&png) == DV_E_FORMATETC);

    // Saved into the container's storage, as OleSave does, and loaded again: the same pages.
    CHECK(ole::embedSource({}) == nullptr);
    IStorage* saved = ole::embedSource({{"Page 1", Document{}}});
    REQUIRE(saved);
    REQUIRE(persist->Save(saved, FALSE) == S_OK);
    CHECK(persist->SaveCompleted(nullptr) == S_OK);
    MainWindow again;
    IUnknown* second = ole::newObject(again);
    IPersistStorage* persist2 = nullptr;
    REQUIRE(SUCCEEDED(second->QueryInterface(IID_IPersistStorage, reinterpret_cast<void**>(&persist2))));
    second->Release();
    REQUIRE(persist2->Load(saved) == S_OK);
    CHECK(again.sheets() == w.sheets());

    // Insert → Object: an empty drawing, one page.
    CHECK(persist2->InitNew(saved) == S_OK);
    CHECK(again.sheets().size() == 1);
    CHECK(again.sheets()[0].doc.empty());

    for (IUnknown* u : std::initializer_list<IUnknown*>{persist, data, oleObject, persist2, storage, saved}) u->Release();
    OleUninitialize();
}
#endif

TEST_CASE("the predicted shifts view is saved and drawn with its notice (#403)") {
    App app;
    Document doc = *chem::fromSmiles("CCO");
    const QRectF plain = outputBounds(doc);
    doc.showShifts = true;
    CHECK(Document::fromJson(doc.toJson())->showShifts);
    const QRectF shown = outputBounds(doc);
    CHECK(shown.bottom() > plain.bottom() + 12);  // the nmrshiftdb2 notice (three lines), under the drawing
    CHECK(shown.width() > plain.width());
}

TEST_CASE("Open says when a file is missing, not that it isn't a structure (#484)") {
    App app;
    MainWindow w;
    QString shown;
    QTimer::singleShot(0, &w, [&] {
        if (auto* box = qobject_cast<QMessageBox*>(QApplication::activeModalWidget())) {
            shown = box->text();
            box->button(QMessageBox::Ok)->click();
        }
    });
    const QString gone = QDir::tempPath() + "/penzene-gone-488.mol", unmounted = "/no/such/dir/gone.mol";
    QSettings().setValue("recentFiles", QStringList{gone, unmounted});
    CHECK_FALSE(w.openFile(unmounted));
    CHECK(shown.contains("doesn't exist"));
    // Open Recent forgets a file missing from a folder that is there, and keeps one whose folder isn't (#488).
    QTimer::singleShot(0, &w, [] { qobject_cast<QMessageBox*>(QApplication::activeModalWidget())->button(QMessageBox::Ok)->click(); });
    CHECK_FALSE(w.openFile(gone));
    CHECK(w.recentFiles() == QStringList{unmounted});
    QSettings().remove("recentFiles");
}

TEST_CASE("Structure → Variable Attachment: a bond on any one of the selected atoms (#455)") {
    Fixture f;
    f.canvas.setDocumentSilently(*chem::fromSmiles("c1ccccc1"));
    f.canvas.setSelection({0, 1});
    f.canvas.variableAttachment();
    const Document d = f.doc();
    REQUIRE(d.atoms.size() == 8);
    const Atom& point = d.atoms[6];
    CHECK(point.z == 0);
    CHECK(point.attachments == std::vector<int>{0, 1});
    const QLineF bond(point.pos, d.atoms[7].pos), edge(d.atoms[0].pos, d.atoms[1].pos);
    CHECK(bond.intersects(edge) == QLineF::BoundedIntersection);  // drawn across the ring's edge
    CHECK(chem::properties(d)->formula == "C7H8");  // toluene, wherever it sits
    CHECK(edit::moleculeOf(d, 7).size() == 8);  // moves with its ring

    CHECK(*Document::fromJson(d.toJson()) == d);
    const std::string mol = chem::toMolBlock(d, true);
    CHECK(mol.find("ENDPTS=(2 1 2)") != std::string::npos);
    for (const auto& back : {chem::fromMolBlock(mol), chem::fromChemDraw(chem::toCdxml(d))}) {
        REQUIRE(back);
        const auto it = std::find_if(back->atoms.begin(), back->atoms.end(), [](const Atom& a) { return a.z == 0; });
        REQUIRE(it != back->atoms.end());
        CHECK(it->attachments.size() == 2);
    }

    Document less = d;
    less.removeAtoms({0});  // keeps what's left of its atoms
    CHECK(less.atoms[5].attachments == std::vector<int>{0});
    Document more = d;
    more.append(d);
    CHECK(more.atoms[14].attachments == std::vector<int>{8, 9});
}

TEST_CASE("Arrange → Group: grouped objects select, arrange and save as one (#410)") {
    App app;
    MainWindow w;
    w.resize(900, 600);
    w.show();
    auto* canvas = w.findChild<Canvas*>();
    canvas->setTool(Canvas::Tool::Select);
    Document d = *chem::fromSmiles("CCO");  // ethanol, a benzene well to its right, and a note under ethanol
    d.append(*chem::fromSmiles("c1ccccc1"), {10 * kBondLength, 0});
    d.texts.push_back({{0, 3 * kBondLength}, "crude"});
    canvas->setDocumentSilently(d);
    QAction *group = nullptr, *ungroup = nullptr;
    for (auto* a : w.findChildren<QAction*>()) {
        if (a->text() == "&Group") group = a;
        if (a->text() == "&Ungroup") ungroup = a;
    }
    REQUIRE(group);
    REQUIRE(ungroup);
    CHECK(group->shortcut() == QKeySequence(Qt::CTRL | Qt::Key_G));
    canvas->setSelection({0}, {}, {0});
    group->trigger();
    const Document& doc = canvas->document();
    const int id = doc.texts[0].group;
    CHECK(id >= 0);
    for (int i : {0, 1, 2}) CHECK(doc.atoms[i].group == id);  // the whole molecule, not just the picked atom
    CHECK(doc.atoms[3].group == -1);

    auto click = [&](QPointF scene, Qt::KeyboardModifiers keys = {}) {
        QTest::mouseClick(canvas->viewport(), Qt::LeftButton, keys, canvas->mapFromScene(scene));
    };
    canvas->centerOn(doc.atoms[1].pos);
    {  // a selected group is one even tint: beside an atom, where dot and bond marks used to overlap, as mid-bond
        canvas->setSelection({0, 1, 2}, {}, {0});
        const QImage shot = canvas->viewport()->grab().toImage();
        const QPointF a = doc.atoms[0].pos, b = doc.atoms[1].pos, along = (b - a) / len(b - a), across(-along.y(), along.x());
        const QPoint nearAtom = canvas->mapFromScene(a + 2 * along + 1.5 * across), midBond = canvas->mapFromScene((a + b) / 2 + 1.5 * across);
        const qreal ratio = shot.devicePixelRatio();
        const QColor x = shot.pixelColor(nearAtom * ratio), y = shot.pixelColor(midBond * ratio);
        CHECK(std::abs(x.red() - y.red()) + std::abs(x.green() - y.green()) + std::abs(x.blue() - y.blue()) < 25);  // a grid line may cross
        CHECK(shot.pixelColor(midBond * ratio) != shot.pixelColor(canvas->mapFromScene(a + 20 * across) * ratio));  // tinted
    }
    canvas->setSelection({});
    click(doc.atoms[1].pos);
    CHECK(canvas->selection() == QSet<int>{0, 1, 2});
    CHECK(canvas->selectedTexts() == QSet<int>{0});  // one member takes the group
    click(doc.atoms[1].pos, Qt::ControlModifier);
    CHECK(canvas->selection() == QSet<int>{1});  // Ctrl/Cmd: into the group, just the one
    {  // and it shows as that part: the rest of the group isn't tinted
        const QImage shot = canvas->viewport()->grab().toImage();
        const QPointF a = doc.atoms[0].pos, b = doc.atoms[1].pos, along = (b - a) / len(b - a), across(-along.y(), along.x());
        const qreal ratio = shot.devicePixelRatio();
        CHECK(shot.pixelColor(canvas->mapFromScene((a + b) / 2 - 1.5 * across) * ratio) ==
              shot.pixelColor(canvas->mapFromScene(a + 20 * across) * ratio));
    }
    canvas->setSelection({0, 1, 2});
    const QPoint on = canvas->mapFromScene(doc.atoms[1].pos);
    QTest::mousePress(canvas->viewport(), Qt::LeftButton, Qt::ControlModifier | Qt::AltModifier, on);
    CHECK(canvas->selection() == QSet<int>{0, 1, 2});  // Ctrl+Alt: a snapped rotation of the whole selection
    QTest::mouseRelease(canvas->viewport(), Qt::LeftButton, Qt::ControlModifier | Qt::AltModifier, on);

    canvas->setSelection({});
    const QPointF offset = doc.texts[0].pos - doc.atoms[0].pos;
    canvas->arrangeScheme();
    CHECK(len(doc.texts[0].pos - doc.atoms[0].pos - offset) < 1e-6);  // one piece: the note kept its place

    const auto saved = Document::fromJson(doc.toJson());
    REQUIRE(saved);
    CHECK(saved->texts[0].group == id);
    CHECK(saved->atoms[0].group == id);
    Document twice = *saved;
    twice.append(*saved);
    CHECK(twice.atoms[twice.atoms.size() - 9].group != id);  // a second copy is a group of its own

    canvas->setSelection({0});
    ungroup->trigger();
    CHECK(doc.texts[0].group == -1);
    for (const Atom& a : doc.atoms) CHECK(a.group == -1);
    canvas->arrangeScheme();
    CHECK(len(doc.texts[0].pos - doc.atoms[0].pos - offset) > 1);  // apart again
    canvas->setSelection({0, 3});  // two molecules alone, no text to carry the change
    group->trigger();
    CHECK(doc.atoms[3].group >= 0);
    CHECK(doc.atoms[0].group == doc.atoms[3].group);
}

TEST_CASE("atoms added to a grouped molecule join its group: explicit hydrogens go with it (#590)") {
    App app;
    MainWindow w;
    w.resize(900, 600);
    w.show();
    auto* canvas = w.findChild<Canvas*>();
    canvas->setTool(Canvas::Tool::Select);
    Document d = *chem::fromSmiles("CCO");
    d.append(*chem::fromSmiles("C"), {10 * kBondLength, 0});  // methane, ungrouped
    d.texts.push_back({{0, 3 * kBondLength}, "crude"});
    canvas->setDocumentSilently(d);
    QAction *group = nullptr, *hydrogens = nullptr;
    for (auto* a : w.findChildren<QAction*>()) {
        if (a->text() == "&Group") group = a;
        if (a->text() == "Add Explicit &Hydrogens") hydrogens = a;
    }
    REQUIRE(group);
    REQUIRE(hydrogens);
    canvas->setSelection({0}, {}, {0});
    group->trigger();
    const int id = canvas->document().texts[0].group;
    canvas->setSelection({}, {}, {});
    hydrogens->trigger();
    const Document& doc = canvas->document();
    REQUIRE(doc.atoms.size() == 4 + 6 + 4);
    for (int a : edit::moleculeOf(doc, 0)) CHECK(doc.atoms[a].group == id);  // ethanol's six H too
    for (int a : edit::moleculeOf(doc, 3)) CHECK(doc.atoms[a].group == -1);  // methane's stay out
    canvas->centerOn(doc.atoms[1].pos);
    QTest::mouseClick(canvas->viewport(), Qt::LeftButton, {}, canvas->mapFromScene(doc.atoms[1].pos));
    canvas->deleteSelection();
    CHECK(canvas->document().atoms.size() == 5);  // the group went whole: only methane is left
}

TEST_CASE("drag anywhere in a selected molecule's box moves it; a group has one box (#410)") {
    App app;
    MainWindow w;
    w.resize(900, 600);
    w.show();
    auto* canvas = w.findChild<Canvas*>();
    canvas->setTool(Canvas::Tool::Select);
    Document d = *chem::fromSmiles("c1ccccc1");
    d.append(*chem::fromSmiles("c1ccccc1"), {5 * kBondLength, 0});
    canvas->setDocumentSilently(d);
    const Document& doc = canvas->document();
    QPointF left, right;
    for (int i = 0; i < 6; ++i) left += doc.atoms[i].pos / 6, right += doc.atoms[i + 6].pos / 6;
    canvas->centerOn((left + right) / 2);
    auto drag = [&](QPointF from, QPointF by) {
        const QPoint a = canvas->mapFromScene(from), b = canvas->mapFromScene(from + by);
        QTest::mousePress(canvas->viewport(), Qt::LeftButton, {}, a);
        QTest::mouseMove(canvas->viewport(), (a + b) / 2);
        QTest::mouseMove(canvas->viewport(), b);
        QTest::mouseRelease(canvas->viewport(), Qt::LeftButton, {}, b);
    };
    const QSet<int> both = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
    canvas->setSelection(both);
    const QPointF was = doc.atoms[0].pos;
    drag(left, {0, 2 * kBondLength});  // the ring's empty middle: inside its box, on no atom or bond
    CHECK(len(doc.atoms[0].pos - was - QPointF(0, 2 * kBondLength)) < 1);
    CHECK(len(doc.atoms[6].pos - (right + QPointF(0, 2 * kBondLength))) > 1e-6);  // both moved together
    const QPointF gap = (left + right) / 2 + QPointF(0, 2 * kBondLength);
    canvas->setSelection(both);
    const QPointF before = doc.atoms[0].pos;
    drag(gap, {0, 3 * kBondLength});  // between the two boxes: a marquee, nothing moves
    CHECK(doc.atoms[0].pos == before);

    canvas->setSelection(both);
    for (auto* a : w.findChildren<QAction*>())
        if (a->text() == "&Group") a->trigger();
    canvas->setSelection(both);
    drag(gap, {0, 3 * kBondLength});  // grouped: one box, so the gap is inside it
    CHECK(len(doc.atoms[0].pos - before - QPointF(0, 3 * kBondLength)) < 1);
}

TEST_CASE("the selection tint covers whole labels and the insides of whole rings (#410)") {
    App app;
    MainWindow w;
    w.resize(900, 600);
    w.show();
    auto* canvas = w.findChild<Canvas*>();
    Document d = *chem::fromSmiles("O");  // water: OH₂, its H₂ to the right
    d.append(*chem::fromSmiles("c1ccccc1"), {5 * kBondLength, 0});
    canvas->setDocumentSilently(d);
    const Document& doc = canvas->document();
    QPointF ring;
    for (int i = 1; i < 7; ++i) ring += doc.atoms[i].pos / 6;
    canvas->centerOn((doc.atoms[0].pos + ring) / 2);
    auto shot = [&] { return canvas->viewport()->grab().toImage(); };
    auto at = [&](const QImage& img, QPointF p) { return img.pixelColor(canvas->mapFromScene(p) * img.devicePixelRatio()); };
    const QPointF h = doc.atoms[0].pos + QPointF(12, 0), clear = doc.atoms[0].pos + QPointF(0, 40);
    const QImage before = shot();
    canvas->setSelection({0, 1, 2, 3, 4, 5, 6});
    const QImage after = shot();
    CHECK(at(after, h) != at(before, h));  // the H, well past the O
    CHECK(at(after, ring) != at(before, ring));  // the ring's middle
    CHECK(at(after, clear) == at(before, clear));
    canvas->setSelection({1, 2, 3, 4, 5});  // not the whole ring: its middle stays clear
    CHECK(at(shot(), ring) == at(before, ring));
}

TEST_CASE("a flyout opens beside the open ones, and the tool layout is remembered until reset (#408)") {
    App app;
    QSettings().remove("toolLayout");
    auto rail = [](MainWindow& w, const QString& g) {
        for (auto* b : w.findChildren<QToolButton*>("railButton"))
            if (b->text() == g) return b;
        return static_cast<QToolButton*>(nullptr);
    };
    auto flyout = [](MainWindow& w, const QString& g) {
        for (auto* f : w.findChildren<QFrame*>("toolFlyout"))
            if (f->findChild<QLabel*>("flyoutTitle")->text() == g.toUpper()) return f;
        return static_cast<QFrame*>(nullptr);
    };
    QPoint ringsAt;
    int ringsWidth = 0;
    {
        MainWindow w;
        w.resize(1400, 900);
        w.show();
        rail(w, "Bonds")->click();
        rail(w, "Rings")->click();
        auto *bonds = flyout(w, "Bonds"), *rings = flyout(w, "Rings");
        CHECK_FALSE(rings->geometry().intersects(bonds->geometry()));  // beside Bonds, not over it
        rings->move(500, 400);
        rail(w, "Atoms")->click();
        for (auto* other : {bonds, rings}) CHECK_FALSE(flyout(w, "Atoms")->geometry().intersects(other->geometry()));
        for (const char* g : {"Select", "Arrows", "Shapes"}) rail(w, g)->click();
        w.resize(700, 500);  // no room for them all: each still inside the window
        for (auto* f : w.findChildren<QFrame*>("toolFlyout"))
            if (f->isVisible()) CHECK(w.rect().contains(f->geometry()));
        w.resize(1400, 900);
        for (const char* g : {"Bonds", "Atoms", "Select", "Arrows", "Shapes"}) flyout(w, g)->hide();
        rings->move(500, 400);
        const QPoint corner(rings->width() - 3, rings->height() - 3);  // wider: more columns
        QTest::mousePress(rings, Qt::LeftButton, {}, corner);
        QTest::mouseMove(rings, corner + QPoint(200, 0));
        QTest::mouseRelease(rings, Qt::LeftButton, {}, corner + QPoint(200, 0));
        ringsAt = rings->pos(), ringsWidth = rings->width();
        CHECK(w.close());
    }
    {
        MainWindow w;  // the next launch
        w.resize(1400, 900);
        w.show();
        auto* rings = flyout(w, "Rings");
        CHECK(rings->isVisible());
        CHECK(rings->pos() == ringsAt);
        CHECK(rings->width() == ringsWidth);
        CHECK_FALSE(flyout(w, "Bonds")->isVisible());
        QAction* reset = nullptr;
        for (auto* a : w.findChildren<QAction*>())
            if (a->text() == "Reset Tool &Layout") reset = a;
        REQUIRE(reset);
        reset->trigger();
        CHECK_FALSE(rings->isVisible());
        CHECK_FALSE(QSettings().contains("toolLayout/2"));
        rail(w, "Rings")->click();
        CHECK(rings->width() < ringsWidth);  // its usual size again
    }
    QSettings().remove("toolLayout");
}

TEST_CASE("TLC plate: Rf beside each spot, kept up to date as spots and the front move (#506)") {
    CHECK(std::abs(edit::rf(70, 100, 0) - 0.3) < 1e-12);
    CHECK(std::isnan(edit::rf(70, 100, 100)));  // no distance run
    const auto& all = builtinTemplates();
    const auto tlc = std::find_if(all.begin(), all.end(), [](const Template& t) { return t.name == "TLC plate"; });
    REQUIRE(tlc != all.end());
    const Document plate = templateDocument(*tlc);
    auto rfs = [](const Document& d) {
        std::vector<QString> out;
        for (const auto& t : d.texts)
            if (t.rf) out.push_back(t.text);
        return out;
    };
    CHECK(rfs(plate) == std::vector<QString>{"0.62", "0.62", "0.31", "0.31"});

    Fixture f;
    f.canvas.insert(plate, "Insert");
    const Document placed = f.doc();
    CHECK(rfs(placed) == rfs(plate));
    std::vector<int> spots;
    int frontLine = -1;
    for (int i = 0; i < int(placed.arrows.size()); ++i) {
        if (placed.arrows[i].kind == ArrowKind::Ellipse) spots.push_back(i);
        if (placed.arrows[i].dashed) frontLine = i;
    }
    REQUIRE(spots.size() == 4);
    REQUIRE(frontLine >= 0);
    CHECK(std::all_of(placed.texts.begin(), placed.texts.end(), [&](const Text& t) { return t.group == placed.arrows[0].group; }));

    Document next = placed;  // the first spot moved up to the front
    const double front = next.arrows[frontLine].from.y();
    const QPointF up(0, front - QRectF(next.arrows[spots[0]].from, next.arrows[spots[0]].to).center().y());
    next.arrows[spots[0]].from += up, next.arrows[spots[0]].to += up;
    f.canvas.commit(next, "Move");
    CHECK(rfs(f.doc())[0] == "1.00");
    CHECK(f.doc().texts.size() == placed.texts.size());  // rewritten in place

    next = f.doc();  // the front moved down to halfway, so the other spots run further
    for (QPointF* p : {&next.arrows[frontLine].from, &next.arrows[frontLine].to})
        p->setY((front + next.arrows[1].from.y()) / 2);
    f.canvas.commit(next, "Move");
    CHECK(rfs(f.doc())[3] == "0.62");

    next = f.doc();  // without a front there's no Rf
    next.arrows.erase(next.arrows.begin() + frontLine);
    f.canvas.commit(next, "Delete");
    CHECK(rfs(f.doc()).empty());
    CHECK(f.doc().texts.size() == 3);  // the lane labels stay

    auto back = Document::fromJson(placed.toJson());
    REQUIRE(back);
    CHECK(*back == placed);
    auto cdx = chem::fromChemDraw(chem::toCdxml(placed));
    REQUIRE(cdx);
    CHECK(std::count_if(cdx->arrows.begin(), cdx->arrows.end(), [](const Arrow& a) { return a.filled; }) == 4);
}

// The keyboard tutorial's caffeine (#433): one click for a ring, then only hotkeys. With
// PENZENE_TUTORIAL_FRAMES set, each step is saved there as a frame for docs/tutorials/make-gif.py.
TEST_CASE("caffeine drawn with hotkeys, as the keyboard tutorial shows") {
    Fixture f;
    f.canvas.resize(420, 340);
    const QString dir = qEnvironmentVariable("PENZENE_TUTORIAL_FRAMES");
    QStringList captions;
    auto frame = [&](const QString& caption) {
        if (dir.isEmpty()) return;
        const QString name = QString("frame-%1.png").arg(captions.size(), 2, 10, QChar('0'));
        f.canvas.grab().save(dir + "/" + name);
        captions << name + "\t" + caption;
    };
    f.canvas.setTool(Canvas::Tool::Ring);
    f.canvas.setRing(6, false);
    f.click({0, 0});
    REQUIRE(f.doc().atoms.size() == 6);
    f.canvas.fitToDocument();
    f.canvas.zoomBy(0.8);  // room for the five-membered ring and the methyls
    f.canvas.setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    f.canvas.setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    frame("Click once with the Ring tool (cyclohexane). The rest is keys: point, then type.");
    auto mid = [&](int a, int b) { return (f.doc().atoms[a].pos + f.doc().atoms[b].pos) / 2; };
    auto step = [&](QPointF at, const QString& keys, const QString& what) {
        f.hover(at);
        frame(what);
        f.key(keys);
        f.hover(at);
        frame(what);
    };
    // The five-membered ring goes on the upper-right bond; walk the six-ring from its far end.
    int c4 = 0, c5 = 1;
    double best = -1e9;
    for (const auto& b : f.doc().bonds)
        if (const QPointF m = mid(b.a, b.b); m.x() - m.y() > best) best = m.x() - m.y(), c4 = b.a, c5 = b.b;
    if (f.doc().atoms[c4].pos.y() > f.doc().atoms[c5].pos.y()) std::swap(c4, c5);  // c5 the lower
    auto next = [&](int from, int prev) {
        for (int n : f.doc().neighbors(from))
            if (n != prev) return n;
        return -1;
    };
    const int c6 = next(c5, c4), n1 = next(c6, c5), c2 = next(n1, c6), n3 = next(c2, n1);
    step(mid(c4, c5), "5", "Bond: 5 fuses a five-membered ring");
    const int n7 = int(f.doc().atoms.size()) - 1, n9 = n7 - 2, c8 = n7 - 1;  // the new ring's atoms, from c4's side
    REQUIRE(f.doc().bondBetween(n9, c4) >= 0);
    REQUIRE(f.doc().bondBetween(n7, c5) >= 0);
    step(f.doc().atoms[c6].pos, "2", "Atom: 2 makes a C=O");
    step(f.doc().atoms[c2].pos, "2", "Atom: 2 again for the second C=O");
    step(f.doc().atoms[n1].pos, "n", "Atom: n makes it nitrogen");
    step(f.doc().atoms[n1].pos, "1", "Atom: 1 sprouts a methyl");
    step(f.doc().atoms[n3].pos, "n", "Atom: n");
    step(f.doc().atoms[n3].pos, "1", "Atom: 1");
    step(f.doc().atoms[n7].pos, "n", "Atom: n");
    step(f.doc().atoms[n7].pos, "1", "Atom: 1");
    step(f.doc().atoms[n9].pos, "n", "Atom: n");
    step(mid(c4, c5), "2", "Bond: 2 makes it double");
    step(mid(c8, n9), "2", "Bond: 2");
    QTest::keyClick(f.canvas.viewport(), Qt::Key_Escape);  // no hotspot on the finished drawing
    frame("Caffeine: one click and 12 keys");
    CHECK(chem::toSmiles(f.doc()) == chem::toSmiles(*chem::fromSmiles("Cn1c(=O)c2c(ncn2C)n(C)c1=O")));
    if (!dir.isEmpty()) {
        QFile list(dir + "/frames.txt");
        REQUIRE(list.open(QIODevice::WriteOnly | QIODevice::Text));
        list.write((captions.join('\n') + '\n').toUtf8());
    }
}

TEST_CASE("Open says which records of a library it couldn't read (#584)") {
    QStringList unreadable;
    const auto doc = chem::readFile(QString(PENZENE_TEST_DATA) + "/broken.smi", &unreadable);
    REQUIRE(doc);
    CHECK(unreadable == QStringList{"broken"});

    MainWindow w;
    QString shown;
    QTimer::singleShot(0, &w, [&] {
        if (auto* box = qobject_cast<QMessageBox*>(QApplication::activeModalWidget())) shown = box->text(), box->accept();
    });
    REQUIRE(w.openFile(QString(PENZENE_TEST_DATA) + "/broken.smi"));
    CHECK(shown.contains("broken"));
    CHECK(chem::toSmiles(w.findChild<Canvas*>()->document()) == "CCO");
}

TEST_CASE("metal ligand labels and dative hotkeys work on the canvas (#616)") {
    Fixture f;
    Document d;
    d.addAtom({0, 0}, 78), d.addAtom({kBondLength, 0}, 7);
    d.bonds.push_back({0, 1});
    REQUIRE(edit::applyLabel(d, 1, "NH3"));
    f.canvas.setDocumentSilently(d);
    CHECK(chem::atomInfo(f.doc())[1].hydrogens == 3);
    const auto preview = qEnvironmentVariable("PENZENE_METAL_PREVIEW");
    if (!preview.isEmpty()) {
        f.canvas.zoomBy(4);
        f.canvas.centerOn({kBondLength / 2, 0});
        f.canvas.grab().save(preview + "/before-arrow.png");
    }
    f.hover({kBondLength / 2, 0});
    f.key("k");
    CHECK(f.doc().bonds[0].stereo == BondStereo::Dative);
    CHECK(f.doc().bonds[0].a == 1);
    if (!preview.isEmpty()) f.canvas.grab().save(preview + "/after-arrow.png");
    f.undo.undo();
    CHECK(f.doc().bonds[0].stereo == BondStereo::None);
    CHECK(chem::atomInfo(f.doc())[1].hydrogens == 3);
}


TEST_CASE("format toolbar: mixed selections, undo and enabled controls (#621)") {
    MainWindow w;
    w.resize(1100, 750);
    w.show();
    auto* canvas = w.findChild<Canvas*>();
    auto* undo = w.findChild<QUndoStack*>();
    auto* colour = w.findChild<QToolButton*>("formatColour");
    auto* fill = w.findChild<QToolButton*>("formatFill");
    auto* font = w.findChild<QComboBox*>("formatFont");
    auto* size = w.findChild<QDoubleSpinBox*>("formatSize");
    auto* bold = w.findChild<QCheckBox*>("formatBold");
    auto* italic = w.findChild<QCheckBox*>("formatItalic");
    REQUIRE(canvas); REQUIRE(undo); REQUIRE(colour); REQUIRE(fill); REQUIRE(font); REQUIRE(size); REQUIRE(bold); REQUIRE(italic);
    CHECK(!canvas->colour().isValid());
    CHECK(!fill->isEnabled());
    CHECK(!font->isEnabled());
    Document d;
    d.atoms = {{{0, 0}}, {{kBondLength, 0}, 8}};
    d.bonds = {{0, 1}};
    d.arrows = {{{-20, -20}, {50, 20}, ArrowKind::Box}, {{60, -20}, {100, 20}, ArrowKind::Ellipse}};
    d.arrows[1].filled = true, d.arrows[1].fill = Qt::red;
    d.texts = {{{0, 40}, "heat"}, {{0, 65}, "H2O"}};
    d.texts[1].bold = true, d.texts[1].italic = true, d.texts[1].scale = 2, d.texts[1].font = "Times New Roman";
    canvas->setDocumentSilently(d);
    canvas->setTool(Canvas::Tool::Select);
    canvas->setSelection({0});
    CHECK(!fill->isEnabled()); CHECK(!font->isEnabled()); CHECK(!bold->isEnabled());
    canvas->setSelection({}, {0, 1});
    CHECK(fill->isEnabled()); CHECK(!font->isEnabled());
    CHECK(fill->text().contains("Mixed"));
    canvas->setSelection({}, {}, {0, 1});
    CHECK(!fill->isEnabled()); CHECK(font->isEnabled()); CHECK(size->isEnabled()); CHECK(bold->isEnabled());
    CHECK(font->currentIndex() == -1);
    CHECK(size->value() == 0);
    CHECK(bold->checkState() == Qt::PartiallyChecked);
    CHECK(italic->checkState() == Qt::PartiallyChecked);
    const int before = undo->count();
    bold->click();
    CHECK(undo->count() == before + 1);
    CHECK(canvas->document().texts[0].bold); CHECK(canvas->document().texts[1].bold);
    undo->undo();
    CHECK(canvas->document() == d);
    CHECK(bold->checkState() == Qt::PartiallyChecked);
    canvas->setSelection({0, 1}, {0}, {0, 1});
    const int colourBefore = undo->index();
    auto* swatch = colour->menu()->findChild<QWidgetAction*>()->defaultWidget()->findChildren<QToolButton*>()[1];
    swatch->click();
    CHECK(undo->index() == colourBefore + 1);
    const QColor blue(0x30, 0x50, 0xF8);
    CHECK(canvas->colour() == blue);
    CHECK(canvas->document().atoms[0].color == blue);
    CHECK(canvas->document().bonds[0].color == blue);
    CHECK(canvas->document().arrows[0].color == blue);
    CHECK(canvas->document().texts[1].color == blue);
    colour->menu()->actions()[0]->trigger();  // Automatic
    CHECK(!canvas->colour().isValid());
    CHECK(!canvas->document().texts[0].color.isValid());
    canvas->setSelection({}, {0});
    fill->menu()->findChild<QWidgetAction*>()->defaultWidget()->findChildren<QToolButton*>()[0]->click();
    CHECK(canvas->document().arrows[0].filled);
    CHECK(canvas->document().arrows[0].fill == QColor(207, 227, 255));
    fill->menu()->actions()[0]->trigger();  // None
    CHECK(!canvas->document().arrows[0].filled);
    canvas->setSelection({}, {}, {0, 1});
    size->setValue(12);
    const int fontBefore = undo->index();
    font->lineEdit()->selectAll();
    QTest::keyClicks(font->lineEdit(), "Times New Roman");
    QTest::keyClick(font->lineEdit(), Qt::Key_Return);
    CHECK(undo->index() == fontBefore + 1);
    CHECK(canvas->document().texts[0].font == "Times New Roman");
    CHECK(canvas->document().texts[1].font == "Times New Roman");
    canvas->setTextBold(true);
    canvas->setTextItalic(true);
    CHECK(size->value() == 12);
    CHECK(font->currentText() == "Times New Roman");
    CHECK(bold->isChecked()); CHECK(italic->isChecked());
    const QString previews = qEnvironmentVariable("PENZENE_FORMAT_PREVIEW");
    if (!previews.isEmpty()) {
        canvas->setSelection({}, {0}); canvas->setFillColor(QColor(207, 227, 255)); canvas->fillSelection();
        canvas->setSelection({}, {}, {0, 1}); canvas->fitToDocument();
        QApplication::processEvents();
        CHECK(w.grab().save(previews + "/formatted.png"));
        QFile sample(previews + "/formatted.penz");
        REQUIRE(sample.open(QIODevice::WriteOnly)); sample.write(canvas->document().toJson());
        canvas->setDocumentSilently(d); canvas->setSelection({}, {}, {0, 1});
        QApplication::processEvents(); CHECK(w.grab().save(previews + "/mixed.png"));
    }
    canvas->setSelection({});
    canvas->setTool(Canvas::Tool::Text);
    CHECK(font->isEnabled()); CHECK(!fill->isEnabled());
    canvas->setArrow(ArrowKind::Box, false);
    canvas->setTool(Canvas::Tool::Arrow);
    CHECK(fill->isEnabled()); CHECK(!font->isEnabled());
    Document old = d;
    old.arrows[0].filled = true;
    canvas->setDocumentSilently(old); canvas->setSelection({}, {0});
    CHECK(fill->text().contains("Automatic"));
    if (qEnvironmentVariableIsSet("PENZENE_FORMAT_PREVIEW_LIVE")) {
        canvas->setTool(Canvas::Tool::Select);
        canvas->setSelection({}, {}, {0, 1}); canvas->fitToDocument();
        QEventLoop preview;
        QObject::connect(qApp, &QApplication::lastWindowClosed, &preview, &QEventLoop::quit);
        preview.exec();
    }
}

TEST_CASE("colour and fill menus draw their swatches (#674)") {
    MainWindow w;
    w.resize(1100, 750);
    w.show();
    auto* colour = w.findChild<QToolButton*>("formatColour");
    auto* fill = w.findChild<QToolButton*>("formatFill");
    REQUIRE(colour);
    REQUIRE(fill);
    colour->menu()->show();
    QApplication::processEvents();
    auto* holder = colour->menu()->findChild<QWidgetAction*>()->defaultWidget();
    REQUIRE(holder);
    const QImage image = holder->grab().toImage();
    const QList<QColor> expected{Qt::black, QColor(0x30, 0x50, 0xF8), QColor(0xFF, 0x0D, 0x0D), QColor(0, 0x90, 0x30)};
    const auto buttons = holder->findChildren<QToolButton*>();
    REQUIRE(buttons.size() == expected.size() + 1);  // swatches and Custom…
    for (int i = 0; i < expected.size(); ++i)
        CHECK(image.pixelColor(buttons[i]->geometry().center()) == expected[i]);
    colour->menu()->close();
    fill->menu()->show();
    QApplication::processEvents();
    holder = fill->menu()->findChild<QWidgetAction*>()->defaultWidget();
    REQUIRE(holder);
    const QImage fillImage = holder->grab().toImage();
    const auto fillButtons = holder->findChildren<QToolButton*>();
    REQUIRE(fillButtons.size() == expected.size() + 1);
    CHECK(fillImage.pixelColor(fillButtons[0]->geometry().center()) == QColor(207, 227, 255));
    fill->menu()->close();
}

TEST_CASE("format defaults reach new shapes, text, bonds and ring fills (#621)") {
    Fixture f;
    f.canvas.setColour(Qt::blue);
    f.canvas.setFillColor(Qt::red);
    f.canvas.setTextFont("Times New Roman");
    f.canvas.setTextSize(12);
    f.canvas.setTextBold(true);
    f.canvas.setTextItalic(true);
    CHECK(f.undo.count() == 0);
    f.canvas.setArrow(ArrowKind::Box, false);
    f.canvas.setTool(Canvas::Tool::Arrow);
    f.drag({0, 0}, {40, 25});
    REQUIRE(f.doc().arrows.size() == 1);
    CHECK(f.doc().arrows[0].color == Qt::blue);
    CHECK(f.doc().arrows[0].filled);
    CHECK(f.doc().arrows[0].fill == Qt::red);
    QTimer::singleShot(0, [] { auto* dialog = qobject_cast<QInputDialog*>(QApplication::activeModalWidget()); REQUIRE(dialog); dialog->setTextValue("H2O"); dialog->accept(); });
    f.canvas.editText(-1, {0, 50});
    REQUIRE(f.doc().texts.size() == 1);
    CHECK(f.doc().texts[0].font == "Times New Roman");
    CHECK(f.doc().texts[0].bold); CHECK(f.doc().texts[0].italic);
    CHECK(f.doc().texts[0].color == Qt::blue);
    CHECK(std::abs(f.doc().texts[0].scale * documentStyle(f.doc()).fontSize * exportScale(f.doc()) - 12) < 1e-6);
    f.canvas.setTool(Canvas::Tool::Bond);
    f.click({70, 50});
    REQUIRE(f.doc().bonds.size() == 1);
    CHECK(f.doc().bonds[0].color == Qt::blue);
    CHECK(f.doc().atoms[0].color == Qt::blue);
    f.canvas.setTool(Canvas::Tool::Ring);
    f.click({70, 90});
    const auto rings = chem::rings(f.doc()); REQUIRE(rings.size() == 1);
    f.canvas.setSelection(QSet<int>(rings[0].begin(), rings[0].end()));
    const int index = f.undo.index();
    f.canvas.fillSelection();
    CHECK(f.undo.index() == index + 1);
    REQUIRE(f.doc().fills.size() == 1); CHECK(f.doc().fills[0].color == Qt::red);
    f.canvas.setFillColor({}); f.canvas.fillSelection(); CHECK(f.doc().fills.empty());
    f.undo.undo(); REQUIRE(f.doc().fills.size() == 1);
}

TEST_CASE("measurement fields keep precision and localized keyboard entry without trailing zeroes") {
    App app;
    for (const QLocale& locale : {QLocale::c(), QLocale(QLocale::German), QLocale(QLocale::Arabic)}) {
        NumberBox box;
        box.setLocale(locale);
        box.setDecimals(6);
        box.setRange(0, 1000);
        box.setValue(20);
        CHECK(box.cleanText() == locale.toString(20));
        box.setValue(0.595276);
        CHECK(box.value() == Catch::Approx(0.595276));
        CHECK(box.cleanText() == locale.toString(0.595276, 'f', 6));
        box.show();
        box.setFocus();
        box.selectAll();
        if (locale.language() == QLocale::Arabic)
            box.findChild<QLineEdit*>()->setText(locale.toString(1.125, 'f', 3));
        else QTest::keyClicks(&box, locale.language() == QLocale::German ? "1,125" : "1.125");
        QTest::keyClick(&box, Qt::Key_Return);
        CHECK(box.value() == Catch::Approx(1.125));
        box.stepUp();
        CHECK(box.value() == Catch::Approx(2.125));
    }
}

TEST_CASE("embedded custom styles can be edited without the original style library") {
    App app;
    MainWindow window;
    Document doc = *chem::fromSmiles("CCO");
    auto style = styleSettings(drawingStyle("RSC"));
    style.name = "Portable embedded style regression";
    style.bondLength = 18.75;
    style.captionSize = 6.5;
    doc.styleSettings = style;
    QTimer::singleShot(0, [&] {
        auto* dialog = qobject_cast<QDialog*>(QApplication::activeModalWidget());
        REQUIRE(dialog);
        CHECK(dialog->findChild<QLineEdit*>("styleName")->text() == style.name);
        CHECK(dialog->findChild<QDoubleSpinBox*>("bondLength")->value() == Catch::Approx(18.75));
        CHECK(dialog->findChild<QDoubleSpinBox*>("captionSize")->value() == Catch::Approx(6.5));
        CHECK(dialog->findChild<QPushButton*>("saveStyle")->isEnabled());
        CHECK_FALSE(dialog->findChild<QPushButton*>("deleteStyle")->isEnabled());
        dialog->reject();
    });
    editDrawingStyles(&window, doc);
    CHECK(doc.styleSettings == style);
}
