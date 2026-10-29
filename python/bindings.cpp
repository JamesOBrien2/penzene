// Python bindings over penzene_core: the same engine and renderer as the app.
#include "Chem.h"
#include "Edit.h"
#include "Render.h"
#include "Styles.h"

#include <QFileInfo>
#include <QGuiApplication>
#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace {

// Fonts and painting need a QGuiApplication. Reuse the host's if there is one
// (e.g. PySide6); otherwise make an offscreen one that lives as long as the process.
void ensureApp() {
    if (QCoreApplication::instance()) return;
    if (qEnvironmentVariableIsEmpty("QT_QPA_PLATFORM")) qputenv("QT_QPA_PLATFORM", "offscreen");
    static int argc = 1;
    static char name[] = "penzene";
    static char* argv[] = {name, nullptr};
    new QGuiApplication(argc, argv);  // never freed; the process owns it until exit
}

QString qs(const std::string& s) { return QString::fromStdString(s); }

Document fromSmiles(const std::string& smiles) {
    auto doc = chem::fromSmiles(smiles);
    if (!doc) throw nb::value_error(("not a valid SMILES: " + smiles).c_str());
    return *doc;
}

// A missing file is a FileNotFoundError, and still the ValueError it used to raise (the stability policy).
[[noreturn]] void missingFile(const std::string& path) {
    static PyObject* type = [] {
        PyObject* bases = PyTuple_Pack(2, PyExc_FileNotFoundError, PyExc_ValueError);
        PyObject* t = PyErr_NewException("penzene.MissingFile", bases, nullptr);
        Py_DECREF(bases);
        return t;
    }();
    PyErr_SetString(type, ("no such file: " + path).c_str());
    throw nb::python_error();
}

Document readPath(const std::string& path) {
    if (!QFileInfo::exists(qs(path))) missingFile(path);
    QStringList unreadable;
    auto doc = chem::readFile(qs(path), &unreadable);
    if (!doc) throw nb::value_error(("cannot read " + path).c_str());
    if (!unreadable.isEmpty()) {  // a partial library is never a plain success (#584)
        const std::string message = "left out unreadable records of " + path + ": " + unreadable.join(", ").toStdString();
        if (PyErr_WarnEx(PyExc_UserWarning, message.c_str(), 1) < 0) throw nb::python_error();
    }
    return *doc;
}

void save(const Document& doc, const std::string& path) {
    const QString ext = QFileInfo(qs(path)).suffix().toLower();
    QByteArray data;
    if (ext == "penz") data = doc.toJson();
    else if (ext == "cdxml") data = chem::toCdxml(doc);
    else if (ext == "cdx") data = chem::toCdx(doc);
    else if (ext == "mol" || ext.isEmpty()) data = QByteArray::fromStdString(chem::toMolBlock(doc));
    else  // an image or library path would get MOL text (#327)
        throw nb::value_error(("save writes .penz, .mol, .cdxml or .cdx; use export() for ." + ext.toStdString()).c_str());
    // Written aside and swapped in, so a failed save never truncates the file already there.
    if (data.isEmpty() || !writeWhole(qs(path), data)) throw nb::value_error(("cannot write " + path).c_str());
}

std::string symbol(const Atom& a) { return a.label.isEmpty() ? chem::symbol(a.z) : a.label.toStdString(); }

Atom& atomAt(Document& d, int atom) {
    if (atom < 0 || atom >= int(d.atoms.size())) throw nb::index_error("no such atom");
    return d.atoms[atom];
}

void inRange(int value, int lo, int hi, const char* what) {
    if (value < lo || value > hi)
        throw nb::value_error((std::string(what) + " must be " + std::to_string(lo) + " to " + std::to_string(hi)).c_str());
}

[[noreturn]] void readOnly(const char* type, const std::string& field, const char* hint) {
    throw nb::attribute_error((std::string(type) + "." + field + " is read-only; " + hint).c_str());
}

void requireFinite(double value) {
    if (!std::isfinite(value)) throw nb::value_error("coordinates and angles must be finite");
}

ExportOptions pythonExportOptions(const Document& d, double margin, const std::optional<std::string>& background,
                                  double scale, double dpi = 300) {
    for (double value : {margin, scale, dpi}) requireFinite(value);
    if (margin < 0 || scale <= 0 || dpi <= 0) throw nb::value_error("margin must be non-negative; scale and dpi must be positive");
    QColor color = Qt::transparent;
    if (background) {
        color = QColor(qs(*background));
        if (!color.isValid()) throw nb::value_error("background must be None or a colour name or hex string");
    }
    const double s = exportScale(d) * scale;
    const QRectF page = pageRect(d);
    const QSizeF size = page.isEmpty() ? outputBounds(d).size() * s + QSizeF(2 * margin, 2 * margin) : page.size() * s;
    const double pixels = std::max(size.width(), size.height()) * dpi / 72;
    if (!std::isfinite(s) || s <= 0 || !std::isfinite(pixels) || pixels >= std::numeric_limits<int>::max()
        || dpi / 0.0254 >= std::numeric_limits<int>::max())
        throw nb::value_error("export dimensions or resolution are too large");
    return {dpi, color, scale, margin};
}

std::vector<int> selectedAtoms(Document& d, std::vector<int> atoms) {
    if (atoms.empty())
        for (int i = 0; i < int(d.atoms.size()); ++i) atoms.push_back(i);
    for (int i : atoms) atomAt(d, i);
    std::sort(atoms.begin(), atoms.end());
    atoms.erase(std::unique(atoms.begin(), atoms.end()), atoms.end());
    return atoms;
}

void transformAtoms(Document& d, const std::vector<int>& atoms, const QTransform& transform) {
    Document next = d;
    edit::transformAtoms(next, atoms, transform);
    for (int i : atoms) requireFinite(next.atoms[i].pos.x()), requireFinite(next.atoms[i].pos.y());
    edit::followAnchors(d, next);
    edit::followNumbers(d, next);
    d = std::move(next);
}

}  // namespace

NB_MODULE(_penzene, m) {
    ensureApp();
    m.doc() = "Penzene: 2D chemical structures, drawn by the same engine as the app.";
    m.attr("__version__") = PENZENE_VERSION;

    nb::class_<Atom>(m, "Atom", "An atom of a Document (read-only; edit through Document: set_label, remove_atoms). (since 0.4)")
        .def("__setattr__", [](const Atom&, const std::string& field, nb::handle) {
            const std::map<std::string, const char*> hints = {
                {"x", "use Document.set_position(atom, x, y)"},
                {"y", "use Document.set_position(atom, x, y)"},
                {"symbol", "use Document.set_label(atom, label)"},
                {"label", "use Document.set_label(atom, label)"},
                {"charge", "use Document.set_charge(atom, charge)"},
                {"radicals", "use Document.set_radicals(atom, radicals)"},
                {"partial", "use Document.set_partial(atom, partial)"}};
            auto hint = hints.find(field);
            readOnly("Atom", field, hint == hints.end() ? "edit through Document; atoms are snapshots" : hint->second);
        }, "name"_a, "value"_a.none())
        .def_prop_ro("symbol", &symbol, "Element symbol, or the abbreviation (Boc, OMe…) if it has one (since 0.4)")
        .def_prop_ro("x", [](const Atom& a) { return a.pos.x(); }, "Points, x right (since 0.4)")
        .def_prop_ro("y", [](const Atom& a) { return a.pos.y(); }, "Points, y down (since 0.4)")
        .def_ro("charge", &Atom::charge, "Formal charge (since 0.4)")
        .def_ro("isotope", &Atom::isotope, "Mass number (13 for carbon-13, 2 for deuterium); 0 for natural abundance (since 1.3)")
        .def_prop_ro("label", [](const Atom& a) { return a.label.toStdString(); },
                     "The label drawn instead of the element (Boc, OMe, R, or its own symbol when shown); \"\" for none (since 2.0)")
        .def_ro("radicals", &Atom::radicals, "Radical electrons, 0 to 2 (since 2.0)")
        .def_ro("partial", &Atom::partial, "Partial charge mark: 1 for δ+, -1 for δ−, 0 for none (since 2.0)")
        .def("__repr__", [](const Atom& a) { return "<Atom " + symbol(a) + ">"; });

    nb::class_<Bond>(m, "Bond", "A bond of a Document (read-only). (since 0.4)")
        .def("__setattr__", [](const Bond&, const std::string& field, nb::handle) {
            const char* hint = field == "order" ? "use Document.set_bond_order(bond, order)"
                             : field == "style" ? "use Document.set_bond_style(bond, style)"
                                                : "use Document.add_bond(a, b) to connect atoms";
            readOnly("Bond", field, hint);
        }, "name"_a, "value"_a.none())
        .def_ro("a", &Bond::a, "Index of the first atom (a wedge starts here) (since 0.4)")
        .def_ro("b", &Bond::b, "Index of the second atom (since 0.4)")
        .def_ro("order", &Bond::order, "1, 2, 3 or 4 (since 0.4)")
        .def_prop_ro("style", [](const Bond& b) { return std::string(kBondStereoNames[int(b.stereo)]); },
                     "How it's drawn, as set_bond_style takes it: 'none', 'wedge', 'partial'… (since 2.0)")
        .def("__repr__", [](const Bond& b) {
            return "<Bond " + std::to_string(b.a) + "-" + std::to_string(b.b) + " order " + std::to_string(b.order) + ">";
        });

    nb::class_<Document>(m, "Document", "A drawing: molecules, arrows and text, as the app holds it. (since 0.4)")
        .def(nb::init<>(), "An empty drawing. (since 0.4)")
        .def_ro("atoms", &Document::atoms, "The atoms, in order (their index is what add_bond and hotkeys take) (since 0.4)")
        .def_ro("bonds", &Document::bonds, "The bonds (since 0.4)")
        .def_prop_rw(
            "style", [](const Document& d) { return documentStyle(d).name.toStdString(); },
            [](Document& d, const std::string& name) {
                const DrawingStyle& s = drawingStyle(qs(name));
                if (s.name != qs(name)) {
                    QStringList names;
                    for (const auto& known : drawingStyles()) names << "'" + known.name + "'";
                    throw nb::value_error(("unknown drawing style: " + name + " (use " + names.join(", ").toStdString() + ")").c_str());
                }
                applyDrawingStyle(d, qs(name));
            },
            "Drawing style name, from drawing_styles(); includes saved user styles (since 0.4)")
        .def(
            "add_atom",
            [](Document& d, const std::string& label, double x, double y, bool generic) {
                int i = d.addAtom({x, y});
                if (!edit::applyLabel(d, i, qs(label), generic)) {
                    d.atoms.pop_back();
                    throw nb::value_error(("not an element, abbreviation or SMILES: " + label).c_str());
                }
                return i;
            },
            "label"_a = "C", "x"_a = 0.0, "y"_a = 0.0, "generic"_a = false,
            "Add an atom: an element, an abbreviation (OMe, Boc…) or a SMILES fragment. Returns its index. "
            "With generic=True (since 2.0), any other text is a generic atom as the app takes a typed label "
            "(R, R1, X, Ar). (since 0.4)")
        .def(
            "add_bond",
            [](Document& d, int a, int b, int order) {
                const int n = int(d.atoms.size());
                if (a < 0 || b < 0 || a >= n || b >= n || a == b) throw nb::index_error("no such atoms");
                if (order < 1 || order > 4) throw nb::value_error("order must be 1, 2, 3 or 4");
                edit::link(d, a, b, order);
            },
            "a"_a, "b"_a, "order"_a = 1, "Bond atoms a and b (order 1, 2, 3 or 4). (since 0.4)")
        .def(
            "bond_between",
            [](Document& d, int a, int b) -> std::optional<int> {
                atomAt(d, a), atomAt(d, b);
                const int bond = d.bondBetween(a, b);
                return bond < 0 ? std::nullopt : std::optional<int>(bond);
            },
            "a"_a, "b"_a,
            "Index of the bond joining atoms a and b in either order, or None if unbonded (including a == b). "
            "Invalid atom indices raise IndexError. Look up again after removing atoms or bonds, which reindexes them. (since 2.1)")
        .def(
            "set_label",
            [](Document& d, int atom, const std::string& label, bool generic) {
                if (atom < 0 || atom >= int(d.atoms.size())) throw nb::index_error("no such atom");
                Document next = d;
                if (!edit::applyLabel(next, atom, qs(label), generic))
                    throw nb::value_error(("not an element, abbreviation or SMILES: " + label).c_str());
                d = next;
            },
            "atom"_a, "label"_a, "generic"_a = false,
            "Relabel an atom, as add_atom takes a label: an element (\"C\" undoes a label; show_symbol shows it), an "
            "abbreviation or a SMILES fragment (drawn out from it, its new atoms added at the end), or with "
            "generic=True any other text (R, X). (since 2.0)")
        .def(
            "set_bond_order",
            [](Document& d, int bond, int order) {
                if (bond < 0 || bond >= int(d.bonds.size())) throw nb::index_error("no such bond");
                if (order < 1 || order > 4) throw nb::value_error("order must be 1, 2, 3 or 4");
                Bond& b = d.bonds[bond];
                b.order = order;
                if (order != 1) b.stereo = BondStereo::None;  // a wedge or hash is a single bond's
            },
            "bond"_a, "order"_a, "Make a bond single, double, triple or quadruple (a wedge or hash goes with a single bond). (since 2.0)")
        .def(
            "set_charge", [](Document& d, int atom, int charge) { inRange(charge, -8, 8, "charge"), atomAt(d, atom).charge = charge; },
            "atom"_a, "charge"_a, "Set an atom's formal charge, -8 to 8. (since 2.0)")
        .def(
            "set_radicals", [](Document& d, int atom, int n) { inRange(n, 0, 2, "radicals"), atomAt(d, atom).radicals = n; },
            "atom"_a, "radicals"_a, "Set an atom's radical electrons, 0 to 2 (each takes the place of an H). (since 2.0)")
        .def(
            "set_partial", [](Document& d, int atom, int partial) { inRange(partial, -1, 1, "partial"), atomAt(d, atom).partial = partial; },
            "atom"_a, "partial"_a, "Mark an atom δ+ (1) or δ− (-1), or neither (0); drawn only, not chemistry. (since 2.0)")
        .def(
            "show_symbol",
            [](Document& d, int atom, bool shown) {
                Atom& a = atomAt(d, atom);
                const QString own = QString::fromStdString(chem::symbol(a.z));
                if (a.z == 0 || (!a.label.isEmpty() && a.label != own))
                    throw nb::value_error("show_symbol is for an element; this atom has a label of its own");
                a.label = shown ? own : QString();
            },
            "atom"_a, "shown"_a = true,
            "Draw an atom's element symbol even where skeletal style leaves it out (the C of R-N≡C), or "
            "(shown=False) hide it again. Its chemistry doesn't change. (since 2.0)")
        .def(
            "set_bond_style",
            [](Document& d, int bond, const std::string& style) {
                if (bond < 0 || bond >= int(d.bonds.size())) throw nb::index_error("no such bond");
                const auto* name = std::find_if(std::begin(kBondStereoNames), std::end(kBondStereoNames),
                                                [&](const char* n) { return style == n; });
                if (name == std::end(kBondStereoNames)) {
                    std::string names;
                    for (const char* n : kBondStereoNames) names += std::string(names.empty() ? "'" : ", '") + n + "'";
                    throw nb::value_error(("unknown bond style: " + style + " (use " + names + ")").c_str());
                }
                const auto s = BondStereo(name - std::begin(kBondStereoNames));
                Bond& b = d.bonds[bond];
                if ((s == BondStereo::Wedge || s == BondStereo::Hash || s == BondStereo::Wavy) && b.order != 1)
                    throw nb::value_error(("a " + style + " bond is single").c_str());
                b.stereo = s;
            },
            "bond"_a, "style"_a,
            "Draw a bond as 'none' (plain), 'wedge', 'hash', 'bold', 'dashed', 'wavy', 'interaction' "
            "(dotted, not chemistry), 'partial' (forming or breaking in a transition state), 'dative' or "
            "'double-dashed' (two dashed lines). "
            "A wedge, hash or wavy bond must be single. (since 2.0)")
        .def(
            "add_bracket",
            [](Document& d, std::vector<int> atoms, bool round, const std::string& label, const std::string& superscript) {
                if (atoms.empty()) throw nb::value_error("a bracket needs atoms");
                for (int i : atoms) atomAt(d, i);
                std::sort(atoms.begin(), atoms.end());
                atoms.erase(std::unique(atoms.begin(), atoms.end()), atoms.end());
                d.brackets.push_back({atoms, !round, qs(label), qs(superscript)});
                return int(d.brackets.size()) - 1;
            },
            "atoms"_a, "round"_a = false, "label"_a = "", "superscript"_a = "",
            "Brackets around atoms, square or round, with a subscript label at the bottom right (n) and a "
            "superscript at the top right (‡, •+). Returns its index. (since 2.0)")
        .def(
            "add_text",
            [](Document& d, double x, double y, const std::string& text) {
                if (text.empty()) throw nb::value_error("no text");
                d.texts.push_back({{x, y}, qs(text)});
                return int(d.texts.size()) - 1;
            },
            "x"_a, "y"_a, "text"_a,
            "Add free text, (x, y) the left end of its first baseline; formulas get subscripts (H2O). "
            "Returns its index. (since 2.0)")
        .def(
            "add_arrow",
            [](Document& d, double x0, double y0, double x1, double y1, const std::string& kind,
               const std::string& above, const std::string& below) {
                for (double value : {x0, y0, x1, y1}) requireFinite(value);
                const std::map<std::string, ArrowKind> kinds = {
                    {"reaction", ArrowKind::Reaction}, {"resonance", ArrowKind::Resonance},
                    {"equilibrium", ArrowKind::Equilibrium}, {"retro", ArrowKind::Retro}};
                const auto found = kinds.find(kind);
                if (found == kinds.end()) throw nb::value_error("kind must be reaction, resonance, equilibrium or retro");
                if (x0 == x1 && y0 == y1) throw nb::value_error("an arrow needs different endpoints");
                Arrow arrow{{x0, y0}, {x1, y1}, found->second};
                const QRectF box = arrowPath(arrow).boundingRect();
                std::vector<Text> texts;
                for (const auto& [label, top] : {std::pair{above, true}, std::pair{below, false}}) {
                    if (label.empty()) continue;
                    Text text{{0, 0}, qs(label)};
                    const QRectF ink = textPath(text, documentStyle(d)).boundingRect();
                    text.pos = {box.center().x() - ink.center().x(),
                                top ? box.top() - 0.3 * kBondLength - ink.bottom()
                                    : box.bottom() + 0.3 * kBondLength - ink.top()};
                    requireFinite(text.pos.x()), requireFinite(text.pos.y());
                    texts.push_back(text);
                }
                const int index = int(d.arrows.size());
                d.arrows.push_back(arrow);
                d.texts.insert(d.texts.end(), texts.begin(), texts.end());
                return index;
            },
            "x0"_a, "y0"_a, "x1"_a, "y1"_a, "kind"_a = "reaction", "text_above"_a = "", "text_below"_a = "",
            "Add an arrow from (x0, y0) to (x1, y1) in points: reaction, resonance, equilibrium or retro. "
            "Returns its arrow index. Condition text is placed horizontally, centred above/below the arrow as "
            "ordinary editable text; it is not attached to later endpoint edits. (since 2.1)")
        .def(
            "move_atoms",
            [](Document& d, std::vector<int> atoms, double dx, double dy) {
                requireFinite(dx), requireFinite(dy);
                transformAtoms(d, selectedAtoms(d, std::move(atoms)), QTransform::fromTranslate(dx, dy));
            },
            "atoms"_a, "dx"_a, "dy"_a,
            "Move atoms by (dx, dy) points; [] moves all atoms. Duplicate indices move once. (since 2.1)")
        .def(
            "set_position",
            [](Document& d, int atom, double x, double y) {
                requireFinite(x), requireFinite(y);
                const QPointF by = QPointF(x, y) - atomAt(d, atom).pos;
                transformAtoms(d, {atom}, QTransform::fromTranslate(by.x(), by.y()));
            },
            "atom"_a, "x"_a, "y"_a, "Place an atom at (x, y) in points. (since 2.1)")
        .def(
            "flip",
            [](Document& d, std::vector<int> atoms, bool horizontal) {
                atoms = selectedAtoms(d, std::move(atoms));
                if (atoms.empty()) return;
                Document part;
                for (int i : atoms) part.atoms.push_back(d.atoms[i]);
                part.style = d.style, part.styleSettings = d.styleSettings;
                part.labelRatio = d.labelRatio, part.circledCharges = d.circledCharges;
                const QPointF c = documentBounds(part).center();
                transformAtoms(d, atoms, QTransform::fromTranslate(-c.x(), -c.y())
                               * QTransform::fromScale(horizontal ? -1 : 1, horizontal ? 1 : -1)
                               * QTransform::fromTranslate(c.x(), c.y()));
            },
            "atoms"_a = std::vector<int>{}, "horizontal"_a = true,
            "Mirror atoms about their drawn centre, horizontally or vertically; [] means all atoms. "
            "As in the app, wedges stay wedges (a mirror reverses chirality) and explicit double-bond sides flip. (since 2.1)")
        .def(
            "rotate",
            [](Document& d, std::vector<int> atoms, double degrees, std::optional<std::pair<double, double>> about) {
                requireFinite(degrees);
                atoms = selectedAtoms(d, std::move(atoms));
                QPointF c;
                if (about) {
                    requireFinite(about->first), requireFinite(about->second);
                    c = {about->first, about->second};
                } else {
                    for (int i : atoms) c += d.atoms[i].pos / double(atoms.size());
                }
                transformAtoms(d, atoms, QTransform::fromTranslate(-c.x(), -c.y())
                               * QTransform().rotate(degrees) * QTransform::fromTranslate(c.x(), c.y()));
            },
            "atoms"_a = std::vector<int>{}, nb::kw_only(), "degrees"_a, "about"_a = nb::none(),
            "Rotate atoms clockwise by degrees (y points down), about (x, y) in points or their mean position. "
            "[] means all atoms; bonds and wedges stay attached. (since 2.1)")
        .def(
            "remove_atoms",
            [](Document& d, std::vector<int> atoms) {
                for (int i : atoms)
                    if (i < 0 || i >= int(d.atoms.size())) throw nb::index_error("no such atom");
                std::sort(atoms.begin(), atoms.end());
                atoms.erase(std::unique(atoms.begin(), atoms.end()), atoms.end());
                d.removeAtoms(atoms);
            },
            "atoms"_a,
            "Remove atoms and their bonds. The atoms after each one removed move down to fill its place, so "
            "indices taken before are stale: remove everything in one call. (since 2.0)")
        .def(
            "remove_bonds",
            [](Document& d, std::vector<int> bonds) {
                for (int i : bonds)
                    if (i < 0 || i >= int(d.bonds.size())) throw nb::index_error("no such bond");
                std::sort(bonds.begin(), bonds.end());
                bonds.erase(std::unique(bonds.begin(), bonds.end()), bonds.end());
                for (auto k = bonds.rbegin(); k != bonds.rend(); ++k) d.bonds.erase(d.bonds.begin() + *k);
            },
            "bonds"_a,
            "Remove bonds; their atoms stay. Later bonds move down as in remove_atoms, and atom indices "
            "don't change. (since 2.0)")
        .def_prop_ro("r_sites", &chem::rSites,
                     "The numbers of the R sites (atoms labelled R1, R2…, or [1*] and [*:1] in SMILES), sorted. (since 2.0)")
        .def(
            "substitute",
            [](const Document& d, const nb::dict& fragments) {
                std::map<int, std::string> sites;
                for (auto [key, value] : fragments) {
                    const std::string k = nb::isinstance<nb::int_>(key) ? "R" + std::to_string(nb::cast<int>(key)) : nb::cast<std::string>(key);
                    if (k.size() < 2 || k[0] != 'R' || !std::all_of(k.begin() + 1, k.end(), ::isdigit))
                        throw nb::value_error(("not an R site: " + k + " (use 1 or \"R1\")").c_str());
                    sites[std::stoi(k.substr(1))] = nb::cast<std::string>(value);
                }
                return chem::substitute(d, sites);  // std::invalid_argument arrives as ValueError
            },
            "fragments"_a,
            "A copy with each R site replaced, from a dict keyed by site (1 or \"R1\"): an abbreviation (OMe, Ph) "
            "or SMILES whose first atom, or the atom bonded to its lone * (\"*OC\"), takes the site's place; \"H\" "
            "removes it. The scaffold keeps its drawing and wedges. ValueError names a missing or unknown site, an R "
            "site with more than one bond, or a fragment that doesn't parse. (since 2.0)")
        .def(
            "hotkeys",
            [](Document& d, int atom, const std::string& keys, std::optional<int> bond) {
                edit::Hotspot h{bond ? -1 : atom, bond.value_or(-1)};
                if (h.atom >= int(d.atoms.size()) || h.bond >= int(d.bonds.size()))
                    throw nb::index_error("no such atom or bond");
                for (QChar k : qs(keys)) {
                    h = edit::hotkey(d, h, QString(k));
                    if (!h.valid()) throw nb::value_error(("not a hotkey here: " + QString(k).toStdString()).c_str());
                }
                return std::make_tuple(h.atom, h.bond);
            },
            "atom"_a, "keys"_a, "bond"_a = nb::none(),
            "Type hotkeys with this atom as the hotspot (or atom=-1, bond=i for a bond). "
            "Returns the final (atom, bond) hotspot. (since 0.4)")
        .def(
            "clean", [](Document& d, std::vector<int> atoms) { d = chem::clean2D(d, atoms); }, "atoms"_a = std::vector<int>{},
            "Lay out afresh with RDKit; with atoms, only the molecules containing them. (since 0.4)")
        .def("to_smiles", [](const Document& d) { return chem::toSmiles(d); }, "Canonical SMILES (\"\" if the drawing isn't valid chemistry), as CXSMILES carrying only the stereo groups (&1, or1) when there are any. (since 0.4)")
        .def("to_molblock", [](const Document& d) { return chem::toMolBlock(d); }, "MDL MOL, with the drawing's wedges; V3000 when coordinate bonds require it. (since 0.4)")
        .def("to_inchi", [](const Document& d) { return chem::toInchi(d); }, "Standard InChI. (since 0.4)")
        .def("to_inchikey", [](const Document& d) { return chem::toInchiKey(d); }, "Standard InChIKey. (since 0.4)")
        .def("to_json", [](const Document& d) { return d.toJson().toStdString(); }, "The .penz document (JSON). (since 0.4)")
        .def(
            "bounds",
            [](Document& d, std::vector<int> atoms) {
                Document part = d;
                if (!atoms.empty()) {
                    atoms = selectedAtoms(d, std::move(atoms));
                    std::vector<int> drop;
                    for (int i = 0; i < int(d.atoms.size()); ++i)
                        if (!std::binary_search(atoms.begin(), atoms.end(), i)) drop.push_back(i);
                    part.arrows.clear(), part.texts.clear();
                    part.removeAtoms(drop);
                }
                const QRectF r = outputBounds(part);
                const double s = exportScale(d);
                return std::make_tuple(r.left() * s, r.top() * s, r.right() * s, r.bottom() * s);
            },
            "atoms"_a = std::vector<int>{},
            "Painted bounds (x0, y0, x1, y1) in drawing-style points, without export margin or scale. "
            "Includes labels, charges, radicals and the renderer's small antialiasing allowance. [] measures the "
            "whole drawing; otherwise measures the selected atoms as a standalone subdrawing, without free text "
            "or arrows. Empty drawings return (0, 0, 0, 0). (since 2.1)")
        .def(
            "to_svg",
            [](const Document& d, double margin, std::optional<std::string> background, double scale) {
                return renderSvg(d, pythonExportOptions(d, margin, background, scale)).toStdString();
            },
            nb::kw_only(), "margin"_a = 0.0, "background"_a = nb::none(), "scale"_a = 1.0,
            "SVG with the drawing embedded. Options (since 2.1): margin in points (0), background None for transparent "
            "or a colour name/hex string, and positive scale (1). As in the app, a set page keeps its full frame. (since 0.4)")
        .def(
            "to_png",
            [](const Document& d, double dpi, double margin, std::optional<std::string> background, double scale) {
                const QByteArray png = renderPng(d, pythonExportOptions(d, margin, background, scale, dpi));
                if (png.isEmpty()) throw nb::value_error("cannot render PNG at these dimensions");
                return nb::bytes(png.constData(), png.size());
            },
            "dpi"_a = 300, nb::kw_only(), "margin"_a = 0.0, "background"_a = nb::none(), "scale"_a = 1.0,
            "PNG bytes at positive dpi (300). Options (since 2.1): margin in points (0), background None for transparent "
            "or a colour name/hex string, and positive scale (1). A set page keeps its full frame. (since 0.4)")
        .def_prop_ro("formula", [](const Document& d) { auto p = chem::properties(d); return p ? p->formula : ""; },
                     "Molecular formula, Hill order, all molecules together (since 0.4)")
        .def_prop_ro("mw", [](const Document& d) { auto p = chem::properties(d); return p ? p->mw : 0.0; }, "Molecular weight (since 0.4)")
        .def_prop_ro("exact_mass", [](const Document& d) { auto p = chem::properties(d); return p ? p->exactMass : 0.0; },
                     "Monoisotopic mass (since 0.4)")
        .def("save", &save, "path"_a, "Write .penz (full fidelity), .mol, or ChemDraw .cdxml or .cdx; images go through export(). (since 0.4)")
        .def(
            "export",
            [](const Document& d, const std::string& path, double margin, std::optional<std::string> background, double scale, double dpi) {
                if (!exportDocument(d, qs(path), pythonExportOptions(d, margin, background, scale, dpi)))
                    throw nb::value_error(("cannot export " + path).c_str());
            },
            "path"_a, nb::kw_only(), "margin"_a = 0.0, "background"_a = nb::none(), "scale"_a = 1.0, "dpi"_a = 300,
            "Write .svg, .png or .pdf as the app exports. Options (since 2.1): margin in points (0), background None "
            "for transparent or a colour name/hex string, positive scale (1), and PNG dpi (300). "
            "A set page keeps its full frame. (since 0.4)")
        .def("_repr_svg_", [](const Document& d) { return renderSvg(d).toStdString(); })
        .def("__repr__", [](const Document& d) {
            auto p = chem::properties(d);
            return "<penzene.Document " + (p ? p->formula : std::to_string(d.atoms.size()) + " atoms") + ">";
        });

    m.def("from_smiles", &fromSmiles, "smiles"_a, "A structure from SMILES, laid out in 2D. (since 0.4)");
    m.def("from_json", [](const std::string& json) {
        auto d = Document::fromJson(QByteArray::fromStdString(json));
        if (!d) throw nb::value_error("not a .penz document");
        return *d;
    }, "json"_a, "A Document from .penz JSON (as to_json writes). (since 0.4)");
    m.def("read", &readPath, "path"_a,
          "Open a file: .penz, MOL, ChemDraw .cdxml/.cdx, .rxn, a Penzene SVG/PNG, or every record of an "
          "SDF, .smi or .inchi file laid out as a grid; records that can't be read are left out with a UserWarning "
          "naming them. (since 0.4)");
    // Private: cmake/nmr-table.py builds the NMR shift table with the app's own HOSE codes (#403).
    m.def("_hose_codes", [](const std::string& molblock, int spheres) { return chem::hoseCodes(molblock, spheres); },
          "molblock"_a, "spheres"_a = 4);
    m.def("drawing_styles", [] {
        std::vector<std::string> out;
        for (const auto& s : drawingStyles()) out.push_back(s.name.toStdString());
        return out;
    }, "The drawing style names Document.style takes. (since 0.4)");
}
