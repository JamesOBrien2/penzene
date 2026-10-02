// Python bindings over penzene_core: the same engine and renderer as the app.
#include "Chem.h"
#include "Edit.h"
#include "Render.h"

#include <QBuffer>
#include <QFileInfo>
#include <QGuiApplication>
#include <QImage>
#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <algorithm>
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
    auto doc = chem::readFile(qs(path));
    if (!doc) throw nb::value_error(("cannot read " + path).c_str());
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

}  // namespace

NB_MODULE(_penzene, m) {
    ensureApp();
    m.doc() = "Penzene: 2D chemical structures, drawn by the same engine as the app.";
    m.attr("__version__") = PENZENE_VERSION;

    nb::class_<Atom>(m, "Atom", "An atom of a Document (read-only; edit through Document: set_label, remove_atoms). (since 0.4)")
        .def_prop_ro("symbol", &symbol, "Element symbol, or the abbreviation (Boc, OMe…) if it has one (since 0.4)")
        .def_prop_ro("x", [](const Atom& a) { return a.pos.x(); }, "Points, x right (since 0.4)")
        .def_prop_ro("y", [](const Atom& a) { return a.pos.y(); }, "Points, y down (since 0.4)")
        .def_ro("charge", &Atom::charge, "Formal charge (since 0.4)")
        .def_ro("isotope", &Atom::isotope, "Mass number (13 for carbon-13, 2 for deuterium); 0 for natural abundance (since 1.3)")
        .def("__repr__", [](const Atom& a) { return "<Atom " + symbol(a) + ">"; });

    nb::class_<Bond>(m, "Bond", "A bond of a Document (read-only). (since 0.4)")
        .def_ro("a", &Bond::a, "Index of the first atom (a wedge starts here) (since 0.4)")
        .def_ro("b", &Bond::b, "Index of the second atom (since 0.4)")
        .def_ro("order", &Bond::order, "1, 2 or 3 (since 0.4)")
        .def("__repr__", [](const Bond& b) {
            return "<Bond " + std::to_string(b.a) + "-" + std::to_string(b.b) + " order " + std::to_string(b.order) + ">";
        });

    nb::class_<Document>(m, "Document", "A drawing: molecules, arrows and text, as the app holds it. (since 0.4)")
        .def(nb::init<>(), "An empty drawing. (since 0.4)")
        .def_ro("atoms", &Document::atoms, "The atoms, in order (their index is what add_bond and hotkeys take) (since 0.4)")
        .def_ro("bonds", &Document::bonds, "The bonds (since 0.4)")
        .def_prop_rw(
            "style", [](const Document& d) { return drawingStyle(d.style).name.toStdString(); },
            [](Document& d, const std::string& name) {
                const DrawingStyle& s = drawingStyle(qs(name));
                if (s.name != qs(name)) {
                    QStringList names;
                    for (const auto& known : drawingStyles()) names << "'" + known.name + "'";
                    throw nb::value_error(("unknown drawing style: " + name + " (use " + names.join(", ").toStdString() + ")").c_str());
                }
                d.style = s.name == drawingStyles()[0].name ? QString() : s.name;
            },
            "Drawing style: 'ACS 1996' (default), 'JDP' or 'RSC' (since 0.4)")
        .def(
            "add_atom",
            [](Document& d, const std::string& label, double x, double y) {
                int i = d.addAtom({x, y});
                if (!edit::applyLabel(d, i, qs(label))) {
                    d.atoms.pop_back();
                    throw nb::value_error(("not an element, abbreviation or SMILES: " + label).c_str());
                }
                return i;
            },
            "label"_a = "C", "x"_a = 0.0, "y"_a = 0.0,
            "Add an atom: an element, an abbreviation (OMe, Boc…) or a SMILES fragment. Returns its index. (since 0.4)")
        .def(
            "add_bond",
            [](Document& d, int a, int b, int order) {
                const int n = int(d.atoms.size());
                if (a < 0 || b < 0 || a >= n || b >= n || a == b) throw nb::index_error("no such atoms");
                if (order < 1 || order > 3) throw nb::value_error("order must be 1, 2 or 3");
                edit::link(d, a, b, order);
            },
            "a"_a, "b"_a, "order"_a = 1, "Bond atoms a and b (order 1, 2 or 3). (since 0.4)")
        .def(
            "set_label",
            [](Document& d, int atom, const std::string& label) {
                if (atom < 0 || atom >= int(d.atoms.size())) throw nb::index_error("no such atom");
                Document next = d;
                if (!edit::applyLabel(next, atom, qs(label)))
                    throw nb::value_error(("not an element, abbreviation or SMILES: " + label).c_str());
                d = next;
            },
            "atom"_a, "label"_a,
            "Relabel an atom, as add_atom takes a label: an element (\"C\" undoes a label), an abbreviation or a "
            "SMILES fragment (drawn out from it, its new atoms added at the end). (since 2.0)")
        .def(
            "set_bond_order",
            [](Document& d, int bond, int order) {
                if (bond < 0 || bond >= int(d.bonds.size())) throw nb::index_error("no such bond");
                if (order < 1 || order > 3) throw nb::value_error("order must be 1, 2 or 3");
                Bond& b = d.bonds[bond];
                b.order = order;
                if (order != 1) b.stereo = BondStereo::None;  // a wedge or hash is a single bond's
            },
            "bond"_a, "order"_a, "Make a bond single, double or triple (a wedge or hash goes with a single bond). (since 2.0)")
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
        .def("to_molblock", [](const Document& d) { return chem::toMolBlock(d); }, "MDL MOL (V2000), with the drawing's wedges. (since 0.4)")
        .def("to_inchi", [](const Document& d) { return chem::toInchi(d); }, "Standard InChI. (since 0.4)")
        .def("to_inchikey", [](const Document& d) { return chem::toInchiKey(d); }, "Standard InChIKey. (since 0.4)")
        .def("to_json", [](const Document& d) { return d.toJson().toStdString(); }, "The .penz document (JSON). (since 0.4)")
        .def("to_svg", [](const Document& d) { return renderSvg(d).toStdString(); }, "SVG, as the app exports it (the drawing embedded, so it reopens editable). (since 0.4)")
        .def(
            "to_png",
            [](const Document& d, double dpi) {
                QByteArray png;
                QBuffer buf(&png);
                buf.open(QIODevice::WriteOnly);
                renderImage(d, {dpi}).save(&buf, "PNG");
                return nb::bytes(png.constData(), png.size());
            },
            "dpi"_a = 300, "PNG bytes at this resolution. (since 0.4)")
        .def_prop_ro("formula", [](const Document& d) { auto p = chem::properties(d); return p ? p->formula : ""; },
                     "Molecular formula, Hill order, all molecules together (since 0.4)")
        .def_prop_ro("mw", [](const Document& d) { auto p = chem::properties(d); return p ? p->mw : 0.0; }, "Molecular weight (since 0.4)")
        .def_prop_ro("exact_mass", [](const Document& d) { auto p = chem::properties(d); return p ? p->exactMass : 0.0; },
                     "Monoisotopic mass (since 0.4)")
        .def("save", &save, "path"_a, "Write .penz (full fidelity), .mol, or ChemDraw .cdxml or .cdx; images go through export(). (since 0.4)")
        .def(
            "export",
            [](const Document& d, const std::string& path) {
                if (!exportDocument(d, qs(path))) throw nb::value_error(("cannot export " + path).c_str());
            },
            "path"_a, "Write .svg, .png or .pdf, exactly as the app exports. (since 0.4)")
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
          "SDF, .smi or .inchi file laid out as a grid. (since 0.4)");
    // Private: cmake/nmr-table.py builds the NMR shift table with the app's own HOSE codes (#403).
    m.def("_hose_codes", [](const std::string& molblock, int spheres) { return chem::hoseCodes(molblock, spheres); },
          "molblock"_a, "spheres"_a = 4);
    m.def("drawing_styles", [] {
        std::vector<std::string> out;
        for (const auto& s : drawingStyles()) out.push_back(s.name.toStdString());
        return out;
    }, "The drawing style names Document.style takes. (since 0.4)");
}
