#pragma once
// The only file that talks to RDKit. Everything else sees a Document.
#include "Document.h"
#include <QStringList>
#include <array>
#include <optional>
#include <string>
#include <vector>

namespace chem {
std::optional<Document> fromSmiles(const std::string& smiles);
std::optional<Document> fromMolBlock(const std::string& block);
std::optional<Document> fromInchi(const std::string& inchi);
// A peptide from its sequence: one-letter (GFLS; lower case for D) or three-letter (Gly-Phe-Leu-Ser).
std::optional<Document> fromSequence(const QString& sequence);
// ChemDraw .cdxml (molecules, arrows, text) or binary .cdx (molecules only,
// where RDKit was built with ChemDraw support).
std::optional<Document> fromChemDraw(const QByteArray& data);
// .penz, .cdxml/.cdx, a Penzene SVG/PNG, MOL, or every record of an SDF,
// .smi or .inchi file laid out as a grid, by extension.
std::optional<Document> readFile(const QString& path);
// Each record of a multi-record file (SDF, .smi "SMILES name" lines, .inchi
// lines), named by the file (its title or name column) or as base-N.
struct Record {
    QString name;
    std::optional<Document> doc;
};
std::vector<Record> readRecords(const QString& path);
std::string toMolBlock(const Document& doc, bool v3000 = false);
std::string toSdf(const Document& doc, bool v3000 = false);  // one record per molecule

// Rotation out of the page: the molecules of `atoms` given a 3D shape (an RDKit
// conformer posed to match the drawing), which project3D turns about the page's
// x then y axis (degrees) and projects back, re-choosing wedges so stereo is kept.
// Unturned it's the drawing exactly; over the first 30° it becomes the conformer.
struct Pose3D {
    std::vector<int> atoms;
    std::vector<std::array<double, 3>> drawn;      // Å, y up, about `centre`: the drawing, with the conformer's depth
    std::vector<std::array<double, 3>> conformer;  // the conformer itself, posed like the drawing
    std::array<double, 3> centre{};
};
std::optional<Pose3D> pose3D(const Document& doc, const std::vector<int>& atoms);
Document project3D(const Document& doc, const Pose3D& pose, double aboutX, double aboutY);
QByteArray toCdxml(const Document& doc);  // molecules, arrows and text
QByteArray toCdx(const Document& doc);    // binary CDX: all that toCdxml writes
// Binary CDX ⇄ CDXML (Revvity's ChemDraw file library); empty if it can't be read.
QByteArray cdxToCdxml(const QByteArray& cdx);
QByteArray cdxmlToCdx(const QByteArray& cdxml);
std::string toSmiles(const Document& doc);  // "" if the structure isn't valid
std::string toSmarts(const Document& doc);  // query atoms (A, Q, X, M, [N,O,S]) as queries
QString queryMeaning(const Atom& a);  // "a metal" for M, "one of N, O, S" for X standing for them; "" for no query

// A drawn reaction: the molecules before, alongside and after its arrow.
struct Reaction {
    std::vector<Document> reactants, agents, products;
};
std::vector<Reaction> reactionsOf(const Document& doc);    // one per reaction arrow, in reading order
std::optional<Reaction> reactionOf(const Document& doc);  // the first; nullopt without a reaction arrow
std::string toReactionSmiles(const Reaction& r);           // reactants>agents>products; empty if any molecule is invalid
std::string toReactionSmiles(const std::vector<Reaction>& steps);  // one line per step; empty if any molecule is invalid
std::string toRxn(const Reaction& r);                      // MDL Rxnfile (V2000)
std::string toRdf(const std::vector<Reaction>& steps);     // MDL RD file: one Rxnfile per step
Document layoutReaction(const std::vector<Reaction>& steps);  // the steps left to right, as one scheme
std::optional<Document> fromReactionSmiles(const std::string& smiles);
std::optional<Document> fromRxn(const std::string& text);
std::optional<Document> fromRdf(const std::string& text);
// New layout, same atom order; each molecule keeps its centroid. With `only`,
// just the molecules containing those atoms are touched.
Document clean2D(const Document& doc, const std::vector<int>& only = {});

struct Properties {
    std::string formula;  // Hill order, all fragments together
    double mw = 0, exactMass = 0;
};
std::optional<Properties> properties(const Document& doc);  // nullopt if empty or invalid
// The properties panel: everything an SI table or a med-chem glance needs.
struct Profile {
    Properties basic;
    double logP = 0, tpsa = 0;
    int hbd = 0, hba = 0, rotatable = 0, heavyAtoms = 0;
    std::vector<std::pair<std::string, double>> elemental;  // symbol, mass %, in Hill order
    int lipinskiViolations = 0;  // Ro5: MW > 500, logP > 5, HBD > 5, HBA > 10
    bool veber = true;           // rotatable bonds <= 10 and TPSA <= 140
};
std::optional<Profile> profile(const Document& doc);  // nullopt if empty or invalid
// The mass spectrum's isotope pattern for an ion of everything in `doc` (all fragments, like the
// formula), from natural abundances; drawn isotopes (13C) count as that isotope only.
enum class Ion { M, MplusH, MplusNa, MminusH };  // M: M+• of a neutral molecule, else the drawn ion
struct Peak {
    double mz, intensity;  // intensity: the tallest is 100
};
std::vector<Peak> isotopePattern(const Document& doc, Ion ion);  // by m/z; empty if it can't be worked out
// The ion's calculated mass as a supporting-information line, e.g. "HRMS (ESI) m/z: [M+H]+ calcd
// for C9H9O4 181.0495"; [M] of a neutral molecule is EI's M+. Empty if it can't be worked out.
QString hrmsLine(const Document& doc, Ion ion);
// The textbook EI ions a chemist looks for, from the groups present (M − 15 for a methyl, m/z 91 for a
// benzyl, McLafferty…), M⁺• first: candidates to check, not a predicted spectrum.
struct EiIon {
    double mz;
    QString formula, from;  // "C7H7+", "benzyl (tropylium)"
};
std::vector<EiIon> eiIons(const Document& doc);
std::string toInchi(const Document& doc);                   // "" if invalid
// One CSV row per record: identifiers and descriptors (docs/cli.md defines them). A record that
// isn't valid chemistry keeps its row, with the reason under "error". columns: a subset, in order.
QStringList descriptorColumns();
std::string descriptorsCsv(const std::vector<Record>& records, const QStringList& columns = {});
std::vector<Document> molecules(const Document& doc);  // each connected piece of a drawing
std::string toInchiKey(const Document& doc);

// Smallest set of smallest rings, each in ring order (abbreviations excluded).
std::vector<std::vector<int>> rings(const Document& doc);

// Explicit hydrogens: add them where atoms have implicit ones (placed by RDKit),
// or remove plain terminal H atoms again (wedged/hashed ones carry stereo, so stay).
Document addHydrogens(const Document& doc);
Document removeHydrogens(const Document& doc);
// CIP descriptors: (R)/(S) (or r/s) on atoms, (E)/(Z) on double bonds.
struct StereoLabel {
    int atom = -1, bond = -1;
    QString text;
};
std::vector<StereoLabel> stereoLabels(const Document& doc);
// Things a chemist would want fixed before publishing: valence errors,
// stereocentres without a wedge, wedges on non-stereocentres, unknown labels,
// overlapping atoms. Each names the atoms involved, to select them.
struct Problem {
    QString message;
    std::vector<int> atoms;
};
std::vector<Problem> checkStructure(const Document& doc);
// Rings RDKit perceives as aromatic (every atom aromatic), in ring order.
std::vector<std::vector<int>> aromaticRings(const Document& doc);

struct AtomInfo {
    int hydrogens = 0;
    bool valenceError = false;
};
std::vector<AtomInfo> atomInfo(const Document& doc);

// Abbreviations (Me, OMe, Boc…): drawn as a label, expanded for chemistry.
std::optional<Atom> abbreviationHead(const QString& label);  // attaching atom; nullopt if unknown
QStringList abbreviations();
// Replaces `atom` with the first atom of `smiles` (or an abbreviation) and lays
// the rest out away from its bonds. New atoms are appended, so indices stay valid.
bool attach(Document& doc, int atom, const std::string& smilesOrAbbreviation);
bool expandLabel(Document& doc, int atom);  // one label drawn out in full, a charged group keeping its charge (N3-)
Document expanded(const Document& doc);  // abbreviations drawn out in full
// Fischer crossings and Haworth rings redrawn with the wedges they mean (chemistry uses this).
Document projectionsAsWedges(const Document& doc);

// NMR (#403). HOSE codes (Bremser, written as CDK writes them), H-suppressed: codes[atom][s - 1] is
// the atom's code to s spheres, s = 1…maxSpheres; hydrogens get none. Empty if the MOL block can't be read.
std::vector<std::vector<std::string>> hoseCodes(const std::string& molBlock, int maxSpheres = 4);
std::vector<std::vector<std::string>> hoseCodes(const Document& doc, int maxSpheres = 4);
// Predicted 13C (carbons) and 1H (H-bearing atoms) shifts in ppm, looked up by HOSE code in a table built
// from nmrshiftdb2 (resources/nmr). spheres: how many spheres matched (4 best, 1 worst); 0 = no prediction.
struct Shift {
    int atom = -1;
    double carbon = 0, proton = 0;
    int carbonSpheres = 0, protonSpheres = 0;
    int hydrogens = 0;  // on the atom, drawn or implicit
    int symmetry = -1;  // equal for atoms the molecule can't tell apart
    int coupled = 0;    // H on neighbouring carbons outside the atom's own set: n in the n + 1 rule (#552)
};
std::vector<Shift> predictShifts(const Document& doc);  // one per atom with a prediction, in atom order
// A predicted 13C or 1H spectrum as sticks, one per set of equivalent atoms, highest ppm first; count is the
// carbons or hydrogens under it. `only`: just these atoms (predicted in the whole drawing, as bonded there).
struct NmrStick {
    double ppm = 0;
    int count = 0;
    std::vector<int> atoms;
    bool weak = false;  // fewer than three spheres matched
    int coupled = 0;    // 1H: first order, the stick is a multiplet of coupled + 1 lines
    QString multiplicity() const { return coupled < 4 ? QString("sdtq"[coupled]) : "m"; }
};
std::vector<NmrStick> nmrSticks(const Document& doc, bool proton, const std::vector<int>& only = {});
// The sticks as a supporting-information line to fill in, e.g. "1H NMR (predicted) δ 3.69 (q, 2H), 1.22 (t, 3H)."
// Empty without a prediction.
QString nmrLine(const Document& doc, bool proton, const std::vector<int>& only = {});

std::string symbol(int z);
std::string elementName(int z);  // "Carbon"
int atomicNumber(const std::string& symbol);  // 0 if unknown
}  // namespace chem
