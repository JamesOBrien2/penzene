#pragma once
// Structure editing without a UI: the geometry behind the drawing tools and
// every atom/bond hotkey. The canvas and (later) the Python API both use it.
#include "Document.h"

#include <QString>
#include <map>
#include <vector>

namespace edit {
constexpr double kMergeRadius = 0.3 * kBondLength;  // closer than this is the same atom

// Direction for a new bond from `atom` that avoids existing bonds; `newOrder`
// makes it straight on at an sp centre.
QPointF freeDirection(const Document& doc, int atom, int newOrder = 1);
QPointF snapped(QPointF from, QPointF to);  // unit vector rounded to 30°
int atomNear(const Document& doc, QPointF p, double r, int skip = -1);
// What a curved arrow's end at p rests on (Arrow::fromAt): an atom or its electron pair, else a bond.
std::array<int, 2> anchorAt(const Document& doc, QPointF p);
// Where that end settles: a bond's middle, or just off an atom on the side it was drawn (toward
// `other` when drawn on the atom itself).
QPointF snapToAnchor(const Document& doc, std::array<int, 2> at, QPointF p, QPointF other);
// Arrow ends in `after` follow the atoms they rest on as they moved from `before`; an arrow
// moved off its atoms takes whatever it now rests on.
void followAnchors(const Document& before, Document& after);
std::vector<int> moleculeOf(const Document& doc, int atom);  // the atoms bonded to it, directly or not, itself included
// Compound numbers (#504) in scheme order: rows top to bottom, left to right along each. A number's
// suffix is kept, and suffixed numbers that shared a number (2a, 2b) still share one. `from`: the
// numbers taken on the pages before (#569); returns them with this page's added.
struct CompoundCount {
    int last = 0;
    std::map<QString, int> series;  // a suffixed number's old number: its new one
};
CompoundCount renumberCompounds(Document& doc, CompoundCount from = {});
// A compound number in `after` that wasn't moved itself follows the foot of its molecule from `before`.
void followNumbers(const Document& before, Document& after);
int atomAtOrNew(Document& doc, QPointF p, int z = 6);
void link(Document& doc, int a, int b, int order = 1, BondStereo stereo = BondStereo::None);
std::vector<int> addRing(Document& doc, const std::vector<QPointF>& verts, bool aromatic);
std::vector<QPointF> polygon(QPointF centre, QPointF firstVertex, int n);
double circumradius(int n);  // of a ring with standard bonds
std::vector<int> ringAt(Document& doc, QPointF centre, int n, bool aromatic);
std::vector<int> ringOnAtom(Document& doc, int atom, int n, bool aromatic);
void ringOnBond(Document& doc, int bond, int n, bool aromatic);
void chairOnBond(Document& doc, int bond, int edge);
void straightenSp(Document& doc, int bond);
// Fuses each (keep, drop) pair: drop's bonds and ring fills move to keep, then
// drop is removed. Duplicate and self bonds are dropped (the higher order wins).
void mergeAtoms(Document& doc, const std::vector<std::pair<int, int>>& keepDrop);

// Element symbol, abbreviation (drawn as its label) or SMILES (drawn out).
bool applyLabel(Document& doc, int atom, const QString& label, bool anyText = false);
// What applyLabel reads back to the same atom: "C", "13C", "D", "Boc".
QString atomText(const Atom& a);

// ChemDraw's hotkeys, typed with `h` as the hotspot. Returns the new hotspot,
// or an empty one ({-1, -1}) if the key means nothing there.
struct Hotspot {
    int atom = -1, bond = -1;
    bool valid() const { return atom >= 0 || bond >= 0; }
};
Hotspot hotkey(Document& doc, Hotspot h, const QString& key);

// Layers, as in ChemDraw: arrows, shapes and orbitals stack behind or in front of the
// molecule, which counts as one layer. Returns where the moved arrows ended up.
enum class Restack { Front, Forward, Backward, Back };
std::vector<int> restack(Document& doc, const std::vector<int>& arrows, Restack how);
}  // namespace edit
