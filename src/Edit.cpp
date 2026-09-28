#include "Edit.h"
#include "Chem.h"
#include "Geometry.h"

#include <QHash>
#include <QRegularExpression>
#include <algorithm>

namespace edit {


// Direction for a new bond from `atom` that avoids existing bonds.
// `newOrder` is the order of the bond about to be added: it makes the atom sp
// (straight on) after a triple bond, or when it cumulates two double bonds.
QPointF freeDirection(const Document& doc, int atom, int newOrder) {
    auto nbs = doc.neighbors(atom);
    QPointF p = doc.atoms[atom].pos;
    if (nbs.empty()) return dirAt(-30);
    if (nbs.size() == 1) {
        QPointF back = unit(doc.atoms[nbs[0]].pos - p);
        int have = doc.bonds[doc.bondBetween(atom, nbs[0])].order;
        if (have == 3 || newOrder == 3 || (have == 2 && newOrder == 2)) return -back;
        double base = qRadiansToDegrees(std::atan2(back.y(), back.x()));
        // Zig-zag: of the two 120° options, take the one farther from everything else.
        QPointF best;
        double bestScore = -1;
        for (double off : {120.0, -120.0}) {
            QPointF cand = p + dirAt(base + off) * kBondLength;
            double score = 1e9;
            for (size_t j = 0; j < doc.atoms.size(); ++j)
                if (int(j) != atom) score = std::min(score, len(doc.atoms[j].pos - cand));
            if (score > bestScore) bestScore = score, best = dirAt(base + off);
        }
        return best;
    }
    return doc.awayDirection(atom);
}

QPointF snapped(QPointF from, QPointF to) {
    QPointF v = to - from;
    double deg = qRadiansToDegrees(std::atan2(v.y(), v.x()));
    return dirAt(std::round(deg / 30) * 30);
}

int atomNear(const Document& doc, QPointF p, double r, int skip) {
    int best = -1;
    for (size_t i = 0; i < doc.atoms.size(); ++i) {
        double d = len(doc.atoms[i].pos - p);
        if (int(i) != skip && d < r) r = d, best = int(i);
    }
    return best;
}

// Returns the atom at `p`, creating one if nothing is close enough.
int atomAtOrNew(Document& doc, QPointF p, int z) {
    int i = atomNear(doc, p, kMergeRadius);
    return i >= 0 ? i : doc.addAtom(p, z);
}

void link(Document& doc, int a, int b, int order, BondStereo stereo) {
    if (a == b || doc.bondBetween(a, b) >= 0) return;
    doc.bonds.push_back({a, b, order, stereo});
}

bool hasDouble(const Document& doc, int atom) {
    return std::any_of(doc.bonds.begin(), doc.bonds.end(),
                       [&](const Bond& b) { return b.order > 1 && (b.a == atom || b.b == atom); });
}

// Adds a ring through `verts` (merging with existing atoms) and optionally
// alternates double bonds where the atoms are still free for one.
std::vector<int> addRing(Document& doc, const std::vector<QPointF>& verts, bool aromatic) {
    std::vector<int> ids;
    for (QPointF v : verts) ids.push_back(atomAtOrNew(doc, v));
    const int n = int(ids.size());
    for (int k = 0; k < n; ++k) link(doc, ids[k], ids[(k + 1) % n]);
    if (!aromatic) return ids;
    for (int k = 1; k <= n; ++k) {  // start after the (possibly shared) first edge
        int a = ids[k % n], b = ids[(k + 1) % n];
        int bi = doc.bondBetween(a, b);
        if (bi >= 0 && doc.bonds[bi].order == 1 && !hasDouble(doc, a) && !hasDouble(doc, b))
            doc.bonds[bi].order = 2;
    }
    return ids;
}

std::vector<QPointF> polygon(QPointF centre, QPointF firstVertex, int n) {
    std::vector<QPointF> out;
    QPointF r = firstVertex - centre;
    for (int k = 0; k < n; ++k) {
        double t = 2 * M_PI * k / n;
        out.push_back(centre + QPointF(r.x() * std::cos(t) - r.y() * std::sin(t),
                                       r.x() * std::sin(t) + r.y() * std::cos(t)));
    }
    return out;
}

double circumradius(int n) { return kBondLength / (2 * std::sin(M_PI / n)); }

std::vector<int> ringAt(Document& doc, QPointF centre, int n, bool aromatic) {
    return addRing(doc, polygon(centre, centre + QPointF(0, -circumradius(n)), n), aromatic);
}

// Ring through the atom, pointing away from its bonds so they bisect the ring's outside angle.
std::vector<int> ringOnAtom(Document& doc, int atom, int n, bool aromatic) {
    QPointF p = doc.atoms[atom].pos;
    return addRing(doc, polygon(p + doc.awayDirection(atom) * circumradius(n), p, n), aromatic);
}

// Ring fused onto the bond, on the side away from the other neighbours.
void ringOnBond(Document& doc, int bond, int n, bool aromatic) {
    const Bond& b = doc.bonds[bond];
    QPointF pa = doc.atoms[b.a].pos, pb = doc.atoms[b.b].pos, d = unit(pb - pa);
    double side = 0;
    for (int end : {b.a, b.b})
        for (int nb : doc.neighbors(end))
            if (nb != b.a && nb != b.b) side += cross(d, doc.atoms[nb].pos - pa);
    double apothem = kBondLength / (2 * std::tan(M_PI / n));
    QPointF centre = (pa + pb) / 2 + perp(d) * (side > 0 ? -apothem : apothem);
    auto verts = polygon(centre, pa, n);
    if (len(verts[1] - pb) > 1) verts = polygon(centre, pb, n);  // wind the right way
    addRing(doc, verts, aromatic);
}

// Chair cyclohexane fused onto the bond, on the side away from the other neighbours. The bond is
// one end of the chair, whose other atoms all lie on one side of it: a long edge has atoms on both
// sides, so they reached back over a ring already there (#417). `edge` (0 or 1, the ChemDraw 9 / 0
// keys) picks the mirror image: the chair's pointed end at one atom of the bond or the other.
void chairOnBond(Document& doc, int bond, int edge) {
    // Opposite edges parallel; roughly unit bonds. The other atoms are all left of 2 -> 3.
    static const QPointF chair[6] = {{0, 0}, {0.95, 0.35}, {1.95, 0.05}, {2.55, 0.75}, {1.6, 0.4}, {0.6, 0.7}};
    const Bond& b = doc.bonds[bond];
    QPointF pa = doc.atoms[b.a].pos, pb = doc.atoms[b.b].pos;
    double side = 0;
    for (int end : {b.a, b.b})
        for (int nb : doc.neighbors(end))
            if (nb != b.a && nb != b.b) side += cross(pb - pa, doc.atoms[nb].pos - pa);
    if (side > 0) std::swap(pa, pb);  // the neighbours on the right, the chair on the left
    const int step = edge == 0 ? 1 : -1, first = edge == 0 ? 2 : 3;  // 0: mirrored, 3 -> 2
    auto t = [&](int k) {
        QPointF p = chair[(first + step * k + 6) % 6];
        return QPointF(p.x(), p.y() * step);
    };
    const QPointF d = pb - pa, td = t(1) - t(0);
    const double turn = qRadiansToDegrees(std::atan2(d.y(), d.x()) - std::atan2(td.y(), td.x()));
    std::vector<QPointF> verts;
    for (int k = 0; k < 6; ++k) verts.push_back(pa + rotated(t(k) - t(0), turn) * (len(d) / len(td)));
    addRing(doc, verts, false);
}

// After a bond order change: if an end became an sp centre with two neighbours,
// swing a terminal neighbour into line (Clean handles the general case).
void straightenSp(Document& doc, int bond) {
    for (int e : {doc.bonds[bond].a, doc.bonds[bond].b}) {
        auto nbs = doc.neighbors(e);
        if (nbs.size() != 2 || !isSp(doc, e)) continue;
        for (int k : {0, 1}) {
            int mover = nbs[k], anchor = nbs[1 - k];
            if (doc.neighbors(mover).size() != 1) continue;
            QPointF c = doc.atoms[e].pos;
            doc.atoms[mover].pos = c + unit(c - doc.atoms[anchor].pos) * len(doc.atoms[mover].pos - c);
            break;
        }
    }
}

void mergeAtoms(Document& doc, const std::vector<std::pair<int, int>>& keepDrop) {
    if (keepDrop.empty()) return;
    std::vector<int> target(doc.atoms.size());
    for (int i = 0; i < int(target.size()); ++i) target[i] = i;
    for (auto [keep, drop] : keepDrop) target[drop] = keep;
    std::vector<Bond> bonds;
    for (Bond b : doc.bonds) {
        b.a = target[b.a], b.b = target[b.b];
        if (b.a == b.b) continue;
        auto same = std::find_if(bonds.begin(), bonds.end(), [&](const Bond& o) {
            return (o.a == b.a && o.b == b.b) || (o.a == b.b && o.b == b.a);
        });
        if (same == bonds.end()) bonds.push_back(b);
        else if (b.order > same->order) *same = b;
    }
    doc.bonds = std::move(bonds);
    for (auto& f : doc.fills)
        for (int& i : f.atoms) i = target[i];
    auto moveOnto = [&](std::vector<int>& ids) {  // the kept atom takes the dropped one's place
        for (int& i : ids) i = target[i];
        std::sort(ids.begin(), ids.end());
        ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    };
    for (auto& b : doc.brackets) moveOnto(b.atoms);
    for (auto& ring : doc.aromaticCircleOverrides) moveOnto(ring);
    std::vector<int> drops;
    for (auto [keep, drop] : keepDrop) drops.push_back(drop);
    doc.removeAtoms(drops);
}

namespace {
enum class Site { Primary, Secondary, Tertiary, Aromatic };

bool inRing(const Document& doc, int at) {
    auto nbs = doc.neighbors(at);
    for (int start : nbs) {  // can another neighbour be reached without passing through `at`?
        std::vector<bool> seen(doc.atoms.size());
        seen[at] = seen[start] = true;
        std::vector<int> stack{start};
        while (!stack.empty()) {
            int i = stack.back();
            stack.pop_back();
            for (int nb : doc.neighbors(i)) {
                if (nb == at && i != start) return true;
                if (!seen[nb]) seen[nb] = true, stack.push_back(nb);
            }
        }
    }
    return false;
}

Site site(const Document& doc, int at) {
    size_t deg = doc.neighbors(at).size();
    if (deg <= 1) return Site::Primary;
    if (deg >= 3) return Site::Tertiary;
    return inRing(doc, at) && hasDouble(doc, at) ? Site::Aromatic : Site::Secondary;
}

// The two 120° directions open to an atom with at most one bond: the zig-zag
// ("linear mode") one and the other ("cyclic mode").
std::pair<QPointF, QPointF> openDirections(const Document& doc, int at) {
    QPointF lin = freeDirection(doc, at);
    auto nbs = doc.neighbors(at);
    if (nbs.empty()) return {lin, dirAt(90)};
    QPointF back = unit(doc.atoms[nbs[0]].pos - doc.atoms[at].pos);
    // Reflect `lin` across the bond axis.
    QPointF other = 2 * QPointF::dotProduct(lin, back) * back - lin;
    return {lin, other};
}


int sprout(Document& doc, int from, QPointF dir, int order = 1, int z = 6,
           BondStereo stereo = BondStereo::None, double length = kBondLength) {
    int n = atomAtOrNew(doc, doc.atoms[from].pos + dir * length, z);
    link(doc, from, n, order, stereo);
    return n;
}

// "Adds a C-C bond and then ...": tertiary and aromatic sites grow from a new carbon.
int linker(Document& doc, int at) { return sprout(doc, at, doc.awayDirection(at)); }

int farthestFrom(const Document& doc, const std::vector<int>& ids, int from) {
    int best = from;
    double d = -1;
    for (int i : ids)
        if (double l = len(doc.atoms[i].pos - doc.atoms[from].pos); l > d) d = l, best = i;
    return best;
}

constexpr int kUnhandled = -2;

// Atom "sprout" hotkeys. Returns the new hotspot atom, or kUnhandled.
int sproutHotkey(Document& doc, int at, const QString& key) {
    const Site s = site(doc, at);
    const bool needsLinker = s == Site::Tertiary || s == Site::Aromatic;

    if (key == "1") return sprout(doc, at, freeDirection(doc, at));
    if (key == "0") {
        if (s == Site::Primary) return sprout(doc, at, openDirections(doc, at).second);
        return sprout(doc, at, doc.awayDirection(at), 1, 6, BondStereo::None, 1.5 * kBondLength);
    }
    if (key == "2") {  // carbonyl / acetyl
        if (s == Site::Secondary) {
            sprout(doc, at, doc.awayDirection(at), 2, 8);
            return at;
        }
        int c = needsLinker ? linker(doc, at) : at;
        auto [lin, other] = openDirections(doc, c);
        sprout(doc, c, other, 2, 8);
        return sprout(doc, c, lin);
    }
    if (key == "3" || key == "a") {  // phenyl
        int c = s == Site::Primary ? at : linker(doc, at);
        return farthestFrom(doc, ringOnAtom(doc, c, 6, true), c);
    }
    if (key == "4" || key == "5") {  // wedged / hashed methyl
        const BondStereo st = key == "4" ? BondStereo::Wedge : BondStereo::Hash;
        if (s == Site::Secondary || s == Site::Tertiary) {
            sprout(doc, at, doc.awayDirection(at), 1, 6, st);
            return at;
        }
        int c = s == Site::Aromatic ? linker(doc, at) : at;
        auto [lin, other] = openDirections(doc, c);
        sprout(doc, c, other, 1, 6, st);
        return sprout(doc, c, lin);
    }
    static const QHash<QString, int> rings{{"6", 6}, {"7", 5}, {"u", 4}, {"v", 3}};
    if (rings.contains(key)) {  // cycloalkyl; spiro on a secondary carbon
        int c = needsLinker ? linker(doc, at) : at;
        return farthestFrom(doc, ringOnAtom(doc, c, rings[key], false), c);
    }
    if (key == "8") {  // methylidene
        int c = needsLinker ? linker(doc, at) : at;
        QPointF dir = site(doc, c) == Site::Primary ? freeDirection(doc, c, 2) : doc.awayDirection(c);
        return sprout(doc, c, dir, 2);
    }
    if (key == "9") {  // dimethyl / gem-dimethyl / isopropyl
        if (s == Site::Secondary) {
            QPointF away = doc.awayDirection(at);
            sprout(doc, at, rotated(away, 60));
            sprout(doc, at, rotated(away, -60));
            return at;
        }
        int c = needsLinker ? linker(doc, at) : at;
        auto [lin, other] = openDirections(doc, c);
        sprout(doc, c, lin);
        sprout(doc, c, other);
        return c;
    }
    if (key == "z") {  // alkyne, linear
        QPointF dir = freeDirection(doc, at);
        int c1 = sprout(doc, at, dir);
        return sprout(doc, c1, dir, 3);
    }
    if (key == "k") {  // sulfonyl
        QPointF dir = freeDirection(doc, at);
        int sulfur = sprout(doc, at, dir, 1, 16);
        sprout(doc, sulfur, perp(dir), 2, 8);
        sprout(doc, sulfur, -perp(dir), 2, 8);
        return sulfur;
    }
    if (key == "K") {  // t-Bu at 90°
        QPointF dir = freeDirection(doc, at);
        int c = sprout(doc, at, dir);
        sprout(doc, c, dir);
        sprout(doc, c, perp(dir));
        sprout(doc, c, -perp(dir));
        return at;
    }
    return kUnhandled;
}

}  // namespace

// Atom label hotkeys (the hotspot atom becomes this element or group).
QString labelHotkey(const QString& key) {
    static const QHash<QString, QString> k{
        {"c", "C"},   {"n", "N"},    {"w", "N"},     {"o", "O"},    {"q", "O"},   {"s", "S"},
        {"p", "P"},   {"f", "F"},    {"l", "Cl"},    {"C", "Cl"},   {"b", "Br"},  {"i", "I"},
        {"h", "H"},   {"d", "D"},    {"B", "B"},     {"S", "Si"},   {"L", "Li"},  {"m", "Me"},
        {"e", "Et"},  {"A", "Ac"},   {"P", "Ph"},    {"F", "CF3"},  {"N", "NO2"}, {"O", "OMe"},
        {"E", "CO2Me"}, {"Z", "N3"}, {"M", "MgBr"},  {"Q", "Fmoc"}, {"H", "Cbz"}, {"Y", "Boc"}, {"y", "Boc"},
        {"x", "X"},   {"r", "R"},
    };
    return k.value(key);
}

// Element symbol, abbreviation (drawn as its label) or SMILES (drawn out); with
// anyText, other text (R, X, MgEt) too: shown as written, its chemistry unspecified (z 0).
bool applyLabel(Document& doc, int at, const QString& label, bool anyText) {
    Atom& a = doc.atoms[at];
    // Abbreviations first: in a drawing, Ac, Pr and Ts mean acetyl, propyl and
    // tosyl, not actinium, praseodymium and tennessine (#125).
    if (auto head = chem::abbreviationHead(label)) {
        a.z = head->z, a.charge = head->charge, a.label = label, a.isotope = 0;
        return true;
    }
    // In a drawing "Ar" is an aryl group, not argon (the strict API keeps the element).
    if (anyText && label.trimmed() == "Ar") {
        a.z = 0, a.charge = 0, a.label = "Ar", a.isotope = 0;
        return true;
    }
    // "OH", "NH2": the element; hydrogens are implicit. A mass number before it is an
    // isotope ("13C", "18OH"); D and T are hydrogen-2 and -3. Gives {z, mass}, z 0 if not one.
    static const QRegularExpression element("^(\\d{0,3})([A-Z][a-z]?)(H\\d*)?$");
    auto elementOf = [](const QString& s) -> std::pair<int, int> {
        const auto m = element.match(s.trimmed());
        if (!m.hasMatch()) return {0, 0};
        const QString sym = m.captured(2);
        const bool heavyH = (sym == "D" || sym == "T") && m.captured(1).isEmpty();
        const int z = heavyH ? 1 : chem::atomicNumber(sym.toStdString()), mass = heavyH ? (sym == "D" ? 2 : 3) : m.captured(1).toInt();
        return z > 0 && (mass == 0 || mass >= z) ? std::pair{z, mass} : std::pair{0, 0};  // no mass below the proton count
    };
    if (const auto [z, mass] = elementOf(label); z > 0) {
        a.z = z, a.isotope = mass, a.label.clear();
        return true;
    }
    // With a charge: NH3+, O-, Na+, O2- (#324). Digits before the sign are the
    // charge only on a bare element that isn't an abbreviation (Fe3+); NH3+ is NH3 with one plus.
    static const QRegularExpression metalIon("^(?<base>[A-Z][a-z]?)(?<n>\\d+)(?<sign>[+\\-\\x{2212}])$"),
        ion("^(?<base>.+?)(?<sign>[+\\-\\x{2212}])(?<n>\\d*)$");
    auto charge = [](const QRegularExpressionMatch& m) {
        return (m.captured("sign") == "+" ? 1 : -1) * (m.captured("n").isEmpty() ? 1 : m.captured("n").toInt());
    };
    // An abbreviation first: N3- is the azide anion, the charge on the group's attaching atom (#370).
    if (const auto m = ion.match(label); m.hasMatch())
        if (const auto head = chem::abbreviationHead(m.captured("base"))) {
            a.z = head->z, a.charge = head->charge + charge(m), a.label = m.captured("base"), a.isotope = 0;
            return true;
        }
    auto m = metalIon.match(label);
    if (!m.hasMatch() || !elementOf(m.captured("base")).first) m = ion.match(label);
    if (const auto [z, mass] = m.hasMatch() ? elementOf(m.captured("base")) : std::pair{0, 0}; z > 0) {
        a.z = z, a.isotope = mass, a.label.clear();
        a.charge = charge(m);
        return true;
    }
    if (chem::attach(doc, at, label.toStdString())) return true;
    if (!anyText || label.trimmed().isEmpty()) return false;
    a.z = 0, a.charge = 0, a.label = label.trimmed(), a.isotope = 0;
    return true;
}

QString atomText(const Atom& a) {
    if (!a.label.isEmpty()) return a.label;
    if (a.z == 1 && (a.isotope == 2 || a.isotope == 3)) return a.isotope == 2 ? "D" : "T";
    const QString sym = QString::fromStdString(chem::symbol(a.z));
    return a.isotope ? QString::number(a.isotope) + sym : sym;
}


Hotspot hotkey(Document& doc, Hotspot h, const QString& t) {
    if (h.atom >= 0) {
        const int at = h.atom;
        if (t == "+" || t == "-") {
            doc.atoms[at].charge += t == "+" ? 1 : -1;
            return h;
        }
        if (t == ".") {  // attachment point: a wavy bond to a bare point (* in SMILES)
            sprout(doc, at, freeDirection(doc, at), 1, 0, BondStereo::Wavy);
            return h;
        }
        if (t == "j" || t == "J") {  // η5-cyclopentadienyl / η6-benzene, bonded through the ring's centre
            // ponytail: η-bonds have no SMILES or MOL form; the centroid is a bare dummy (*) there.
            const int n = t == "j" ? 5 : 6;
            const QPointF centre = doc.atoms[at].pos + doc.awayDirection(at) * (1.6 * kBondLength);
            const auto ring = ringAt(doc, centre, n, n == 6);  // may reuse atoms already there
            if (n == 5) {  // Cp⁻: two double bonds and the charge
                doc.bonds[doc.bondBetween(ring[0], ring[1])].order = 2;
                doc.bonds[doc.bondBetween(ring[2], ring[3])].order = 2;
                doc.atoms[ring[4]].charge = -1;
            }
            link(doc, at, doc.addAtom(centre, 0));
            return h;
        }
        if (t == ":") {  // lone pairs: none, 1, 2, 3, none
            doc.atoms[at].lonePairs = (doc.atoms[at].lonePairs + 1) % 4;
            return h;
        }
        if (t == "*") {  // radical dot on / off
            doc.atoms[at].radicals = doc.atoms[at].radicals ? 0 : 1;
            return h;
        }
        if (t == "'") {  // atom-map number: the next free one, or off again
            int next = 0;
            for (const auto& a : doc.atoms) next = std::max(next, a.map);
            doc.atoms[at].map = doc.atoms[at].map ? 0 : next + 1;
            return h;
        }
        int hot = sproutHotkey(doc, at, t);
        if (hot == kUnhandled) {
            QString label = labelHotkey(t);
            if (label.isEmpty() || !applyLabel(doc, at, label, true)) return {};
            hot = at;
        }
        return {hot, -1};
    }
    if (h.bond < 0) return {};
    Bond& b = doc.bonds[h.bond];
    static const QHash<QString, std::pair<int, bool>> fuse{
        {"a", {6, true}}, {"z", {5, true}}, {"v", {3, false}}, {"4", {4, false}},
        {"5", {5, false}}, {"6", {6, false}}, {"7", {7, false}}, {"8", {8, false}}};
    if (t == "2" && b.order == 2 && b.stereo == BondStereo::None) {
        // Already double: move the second line to the other side (centred goes to one side).
        b.position = doubleBondSide(doc, b) > 0 ? BondPosition::Left : BondPosition::Right;
    } else if (t == "1" || t == "2" || t == "3") {
        b.order = t.toInt(), b.stereo = BondStereo::None, b.position = BondPosition::Auto;
        straightenSp(doc, h.bond);
    } else if (t == "w" || t == "h" || t == "H" || t == "W") {
        BondStereo s = t == "w" ? BondStereo::Wedge : BondStereo::Hash;
        if (b.stereo == s) std::swap(b.a, b.b);  // again: flip which end is narrow
        b.stereo = s, b.order = 1;
    } else if (fuse.contains(t)) {
        ringOnBond(doc, h.bond, fuse[t].first, fuse[t].second);
    } else if (t == "9" || t == "0") {
        chairOnBond(doc, h.bond, t == "9" ? 0 : 1);
    } else if (t == "d" || t == "b" || t == "y" || t == "D" || t == "B") {
        static const QHash<QString, BondStereo> styles{{"d", BondStereo::Dashed}, {"b", BondStereo::Bold},
                                                       {"y", BondStereo::Wavy},   {"D", BondStereo::Dashed},
                                                       {"B", BondStereo::Bold}};
        b.stereo = styles[t];
        b.order = t == "D" || t == "B" ? 2 : 1;
    } else if (t == "f") {  // bring to front: drawn last, over bonds it crosses
        const Bond front = b;
        doc.bonds.erase(doc.bonds.begin() + h.bond);
        doc.bonds.push_back(front);
        return {-1, int(doc.bonds.size()) - 1};
    } else if (t == "i" || t == "p" || t == "P") {  // interaction; partial (forming/breaking) single or double
        b.stereo = t == "i" ? BondStereo::Interaction : BondStereo::Partial;
        b.order = t == "P" ? 2 : 1;
        b.position = BondPosition::Auto;
    } else if (t == "l" || t == "c" || t == "r") {
        if (b.order != 2) b.order = 2, b.stereo = BondStereo::None;
        b.position = t == "l" ? BondPosition::Left : t == "c" ? BondPosition::Centre : BondPosition::Right;
    } else {
        return {};
    }
    return h;
}


std::vector<int> restack(Document& doc, const std::vector<int>& arrows, Restack how) {
    constexpr int kMolecule = -1;
    std::vector<int> order;  // bottom to top
    for (int i = 0; i < int(doc.arrows.size()); ++i)
        if (doc.arrows[i].behind) order.push_back(i);
    order.push_back(kMolecule);
    for (int i = 0; i < int(doc.arrows.size()); ++i)
        if (!doc.arrows[i].behind) order.push_back(i);
    auto moving = [&](int t) { return t != kMolecule && std::find(arrows.begin(), arrows.end(), t) != arrows.end(); };
    if (how == Restack::Front || how == Restack::Back) {
        std::stable_partition(order.begin(), order.end(), [&](int t) { return moving(t) == (how == Restack::Back); });
    } else if (how == Restack::Forward) {  // each one up past the next unselected layer
        for (int k = int(order.size()) - 2; k >= 0; --k)
            if (moving(order[k]) && !moving(order[k + 1])) std::swap(order[k], order[k + 1]);
    } else {
        for (int k = 1; k < int(order.size()); ++k)
            if (moving(order[k]) && !moving(order[k - 1])) std::swap(order[k], order[k - 1]);
    }
    std::vector<Arrow> stacked;
    std::vector<int> moved;
    bool behind = true;
    for (int t : order) {
        if (t == kMolecule) {
            behind = false;
            continue;
        }
        if (moving(t)) moved.push_back(int(stacked.size()));
        stacked.push_back(doc.arrows[t]);
        stacked.back().behind = behind;
    }
    doc.arrows = std::move(stacked);
    return moved;
}
}  // namespace edit
