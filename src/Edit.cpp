#include "Edit.h"
#include "Chem.h"
#include "Geometry.h"

#include <QHash>
#include <QLineF>
#include <QRectF>
#include <QRegularExpression>
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <tuple>

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

// In bond lengths: an arrow end rests on an atom within reach of it (an end drawn beside a label
// counts), and settles this far from a labelled one.
constexpr double kReach = 0.75, kLabelGap = 0.55;

// Whichever is nearer, an atom or a bond's middle: a lone pair sits beside its atom, and an
// arrow from a π bond starts along it.
std::array<int, 2> anchorAt(const Document& doc, QPointF p) {
    std::array<int, 2> best{-1, -1};
    double r = kReach * kBondLength;
    if (int i = atomNear(doc, p, r); i >= 0) best = {i, -1}, r = len(doc.atoms[i].pos - p);
    for (const Bond& b : doc.bonds) {
        const QPointF a = doc.atoms[b.a].pos, ab = doc.atoms[b.b].pos - a;
        const double t = std::clamp(QPointF::dotProduct(p - a, ab) / std::max(1e-9, QPointF::dotProduct(ab, ab)), 0.0, 1.0);
        if (const double mid = len(a + ab / 2 - p); len(a + ab * t - p) < kBondLength / 3 && mid < r) best = {b.a, b.b}, r = mid;
    }
    return best;
}

// The atom, or the middle of the bond.
static QPointF anchorPos(const Document& doc, std::array<int, 2> at) {
    return at[1] < 0 ? doc.atoms[at[0]].pos : (doc.atoms[at[0]].pos + doc.atoms[at[1]].pos) / 2;
}


QPointF snapToAnchor(const Document& doc, std::array<int, 2> at, QPointF p, QPointF other) {
    if (at[0] < 0) return p;
    const QPointF c = anchorPos(doc, at);
    if (at[1] >= 0) return c;
    const QPointF d = len(p - c) > 0.15 * kBondLength ? p - c : other - c;
    const Atom& atom = doc.atoms[at[0]];
    const double gap = (atom.z != 6 || !atom.label.isEmpty() ? kLabelGap : 0.35) * kBondLength;  // clear of a label
    return len(d) < 1e-6 ? c : c + unit(d) * gap;
}

// An end that follows is shifted, not turned, with its atoms: an electron pair drawn
// beside an atom stays on the same side of it when only the structure is rotated.
void followAnchors(const Document& before, Document& after) {
    if (before.atoms.size() != after.atoms.size() || before.arrows.size() != after.arrows.size()) return;
    for (size_t k = 0; k < after.arrows.size(); ++k) {
        Arrow& a = after.arrows[k];
        const Arrow& was = before.arrows[k];
        if (a.fromAt != was.fromAt || a.toAt != was.toAt) continue;  // not the same arrow (restacked), or just anchored
        const double chord = len(a.to - a.from);
        const bool slid = len((a.from - was.from) - (a.to - was.to)) < 1e-6 && a.bend == was.bend;  // moved, not turned or flipped
        bool followed = false;
        for (auto [end, from, at] : {std::tuple{&a.from, was.from, &a.fromAt}, std::tuple{&a.to, was.to, &a.toAt}}) {
            if ((*at)[0] < 0) continue;
            const QPointF moved = anchorPos(after, *at) - anchorPos(before, *at);
            if (len(moved) < 1e-6) {
                if (len(*end - from) > 1e-6 && len(*end - anchorPos(after, *at)) > kReach * kBondLength)
                    *at = anchorAt(after, *end);  // the arrow moved off: onto what's there now
            } else if (slid) {
                *end = from + moved, followed = true;  // the structure moved under it (Arrange: each its own way)
            }  // else turned, flipped or scaled with its atoms: already in place
        }
        if (followed && chord > 1e-9) a.bend *= len(a.to - a.from) / chord;  // the curve keeps its shape
    }
}

std::vector<int> moleculeOf(const Document& doc, int atom) {
    const auto joined = doc.joined();
    std::vector<bool> seen(doc.atoms.size());
    std::vector<int> out{atom};
    seen[atom] = true;
    for (size_t k = 0; k < out.size(); ++k)
        for (int nb : joined[out[k]])
            if (!seen[nb]) seen[nb] = true, out.push_back(nb);
    return out;
}

static QRectF atomBox(const Document& doc, const std::vector<int>& atoms) {
    QPointF lo = doc.atoms[atoms[0]].pos, hi = lo;
    for (int a : atoms) {
        const QPointF p = doc.atoms[a].pos;
        lo = {std::min(lo.x(), p.x()), std::min(lo.y(), p.y())}, hi = {std::max(hi.x(), p.x()), std::max(hi.y(), p.y())};
    }
    return QRectF(lo, hi);
}

CompoundCount renumberCompounds(Document& doc, CompoundCount count) {
    struct Item {
        int text;
        double x, y;
    };
    std::vector<Item> items;
    for (int i = 0; i < int(doc.texts.size()); ++i)
        if (const Text& t = doc.texts[i]; t.compound)  // a row is where the molecules sit, whatever their height
            items.push_back({i, t.pos.x(), t.anchor >= 0 ? atomBox(doc, moleculeOf(doc, t.anchor)).center().y() : t.pos.y()});
    std::stable_sort(items.begin(), items.end(), [](const Item& a, const Item& b) { return a.y < b.y; });
    for (size_t row = 0; row < items.size();) {  // from its highest item to two bond lengths below it
        size_t end = row;
        while (end < items.size() && items[end].y - items[row].y < 2 * kBondLength) ++end;
        std::stable_sort(items.begin() + row, items.begin() + end, [](const Item& a, const Item& b) { return a.x < b.x; });
        row = end;
    }
    noteSeries(doc, count);
    int& next = count.last;
    for (const Item& item : items) {
        Text& t = doc.texts[item.text];
        int digits = 0;
        while (digits < t.text.size() && t.text[digits].isDigit()) ++digits;
        const QString old = t.text.left(digits), suffix = t.text.mid(digits);
        if (suffix.isEmpty()) {
            t.series = 0, t.text = QString::number(++next);
            continue;
        }
        if (!t.series) t.series = count.series.count(old) ? count.series[old] : (count.series[old] = ++count.lastSeries);
        const auto [n, fresh] = count.numbers.try_emplace(t.series, next + 1);
        if (fresh) ++next;
        t.text = QString::number(n->second) + suffix;
    }
    return count;
}

void noteSeries(const Document& doc, CompoundCount& count) {
    for (const Text& t : doc.texts) {
        if (!t.compound || !t.series) continue;
        int digits = 0;
        while (digits < t.text.size() && t.text[digits].isDigit()) ++digits;
        count.series.try_emplace(t.text.left(digits), t.series);
        count.lastSeries = std::max(count.lastSeries, t.series);
    }
}

std::vector<int> syncLegends(Document& doc) {
    if (std::none_of(doc.atoms.begin(), doc.atoms.end(), [](const Atom& a) { return !a.standsFor.isEmpty(); }) &&
        std::none_of(doc.texts.begin(), doc.texts.end(), [](const Text& t) { return t.legend; }))
        return {};  // nothing to write or erase: skip the molecule scan on every commit
    const int n = int(doc.atoms.size());
    std::vector<int> mol(n, -1);  // each atom's molecule, numbered by its first atom
    for (int i = 0; i < n; ++i)
        if (mol[i] < 0)
            for (int a : moleculeOf(doc, i)) mol[a] = i;
    auto lines = [&](int m) {
        QStringList out;
        for (int a = 0; a < n; ++a)
            if (const Atom& at = doc.atoms[a]; mol[a] == m && !at.label.isEmpty() && !at.standsFor.trimmed().isEmpty())
                out << at.label + " = " + at.standsFor.trimmed();
        out.removeDuplicates();
        out.sort();
        return out.join('\n');
    };
    std::set<int> hasLegend;
    std::erase_if(doc.texts, [&](Text& t) {  // rewritten from the definitions; gone with the last one
        if (!t.legend) return false;
        const int m = t.anchor >= 0 && t.anchor < n ? mol[t.anchor] : -1;
        if (m < 0 || hasLegend.count(m)) return true;
        t.text = lines(m);
        hasLegend.insert(m);
        return t.text.isEmpty();
    });
    std::vector<int> placed;
    for (int m = 0; m < n; ++m) {
        if (mol[m] != m || hasLegend.count(m)) continue;
        const QString text = lines(m);
        if (text.isEmpty()) continue;
        std::vector<int> atoms;
        for (int a = 0; a < n; ++a)
            if (mol[a] == m) atoms.push_back(a);
        const QRectF box = atomBox(doc, atoms);
        double y = box.bottom() + 1.3 * kBondLength;  // under the molecule, and under its compound number
        for (const Text& t : doc.texts)
            if (t.compound && t.anchor >= 0 && t.anchor < n && mol[t.anchor] == m) y = std::max(y, t.pos.y() + 1.2 * kBondLength);
        Text t{{box.center().x(), y}, text};
        t.legend = true, t.anchor = m;
        doc.texts.push_back(t);
        placed.push_back(int(doc.texts.size()) - 1);
    }
    return placed;
}

// A group holds whole molecules (#410), so an atom bonded into a grouped one (explicit H, an
// expanded label, a sprouted bond) but in no group itself joins it, as does everything joined to that.
void joinGroups(Document& doc) {
    if (std::none_of(doc.atoms.begin(), doc.atoms.end(), [](const Atom& a) { return a.group >= 0; })) return;
    const auto joined = doc.joined();
    std::vector<int> todo;
    for (int i = 0; i < int(doc.atoms.size()); ++i)
        if (doc.atoms[i].group >= 0) todo.push_back(i);
    while (!todo.empty()) {
        const int a = todo.back();
        todo.pop_back();
        for (int nb : joined[a])
            if (doc.atoms[nb].group < 0) doc.atoms[nb].group = doc.atoms[a].group, todo.push_back(nb);
    }
}

double rf(double spot, double baseline, double front) {
    return baseline - front > 1e-6 ? (baseline - spot) / (baseline - front) : std::nan("");
}

void syncPlates(Document& doc) {
    const bool any = std::any_of(doc.arrows.begin(), doc.arrows.end(), [](const Arrow& a) { return a.plate; });
    if (!any && std::none_of(doc.texts.begin(), doc.texts.end(), [](const Text& t) { return t.rf; })) return;
    std::vector<Text> labels;  // what the Rf texts should be, plate by plate
    for (const Arrow& plate : doc.arrows) {
        if (!plate.plate) continue;
        const QRectF r = QRectF(plate.from, plate.to).normalized();
        const double slack = 0.05 * r.width();
        std::vector<double> lines;  // the y of each line drawn across the plate
        for (const Arrow& a : doc.arrows) {
            const QLineF l(a.from, a.to);
            if (a.kind == ArrowKind::Line && std::abs(l.dy()) < 0.1 * std::abs(l.dx()) && std::abs(l.dx()) > 0.5 * r.width() &&
                r.adjusted(-slack, -slack, slack, slack).contains(a.from) && r.adjusted(-slack, -slack, slack, slack).contains(a.to))
                lines.push_back(l.center().y());
        }
        if (lines.size() < 2) continue;
        const auto [front, baseline] = std::minmax_element(lines.begin(), lines.end());
        for (const Arrow& a : doc.arrows) {
            const QRectF s = QRectF(a.from, a.to).normalized();
            const bool spot = (a.kind == ArrowKind::Ellipse || a.kind == ArrowKind::Box || a.kind == ArrowKind::RoundedBox) &&
                              !a.plate && r.contains(s.center()) && s.width() < 0.5 * r.width();
            if (!spot) continue;
            const double v = rf(s.center().y(), *baseline, *front);
            if (std::isnan(v)) continue;
            Text t{{s.right() + 0.15 * kBondLength, s.center().y() + 0.17 * kBondLength}, QString::number(v, 'f', 2), 0.7};
            t.rf = true, t.group = plate.group;
            labels.push_back(t);
        }
    }
    size_t k = 0;  // rewritten in place, so the other texts keep their indices
    std::erase_if(doc.texts, [&](Text& t) {
        if (!t.rf) return false;
        if (k == labels.size()) return true;
        const QColor color = t.color;  // a recoloured Rf keeps its colour
        t = labels[k++], t.color = color;
        return false;
    });
    doc.texts.insert(doc.texts.end(), labels.begin() + k, labels.end());
}

void followNumbers(const Document& before, Document& after) {
    if (before.atoms.size() != after.atoms.size() || before.texts.size() != after.texts.size()) return;
    for (size_t k = 0; k < after.texts.size(); ++k) {
        Text& t = after.texts[k];
        const Text& was = before.texts[k];
        if (t.anchor < 0 || t.anchor != was.anchor || t.pos != was.pos) continue;  // free, re-anchored or moved itself
        const auto mol = moleculeOf(after, t.anchor);
        auto foot = [&](const Document& d) {
            const QRectF box = atomBox(d, mol);
            return QPointF(box.center().x(), box.bottom());
        };
        t.pos += foot(after) - foot(before);
    }
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
    auto build = [&](QPointF from, QPointF to) {
        const QPointF d = to - from, td = t(1) - t(0);
        const double turn = qRadiansToDegrees(std::atan2(d.y(), d.x()) - std::atan2(td.y(), td.x()));
        std::vector<QPointF> verts;
        for (int k = 0; k < 6; ++k) verts.push_back(from + rotated(t(k) - t(0), turn) * (len(d) / len(td)));
        return verts;
    };
    // The far side of the bond too, for when the neighbours don't say which side is free (a bond of
    // a chair has atoms on both): keep the first that crosses or crowds nothing already drawn (#449).
    auto clashes = [&](const std::vector<QPointF>& verts) {
        auto same = [](QPointF x, QPointF y) { return len(x - y) < 1e-6; };
        int n = 0;
        for (int k = 2; k < 6; ++k)  // 0 and 1 are the bond's own atoms
            for (const Atom& a : doc.atoms)
                n += len(a.pos - verts[k]) < 0.45 * kBondLength;
        for (int k = 1; k < 6; ++k) {  // the edges after the bond
            const QPointF u = verts[k], v = verts[(k + 1) % 6];
            for (const Bond& o : doc.bonds) {
                const QPointF x = doc.atoms[o.a].pos, y = doc.atoms[o.b].pos;
                n += !same(x, u) && !same(x, v) && !same(y, u) && !same(y, v) &&
                     QLineF(u, v).intersects(QLineF(x, y), nullptr) == QLineF::BoundedIntersection;
            }
        }
        return n;
    };
    auto verts = build(pa, pb), other = build(pb, pa);
    if (clashes(other) < clashes(verts)) verts = other;
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
    for (auto& a : doc.arrows)  // a curved arrow on the dropped atom rests on the kept one
        for (int* i : {&a.fromAt[0], &a.fromAt[1], &a.toAt[0], &a.toAt[1]})
            if (*i >= 0) *i = target[*i];
    for (auto& t : doc.texts)
        if (t.anchor >= 0) t.anchor = target[t.anchor];
    auto moveOnto = [&](std::vector<int>& ids) {  // the kept atom takes the dropped one's place
        for (int& i : ids) i = target[i];
        std::sort(ids.begin(), ids.end());
        ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    };
    for (auto& b : doc.brackets) moveOnto(b.atoms);
    for (auto& a : doc.atoms) moveOnto(a.attachments);
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
    if (label.trimmed() != a.label) a.standsFor.clear();  // a definition belongs to its label (#505)
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
    // Alternatives typed as one label (N,O,S, N/O/S, [N,O,S]): one variable atom standing for them (#587), not a chain.
    static const QRegularExpression brackets("^\\[|\\]$"), separator("\\s*[,/]\\s*");
    if (const QStringList parts = label.trimmed().remove(brackets).split(separator);
        anyText && parts.size() > 1 && std::all_of(parts.begin(), parts.end(), [](const QString& p) {
            const int z = chem::atomicNumber(p.toStdString());
            return z > 0 && QString::fromStdString(chem::symbol(z)) == p;
        })) {
        a.z = 0, a.charge = 0, a.isotope = 0, a.label.clear();
        a.label = chem::freeVariableName(doc), a.standsFor = parts.join(", ");
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
            // η-bonds have no SMILES or MOL form; the centroid is a bare dummy (*) there.
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
