#include "Render.h"
#include "Chem.h"
#include "Geometry.h"

#include <QBuffer>
#include <QFile>
#include <QFileInfo>
#include <QImage>
#include <QPaintEngine>
#include <QPainter>
#include <QPicture>
#include <QSet>
#include <QPainterPath>
#include <QRegularExpression>
#include <QPdfWriter>
#include <QSaveFile>
#include <QSvgGenerator>
#include <algorithm>
#include <limits>
#include <numbers>
#include <optional>
#include <tuple>
#include <utility>
#include <numeric>

// Presets. Values come from the ChemDraw stationery (.cds) of the same name;
// wedge width, hash spacing and label gap keep ACS's proportions to ours.
// Model units per point for a style's native values (bond length / 14.4 undone).
constexpr double kJdp = kBondLength / 14.17, kRsc = kBondLength / 12.2;

const std::vector<DrawingStyle>& drawingStyles() {
    // Wedge width and hash spacing keep ACS's proportion to the bold width and
    // hash setting (4.5/2.0 and 2.2/2.5); label radius is half a glyph plus the margin.
    static const std::vector<DrawingStyle> styles{
        {"ACS 1996", 14.4, 0.6, 2.0, 4.5, 2.2, 0.18, 5.5, "Arial", QFont::Normal, 10},
        // JDPReport.cds: bond 14.17 (0.5 cm), line 0.879, bold 1.814, hash 1.814, margin 1.162,
        // IBM Plex Sans Light 10 pt.
        {"JDP", 14.17, 0.879 * kJdp, 1.814 * kJdp, 4.5 * 1.814 / 2.0 * kJdp, 2.2 * 1.814 / 2.5 * kJdp, 0.18,
         (3.9 + 1.162) * kJdp, "IBM Plex Sans", QFont::Light, 10 * kJdp},
        // RSC (1 Column).cds (2 Column only differs in page size): bond 12.2, line 0.449, bold 1.602,
        // hash 1.75, margin 1.25, spacing 20%, Helvetica 7 pt.
        {"RSC", 12.2, 0.449 * kRsc, 1.602 * kRsc, 4.5 * 1.602 / 2.0 * kRsc, 2.2 * 1.75 / 2.5 * kRsc, 0.20,
         (3.9 * 0.7 + 1.25) * kRsc, "Helvetica", QFont::Normal, 7 * kRsc},
    };
    return styles;
}

const DrawingStyle& drawingStyle(const QString& name) {
    for (const auto& s : drawingStyles())
        if (s.name == name) return s;
    return drawingStyles()[0];
}
// Carbons are skeletal unless alone or the document asks for their labels.
static bool hasLabel(const Document& doc, int i, const std::vector<int>& degree) {
    if (doc.atoms[i].z == 0 && doc.atoms[i].label.isEmpty()) return false;  // attachment point, η centroid: a bare point
    if (doc.atoms[i].z != 6 || degree[i] == 0 || !doc.atoms[i].label.isEmpty() || doc.atoms[i].isotope) return true;
    return doc.carbonLabels == Document::CarbonLabels::All ||
           (doc.carbonLabels == Document::CarbonLabels::Terminal && degree[i] == 1);
}

// ---------------------------------------------------------------- rendering

DrawingStyle documentStyle(const Document& doc) {
    DrawingStyle st = drawingStyle(doc.style);
    if (doc.labelRatio > 0) {  // text and the bonds' clearance around labels scale together
        const double k = doc.labelRatio * kBondLength / st.fontSize;
        st.fontSize *= k, st.labelRadius *= k;
    }
    return st;
}

QFont labelFont(const DrawingStyle& s, double scale) {
    QFont f(s.font);
    f.setWeight(s.weight);
    f.setPixelSize(int(s.fontSize * scale));  // 1 px == 1 pt in our scene units
    return f;
}

// Text as outlines, so it scales identically on screen, SVG, PDF and PNG.
static void drawText(QPainter& p, const QString& s, QPointF baselineLeft, const QFont& f) {
    QPainterPath path;
    path.addText(baselineLeft, f, s);
    p.fillPath(path, p.pen().color());
}

// An isotope's mass number, written before the symbol (¹³C); none for D and T, which have letters.
static QString massNumber(const Atom& a) {
    return a.isotope && a.label.isEmpty() && !(a.z == 1 && (a.isotope == 2 || a.isotope == 3)) ? QString::number(a.isotope) : QString();
}

// Where the mass number sits, relative to the atom (as drawLabel puts it), with a point to spare.
static QRectF massBox(const Atom& a, const DrawingStyle& st) {
    QFontMetricsF fm(labelFont(st)), sm(labelFont(st, 0.7));
    const double symbol = fm.horizontalAdvance(QString::fromStdString(chem::symbol(a.z))), cap = fm.capHeight();
    const double right = -symbol / 2, baseline = cap / 2 - 0.7 * cap;
    return QRectF(QPointF(right - sm.horizontalAdvance(massNumber(a)), baseline - sm.capHeight()), QPointF(right, baseline))
        .adjusted(-1, -1, 0.5, 0.5);
}

// How far along the unit vector `u` a ray from the origin has left `box`; 0 if it misses it.
static double exitAlong(const QRectF& box, QPointF u) {
    double lo = 0, hi = std::numeric_limits<double>::infinity();
    for (auto [o, a0, a1] : {std::tuple{u.x(), box.left(), box.right()}, std::tuple{u.y(), box.top(), box.bottom()}}) {
        if (std::abs(o) < 1e-9) {
            if (a0 > 0 || a1 < 0) return 0;
            continue;
        }
        const double t0 = std::min(a0 / o, a1 / o), t1 = std::max(a0 / o, a1 / o);
        lo = std::max(lo, t0), hi = std::min(hi, t1);
    }
    return lo <= hi ? hi : 0;
}

// Abbreviation written from the right, bond side last: OMe -> MeO.
static QString reversedLabel(const QString& s) {
    static const QHash<QString, QString> r{
        {"OMe", "MeO"}, {"CO2Me", "MeO2C"}, {"CO2Et", "EtO2C"}, {"CO2H", "HO2C"}, {"NO2", "O2N"},
        {"CF3", "F3C"}, {"OAc", "AcO"},     {"CHO", "OHC"},     {"SO2Me", "MeO2S"}, {"OTf", "TfO"},
        {"OTs", "TsO"}, {"OTBS", "TBSO"},   {"CN", "NC"},       {"Bpin", "pinB"},
    };
    return r.value(s, s);
}

// Abbreviation centred on its attaching letter (the first, or the last when written from the right).
static void drawAbbreviation(QPainter& p, const Atom& a, bool fromRight, const DrawingStyle& st) {
    QFontMetricsF fm(labelFont(st));
    const QString s = fromRight ? reversedLabel(a.label) : a.label;
    const double base = a.pos.y() + fm.capHeight() / 2;
    QPainterPath path = textPath({{0, base}, s}, st);
    const double x = fromRight ? a.pos.x() + fm.horizontalAdvance(s.back()) / 2 - path.boundingRect().right()
                               : a.pos.x() - fm.horizontalAdvance(s.front()) / 2;
    p.fillPath(path.translated(x, 0), p.pen().color());
    // A charge beyond the group's own, as typed (N3-, #370), raised after the label like an element's.
    const auto head = chem::abbreviationHead(a.label);
    if (const int extra = a.charge - (head ? head->charge : 0)) {
        QString c = QString(extra > 0 ? "+" : "−");
        if (std::abs(extra) > 1) c.prepend(QString::number(std::abs(extra)));
        drawText(p, c, {x + path.boundingRect().right(), base - fm.capHeight() * 0.7}, labelFont(st, 0.7));
    }
}

// Where a label's implicit H goes: after the symbol, before it, or stacked
// below/above when bonds leave neither side free (a middle CH2).
enum class HSide { Right, Left, Below, Above };

// Which way an atom's charge sits from it (y grows down): none for the usual upper right (after
// an H on the right), unless a bond or an H above comes within 50° of it (N=O of a nitro, R4N+,
// =N+=), then the middle of the widest gap that isn't inside a ring, ties to the upper right (#494).
static std::optional<QPointF> chargeDirection(const Document& doc, int i, int hydrogens, HSide side) {
    constexpr double upRight = -std::numbers::pi / 4;
    if (hydrogens && side == HSide::Right) return std::nullopt;
    std::vector<std::pair<double, int>> taken;  // (angle, neighbour or -1)
    for (const Bond& b : doc.bonds)
        if (b.a == i || b.b == i) {
            const int nb = b.a == i ? b.b : b.a;
            const QPointF d = doc.atoms[nb].pos - doc.atoms[i].pos;
            taken.push_back({std::atan2(d.y(), d.x()), nb});
        }
    if (hydrogens) taken.push_back({side == HSide::Left ? std::numbers::pi : side == HSide::Below ? std::numbers::pi / 2 : -std::numbers::pi / 2, -1});
    if (hydrogens > 1 && (side == HSide::Above || side == HSide::Below))  // the H's count, down or up to its right
        taken.push_back({side == HSide::Above ? -std::numbers::pi / 3 : std::numbers::pi / 3, -1});
    const auto away = [](double x, double y) { return std::abs(std::remainder(x - y, 2 * std::numbers::pi)); };
    if (std::all_of(taken.begin(), taken.end(), [&](auto t) { return away(t.first, upRight) > 50 * std::numbers::pi / 180; }))
        return std::nullopt;
    if (!massNumber(doc.atoms[i]).isEmpty()) taken.push_back({-3 * std::numbers::pi / 4, -1});
    // Whether neighbours u and v close a ring of up to 8 through atom i: the gap between them is its inside.
    const auto ring = [&](int u, int v) {
        std::vector<int> seen{i, u}, front{u};
        for (int step = 0; step < 6 && !front.empty(); ++step) {
            std::vector<int> next;
            for (int x : front)
                for (const Bond& b : doc.bonds) {
                    const int y = b.a == x ? b.b : b.b == x ? b.a : -1;
                    if (y == v) return true;
                    if (y >= 0 && std::find(seen.begin(), seen.end(), y) == seen.end()) seen.push_back(y), next.push_back(y);
                }
            front = std::move(next);
        }
        return false;
    };
    std::sort(taken.begin(), taken.end());
    double best = upRight, score = -1e9;
    for (size_t k = 0; k < taken.size(); ++k) {
        const auto [from, u] = taken[k];
        const auto [to, v] = taken[(k + 1) % taken.size()];
        const double width = (k + 1 < taken.size() ? to : to + 2 * std::numbers::pi) - from;
        // As near the upper right as the gap allows, 60° clear of its edges: opposite a lone bond it reads as a bond (−O–N⁺).
        const double clear = std::min(width / 2, std::numbers::pi / 3);
        const double at = std::clamp(from + std::fmod(upRight - from + 4 * std::numbers::pi, 2 * std::numbers::pi), from + clear, from + width - clear);
        const double sc = width - 0.25 * away(at, upRight) - (u >= 0 && v >= 0 && width < std::numbers::pi && ring(u, v) ? 10 : 0);
        if (sc > score) best = at, score = sc;
    }
    return QPointF(std::cos(best), std::sin(best));
}

// Where a charge `c` in font `sub` is drawn from, its glyph centred just past `clear` along `dir`.
static QPointF chargeAt(QPointF pos, QPointF dir, double clear, const QString& c, const QFont& sub) {
    const QFontMetricsF sm(sub);
    const QRectF glyph(-sm.horizontalAdvance(c) / 2, -sm.capHeight() / 2, sm.horizontalAdvance(c), sm.capHeight());
    return pos + dir * (clear + exitAlong(glyph, dir) + 0.8) + QPointF(glyph.left(), -glyph.top());
}

static void drawLabel(QPainter& p, const Document& doc, int i, int hydrogens, HSide side, const DrawingStyle& st) {
    const bool hLeft = side == HSide::Left;
    const auto& a = doc.atoms[i];
    if (!a.label.isEmpty()) return drawAbbreviation(p, a, hLeft, st);
    QFont f = labelFont(st), sub = labelFont(st, 0.7);
    QFontMetricsF fm(f), sm(sub);
    // Deuterium and tritium by their own letters; any other isotope as a mass number before the symbol (¹³C).
    const bool heavyH = a.z == 1 && (a.isotope == 2 || a.isotope == 3);
    QString sym = heavyH ? (a.isotope == 2 ? "D" : "T") : QString::fromStdString(chem::symbol(a.z));
    const QString mass = massNumber(a);
    double w = fm.horizontalAdvance(sym);
    double base = a.pos.y() + fm.capHeight() / 2;
    double x = a.pos.x() - w / 2;
    drawText(p, sym, {x, base}, f);
    double left = x;
    if (!mass.isEmpty()) left -= sm.horizontalAdvance(mass), drawText(p, mass, {left, base - fm.capHeight() * 0.7}, sub);

    double right = x + w;
    if (hydrogens > 0) {
        QString n = hydrogens > 1 ? QString::number(hydrogens) : QString();
        double hw = fm.horizontalAdvance("H"), nw = sm.horizontalAdvance(n);
        if (side == HSide::Below || side == HSide::Above) {
            const double hx = a.pos.x() - (hw + nw) / 2, dy = fm.capHeight() * 1.35 * (side == HSide::Below ? 1 : -1);
            drawText(p, "H", {hx, base + dy}, f);
            if (!n.isEmpty()) drawText(p, n, {hx + hw, base + dy + fm.capHeight() * 0.35}, sub);
        } else {
            double hx = hLeft ? left - hw - nw : right;
            drawText(p, "H", {hx, base}, f);
            if (!n.isEmpty()) drawText(p, n, {hx + hw, base + fm.capHeight() * 0.35}, sub);
            if (!hLeft) right += hw + nw;
        }
    }
    if (a.charge) {
        QString c = QString(a.charge > 0 ? "+" : "−");
        if (std::abs(a.charge) > 1) c.prepend(QString::number(std::abs(a.charge)));
        if (const auto d = chargeDirection(doc, i, hydrogens, side))
            drawText(p, c, chargeAt(a.pos, *d, exitAlong(QRectF(-w / 2, -fm.capHeight() / 2, w, fm.capHeight()), *d), c, sub), sub);
        else drawText(p, c, {right, base - fm.capHeight() * 0.7}, sub);  // after any H
    }
}

constexpr double kWedgeReach = 1.3;  // how far a wedge's wide corner may stretch to meet a neighbouring bond
static void drawBond(QPainter& p, const Document& doc, const Bond& b, const DrawingStyle& st, const std::vector<int>& degree,
                     const std::vector<bool>& labeled, const BondsAt& at, const std::vector<QPointF>& gaps = {}) {
    QPointF pa = doc.atoms[b.a].pos, pb = doc.atoms[b.b].pos;
    QPointF d = unit(pb - pa), n = perp(d);
    const double gap = st.bondSpacing * kBondLength;  // double-bond spacing
    // Trim at labels.
    // A bond that would cross an isotope's raised mass number (¹³C) stops beyond it; a level one passes beneath.
    auto trim = [&](int atom, QPointF toward) {
        const QString mass = massNumber(doc.atoms[atom]);
        return mass.isEmpty() ? st.labelRadius : std::max(st.labelRadius, 1 + exitAlong(massBox(doc.atoms[atom], st), toward));
    };
    QPointF a = labeled[b.a] ? pa + d * trim(b.a, d) : pa;
    QPointF e = labeled[b.b] ? pb - d * trim(b.b, -d) : pb;
    if (QPointF::dotProduct(e - a, d) <= 0) return;  // labels so close that the trimmed ends cross (#496)

    // A wedge's wide end at a bare atom lies along the bonds beside it rather than jutting past them.
    // Each side runs to the drawn line of the bond nearest it on that side (a double bond's nearer
    // line, so it doesn't overlap it); a side with no bond of its own follows the other side's bond, inward only.
    QPointF corner[2] = {e + n * st.wedgeWidth / 2, e - n * st.wedgeWidth / 2};
    if ((b.stereo == BondStereo::Wedge || b.stereo == BondStereo::Hash) && !labeled[b.b]) {
        std::optional<QPointF> dir[2];  // each side's bond direction,
        QPointF on[2];                  // and a point on its nearer line
        for (int s : {0, 1}) {
            double back = -2;
            for (int k : at[b.b]) {
                const Bond& o = doc.bonds[k];
                const int nb = o.a == b.b ? o.b : o.a;
                const QPointF v = unit(doc.atoms[nb].pos - pb);
                if (nb == b.a || QPointF::dotProduct(v, n) * (s ? -1 : 1) <= 1e-6 || QPointF::dotProduct(v, -d) <= back) continue;
                back = QPointF::dotProduct(v, -d);
                QPointF q = perp(v);
                if (QPointF::dotProduct(q, corner[s] - pb) < 0) q = -q;  // toward this corner
                double off = o.order == 3 ? gap : 0;
                if (o.order == 2) {  // centred lines sit gap/2 either side; an offset one, gap to one side
                    const int side = doubleBondSide(doc, o, at);
                    const QPointF second = perp(unit(doc.atoms[o.b].pos - doc.atoms[o.a].pos)) * (side >= 0 ? 1 : -1);
                    off = side == 0 ? gap / 2 : QPointF::dotProduct(second, q) > 0 ? gap : 0;
                }
                dir[s] = v, on[s] = pb + q * off;
            }
        }
        for (int s : {0, 1}) {
            const bool own = bool(dir[s]);
            const std::optional<QPointF> u = own ? dir[s] : dir[1 - s];
            const QPointF side = corner[s] - a;
            if (!u || std::abs(cross(side, *u)) < 1e-9) continue;
            const double t = cross((own ? on[s] : pb) - a, *u) / cross(side, *u);
            if (t > 0.5 && t < (own ? kWedgeReach : 1)) corner[s] = a + side * t;
        }
    }
    if (b.stereo == BondStereo::Wedge) {
        p.setBrush(p.pen().color());
        p.drawPolygon(QPolygonF{a, corner[0], e, corner[1]});
        p.setBrush(Qt::NoBrush);
        return;
    }
    if (b.stereo == BondStereo::Hash) {
        double L = len(e - a);
        int count = std::max(3, int(L / st.hashSpacing));
        for (int k = 0; k <= count; ++k) {
            double t = double(k) / count;
            p.drawLine(a + (corner[0] - a) * t, a + (corner[1] - a) * t);
        }
        return;
    }

    if (b.stereo == BondStereo::Interaction) {  // dotted, whatever the order
        QPen dots = p.pen();
        dots.setCapStyle(Qt::RoundCap);
        dots.setDashPattern({0.01, 2.2});
        p.setPen(dots);
        p.drawLine(a, e);
        return;
    }
    if (b.stereo == BondStereo::Wavy) {
        QPainterPath wave(a);
        const double L = len(e - a);
        const int bumps = std::max(2, int(std::round(L / 3)));
        for (int k = 0; k < bumps; ++k) {
            QPointF from = a + (e - a) * (double(k) / bumps), to = a + (e - a) * (double(k + 1) / bumps);
            wave.quadTo((from + to) / 2 + n * (k % 2 ? -2.0 : 2.0) * 1.3, to);
        }
        p.drawPath(wave);
        return;
    }
    // Bold styles the main line, Dashed and Partial the other one (or the only one).
    const QPen pen = p.pen();
    auto line = [&](QPointF x, QPointF y, bool main) {
        QPen q = pen;
        if (b.stereo == BondStereo::Bold && main) q.setWidthF(st.boldWidth), q.setCapStyle(Qt::FlatCap);
        const bool dashed = b.stereo == BondStereo::Dashed || b.stereo == BondStereo::Partial;
        if (dashed && (!main || b.order == 1)) q.setDashPattern({2.5, 2.5});
        p.setPen(q);
        // A bond drawn later (in front) crossing this one leaves a gap in it.
        std::vector<double> cuts;  // along x→y, 0..1
        const double length = len(y - x), half = 0.18 * kBondLength;
        for (QPointF g : gaps) cuts.push_back(QPointF::dotProduct(g - x, y - x) / (length * length));
        std::sort(cuts.begin(), cuts.end());
        double from = 0;
        for (double c : cuts) {
            const double to = c - half / length;
            if (to > from) p.drawLine(x + (y - x) * from, x + (y - x) * to);
            from = std::max(from, c + half / length);
        }
        if (from < 1) p.drawLine(x + (y - x) * from, y);
        p.setPen(pen);
    };

    if (b.order == 1) {
        line(a, e, true);
    } else if (b.order == 3) {
        line(a, e, true);
        line(a + n * gap, e + n * gap, true);
        line(a - n * gap, e - n * gap, true);
    } else {
        // Offset the second line toward the side where the neighbours are
        // (inside the ring); centre it for terminal bonds like C=O.
        const int side = doubleBondSide(doc, b, at);
        if (side == 0) {
            QPointF o = n * gap / 2;
            line(a + o, e + o, true);
            line(a - o, e - o, false);
        } else {
            QPointF o = n * (side >= 0 ? gap : -gap);
            QPointF shrink = d * (0.15 * kBondLength);
            QPointF ia = labeled[b.a] ? a : a + shrink, ie = labeled[b.b] ? e : e - shrink;
            line(a, e, true);
            line(ia + o, ie + o, false);
        }
    }

    // A bold line ends flat, so at a skeletal atom fill the bevel to each other bond there, as a pen
    // joins a polyline: no notch at ring corners (#217). Labels keep their clearance.
    if (b.stereo != BondStereo::Bold || (b.order == 2 && doubleBondSide(doc, b, at) == 0) || b.order == 3) return;
    p.save();
    p.setPen(Qt::NoPen);
    p.setBrush(pen.color());
    for (int end : {b.a, b.b}) {
        if (labeled[end]) continue;
        const QPointF c = doc.atoms[end].pos, out = unit(doc.atoms[end == b.a ? b.b : b.a].pos - c);
        for (int other : at[end]) {
            const Bond& o = doc.bonds[other];
            if (std::minmax(o.a, o.b) == std::minmax(b.a, b.b)) continue;  // this bond (b may be a copy)
            const QPointF dir = unit(doc.atoms[o.a == end ? o.b : o.a].pos - c);
            const double w = st.boldWidth / 2, w2 = o.stereo == BondStereo::Bold ? w : pen.widthF() / 2;
            std::vector<QPointF> quad{c + perp(out) * w, c - perp(out) * w, c + perp(dir) * w2, c - perp(dir) * w2};
            std::sort(quad.begin(), quad.end(), [&](QPointF x, QPointF y) {
                return std::atan2(x.y() - c.y(), x.x() - c.x()) < std::atan2(y.y() - c.y(), y.x() - c.x());
            });
            p.drawPolygon(QPolygonF(QList<QPointF>(quad.begin(), quad.end())));
        }
    }
    p.restore();
}

// ---- arrows and text

constexpr double kHeadWidth = 2.2, kEquilibriumGap = 1.6;

// Quadratic control point: puts the curve's midpoint `bend` to the left of from->to.
// Points along the arrow. A curve is a circular arc through both ends whose
// midpoint sits `bend` off the chord, so arcs past 180° (ChemDraw's cycle
// arrows) are exact too.
static std::vector<QPointF> arrowPoints(const Arrow& a) {
    if (!a.bend || len(a.to - a.from) < 1e-6) return {a.from, a.to};
    const QPointF mid = (a.from + a.to) / 2, u = -perp(unit(a.to - a.from)) * (a.bend > 0 ? 1 : -1);
    const double h = len(a.to - a.from) / 2, sag = std::abs(a.bend), radius = (h * h + sag * sag) / (2 * sag);
    const QPointF top = mid + u * sag, c = top - u * radius;
    auto angle = [&](QPointF p) { return std::atan2(p.y() - c.y(), p.x() - c.x()); };
    const double from = angle(a.from), through = angle(top);
    double sweep = angle(a.to) - from;
    // Go the way that passes through the arc's midpoint.
    auto wrap = [](double x) { return std::remainder(x, 2 * M_PI); };
    const double half = wrap(through - from);
    if (half > 0 && sweep < 0) sweep += 2 * M_PI;
    if (half < 0 && sweep > 0) sweep -= 2 * M_PI;
    // About 1.5 pt per segment, capped so a huge bend from a file can't exhaust memory (#314).
    const int n = std::clamp(int(std::min(std::abs(sweep) * radius / 1.5, 720.0)), 8, 720);
    std::vector<QPointF> pts;
    for (int k = 0; k <= n; ++k) {
        double t = from + sweep * k / n;
        pts.push_back(c + QPointF(std::cos(t), std::sin(t)) * radius);
    }
    pts.front() = a.from, pts.back() = a.to;
    return pts;
}

// An orbital's lobes, each with its phase (true: the main one). A lobe is a
// slim teardrop from the node to a rounded tip.
static std::vector<std::pair<QPainterPath, bool>> orbitalLobes(const Arrow& a) {
    const double L = len(a.to - a.from);
    if (a.kind == ArrowKind::SOrbital) {
        QPainterPath s;
        s.addEllipse(a.from, L, L);
        return {{s, true}};
    }
    auto lobe = [&](double length, double width, bool back) {
        QPainterPath l;
        l.moveTo(0, 0);
        l.cubicTo(0.35 * length, width, length, 0.85 * width, length, 0);
        l.cubicTo(length, -0.85 * width, 0.35 * length, -width, 0, 0);
        QTransform t;
        t.translate(a.from.x(), a.from.y());
        t.rotateRadians(std::atan2(a.to.y() - a.from.y(), a.to.x() - a.from.x()) + (back ? M_PI : 0));
        return t.map(l);
    };
    std::vector<std::pair<QPainterPath, bool>> lobes{{lobe(L, 0.42 * L, false), true}};
    if (a.kind == ArrowKind::POrbital) lobes.push_back({lobe(L, 0.42 * L, true), false});
    if (a.kind == ArrowKind::HybridOrbital) lobes.push_back({lobe(0.4 * L, 0.2 * L, true), false});
    return lobes;
}

QPainterPath arrowPath(const Arrow& a) {
    if (isOrbital(a.kind)) {
        QPainterPath path;
        for (const auto& [lobe, main] : orbitalLobes(a)) path.addPath(lobe);
        return path;
    }
    if (isShape(a.kind) && a.kind != ArrowKind::Line) {
        const QRectF r = QRectF(a.from, a.to).normalized();
        QPainterPath path;
        if (a.kind == ArrowKind::Box) path.addRect(r);
        else if (a.kind == ArrowKind::RoundedBox) path.addRoundedRect(r, 0.35 * kBondLength, 0.35 * kBondLength);
        else path.addEllipse(r);
        return path;
    }
    const auto pts = arrowPoints(a);
    QPainterPath path(pts[0]);
    for (size_t k = 1; k < pts.size(); ++k) path.lineTo(pts[k]);
    return path;
}

// Drops `by` points' worth of length from the end of a polyline.
static void trimEnd(std::vector<QPointF>& pts, double by) {
    while (pts.size() > 2 && len(pts.back() - pts[pts.size() - 2]) <= by)
        by -= len(pts.back() - pts[pts.size() - 2]), pts.pop_back();
    pts.back() -= unit(pts.back() - pts[pts.size() - 2]) * std::min(by, len(pts.back() - pts[pts.size() - 2]) - 0.01);
}

// Where a polyline is `by` back from its end, measured along it.
static QPointF pointBack(const std::vector<QPointF>& pts, double by) {
    for (size_t i = pts.size() - 1; i > 0; --i) {
        const double seg = len(pts[i] - pts[i - 1]);
        if (seg >= by) return pts[i] - unit(pts[i] - pts[i - 1]) * by;
        by -= seg;
    }
    return pts.front();
}

// Filled head at `tip` pointing along `dir`; `sides` +1/-1 for a half head (one barb, and no
// sliver back along the shaft).
// Heads are solid even on a dashed arrow.
static QPen solid(QPen pen) {
    pen.setStyle(Qt::SolidLine);
    return pen;
}

static void drawHead(QPainter& p, QPointF tip, QPointF dir, int sides, double size) {
    const double length = kHeadLength * size, width = kHeadWidth * size;
    QPointF d = unit(dir), n = perp(d), base = tip - d * length, notch = tip - d * (length * 0.8);
    QPolygonF head = sides > 0   ? QPolygonF{tip, base + n * width, notch}
                     : sides < 0 ? QPolygonF{tip, notch, base - n * width}
                                 : QPolygonF{tip, base + n * width, notch, base - n * width};
    const QPen shaft = p.pen();
    p.setPen(solid(shaft));
    p.setBrush(shaft.color());
    p.drawPolygon(head);
    p.setBrush(Qt::NoBrush);
    p.setPen(shaft);
}

static QColor mix(const QColor& c, double white) {
    return QColor::fromRgbF(c.redF() + (1 - c.redF()) * white, c.greenF() + (1 - c.greenF()) * white,
                            c.blueF() + (1 - c.blueF()) * white, c.alphaF());
}

// The Enhanced Metafile engine (Windows copy, #394); GDI has no gradients either.
constexpr auto kEmfEngine = QPaintEngine::Type(QPaintEngine::User + 1);

// Each phase filled for its look, then outlined.
static void drawOrbital(QPainter& p, const Arrow& a, const QColor& color, double lineWidth) {
    QPainterPathStroker outline;
    outline.setWidth(lineWidth);
    outline.setJoinStyle(Qt::RoundJoin);
    const QPointF light(-1, -1);  // highlights toward the top left
    for (const auto& [lobe, main] : orbitalLobes(a)) {
        const QRectF r = lobe.boundingRect();
        const double radius = 0.5 * std::max(r.width(), r.height());
        if (a.look == OrbitalLook::Shaded && main) {
            p.fillPath(lobe, color);
        } else if (a.look == OrbitalLook::Gradient) {
            const QColor inner = main ? mix(color, 0.75) : QColor(Qt::white), outer = main ? color : mix(color, 0.7);
            const QPointF focus = r.center() + light * 0.3 * radius;
            if (p.paintEngine() && (p.paintEngine()->type() == QPaintEngine::Pdf || p.paintEngine()->type() == kEmfEngine)) {
                // Qt writes gradient fills with an uncoloured pattern colour space, which
                // Apple's PDF renderer (Preview, Keynote, Word on macOS) skips; so PDFs get 16 vector
                // bands instead. Drop this once Qt writes a plain /Pattern colour space.
                constexpr int kBands = 16;
                p.fillPath(lobe, outer);
                for (int k = kBands - 1; k > 0; --k) {
                    const double t = double(k) / kBands;
                    QPainterPath band;
                    band.addEllipse(focus + (r.center() - focus) * t, radius * t, radius * t);
                    p.fillPath(lobe.intersected(band), QColor::fromRgbF(inner.redF() + (outer.redF() - inner.redF()) * t,
                                                                         inner.greenF() + (outer.greenF() - inner.greenF()) * t,
                                                                         inner.blueF() + (outer.blueF() - inner.blueF()) * t));
                }
            } else {
                QRadialGradient g(r.center(), radius, focus);
                g.setColorAt(0, inner);
                g.setColorAt(1, outer);
                p.fillPath(lobe, g);
            }
        }
        p.fillPath(outline.createStroke(lobe), color);
    }
}

static void drawArrow(QPainter& p, const Arrow& a) {
    if (a.dashed) {
        QPen dashes = p.pen();
        dashes.setDashPattern({3, 2.5});
        p.setPen(dashes);
    }
    if (isShape(a.kind)) {
        p.drawPath(arrowPath(a));
        return;
    }
    QPointF d = unit(a.to - a.from), n = perp(d);
    const double headLength = kHeadLength * a.head, headWidth = kHeadWidth * a.head;
    if (a.kind == ArrowKind::Equilibrium) {  // ⇌: two half-headed lines
        QPointF o = n * kEquilibriumGap;
        p.drawLine(a.from - o, a.to - o - d * headLength * 0.8);
        drawHead(p, a.to - o, d, -1, a.head);
        p.drawLine(a.to + o, a.from + o + d * headLength * 0.8);
        drawHead(p, a.from + o, -d, -1, a.head);
        return;
    }
    if (a.kind == ArrowKind::Retro) {  // ⇒: open double arrow
        QPointF o = n * kEquilibriumGap, back = a.to - d * headLength;
        p.drawLine(a.from + o, back + o + d * kEquilibriumGap);
        p.drawLine(a.from - o, back - o + d * kEquilibriumGap);
        p.setPen(solid(p.pen()));
        p.drawPolyline(QPolygonF{back + n * (headWidth + kEquilibriumGap), a.to, back - n * (headWidth + kEquilibriumGap)});
        return;
    }
    // Heads follow the tangent at each end; the shaft stops inside them so it
    // doesn't poke through the tip.
    // A head is aimed along the chord it covers, not the tangent at the tip: on a tight curve the
    // tangent turns the head off the shaft, which then leaves through one barb.
    std::vector<QPointF> pts = arrowPoints(a);
    const QPointF endDir = pts.back() - pointBack(pts, headLength * 0.8);
    const QPointF startDir = pts.front() - pointBack(std::vector<QPointF>(pts.rbegin(), pts.rend()), headLength * 0.8);
    trimEnd(pts, headLength * 0.7);
    if (a.kind == ArrowKind::Resonance) {
        std::reverse(pts.begin(), pts.end());
        trimEnd(pts, headLength * 0.7);
    }
    p.drawPolyline(pts.data(), int(pts.size()));
    // Fishhook: the barb sits on the outside of the curve.
    drawHead(p, a.to, endDir, a.kind == ArrowKind::Fishhook ? (a.bend >= 0 ? -1 : 1) : 0, a.head);
    if (a.kind == ArrowKind::Resonance) drawHead(p, a.from, startDir, 0, a.head);
}

// Formula-style text, as chemists type it (and mhchem reads it): digits after a letter or
// bracket are counts (H2O); a + or - ending a formula is its charge (NH4+, OH-), with the
// digits before it too for one element or a bracket (Cu2+, [Fe(CN)6]3-). ^ marks a charge
// outright where it would be ambiguous (SO4^2-); the ^ isn't drawn.
std::vector<Script> scripts(const QString& s) {
    const int n = int(s.size());
    std::vector<Script> out(n, Script::Base);
    for (int i = 1; i < n; ++i)
        if (s[i].isDigit() && (s[i - 1].isLetter() || s[i - 1] == ')' || s[i - 1] == ']' || out[i - 1] == Script::Sub))
            out[i] = Script::Sub;
    auto sign = [&](int i) { return s[i] == '+' || s[i] == '-' || s[i] == QChar(0x2212); };
    for (int i = 0; i < n; ++i) {
        if (s[i] == '^') {
            int j = i + 1;
            while (j < n && s[j].isDigit()) ++j;
            while (j < n && sign(j)) ++j;
            if (j == i + 1) continue;
            out[i] = Script::Hidden;
            std::fill(out.begin() + i + 1, out.begin() + j, Script::Super);
            i = j - 1;
            continue;
        }
        static const QString after(" \t,;:/)"), before(" \t,;:/(");
        if (!sign(i) || i == 0 || (i + 1 < n && !after.contains(s[i + 1]))) continue;
        if (!(s[i - 1].isLetterOrNumber() || s[i - 1] == ')' || s[i - 1] == ']')) continue;
        int w = i;  // the formula's first character
        while (w > 0 && !before.contains(s[w - 1])) --w;
        if (!(s[w].isUpper() || s[w] == '[')) continue;  // "cis-", "-78", "2e-": not formulas
        out[i] = Script::Super;
        int d = i;
        while (d > w && s[d - 1].isDigit()) --d;
        static const QRegularExpression element("^[A-Z][a-z]?$");
        if (d < i && (element.match(s.mid(w, d - w)).hasMatch() || (d > w && s[d - 1] == ']')))
            std::fill(out.begin() + d, out.begin() + i, Script::Super);
    }
    return out;
}

// Text as outlines, formula-style subscripts and charges, one line per '\n'. Laid out
// in runs (not per letter) so kerning and spaces match ordinary text.
QPainterPath textPath(const Text& t, const DrawingStyle& st) {
    QFont f = labelFont(st, t.scale), small = labelFont(st, 0.7 * t.scale);
    f.setBold(t.compound), small.setBold(t.compound);
    QFontMetricsF fm(f), sm(small);
    const double tab = kTabSpaces * fm.horizontalAdvance(' ');
    QPainterPath path;
    const auto lines = t.text.split('\n');
    for (int li = 0; li < lines.size(); ++li) {
        const QString& s = lines[li];
        const auto sc = scripts(s);
        double x = 0, y = t.pos.y() + li * fm.lineSpacing();
        for (int i = 0; i < s.size();) {
            if (s[i] == '\t') {
                x = (std::floor(x / tab + 1e-6) + 1) * tab;
                ++i;
                continue;
            }
            if (s[i] == ' ') {  // by hand: some platforms drop leading spaces from a shaped run
                x += fm.horizontalAdvance(' ');
                ++i;
                continue;
            }
            if (sc[i] == Script::Hidden) {
                ++i;
                continue;
            }
            int j = i + 1;
            while (j < s.size() && s[j] != '\t' && s[j] != ' ' && sc[j] == sc[i]) ++j;
            QString run = s.mid(i, j - i);
            if (sc[i] == Script::Super) run.replace('-', QChar(0x2212));  // a charge's minus, not a hyphen
            const bool base = sc[i] == Script::Base;
            const double dy = sc[i] == Script::Sub ? fm.capHeight() * 0.35 : sc[i] == Script::Super ? -fm.capHeight() * 0.7 : 0;
            path.addText(t.pos.x() + x, y + dy, base ? f : small, run);
            x += (base ? fm : sm).horizontalAdvance(run);
            i = j;
        }
    }
    return path;
}

void paintDocument(QPainter& p, const Document& doc, const RenderStyle& style) {
    const DrawingStyle st = documentStyle(doc);
    const double lineWidth = style.lineWidth > 0 ? style.lineWidth : st.lineWidth;
    p.save();
    p.setRenderHint(QPainter::Antialiasing);
    const BondsAt bondsAt = doc.bondsAt();
    std::vector<int> degree(doc.atoms.size(), 0);
    for (const auto& b : doc.bonds) ++degree[b.a], ++degree[b.b];
    std::vector<bool> labeled(doc.atoms.size());
    for (size_t i = 0; i < doc.atoms.size(); ++i) labeled[i] = hasLabel(doc, int(i), degree);
    auto info = chem::atomInfo(doc);

    p.setPen(Qt::NoPen);
    for (const auto& f : doc.fills) {  // under everything else
        QPolygonF poly;
        for (int i : f.atoms) poly << doc.atoms[i].pos;
        p.setBrush(f.color);
        p.drawPolygon(poly);
    }
    p.setBrush(Qt::NoBrush);
    // A colour the user gave an item wins; otherwise the ink (the theme's on screen, black in exports).
    auto ink = [&](const QColor& own) { return own.isValid() ? own : style.ink; };
    // Arrows, shapes and orbitals sent behind the molecule, then (at the end) those over it.
    auto drawArrows = [&](bool behind) {
        for (const auto& a : doc.arrows) {
            if (a.behind != behind) continue;
            if (isOrbital(a.kind)) {
                drawOrbital(p, a, ink(a.color), lineWidth);
                continue;
            }
            p.setPen(QPen(ink(a.color), lineWidth, Qt::SolidLine, Qt::FlatCap, Qt::MiterJoin));
            drawArrow(p, a);
            if (a.crossed && !isShape(a.kind)) {  // ✕ over the middle of the shaft, solid even on a dashed arrow
                const auto pts = arrowPoints(a);
                const QPointF c = pts.size() == 2 ? (pts[0] + pts[1]) / 2 : pts[pts.size() / 2];
                const QPointF d = unit(pts.size() == 2 ? pts[1] - pts[0] : pts[pts.size() / 2 + 1] - pts[pts.size() / 2 - 1]);
                const QPointF u = rotated(d, 45) * kHeadLength * 0.75, v = rotated(d, -45) * kHeadLength * 0.75;
                p.setPen(QPen(ink(a.color), lineWidth, Qt::SolidLine, Qt::RoundCap));
                p.drawLine(c - u, c + u);
                p.drawLine(c - v, c + v);
            }
        }
    };
    drawArrows(true);
    // Aromatic circles: ring bonds drawn single, with a circle inside each ring.
    std::vector<std::vector<int>> circles;
    QSet<int> inCircle;  // bonds drawn single because of a circle
    if (doc.aromaticCircles || !doc.aromaticCircleOverrides.empty()) {
        for (const auto& ring : chem::aromaticRings(doc)) {
            auto ids = ring;
            std::sort(ids.begin(), ids.end());
            const bool override = std::find(doc.aromaticCircleOverrides.begin(),
                                            doc.aromaticCircleOverrides.end(), ids) != doc.aromaticCircleOverrides.end();
            if (doc.aromaticCircles != override) circles.push_back(ring);
        }
        for (const auto& r : circles)
            for (size_t k = 0; k < r.size(); ++k) inCircle.insert(doc.bondBetween(r[k], r[(k + 1) % r.size()]));
    }
    // Crossings: where a later bond (in front; Bring to Front moves one last) crosses an
    // earlier one that shares no atom with it, the earlier one gets a gap.
    // Sweep and prune: bonds sorted by left edge, so only pairs overlapping in x are tested.
    const int nb = int(doc.bonds.size());
    std::vector<QRectF> box(nb);
    for (int i = 0; i < nb; ++i)
        box[i] = QRectF(doc.atoms[doc.bonds[i].a].pos, doc.atoms[doc.bonds[i].b].pos).normalized();
    std::vector<int> byLeft(nb);
    std::iota(byLeft.begin(), byLeft.end(), 0);
    std::sort(byLeft.begin(), byLeft.end(), [&](int i, int j) { return box[i].left() < box[j].left(); });
    std::vector<std::vector<QPointF>> gaps(nb);
    for (int k = 0; k < nb; ++k)
        for (int m = k + 1; m < nb && box[byLeft[m]].left() <= box[byLeft[k]].right(); ++m) {
            const int i = std::min(byLeft[k], byLeft[m]), j = std::max(byLeft[k], byLeft[m]);
            if (box[i].top() > box[j].bottom() || box[j].top() > box[i].bottom()) continue;
            const Bond &x = doc.bonds[i], &y = doc.bonds[j];
            if (x.a == y.a || x.a == y.b || x.b == y.a || x.b == y.b) continue;
            QPointF at;
            if (QLineF(doc.atoms[x.a].pos, doc.atoms[x.b].pos).intersects(QLineF(doc.atoms[y.a].pos, doc.atoms[y.b].pos), &at) ==
                QLineF::BoundedIntersection)
                gaps[i].push_back(at);  // the earlier bond is behind
        }
    for (int bi = 0; bi < int(doc.bonds.size()); ++bi) {
        const Bond& b = doc.bonds[bi];
        p.setPen(QPen(ink(b.color), lineWidth, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
        if (inCircle.contains(bi) && b.order == 2 && b.stereo == BondStereo::None) {
            Bond single = b;
            single.order = 1;
            drawBond(p, doc, single, st, degree, labeled, bondsAt, gaps[bi]);
        } else {
            drawBond(p, doc, b, st, degree, labeled, bondsAt, gaps[bi]);
        }
    }
    p.setPen(QPen(style.ink, lineWidth));
    for (const auto& r : circles) {
        QPointF c;
        for (int i : r) c += doc.atoms[i].pos / double(r.size());
        const double apothem = kBondLength / (2 * std::tan(std::numbers::pi / r.size()));
        p.drawEllipse(c, 0.62 * apothem, 0.62 * apothem);
    }

    // H goes on the side away from the bonds; if bonds leave from both
    // left and right, above or below, whichever is clear.
    auto hSide = [&](int i) {
        double dx = 0;
        bool left = false, right = false, up = false;
        for (int nb : neighbors(doc, bondsAt, i)) {
            const QPointF d = doc.atoms[nb].pos - doc.atoms[i].pos;
            dx += d.x();
            left |= d.x() < -0.3 * kBondLength, right |= d.x() > 0.3 * kBondLength, up |= d.y() < -0.3 * kBondLength;
        }
        return left && right ? (up ? HSide::Below : HSide::Above) : dx > 0.1 ? HSide::Left : HSide::Right;
    };
    for (size_t i = 0; i < doc.atoms.size(); ++i) {
        const auto& a = doc.atoms[i];
        p.setPen(QPen(info[i].valenceError ? style.error : ink(a.color), lineWidth));
        if (labeled[i]) {
            drawLabel(p, doc, int(i), info[i].hydrogens, hSide(int(i)), st);
            if (style.labels) {
                QPicture label;
                QPainter q(&label);
                drawLabel(q, doc, int(i), info[i].hydrogens, hSide(int(i)), st);
                q.end();
                style.labels->resize(doc.atoms.size());
                (*style.labels)[i] = label.boundingRect();
            }
        } else if (a.charge) {
            QFont sub = labelFont(st, 0.7);
            QString c = QString(a.charge > 0 ? "+" : "−");
            if (std::abs(a.charge) > 1) c.prepend(QString::number(std::abs(a.charge)));
            const auto d = chargeDirection(doc, int(i), 0, HSide::Right);
            drawText(p, c, d ? chargeAt(a.pos, *d, 2.5, c, sub) : a.pos + QPointF(2, -3), sub);
        }
        if (info[i].valenceError && !labeled[i]) p.drawEllipse(a.pos, 3, 3);
    }
    drawArrows(false);
    for (const auto& t : doc.texts) p.fillPath(textPath(t, st), ink(t.color));
    // Atom-map numbers (:n) and, when shown, indices from 1: small, in the atom's widest gap,
    // ties (a ring atom's equal gaps) going to the side facing out of the drawing.
    QPointF middle;
    for (const auto& a : doc.atoms) middle += a.pos / double(doc.atoms.size());
    auto numberDirection = [&](int i) {
        const QPointF p = doc.atoms[i].pos;
        std::vector<double> ang;
        for (int nb : neighbors(doc, bondsAt, i)) ang.push_back(std::atan2(doc.atoms[nb].pos.y() - p.y(), doc.atoms[nb].pos.x() - p.x()));
        if (ang.size() < 2) return doc.awayDirection(i);
        std::sort(ang.begin(), ang.end());
        QPointF best;
        double bestGap = -1, bestOut = 0;
        for (size_t k = 0; k < ang.size(); ++k) {
            const double next = k + 1 < ang.size() ? ang[k + 1] : ang[0] + 2 * std::numbers::pi, mid = (ang[k] + next) / 2;
            const QPointF d(std::cos(mid), std::sin(mid));
            const double out = QPointF::dotProduct(d, p - middle);
            if (next - ang[k] > bestGap + 0.1 || (next - ang[k] > bestGap - 0.1 && out > bestOut))
                best = d, bestGap = std::max(bestGap, next - ang[k]), bestOut = out;
        }
        return best;
    };
    // Electron marks and δ±, each in the widest gap left by the bonds, the H and the charge.
    for (size_t i = 0; i < doc.atoms.size(); ++i) {
        const Atom& a = doc.atoms[i];
        const int marks = a.lonePairs + a.radicals + (a.partial ? 1 : 0);
        if (!marks) continue;
        std::vector<double> taken;
        for (int nb : neighbors(doc, bondsAt, int(i)))
            taken.push_back(std::atan2(doc.atoms[nb].pos.y() - a.pos.y(), doc.atoms[nb].pos.x() - a.pos.x()));
        if (labeled[i] && info[i].hydrogens) {
            const HSide s = hSide(int(i));
            taken.push_back(s == HSide::Right ? 0 : s == HSide::Left ? std::numbers::pi : s == HSide::Below ? std::numbers::pi / 2 : -std::numbers::pi / 2);
        }
        if (a.charge) {  // its corner
            const QPointF c = chargeDirection(doc, int(i), labeled[i] ? info[i].hydrogens : 0, hSide(int(i))).value_or(QPointF(1, -1));
            taken.push_back(std::atan2(c.y(), c.x()));
        }
        if (doc.showAtomNumbers || a.map) {  // the atom's number keeps its place (#349)
            const QPointF d = numberDirection(int(i));
            taken.push_back(std::atan2(d.y(), d.x()));
        }
        std::vector<QPointF> dirs;
        if (taken.empty()) dirs.push_back({0, -1}), taken.push_back(-std::numbers::pi / 2);  // a lone atom: on top first
        std::vector<std::pair<double, double>> gaps;  // (start, width)
        std::sort(taken.begin(), taken.end());
        for (size_t k = 0; k < taken.size(); ++k) {
            const double next = k + 1 < taken.size() ? taken[k + 1] : taken[0] + 2 * std::numbers::pi;
            gaps.push_back({taken[k], next - taken[k]});
        }
        while (int(dirs.size()) < marks) {
            auto widest = std::max_element(gaps.begin(), gaps.end(), [](auto x, auto y) { return x.second < y.second; });
            const auto [start, width] = *widest;
            const double mid = start + width / 2;
            dirs.push_back({std::cos(mid), std::sin(mid)});
            *widest = {start, mid - start};
            gaps.push_back({mid, start + width - mid});
        }
        const double r = (labeled[i] ? 0.66 : 0.3) * kBondLength, dot = 0.9 * lineWidth + 0.4;
        p.setPen(Qt::NoPen);
        p.setBrush(ink(a.color));
        size_t k = 0;
        for (int n = 0; n < a.radicals; ++n, ++k) p.drawEllipse(a.pos + dirs[k] * r, dot, dot);
        for (int n = 0; n < a.lonePairs; ++n, ++k)
            for (double side : {-1.0, 1.0}) p.drawEllipse(a.pos + dirs[k] * r + perp(dirs[k]) * side * 1.6, dot, dot);
        p.setBrush(Qt::NoBrush);
        if (a.partial) {
            const QFont f = labelFont(st, 0.7);
            const QString s = a.partial > 0 ? "δ+" : "δ−";
            QFontMetricsF fm(f);
            const QPointF at = a.pos + dirs[k] * (r + 2.5);
            p.setPen(QPen(ink(a.color), lineWidth));
            drawText(p, s, at - QPointF(fm.horizontalAdvance(s) / 2, -fm.capHeight() / 2), f);
        }
    }
    // Brackets: a pair of [ ] or ( ) just outside their atoms, the label at the bottom right.
    for (const Bracket& b : doc.brackets) {
        QPolygonF pts;
        for (int i : b.atoms) pts << doc.atoms[i].pos;
        const QRectF box = pts.boundingRect().adjusted(-0.55 * kBondLength, -0.6 * kBondLength, 0.55 * kBondLength, 0.6 * kBondLength);
        const double lip = 0.2 * kBondLength;
        p.setPen(QPen(style.ink, lineWidth));
        p.setBrush(Qt::NoBrush);
        for (double side : {-1.0, 1.0}) {
            const double x = side < 0 ? box.left() : box.right();
            QPainterPath path;
            if (b.square) {
                path.moveTo(x - side * lip, box.top());
                path.lineTo(x, box.top());
                path.lineTo(x, box.bottom());
                path.lineTo(x - side * lip, box.bottom());
            } else {
                path.moveTo(x - side * lip, box.top());
                path.quadTo(QPointF(x + side * lip * 0.6, box.center().y()), QPointF(x - side * lip, box.bottom()));
            }
            p.drawPath(path);
        }
        if (!b.label.isEmpty()) drawText(p, b.label, QPointF(box.right() + 1.5, box.bottom() + 1), labelFont(st, 0.7));
    }
    for (size_t i = 0; i < doc.atoms.size(); ++i) {
        const Atom& a = doc.atoms[i];
        QStringList parts;
        if (doc.showAtomNumbers) parts << QString::number(i + 1);
        if (a.map) parts << ":" + QString::number(a.map);
        if (parts.isEmpty()) continue;
        const QFont f = labelFont(st, 0.6);
        const QString s = parts.join(' ');
        QFontMetricsF fm(f);
        const QPointF at = a.pos + numberDirection(int(i)) * ((labeled[i] ? 0.8 : 0.4) * kBondLength);
        p.setPen(QPen(ink(a.color), lineWidth));
        drawText(p, s, at - QPointF(fm.horizontalAdvance(s) / 2, -fm.capHeight() / 2), f);
    }
    // Beside a stereocentre, clear of its bonds, and past the atom's number (#325); `further` stacks a second mark.
    auto besideCentre = [&](int i, double further) {
        if (doc.showAtomNumbers || doc.atoms[i].map) {
            const double number = (labeled[i] ? 0.8 : 0.4) * kBondLength;
            return doc.atoms[i].pos + numberDirection(i) * (number + (0.6 + further) * kBondLength);
        }
        return doc.atoms[i].pos + doc.awayDirection(i) * ((0.55 + further) * kBondLength);
    };
    for (size_t i = 0; i < doc.atoms.size(); ++i) {  // stereo groups, small and upright as ChemDraw sets them: abs, &1, or1
        const QString s = stereoGroupTag(doc.atoms[i]);
        if (s.isEmpty()) continue;
        const QFont f = labelFont(st, 0.6);
        QFontMetricsF fm(f);
        p.setPen(QPen(ink(doc.atoms[i].color), lineWidth));
        drawText(p, s, besideCentre(int(i), 0) - QPointF(fm.horizontalAdvance(s) / 2, -fm.capHeight() / 2), f);
    }
    if (doc.showStereo) {  // small italic (R)/(E), clear of the atom's bonds or the double bond's second line
        QFont f = labelFont(st, 0.7);
        f.setItalic(true);
        for (const auto& l : chem::stereoLabels(doc)) {
            QPointF at;
            if (l.atom >= 0) {
                at = besideCentre(l.atom, doc.atoms[l.atom].stereoGroup != StereoGroup::None ? 0.7 : 0);
            } else {
                const Bond& b = doc.bonds[l.bond];
                const QPointF a = doc.atoms[b.a].pos, e = doc.atoms[b.b].pos, n = perp(unit(e - a));
                at = (a + e) / 2 - n * (doubleBondSide(doc, b, bondsAt) > 0 ? 1 : -1) * (0.45 * kBondLength);
            }
            const QString s = "(" + l.text + ")";
            QFontMetricsF fm(f);
            p.setPen(QPen(style.ink, lineWidth));
            drawText(p, s, at - QPointF(fm.horizontalAdvance(s) / 2, -fm.capHeight() / 2), f);
        }
    }
    if (doc.showShifts) {  // predicted 13C (1H) beside each atom (#403), with the notice nmrshiftdb2's licence asks for
        const auto shifts = chem::predictShifts(doc);
        QSet<int> labelled;  // stereocentres with an (R)/(S) drawn: past it
        if (doc.showStereo)
            for (const auto& l : chem::stereoLabels(doc)) labelled.insert(l.atom);
        const QFont f = labelFont(st, 0.55);
        QFontMetricsF fm(f);
        p.setPen(QPen(style.ink, lineWidth));
        auto value = [](double v, int spheres, int decimals) { return (spheres < 3 ? "~" : "") + QString::number(v, 'f', decimals); };
        // The widest gap between the atom's bonds; at a branch or ring fusion, the gap facing out of the molecule.
        QPointF middle;
        for (const Atom& a : doc.atoms) middle += a.pos / double(doc.atoms.size());
        auto shiftDirection = [&](int i) {
            const auto nbs = neighbors(doc, bondsAt, i);
            if (nbs.size() < 3) return doc.awayDirection(i);
            std::vector<double> ang;
            for (int nb : nbs) ang.push_back(std::atan2(doc.atoms[nb].pos.y() - doc.atoms[i].pos.y(), doc.atoms[nb].pos.x() - doc.atoms[i].pos.x()));
            std::sort(ang.begin(), ang.end());
            const QPointF out = doc.atoms[i].pos - middle;
            QPointF best = doc.awayDirection(i);
            double bestScore = -1e9;
            for (size_t k = 0; k < ang.size(); ++k) {
                const double next = k + 1 < ang.size() ? ang[k + 1] : ang[0] + 2 * std::numbers::pi, mid = (ang[k] + next) / 2;
                const QPointF d(std::cos(mid), std::sin(mid));
                const double score = (next - ang[k]) * kBondLength + QPointF::dotProduct(d, out);  // wide, and outward
                if (score > bestScore) bestScore = score, best = d;
            }
            return best;
        };
        for (const auto& s : shifts) {
            QStringList parts;
            if (s.carbonSpheres) parts << value(s.carbon, s.carbonSpheres, 1);
            if (s.protonSpheres) parts << "(" + value(s.proton, s.protonSpheres, 2) + ")";
            const QString text = parts.join(' ');
            const double further = (doc.atoms[s.atom].stereoGroup != StereoGroup::None ? 0.7 : 0) + (labelled.contains(s.atom) ? 0.9 : 0) +
                                   (doc.showAtomNumbers || doc.atoms[s.atom].map ? 0.6 : 0);
            // The text's near edge clear of the atom: pushed out by half its size along the way out.
            const QPointF d = shiftDirection(s.atom), w(fm.horizontalAdvance(text) / 2, fm.capHeight() / 2);
            const QPointF centre = doc.atoms[s.atom].pos + d * ((labeled[s.atom] ? 0.6 : 0.3) + further) * kBondLength + QPointF(d.x() * w.x(), d.y() * w.y());
            drawText(p, text, centre - QPointF(w.x(), -w.y()), f);
        }
        if (!shifts.empty()) {
            QPointF at(documentBounds(doc).left(), documentBounds(doc).bottom() + 1.2 * fm.height());
            for (const QString& line : QStringList{QObject::tr("Predicted shifts in ppm: 13C (1H); ~ marks a weaker match.")} + nmrshiftdbNotice()) {
                drawText(p, line, at, f);
                at.ry() += 1.2 * fm.height();
            }
        }
    }
    p.restore();
}

QRectF documentBounds(const Document& doc) {
    const DrawingStyle st = documentStyle(doc);
    const double fs = st.fontSize;
    if (doc.empty()) return {};
    // Not QRectF::united: it ignores zero-size rects.
    double inf = std::numeric_limits<double>::infinity();
    QPointF lo(inf, inf), hi(-inf, -inf);
    auto grow = [&](QRectF r) {
        lo = {std::min(lo.x(), r.left()), std::min(lo.y(), r.top())};
        hi = {std::max(hi.x(), r.right()), std::max(hi.y(), r.bottom())};
    };
    std::vector<int> degree(doc.atoms.size(), 0);
    for (const auto& b : doc.bonds) ++degree[b.a], ++degree[b.b];
    for (size_t i = 0; i < doc.atoms.size(); ++i) {  // room for labels, which can run either way
        const Atom& a = doc.atoms[i];
        const qsizetype chars = a.label.isEmpty() && a.isotope ? 3 + massNumber(a).size() : a.label.size();  // H₃¹³C
        double w = fs * std::max(1.5, 0.7 * chars);
        // A labelled atom with bonds on both sides may stack its H above or below.
        double h = fs * (degree[i] >= 2 && hasLabel(doc, int(i), degree) ? 2.2 : 1);
        grow(QRectF(a.pos, a.pos).adjusted(-w, -h, w, h));
    }
    for (const auto& a : doc.arrows) grow(arrowPath(a).boundingRect().adjusted(-4, -4, 4, 4));
    for (const auto& b : doc.brackets) {  // the brackets and their label
        QPolygonF pts;
        for (int i : b.atoms) pts << doc.atoms[i].pos;
        grow(pts.boundingRect().adjusted(-0.7 * kBondLength, -0.7 * kBondLength, 1.2 * kBondLength, 0.9 * kBondLength));
    }
    for (const auto& t : doc.texts) grow(textPath(t, st).boundingRect().adjusted(-2, -2, 2, 2));
    return QRectF(lo, hi);
}

QRectF outputBounds(const Document& doc) {
    // Atom numbers, stereo labels and electron marks go where the painter finds
    // room: measure what it paints (#325).
    QPicture painted;
    {
        QPainter p(&painted);
        paintDocument(p, doc);
    }
    const QRect r = painted.boundingRect();
    return r.isEmpty() ? documentBounds(doc) : documentBounds(doc).united(QRectF(r).adjusted(-1, -1, 1, 1));
}

const std::vector<PageSize>& pageSizes() {
    static const std::vector<PageSize> p{
        {"A4", {595.3, 841.9}, 72},
        {"US Letter", {612, 792}, 72},
        {"ACS single column", {240, 684}},  // 3.33 in × 9.5 in
        {"ACS double column", {504, 684}},  // 7 in
        {"RSC single column", {235.3, 660.5}},  // 8.3 cm × 23.3 cm
        {"RSC double column", {484.7, 660.5}},  // 17.1 cm
    };
    return p;
}

QRectF pageRect(const Document& doc) {
    for (const auto& p : pageSizes())
        if (p.name == doc.page) return QRectF(doc.pageOrigin, p.size / exportScale(doc));
    return {};
}

double exportScale(const Document& doc) { return drawingStyle(doc.style).bondLength / kBondLength; }

// The exported area (model units) and the points per model unit.
static std::pair<QRectF, double> exportFrame(const Document& doc, const ExportOptions& o) {
    const double s = exportScale(doc) * o.scale, m = o.margin / s;
    if (const QRectF page = pageRect(doc); !page.isEmpty()) return {page, s};  // the whole page, as laid out
    return {outputBounds(doc).adjusted(-m, -m, m, m), s};
}

// Paints the frame at `s` output units per point, background first.
static void paintFrame(QPainter& p, const Document& doc, const ExportOptions& o, double perPoint) {
    const auto [r, s] = exportFrame(doc, o);
    p.scale(s * perPoint, s * perPoint);
    p.translate(-r.topLeft());
    if (o.background.alpha()) p.fillRect(r, o.background);
    paintDocument(p, doc);
}

Legend placeLegend(const QRectF& plot, QSizeF size, const std::vector<QPointF>& sticks, double clear) {
    Legend best;
    for (bool left : {true, false}) {  // the corner that shrinks the sticks least; left on a tie
        const QRectF r(left ? plot.left() : plot.right() - size.width(), plot.top(), size.width(), size.height());
        double scale = 1;
        for (QPointF s : sticks)
            if (s.x() > r.left() - 2 && s.x() < r.right() + 2 && s.y() > 0)
                scale = std::min(scale, std::max(plot.bottom() - r.bottom() - clear, 0.0) / (s.y() * plot.height()));
        if (left || scale > best.scale) best = {r, scale};
    }
    return best;
}

QStringList nmrshiftdbNotice() {
    return {QObject::tr("Contains information from nmrshiftdb2 (www.nmrshiftdb.org), which is made available here"),
            QObject::tr("under the nmrshiftdb2 Database License (%1).")
                .arg("https://nmrshiftdb.nmr.uni-koeln.de/nmrshiftdbhtml/nmrshiftdb2datalicense.txt")};
}

QImage renderImage(const Document& doc, const ExportOptions& o) {
    const auto [r, s] = exportFrame(doc, o);
    QImage img((r.size() * s * o.dpi / 72.0).toSize().expandedTo({1, 1}), QImage::Format_ARGB32_Premultiplied);
    img.setDotsPerMeterX(int(o.dpi / 0.0254));
    img.setDotsPerMeterY(int(o.dpi / 0.0254));
    img.fill(Qt::transparent);
    QPainter p(&img);
    paintFrame(p, doc, o, o.dpi / 72.0);
    p.end();
    img.setText("penzene", QString::fromUtf8(doc.toJson()));  // reopens as an editable drawing
    return img;
}

QByteArray renderPng(const Document& doc, const ExportOptions& o) {
    QByteArray png;
    QBuffer buf(&png);
    if (!buf.open(QIODevice::WriteOnly) || !renderImage(doc, o).save(&buf, "PNG")) return {};
    return png;
}

QByteArray renderSvg(const Document& doc, const ExportOptions& o) {
    const auto [r, s] = exportFrame(doc, o);
    QBuffer buf;
    QSvgGenerator gen;
    gen.setOutputDevice(&buf);
    gen.setSize((r.size() * s).toSize());
    gen.setViewBox(QRectF(QPointF(), r.size() * s));
    gen.setResolution(72);  // 1 unit == 1 pt
    gen.setTitle("Penzene");
    QPainter p(&gen);
    paintFrame(p, doc, o, 1);
    p.end();
    // The editable drawing rides along (base64: nothing in it can break the XML).
    QByteArray svg = buf.data();
    const qsizetype open = svg.indexOf('>', svg.indexOf("<svg"));
    if (open > 0) svg.insert(open + 1, "\n<metadata id=\"penzene\">" + doc.toJson().toBase64() + "</metadata>");
    return svg;
}

QByteArray renderPdf(const Document& doc, const ExportOptions& o) {
    const auto [r, s] = exportFrame(doc, o);
    QBuffer buf;
    buf.open(QIODevice::WriteOnly);
    {
        QPdfWriter pdf(&buf);
        pdf.setResolution(72);
        pdf.setPageSize(QPageSize(r.size() * s, QPageSize::Point));
        pdf.setPageMargins({});
        pdf.setCreator("Penzene");
        // The editable drawing rides along as attachments: ours in full, and a MOL for anything else.
        pdf.addFileAttachment("drawing.penz", doc.toJson(), "application/x-penzene");
        if (!doc.atoms.empty())
            if (const std::string mol = chem::toMolBlock(doc); !mol.empty())
                pdf.addFileAttachment("structure.mol", QByteArray::fromStdString(mol), "chemical/x-mdl-molfile");
        QPainter p(&pdf);
        paintFrame(p, doc, o, 1);
    }
    return buf.data();
}

bool exportDocument(const Document& doc, const QString& path, const ExportOptions& o) {
    if (doc.empty()) return false;
    const QString ext = QFileInfo(path).suffix().toLower();
    QByteArray data;
    if (ext == "png") {
        data = renderPng(doc, o);
    } else if (ext == "svg") {
        data = renderSvg(doc, o);
    } else if (ext == "pdf") {
        data = renderPdf(doc, o);
    }
    return !data.isEmpty() && writeWhole(path, data);
}

bool writeWhole(const QString& path, const QByteArray& data) {
    QSaveFile f(path);  // written aside and swapped in, so a failed write never truncates the old file
    f.setDirectWriteFallback(true);  // a writable file in a read-only folder (sandbox portals)
    return f.open(QIODevice::WriteOnly) && f.write(data) == data.size() && f.commit();
}

const std::vector<Theme>& themes() {
    auto cat = [](const char* name, bool dark, const char* base, const char* mantle, const char* surface,
                  const char* text, const char* red, const char* blue, const char* green) {
        return Theme{name, dark, QColor(base), QColor(text), QColor(red), QColor(blue), QColor(green),
                     QColor(mantle), QColor(surface), QColor(text)};
    };
    static const std::vector<Theme> t{
        {"System"},
        // White and teal: a white page in a cool grey window.
        {"Light", false, QColor("#FFFFFF"), QColor("#1D2125"), QColor(220, 40, 40),
         QColor("#0F6E56"), QColor("#0F6E56"), QColor("#F3F5F6"), QColor("#FFFFFF"), QColor("#1D2125")},
        {"Dark", true, QColor("#15171A"), QColor("#E8EAED"), QColor(255, 105, 97),
         QColor("#4CC9A0"), QColor("#4CC9A0"), QColor("#1C1F23"), QColor("#24282D"), QColor("#E8EAED")},
        // https://catppuccin.com/palette: base, mantle, surface0, text, red, blue, green (hotspot)
        // Latte's teal, not its green, for the hotspot: the green is under 3:1 on its base.
        cat("Catppuccin Latte", false, "#eff1f5", "#e6e9ef", "#ccd0da", "#4c4f69", "#d20f39", "#1e66f5", "#179299"),
        cat("Catppuccin Frappé", true, "#303446", "#292c3c", "#414559", "#c6d0f5", "#e78284", "#8caaee", "#a6d189"),
        cat("Catppuccin Macchiato", true, "#24273a", "#1e2030", "#363a4f", "#cad3f5", "#ed8796", "#8aadf4", "#a6da95"),
        cat("Catppuccin Mocha", true, "#1e1e2e", "#181825", "#313244", "#cdd6f4", "#f38ba8", "#89b4fa", "#a6e3a1"),
    };
    return t;
}

const Theme& theme(const QString& name) {
    for (const auto& t : themes())
        if (t.name == name) return t;
    return themes()[0];
}

Chrome chrome(const Theme& t) {
    if (t.name == "Light") return {QColor("#DDE1E4"), QColor("#5E6770"), QColor("#E1F5EE")};
    if (t.name == "Dark") return {QColor("#353A40"), QColor("#A3AAB2"), QColor("#0B3B30")};
    // The Catppuccin themes keep their own colours.
    return {t.surface.lighter(125), t.text, t.dark ? t.surface.lighter(145) : t.surface.darker(110)};
}

#ifdef Q_OS_WIN
// Enhanced Metafile, for Office on Windows (#394). Qt has no EMF paint device, so this engine
// replays the export as GDI paths into a metafile DC: vector throughout, text as glyph outlines.
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>

namespace {
constexpr double kTwipsPerPoint = 20;  // MM_TWIPS: 1/1440 inch, y up

// GDI has no alpha: a translucent colour is mixed with white, as it shows on a page.
COLORREF gdiColor(const QColor& c) {
    const double a = c.alphaF();
    auto mixed = [a](int v) { return BYTE(std::lround(v * a + 255 * (1 - a))); };
    return RGB(mixed(c.red()), mixed(c.green()), mixed(c.blue()));
}

class EmfEngine : public QPaintEngine {
public:
    explicit EmfEngine(HDC dc) : QPaintEngine(AllFeatures), dc_(dc) {}
    bool begin(QPaintDevice*) override { return true; }
    bool end() override { return true; }
    Type type() const override { return kEmfEngine; }
    void updateState(const QPaintEngineState& s) override {
        if (s.state() & DirtyPen) pen_ = s.pen();
        if (s.state() & DirtyBrush) brush_ = s.brush();
        if (s.state() & DirtyTransform) transform_ = s.transform();
    }
    void drawPixmap(const QRectF&, const QPixmap&, const QRectF&) override {}  // the renderer draws none
    using QPaintEngine::drawPolygon;
    void drawPolygon(const QPointF* points, int count, PolygonDrawMode mode) override {
        QPainterPath path;
        path.setFillRule(mode == WindingMode ? Qt::WindingFill : Qt::OddEvenFill);
        path.addPolygon(QPolygonF(QList<QPointF>(points, points + count)));
        if (mode != PolylineMode) {
            path.closeSubpath();
            return drawPath(path);
        }
        const QBrush brush = std::exchange(brush_, QBrush());
        drawPath(path);
        brush_ = brush;
    }
    void drawPath(const QPainterPath& path) override {
        const bool fill = brush_.style() != Qt::NoBrush && brush_.color().alpha() > 0;
        const bool stroke = pen_.style() != Qt::NoPen && pen_.color().alpha() > 0;
        if (!fill && !stroke) return;
        const QPainterPath mapped = transform_.map(path);
        auto at = [&](int i) {
            const QPainterPath::Element e = mapped.elementAt(i);
            return POINT{LONG(std::lround(e.x)), LONG(-std::lround(e.y))};
        };
        BeginPath(dc_);
        POINT start{}, last{};
        auto closeIfClosed = [&] {
            if (start.x == last.x && start.y == last.y) CloseFigure(dc_);
        };
        for (int i = 0; i < mapped.elementCount(); ++i) {
            const QPainterPath::Element e = mapped.elementAt(i);
            if (e.isMoveTo()) {
                if (i > 0) closeIfClosed();
                start = last = at(i);
                MoveToEx(dc_, start.x, start.y, nullptr);
            } else if (e.isLineTo()) {
                last = at(i);
                LineTo(dc_, last.x, last.y);
            } else if (e.isCurveTo() && i + 2 < mapped.elementCount()) {
                const POINT curve[3] = {at(i), at(i + 1), at(i + 2)};
                PolyBezierTo(dc_, curve, 3);
                last = curve[2];
                i += 2;
            }
        }
        closeIfClosed();
        EndPath(dc_);
        SetPolyFillMode(dc_, path.fillRule() == Qt::WindingFill ? WINDING : ALTERNATE);
        HGDIOBJ brush = fill ? HGDIOBJ(CreateSolidBrush(gdiColor(brush_.color()))) : GetStockObject(NULL_BRUSH);
        HGDIOBJ pen = stroke ? HGDIOBJ(gdiPen()) : GetStockObject(NULL_PEN);
        HGDIOBJ oldBrush = SelectObject(dc_, brush), oldPen = SelectObject(dc_, pen);
        if (fill && stroke) StrokeAndFillPath(dc_);
        else if (fill) FillPath(dc_);
        else StrokePath(dc_);
        SelectObject(dc_, oldBrush);
        SelectObject(dc_, oldPen);
        if (fill) DeleteObject(brush);
        if (stroke) DeleteObject(pen);
    }

private:
    HPEN gdiPen() {
        const double scale = pen_.isCosmetic() ? 1 : std::sqrt(std::abs(transform_.determinant()));
        const double width = std::max(1.0, pen_.widthF() * scale);
        DWORD style = PS_GEOMETRIC | PS_SOLID;
        std::vector<DWORD> dashes;
        if (pen_.style() != Qt::SolidLine) {  // Qt's dashes are in pen widths
            style = PS_GEOMETRIC | PS_USERSTYLE;
            for (qreal d : pen_.dashPattern()) dashes.push_back(DWORD(std::max(1.0, std::round(d * width))));
        }
        style |= pen_.capStyle() == Qt::RoundCap ? PS_ENDCAP_ROUND : pen_.capStyle() == Qt::SquareCap ? PS_ENDCAP_SQUARE : PS_ENDCAP_FLAT;
        style |= pen_.joinStyle() == Qt::RoundJoin ? PS_JOIN_ROUND : pen_.joinStyle() == Qt::BevelJoin ? PS_JOIN_BEVEL : PS_JOIN_MITER;
        SetMiterLimit(dc_, FLOAT(pen_.miterLimit()), nullptr);
        const LOGBRUSH brush{BS_SOLID, gdiColor(pen_.color()), 0};
        return ExtCreatePen(style, DWORD(std::lround(width)), &brush, DWORD(dashes.size()), dashes.empty() ? nullptr : dashes.data());
    }
    HDC dc_;
    QPen pen_;
    QBrush brush_;
    QTransform transform_;
};

class EmfDevice : public QPaintDevice {
public:
    EmfDevice(HDC dc, QSize size) : engine_(dc), size_(size) {}
    QPaintEngine* paintEngine() const override { return &engine_; }

protected:
    int metric(PaintDeviceMetric m) const override {
        switch (m) {
        case PdmWidth: return size_.width();
        case PdmHeight: return size_.height();
        case PdmWidthMM: return qRound(size_.width() * 25.4 / 1440);
        case PdmHeightMM: return qRound(size_.height() * 25.4 / 1440);
        case PdmDpiX: case PdmDpiY: case PdmPhysicalDpiX: case PdmPhysicalDpiY: return 1440;
        case PdmDepth: return 32;
        case PdmNumColors: return std::numeric_limits<int>::max();
        default: return QPaintDevice::metric(m);
        }
    }

private:
    mutable EmfEngine engine_;
    QSize size_;
};
}  // namespace

QByteArray renderEmf(const Document& doc, const ExportOptions& o) {
    const auto [r, s] = exportFrame(doc, o);
    const QSizeF points = r.size() * s;
    const RECT frame{0, 0, LONG(std::lround(points.width() * 2540 / 72)), LONG(std::lround(points.height() * 2540 / 72))};  // 0.01 mm
    HDC dc = CreateEnhMetaFileW(nullptr, nullptr, &frame, L"Penzene\0Structure\0");
    if (!dc) return {};
    SetMapMode(dc, MM_TWIPS);
    SetBkMode(dc, TRANSPARENT);  // dash gaps stay clear
    {
        EmfDevice device(dc, (points * kTwipsPerPoint).toSize().expandedTo({1, 1}));
        QPainter p(&device);
        paintFrame(p, doc, o, kTwipsPerPoint);
    }
    HENHMETAFILE emf = CloseEnhMetaFile(dc);
    if (!emf) return {};
    QByteArray bytes(qsizetype(GetEnhMetaFileBits(emf, 0, nullptr)), Qt::Uninitialized);
    GetEnhMetaFileBits(emf, UINT(bytes.size()), reinterpret_cast<BYTE*>(bytes.data()));
    DeleteEnhMetaFile(emf);
    return bytes;
}
#endif
