#include "Document.h"

#include <QJsonArray>
#include <QImage>
#include <QRegularExpression>
#include <QJsonDocument>
#include <QJsonObject>
#include <QObject>
#include <QtEndian>
#include <algorithm>
#include <cmath>
#include <numbers>

static const char* kStereo[] = {"none", "wedge", "hash", "bold", "dashed", "wavy", "interaction", "partial"};
static const char* kPosition[] = {"auto", "left", "centre", "right"};
static const char* kArrow[] = {"reaction", "equilibrium", "resonance", "retro", "fishhook",
                               "line", "box", "roundedbox", "ellipse",
                               "s-orbital", "p-orbital", "lobe", "hybrid-orbital"};
static const char* kLook[] = {"outline", "shaded", "gradient"};

QString stereoGroupTag(const Atom& a) {
    switch (a.stereoGroup) {
    case StereoGroup::Abs: return "abs";
    case StereoGroup::And: return "&" + QString::number(a.stereoGroupNumber);
    case StereoGroup::Or: return "or" + QString::number(a.stereoGroupNumber);
    default: return {};
    }
}

static void setStereoGroupTag(Atom& a, const QString& tag) {  // anything else: none
    static const QRegularExpression re("^(abs|&|or)([1-9][0-9]{0,3})?$");
    const auto m = re.match(tag);
    const bool numbered = m.captured(1) != "abs";
    if (!m.hasMatch() || numbered == m.captured(2).isEmpty()) return;
    a.stereoGroup = !numbered ? StereoGroup::Abs : m.captured(1) == "&" ? StereoGroup::And : StereoGroup::Or;
    a.stereoGroupNumber = m.captured(2).toInt();
}

QByteArray Document::toJson() const {
    QJsonArray as, bs;
    for (const auto& a : atoms) {
        QJsonObject o{{"x", a.pos.x()}, {"y", a.pos.y()}, {"z", a.z}};
        if (a.charge) o["charge"] = a.charge;
        if (!a.label.isEmpty()) o["label"] = a.label;
        if (a.color.isValid()) o["color"] = a.color.name();
        if (a.map) o["map"] = a.map;
        if (a.lonePairs) o["lonePairs"] = a.lonePairs;
        if (a.radicals) o["radicals"] = a.radicals;
        if (a.partial) o["partial"] = a.partial;
        if (a.isotope) o["isotope"] = a.isotope;
        if (a.stereoGroup != StereoGroup::None) o["stereoGroup"] = stereoGroupTag(a);
        as.append(o);
    }
    for (const auto& b : bonds) {
        QJsonObject o{{"a", b.a}, {"b", b.b}, {"order", b.order}};
        if (b.stereo != BondStereo::None) o["stereo"] = kStereo[int(b.stereo)];
        if (b.position != BondPosition::Auto) o["position"] = kPosition[int(b.position)];
        if (b.color.isValid()) o["color"] = b.color.name();
        bs.append(o);
    }
    QJsonObject root{{"format", "penzene"}, {"version", 1}, {"atoms", as}, {"bonds", bs}};
    QJsonArray ar, ts;
    for (const auto& a : arrows) {
        QJsonObject o{{"x1", a.from.x()}, {"y1", a.from.y()}, {"x2", a.to.x()}, {"y2", a.to.y()},
                      {"kind", kArrow[int(a.kind)]}};
        if (a.bend) o["bend"] = a.bend;
        if (a.color.isValid()) o["color"] = a.color.name();
        if (a.dashed) o["dashed"] = true;
        if (a.look != OrbitalLook::Outline) o["look"] = kLook[int(a.look)];
        if (a.behind) o["behind"] = true;
        if (a.crossed) o["crossed"] = true;
        if (a.head != 1) o["head"] = a.head;
        for (auto [key, at] : {std::pair{"fromAt", a.fromAt}, std::pair{"toAt", a.toAt}})
            if (at[0] >= 0) o[key] = at[1] >= 0 ? QJsonArray{at[0], at[1]} : QJsonArray{at[0]};
        ar.append(o);
    }
    for (const auto& t : texts) {
        QJsonObject o{{"x", t.pos.x()}, {"y", t.pos.y()}, {"text", t.text}};
        if (t.scale != 1) o["scale"] = t.scale;
        if (t.color.isValid()) o["color"] = t.color.name();
        if (t.compound) o["compound"] = true;
        if (t.anchor >= 0) o["anchor"] = t.anchor;
        ts.append(o);
    }
    if (!ar.isEmpty()) root["arrows"] = ar;
    if (!ts.isEmpty()) root["texts"] = ts;
    QJsonArray fs;
    for (const auto& f : fills) {
        QJsonArray ids;
        for (int i : f.atoms) ids.append(i);
        fs.append(QJsonObject{{"atoms", ids}, {"color", f.color.name()}});
    }
    if (!fs.isEmpty()) root["fills"] = fs;
    QJsonArray bs2;
    for (const auto& b : brackets) {
        QJsonArray ids;
        for (int i : b.atoms) ids.append(i);
        QJsonObject o{{"atoms", ids}};
        if (!b.square) o["round"] = true;
        if (!b.label.isEmpty()) o["label"] = b.label;
        bs2.append(o);
    }
    if (!bs2.isEmpty()) root["brackets"] = bs2;
    if (!style.isEmpty()) root["style"] = style;
    if (carbonLabels != CarbonLabels::None) root["carbonLabels"] = carbonLabels == CarbonLabels::All ? "all" : "terminal";
    if (labelRatio > 0) root["labelRatio"] = labelRatio;
    if (showStereo) root["showStereo"] = true;
    if (showAtomNumbers) root["showAtomNumbers"] = true;
    if (showShifts) root["showShifts"] = true;
    if (aromaticCircles) root["aromaticCircles"] = true;
    if (!page.isEmpty()) root["page"] = QJsonObject{{"name", page}, {"x", pageOrigin.x()}, {"y", pageOrigin.y()}};
    QJsonArray circleOverrides;
    for (const auto& ring : aromaticCircleOverrides) {
        QJsonArray ids;
        for (int i : ring) ids.append(i);
        circleOverrides.append(ids);
    }
    if (!circleOverrides.isEmpty()) root["aromaticCircleOverrides"] = circleOverrides;
    return QJsonDocument(root).toJson(QJsonDocument::Indented);
}

std::optional<Document> Document::fromEmbedded(const QByteArray& file) {
    if (file.startsWith("\x89PNG")) return fromJson(QImage::fromData(file, "PNG").text("penzene").toUtf8());
    if (file.startsWith("%PDF")) {
        // The drawing is an attached file: a /Filespec names its stream as /EF <</F n 0 R>>. Only
        // those are inflated, so pasting a big PDF from elsewhere costs a search, not a decode.
        for (qsizetype ef = file.indexOf("/EF"); ef >= 0; ef = file.indexOf("/EF", ef + 3)) {
            static const QRegularExpression ref(R"(^/EF\s*<<\s*/F\s+(\d+)\s+0\s+R)");
            const auto m = ref.match(QString::fromLatin1(file.mid(ef, 40)));
            if (!m.hasMatch()) continue;
            const QByteArray head = m.captured(1).toLatin1() + " 0 obj";
            qsizetype obj = file.indexOf("\n" + head);
            if (obj < 0) continue;
            const qsizetype at = file.indexOf("stream", obj);
            if (at < 0) continue;
            qsizetype begin = at + 6;
            if (file.mid(begin, 2) == "\r\n") begin += 2;
            else if (file.mid(begin, 1) == "\n") begin += 1;
            const qsizetype end = file.indexOf("endstream", begin);
            if (end < 0) continue;
            QByteArray data = file.mid(begin, end - begin);
            if (file.mid(obj, at - obj).contains("/FlateDecode")) {
                QByteArray sized(4, 0);
                qToBigEndian<quint32>(quint32(data.size() * 16), sized.data());  // qUncompress's size hint
                data = qUncompress(sized + data);
            }
            if (auto doc = fromJson(data)) return doc;
        }
        return std::nullopt;
    }
    static const QRegularExpression svg(R"(<metadata id="penzene">([A-Za-z0-9+/=]*)</metadata>)");
    const auto m = svg.match(QString::fromUtf8(file));
    return m.hasMatch() ? fromJson(QByteArray::fromBase64(m.captured(1).toLatin1())) : std::nullopt;
}

std::optional<Document> Document::fromJson(const QByteArray& data) {
    auto root = QJsonDocument::fromJson(data).object();
    if (root["format"].toString() == "penzene" && root["version"].toInt() == 2) {  // its first page
        const auto sheets = sheetsFromJson(data);
        return sheets.empty() ? std::nullopt : std::optional(sheets[0].doc);
    }
    if (root["format"].toString() != "penzene" || root["version"].toInt() != 1)
        return std::nullopt;
    Document doc;
    doc.style = root["style"].toString();
    const QString carbons = root["carbonLabels"].toString();
    doc.carbonLabels = carbons == "all"        ? Document::CarbonLabels::All
                       : carbons == "terminal" ? Document::CarbonLabels::Terminal
                                               : Document::CarbonLabels::None;
    doc.labelRatio = std::max(0.0, root["labelRatio"].toDouble());
    doc.showStereo = root["showStereo"].toBool();
    doc.showAtomNumbers = root["showAtomNumbers"].toBool();
    doc.showShifts = root["showShifts"].toBool();
    doc.aromaticCircles = root["aromaticCircles"].toBool();
    const auto page = root["page"].toObject();
    doc.page = page["name"].toString();
    doc.pageOrigin = {page["x"].toDouble(), page["y"].toDouble()};
    if (!std::isfinite(doc.pageOrigin.x()) || !std::isfinite(doc.pageOrigin.y())) return std::nullopt;
    for (const auto& v : root["atoms"].toArray()) {
        auto o = v.toObject();
        doc.atoms.push_back({QPointF(o["x"].toDouble(), o["y"].toDouble()),
                             o["z"].toInt(6), o["charge"].toInt(), o["label"].toString(),
                             QColor(o["color"].toString()), std::max(0, o["map"].toInt()),
                             std::clamp(o["lonePairs"].toInt(), 0, 4), std::clamp(o["radicals"].toInt(), 0, 2),
                             std::clamp(o["partial"].toInt(), -1, 1), std::clamp(o["isotope"].toInt(), 0, 300)});
        // An element RDKit doesn't know aborts the app wherever the atom is looked up.
        Atom& a = doc.atoms.back();
        if (a.z < 0 || a.z > 118 || !std::isfinite(a.pos.x()) || !std::isfinite(a.pos.y())) return std::nullopt;
        if (a.isotope < a.z) a.isotope = 0;  // lighter than its protons: no such isotope, as a typed label (#368)
        setStereoGroupTag(a, o["stereoGroup"].toString());
    }
    const int n = int(doc.atoms.size());
    for (const auto& v : root["bonds"].toArray()) {
        auto o = v.toObject();
        Bond b{o["a"].toInt(-1), o["b"].toInt(-1), o["order"].toInt(1)};
        // Untrusted file: reject dangling, self or duplicate bonds rather than crash later.
        if (b.a < 0 || b.a >= n || b.b < 0 || b.b >= n || b.a == b.b || doc.bondBetween(b.a, b.b) >= 0) return std::nullopt;
        b.order = std::clamp(b.order, 1, 3);
        auto index = [](const auto& names, const QString& s) {
            auto it = std::find(std::begin(names), std::end(names), s);
            return it == std::end(names) ? 0 : int(it - std::begin(names));  // unknown: default
        };
        b.stereo = BondStereo(index(kStereo, o["stereo"].toString()));
        b.position = BondPosition(index(kPosition, o["position"].toString()));
        b.color = QColor(o["color"].toString());
        doc.bonds.push_back(b);
    }
    auto finite = [](std::initializer_list<double> v) {
        return std::all_of(v.begin(), v.end(), [](double x) { return std::isfinite(x); });
    };
    for (const auto& v : root["arrows"].toArray()) {
        auto o = v.toObject();
        Arrow a{{o["x1"].toDouble(), o["y1"].toDouble()}, {o["x2"].toDouble(), o["y2"].toDouble()}};
        auto k = std::find(std::begin(kArrow), std::end(kArrow), o["kind"].toString("reaction"));
        a.bend = o["bend"].toDouble();
        if (k == std::end(kArrow) || !finite({a.from.x(), a.from.y(), a.to.x(), a.to.y(), a.bend}))
            return std::nullopt;
        a.kind = ArrowKind(k - std::begin(kArrow));
        a.color = QColor(o["color"].toString());
        a.dashed = o["dashed"].toBool();
        auto look = std::find(std::begin(kLook), std::end(kLook), o["look"].toString("outline"));
        if (look == std::end(kLook)) return std::nullopt;
        a.look = OrbitalLook(look - std::begin(kLook));
        a.behind = o["behind"].toBool();
        a.crossed = o["crossed"].toBool();
        if (const double head = o["head"].toDouble(1); std::isfinite(head)) a.head = std::clamp(head, 0.25, 4.0);
        for (auto [key, at] : {std::pair{"fromAt", &a.fromAt}, std::pair{"toAt", &a.toAt}}) {
            const QJsonArray v = o[key].toArray();
            for (int k = 0; k < std::min<int>(2, v.size()); ++k) (*at)[k] = v[k].toInt(-1);
            auto atom = [&](int i) { return i >= 0 && i < int(doc.atoms.size()); };
            if (!atom((*at)[0]) || ((*at)[1] != -1 && !atom((*at)[1]))) *at = {-1, -1};  // not this drawing's: a free end
        }
        doc.arrows.push_back(a);
    }
    for (const auto& v : root["texts"].toArray()) {
        auto o = v.toObject();
        Text t{{o["x"].toDouble(), o["y"].toDouble()}, o["text"].toString(), o["scale"].toDouble(1),
               QColor(o["color"].toString())};
        if (!(t.scale > 0)) t.scale = 1;  // files saved before #316 could hold 0
        t.compound = o["compound"].toBool();
        if (const int a = o["anchor"].toInt(-1); a >= 0 && a < int(doc.atoms.size())) t.anchor = a;
        if (!finite({t.pos.x(), t.pos.y(), t.scale})) return std::nullopt;
        doc.texts.push_back(t);
    }
    for (const auto& v : root["fills"].toArray()) {
        auto o = v.toObject();
        Fill f{{}, QColor(o["color"].toString())};
        for (const auto& i : o["atoms"].toArray()) f.atoms.push_back(i.toInt(-1));
        if (!f.color.isValid() || f.atoms.size() < 3 ||
            std::any_of(f.atoms.begin(), f.atoms.end(), [n](int i) { return i < 0 || i >= n; }))
            return std::nullopt;
        doc.fills.push_back(f);
    }
    for (const auto& v : root["brackets"].toArray()) {
        auto o = v.toObject();
        Bracket b{{}, !o["round"].toBool(), o["label"].toString()};
        for (const auto& i : o["atoms"].toArray()) b.atoms.push_back(i.toInt(-1));
        if (b.atoms.empty() || std::any_of(b.atoms.begin(), b.atoms.end(), [n](int i) { return i < 0 || i >= n; }))
            return std::nullopt;
        doc.brackets.push_back(b);
    }
    for (const auto& v : root["aromaticCircleOverrides"].toArray()) {
        std::vector<int> ring;
        for (const auto& i : v.toArray()) ring.push_back(i.toInt(-1));
        std::sort(ring.begin(), ring.end());
        if (ring.size() < 3 || ring.front() < 0 || ring.back() >= n ||
            std::adjacent_find(ring.begin(), ring.end()) != ring.end()) return std::nullopt;
        doc.aromaticCircleOverrides.push_back(std::move(ring));
    }
    return doc;
}

void Document::append(const Document& o, QPointF shift) {
    const int base = int(atoms.size());
    for (auto a : o.atoms) a.pos += shift, atoms.push_back(a);
    for (auto b : o.bonds) b.a += base, b.b += base, bonds.push_back(b);
    for (auto a : o.arrows) {
        a.from += shift, a.to += shift;
        for (int* i : {&a.fromAt[0], &a.fromAt[1], &a.toAt[0], &a.toAt[1]})
            if (*i >= 0) *i += base;
        arrows.push_back(a);
    }
    for (auto t : o.texts) {
        t.pos += shift;
        if (t.anchor >= 0) t.anchor += base;
        texts.push_back(t);
    }
    for (auto f : o.fills) {
        for (int& i : f.atoms) i += base;
        fills.push_back(f);
    }
    for (auto b : o.brackets) {
        for (int& i : b.atoms) i += base;
        brackets.push_back(b);
    }
    for (auto ring : o.aromaticCircleOverrides) {
        for (int& i : ring) i += base;
        aromaticCircleOverrides.push_back(std::move(ring));
    }
}

int Document::addAtom(QPointF pos, int z) {
    atoms.push_back({pos, z});
    return int(atoms.size()) - 1;
}

int Document::bondBetween(int a, int b) const {
    for (size_t i = 0; i < bonds.size(); ++i)
        if ((bonds[i].a == a && bonds[i].b == b) || (bonds[i].a == b && bonds[i].b == a))
            return int(i);
    return -1;
}

std::vector<int> Document::neighbors(int atom) const {
    std::vector<int> out;
    for (const auto& b : bonds)
        if (b.a == atom) out.push_back(b.b);
        else if (b.b == atom) out.push_back(b.a);
    return out;
}

std::vector<std::vector<int>> Document::bondsAt() const {
    std::vector<std::vector<int>> at(atoms.size());
    for (int i = 0; i < int(bonds.size()); ++i) at[bonds[i].a].push_back(i), at[bonds[i].b].push_back(i);
    return at;
}

void Document::removeBond(int bond) {
    const int a = bonds[bond].a, b = bonds[bond].b;
    const bool dropA = neighbors(a).size() == 1;
    const bool dropB = neighbors(b).size() == 1;
    bonds.erase(bonds.begin() + bond);
    std::vector<int> drop;
    if (dropA) drop.push_back(a);
    if (dropB) drop.push_back(b);
    if (!drop.empty()) removeAtoms(drop);
}

void Document::removeAtom(int atom) {
    std::vector<int> drop{atom};
    for (int n : neighbors(atom))
        if (neighbors(n).size() == 1) drop.push_back(n);
    removeAtoms(drop);
}

void Document::removeAtoms(const std::vector<int>& drop) {
    std::vector<int> remap(atoms.size(), 0);
    for (int i : drop) remap[i] = -1;
    for (auto& t : texts) {  // a compound number whose atom goes takes the nearest of its molecule's that stays
        std::vector<int> todo{t.anchor};
        for (size_t k = 0; t.anchor >= 0 && remap[t.anchor] < 0 && k < todo.size(); ++k)
            for (int nb : neighbors(todo[k]))
                if (std::find(todo.begin(), todo.end(), nb) == todo.end()) {
                    todo.push_back(nb);
                    if (remap[nb] == 0 && remap[t.anchor] < 0) t.anchor = nb;
                }
    }
    std::vector<Atom> kept;
    for (size_t i = 0; i < atoms.size(); ++i)
        if (remap[i] != -1) remap[i] = int(kept.size()), kept.push_back(atoms[i]);
    atoms = std::move(kept);
    std::erase_if(bonds, [&](const Bond& b) { return remap[b.a] < 0 || remap[b.b] < 0; });
    for (auto& b : bonds) b.a = remap[b.a], b.b = remap[b.b];
    std::erase_if(fills, [&](const Fill& f) {
        return std::any_of(f.atoms.begin(), f.atoms.end(), [&](int i) { return remap[i] < 0; });
    });
    for (auto& f : fills)
        for (int& i : f.atoms) i = remap[i];
    for (auto& b : brackets) {  // a bracket keeps around what's left of its atoms
        std::erase_if(b.atoms, [&](int i) { return remap[i] < 0; });
        for (int& i : b.atoms) i = remap[i];
    }
    std::erase_if(brackets, [](const Bracket& b) { return b.atoms.empty(); });
    std::erase_if(aromaticCircleOverrides, [&](const auto& ring) {
        return std::any_of(ring.begin(), ring.end(), [&](int i) { return remap[i] < 0; });
    });
    for (auto& ring : aromaticCircleOverrides)
        for (int& i : ring) i = remap[i];
    for (auto& a : arrows)
        for (auto* at : {&a.fromAt, &a.toAt}) {
            if ((*at)[0] < 0) continue;
            const bool gone = remap[(*at)[0]] < 0 || ((*at)[1] >= 0 && remap[(*at)[1]] < 0);  // its atom, or its bond
            *at = gone ? std::array{-1, -1} : std::array{remap[(*at)[0]], (*at)[1] < 0 ? -1 : remap[(*at)[1]]};
        }
    for (auto& t : texts)
        if (t.anchor >= 0) t.anchor = remap[t.anchor];  // -1 if gone: the number stays, free
}

// Direction pointing away from all of the atom's bonds: the bisector of the
// widest gap between them (straight on for a terminal atom).
QPointF Document::awayDirection(int atom) const {
    auto nbs = neighbors(atom);
    if (nbs.empty()) return {0, -1};
    QPointF p = atoms[atom].pos;
    std::vector<double> ang;
    for (int nb : nbs) ang.push_back(std::atan2(atoms[nb].pos.y() - p.y(), atoms[nb].pos.x() - p.x()));
    std::sort(ang.begin(), ang.end());
    double bestGap = -1, bestMid = 0;
    for (size_t i = 0; i < ang.size(); ++i) {
        double next = i + 1 < ang.size() ? ang[i + 1] : ang[0] + 2 * std::numbers::pi;
        if (next - ang[i] > bestGap + 1e-6) bestGap = next - ang[i], bestMid = (ang[i] + next) / 2;
    }
    return {std::cos(bestMid), std::sin(bestMid)};
}

static QJsonObject sheetObject(const Sheet& s) {
    QJsonObject o = QJsonDocument::fromJson(s.doc.toJson()).object();
    o["name"] = s.name;
    return o;
}

// Document settings: at the top of a version 2 file when every page has the same, else on each page.
static const char* kShared[] = {"style", "carbonLabels", "labelRatio", "showStereo",
                                "showAtomNumbers", "showShifts", "aromaticCircles"};

QByteArray sheetsToJson(const std::vector<Sheet>& sheets) {
    QJsonArray pages;
    for (const Sheet& s : sheets) {
        QJsonObject o = sheetObject(s);
        o.remove("format"), o.remove("version");
        pages.append(o);
    }
    QJsonObject root{{"format", "penzene"}, {"version", 2}};
    for (const char* key : kShared) {
        const QJsonValue first = pages.at(0).toObject().value(key);
        if (!std::all_of(pages.begin(), pages.end(), [&](const QJsonValue& p) { return p.toObject().value(key) == first; }))
            continue;
        if (!first.isUndefined()) root[key] = first;
        for (qsizetype i = 0; i < pages.size(); ++i) {
            QJsonObject o = pages.at(i).toObject();
            o.remove(key);
            pages[i] = o;
        }
    }
    root["pages"] = pages;
    return QJsonDocument(root).toJson(QJsonDocument::Indented);
}

QByteArray sheetsToJsonV1(const std::vector<Sheet>& sheets) {
    QJsonObject root = sheetObject(sheets.at(0));
    QJsonArray rest;
    for (size_t i = 1; i < sheets.size(); ++i) rest.append(sheetObject(sheets[i]));
    if (!rest.isEmpty()) root["pages"] = rest;
    return QJsonDocument(root).toJson(QJsonDocument::Indented);
}

std::vector<Sheet> sheetsFromJson(const QByteArray& data) {
    const QJsonObject root = QJsonDocument::fromJson(data).object();
    auto name = [](const QJsonObject& o, size_t n) {
        const QString s = o["name"].toString().trimmed();
        return s.isEmpty() ? QObject::tr("Page %1").arg(n) : s;
    };
    if (root["format"].toString() == "penzene" && root["version"].toInt() == 2) {
        std::vector<Sheet> sheets;
        for (const auto& v : root["pages"].toArray()) {
            QJsonObject o{{"format", "penzene"}, {"version", 1}};  // each page reads as a version 1 document
            for (const char* key : kShared)
                if (root.contains(key)) o[key] = root[key];
            const QJsonObject page = v.toObject();
            for (auto it = page.begin(); it != page.end(); ++it) o[it.key()] = it.value();
            auto doc = Document::fromJson(QJsonDocument(o).toJson());
            if (!doc) return {};  // a damaged page: refuse the file rather than drop the page on the next save
            sheets.push_back({name(page, sheets.size() + 1), *doc});
        }
        return sheets;  // none: refused, as there's nothing to open
    }
    auto first = Document::fromJson(data);
    if (!first) return {};
    std::vector<Sheet> sheets{{name(root, 1), *first}};
    for (const auto& v : root["pages"].toArray()) {
        auto doc = Document::fromJson(QJsonDocument(v.toObject()).toJson());
        if (!doc) return {};  // a damaged page: refuse the file rather than drop the page on the next save
        sheets.push_back({name(v.toObject(), sheets.size() + 1), *doc});
    }
    return sheets;
}
