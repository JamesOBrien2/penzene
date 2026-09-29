#include "Chem.h"
#include "Geometry.h"
#include "Render.h"

#include <GraphMol/CIPLabeler/CIPLabeler.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/FileParsers/CDXMLParser.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <Eigen/Dense>
#include <GraphMol/Descriptors/Crippen.h>
#include <GraphMol/Descriptors/Lipinski.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Descriptors/MolSurf.h>
#include <algorithm>
#include <map>
#include <set>
#include <type_traits>
#include <GraphMol/inchi.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/FileWriters.h>
#include <GraphMol/FileParsers/SequenceParsers.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/atomic_data.h>
#include <GraphMol/new_canon.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <QDateTime>
#include <QFile>
#include <QFileInfo>
#include <QHash>
#include <QObject>
#include <QPolygonF>
#include <QRegularExpression>
#include <QLineF>
#include <QSet>
#include <QXmlStreamReader>
#include <QXmlStreamWriter>
#include <QStringList>
#include <cmath>
#include <numbers>
#include <memory>

namespace chem {
using RDKit::RWMol;

// MOL files and RDKit use 1.5 Å bonds with y up; scenes use points with y down.
constexpr double kScale = kBondLength / 1.5;

// SMILES whose first atom is the attachment point.
static const QHash<QString, QString>& groups() {
    static const QHash<QString, QString> g{
        {"Me", "C"},        {"Et", "CC"},          {"iPr", "C(C)C"},   {"tBu", "C(C)(C)C"},
        {"Ph", "c1ccccc1"}, {"OMe", "OC"},         {"NO2", "[N+](=O)[O-]"}, {"CF3", "C(F)(F)F"},
        {"CN", "C#N"},      {"CO2Me", "C(=O)OC"},  {"CO2Et", "C(=O)OCC"}, {"CO2H", "C(=O)O"},
        {"CHO", "C=O"},     {"CH2OH", "CO"},       {"Ac", "C(C)=O"},      {"OAc", "OC(C)=O"},  {"N3", "N=[N+]=[N-]"},
        {"Boc", "C(=O)OC(C)(C)C"}, {"Cbz", "C(=O)OCc1ccccc1"}, {"Fmoc", "C(=O)OCC1c2ccccc2-c2ccccc21"},
        {"Bn", "Cc1ccccc1"}, {"Bz", "C(=O)c1ccccc1"}, {"MgBr", "[Mg]Br"}, {"SO2Me", "S(=O)(=O)C"},
        {"Ts", "S(=O)(=O)c1ccc(C)cc1"}, {"Ms", "S(=O)(=O)C"}, {"Tf", "S(=O)(=O)C(F)(F)F"},
        {"OTf", "OS(=O)(=O)C(F)(F)F"}, {"OTs", "OS(=O)(=O)c1ccc(C)cc1"}, {"TMS", "[Si](C)(C)C"},
        {"TBS", "[Si](C)(C)C(C)(C)C"}, {"OTBS", "O[Si](C)(C)C(C)(C)C"}, {"PMB", "Cc1ccc(OC)cc1"},
        {"Bpin", "B1OC(C)(C)C(C)(C)O1"}, {"nBu", "CCCC"}, {"Pr", "CCC"}, {"Cy", "C1CCCCC1"},
        {"COOH", "C(=O)O"}, {"Bu", "CCCC"}, {"iBu", "CC(C)C"}, {"OEt", "OCC"}, {"NHBoc", "NC(=O)OC(C)(C)C"},
        {"NMe2", "N(C)C"}, {"SO3H", "S(=O)(=O)O"}, {"TIPS", "[Si](C(C)C)(C(C)C)C(C)C"}, {"MOM", "COC"},
        {"THP", "C1CCCCO1"},
    };
    return g;
}

QStringList abbreviations() { return groups().keys(); }

std::optional<Atom> abbreviationHead(const QString& label) {
    auto it = groups().find(label);
    if (it == groups().end()) return std::nullopt;
    auto frag = fromSmiles(it->toStdString());
    return frag ? std::optional(frag->atoms[0]) : std::nullopt;
}

bool attach(Document& doc, int at, const std::string& what) {
    auto frag = fromSmiles(groups().value(QString::fromStdString(what), QString::fromStdString(what)).toStdString());
    if (!frag || frag->atoms.empty()) return false;
    doc.atoms[at].z = frag->atoms[0].z;
    doc.atoms[at].charge = frag->atoms[0].charge;
    doc.atoms[at].isotope = frag->atoms[0].isotope;
    doc.atoms[at].label.clear();
    if (frag->atoms.size() == 1) return true;

    // Rotate so the fragment's bulk points away from `at`'s bonds.
    QPointF origin = frag->atoms[0].pos, bulk;
    for (size_t i = 1; i < frag->atoms.size(); ++i) bulk += frag->atoms[i].pos;
    bulk = bulk / double(frag->atoms.size() - 1) - origin;
    QPointF want = doc.neighbors(at).empty() ? QPointF(1, 0) : doc.awayDirection(at);
    double ang = std::atan2(want.y(), want.x()) - std::atan2(bulk.y(), bulk.x());
    const int base = int(doc.atoms.size()) - 1;  // frag atom i -> base + i
    for (size_t i = 1; i < frag->atoms.size(); ++i) {
        QPointF r = frag->atoms[i].pos - origin;
        Atom a = frag->atoms[i];
        a.pos = doc.atoms[at].pos + QPointF(r.x() * std::cos(ang) - r.y() * std::sin(ang),
                                            r.x() * std::sin(ang) + r.y() * std::cos(ang));
        doc.atoms.push_back(a);
    }
    for (auto b : frag->bonds) {
        b.a = b.a == 0 ? at : base + b.a;
        b.b = b.b == 0 ? at : base + b.b;
        doc.bonds.push_back(b);
    }
    return true;
}

Document expanded(const Document& doc) {
    Document out = doc;
    for (int i = 0; i < int(doc.atoms.size()); ++i) {
        if (doc.atoms[i].label.isEmpty()) continue;
        const auto head = abbreviationHead(doc.atoms[i].label);
        if (!attach(out, i, doc.atoms[i].label.toStdString()))
            out.atoms[i].label.clear();  // unknown label (e.g. from a newer file): keep the atom as is
        else if (head)
            out.atoms[i].charge += doc.atoms[i].charge - head->charge;  // a charged group: N3- (#370)
    }
    return out;
}

// Abbreviations are expanded first (unless `expand` is false, for layout);
// atom i of `doc` is atom i of the mol.
// Projections read as the wedge drawings they stand for, for chemistry only:
// - Fischer: an atom whose four drawn bonds go straight up, down, left and right,
//   none wedged: the horizontal ones come toward the viewer.
// - Haworth: a flattened 5- or 6-ring with a bold (front) edge; a substituent
//   drawn straight up is toward the viewer seen from above, straight down away.
Document projectionsAsWedges(const Document& in) {
    Document doc = in;
    const BondsAt at = in.bondsAt();
    const double tol = std::sin(qDegreesToRadians(4.0));
    auto dir = [&](int from, int to) { return unit(doc.atoms[to].pos - doc.atoms[from].pos); };
    for (int i = 0; i < int(doc.atoms.size()); ++i) {
        if (at[i].size() != 4) continue;
        const auto nbs = neighbors(in, at, i);
        if (nbs.size() != 4) continue;
        int horizontal = 0, vertical = 0;
        bool wedged = false;
        for (int nb : nbs) {
            const QPointF d = dir(i, nb);
            horizontal += std::abs(d.y()) < tol, vertical += std::abs(d.x()) < tol;
            const BondStereo s = in.bonds[in.bondBetween(i, nb)].stereo;
            wedged |= s == BondStereo::Wedge || s == BondStereo::Hash;
        }
        if (horizontal != 2 || vertical != 2 || wedged) continue;
        for (int nb : nbs)
            if (std::abs(dir(i, nb).y()) < tol) {
                Bond& b = doc.bonds[in.bondBetween(i, nb)];
                b.a = i, b.b = nb, b.stereo = BondStereo::Wedge;
            }
    }
    // Haworth rings: from each bold bond, the smallest cycle back through it.
    for (const Bond& front : in.bonds) {
        if (front.stereo != BondStereo::Bold) continue;
        std::vector<int> prev(in.atoms.size(), -1);
        std::vector<int> queue{front.b};
        prev[front.b] = front.b;
        for (size_t q = 0; q < queue.size() && prev[front.a] < 0; ++q)
            for (int nb : neighbors(in, at, queue[q]))
                if (prev[nb] < 0 && !(queue[q] == front.b && nb == front.a)) prev[nb] = queue[q], queue.push_back(nb);
        if (prev[front.a] < 0) continue;
        std::vector<int> ring{front.a};
        while (ring.back() != front.b) ring.push_back(prev[ring.back()]);
        if (ring.size() < 5 || ring.size() > 6) continue;
        QPolygonF poly;
        for (int i : ring) poly << in.atoms[i].pos;
        const QRectF box = poly.boundingRect();
        if (box.height() > 0.75 * box.width()) continue;  // not flattened: an ordinary ring
        const QPointF centre = box.center();
        for (int i : ring) {
            std::vector<int> up, down;
            for (int nb : neighbors(in, at, i)) {
                if (std::find(ring.begin(), ring.end(), nb) != ring.end()) continue;
                const QPointF d = dir(i, nb);
                if (std::abs(d.x()) < 0.35 * std::abs(d.y())) (d.y() < 0 ? up : down).push_back(nb);
            }
            if (up.size() > 1 || down.size() > 1 || up.size() + down.size() == 0) continue;
            const QPointF out = unit(in.atoms[i].pos - centre);
            const bool both = up.size() + down.size() == 2;
            for (int nb : up) {
                doc.atoms[nb].pos = in.atoms[i].pos + rotated(out, both ? 25 : 0) * kBondLength;
                Bond& b = doc.bonds[in.bondBetween(i, nb)];
                b.a = i, b.b = nb, b.stereo = BondStereo::Wedge;
            }
            for (int nb : down) {
                doc.atoms[nb].pos = in.atoms[i].pos + rotated(out, both ? -25 : 0) * kBondLength;
                Bond& b = doc.bonds[in.bondBetween(i, nb)];
                b.a = i, b.b = nb, b.stereo = BondStereo::Hash;
            }
        }
    }
    return doc;
}

static std::unique_ptr<RWMol> toRDKit(const Document& in, bool expand = true) {
    const Document doc = projectionsAsWedges(expand ? expanded(in) : in);
    auto mol = std::make_unique<RWMol>();
    auto* conf = new RDKit::Conformer(doc.atoms.size());
    for (size_t i = 0; i < doc.atoms.size(); ++i) {
        const auto& a = doc.atoms[i];
        auto* atom = new RDKit::Atom(a.z);
        atom->setFormalCharge(a.charge);
        atom->setNumRadicalElectrons(a.radicals);
        if (a.isotope && a.z > 0) atom->setIsotope(a.isotope);
        // Not setAtomMapNum: it logs through rdErrorLog, which the Windows DLL doesn't export.
        if (a.map > 0) atom->setProp(RDKit::common_properties::molAtomMapNumber, a.map);
        // Generic atoms (expansion drops their labels, so they come from `in`): R1…Rn as
        // MDL R-groups (R# with RGP, [n*] in SMILES), any other text (X, Ar) as an MDL atom alias.
        // ponytail: V3000 has no alias block, so there X/Ar become plain * atoms.
        if (a.z == 0 && i < in.atoms.size() && !in.atoms[i].label.isEmpty()) {
            static const QRegularExpression rgroup("^R(\\d+)$");
            const QString label = in.atoms[i].label;
            if (const auto m = rgroup.match(label); m.hasMatch()) {
                const unsigned n = m.captured(1).toUInt();
                atom->setProp(RDKit::common_properties::_MolFileRLabel, n);
                atom->setIsotope(n);
            } else {
                atom->setProp(RDKit::common_properties::molFileAlias, label.toStdString());
            }
        }
        mol->addAtom(atom, true, true);
        conf->setAtomPos(i, {a.pos.x() / kScale, -a.pos.y() / kScale, 0});
    }
    for (const auto& b : doc.bonds) {
        const int order = chemicalOrder(b);
        if (order < 1) continue;  // drawn only
        auto type = order == 2 ? RDKit::Bond::DOUBLE
                    : order == 3 ? RDKit::Bond::TRIPLE
                                 : RDKit::Bond::SINGLE;
        mol->addBond(b.a, b.b, type);
        if (b.stereo == BondStereo::Wedge || b.stereo == BondStereo::Hash || b.stereo == BondStereo::Wavy) {
            auto* bond = mol->getBondBetweenAtoms(b.a, b.b);
            bond->setBondDir(b.stereo == BondStereo::Wedge  ? RDKit::Bond::BEGINWEDGE
                             : b.stereo == BondStereo::Hash ? RDKit::Bond::BEGINDASH
                                                            : RDKit::Bond::UNKNOWN);
            bond->setProp(RDKit::common_properties::_MolFileBondStereo,
                          b.stereo == BondStereo::Wedge ? 1u : b.stereo == BondStereo::Hash ? 6u : 4u);
        }
    }
    // Stereo groups: all abs centres in one, each &n and orn its own, n kept as the group's id.
    std::map<std::pair<StereoGroup, int>, std::vector<RDKit::Atom*>> groups;
    for (size_t i = 0; i < in.atoms.size(); ++i)
        if (const Atom& a = in.atoms[i]; a.stereoGroup != StereoGroup::None)
            groups[{a.stereoGroup, a.stereoGroup == StereoGroup::Abs ? 0 : a.stereoGroupNumber}].push_back(mol->getAtomWithIdx(i));
    std::vector<RDKit::StereoGroup> stereoGroups;
    for (auto& [key, atoms] : groups) {
        const auto type = key.first == StereoGroup::Abs  ? RDKit::StereoGroupType::STEREO_ABSOLUTE
                          : key.first == StereoGroup::And ? RDKit::StereoGroupType::STEREO_AND
                                                          : RDKit::StereoGroupType::STEREO_OR;
        stereoGroups.emplace_back(type, atoms, std::vector<RDKit::Bond*>{}, unsigned(key.second));
        stereoGroups.back().setWriteId(key.second);
    }
    mol->setStereoGroups(std::move(stereoGroups));
    conf->set3D(false);
    mol->addConformer(conf, true);
    return mol;
}

// Best-effort sanitize + stereo from wedges. Drawings in progress are often
// invalid (pentavalent carbon), so failure just leaves a partly perceived mol.
static bool perceive(RWMol& mol, bool aromatic = true) {
    try {
        unsigned failed = 0;
        RDKit::MolOps::sanitizeMol(mol, failed,
                                   aromatic ? RDKit::MolOps::SANITIZE_ALL
                                            : RDKit::MolOps::SANITIZE_ALL ^ RDKit::MolOps::SANITIZE_SETAROMATICITY);
    } catch (...) {
        mol.updatePropertyCache(false);
        RDKit::MolOps::fastFindRings(mol);
        return false;
    }
    RDKit::MolOps::assignChiralTypesFromBondDirs(mol);  // stereocentres from wedges
    if (mol.getNumConformers()) RDKit::MolOps::detectBondStereochemistry(mol);  // E/Z from the drawing
    RDKit::MolOps::assignStereochemistry(mol, true, true);
    return true;
}

// What RDKit's own depictor crowds (#502): a macrocycle (a cyclophane's bridge crossed a ring)
// or a peptide (each C=O ran into the next NH). CoordGen lays these out cleanly, but turns and
// crowds everyday molecules, so it's only for them.
static bool forCoordGen(RWMol& mol) {
    if (!mol.getRingInfo()->isInitialized()) RDKit::MolOps::fastFindRings(mol);
    for (const auto& ring : mol.getRingInfo()->atomRings())
        if (ring.size() >= 9) return true;
    static const std::unique_ptr<RWMol> peptideBond(RDKit::SmartsToMol("[CX4][CX3](=O)[NX3][CX4]"));
    return RDKit::SubstructMatch(mol, *peptideBond).size() >= 2;  // a tripeptide or longer
}

// RDKit's depictor with its ring templates, or CoordGen (through it) where that does better.
static void layout(RWMol& mol) {
    RDDepict::Compute2DCoordParameters params;
    params.canonOrient = true;
    params.useRingTemplates = true;
    RDDepict::preferCoordGen = true;  // honoured only where forceRDKit is off
    params.forceRDKit = !forCoordGen(mol);
    if (params.forceRDKit) {
        RDDepict::compute2DCoords(mol, params);
    } else {  // with the H of each N-H, O-H… in a chain as an atom, so side chains keep off where its label goes (#559)
        std::vector<unsigned> hetero;
        for (const auto* a : mol.atoms())
            if (a->getAtomicNum() != 6 && a->getDegree() >= 2 && a->getTotalNumHs()) hetero.push_back(a->getIdx());
        RWMol withH(mol);
        RDKit::MolOps::addHs(withH, false, false, &hetero);
        RDDepict::compute2DCoords(withH, params);
        auto* conf = new RDKit::Conformer(mol.getNumAtoms());
        conf->set3D(false);
        for (unsigned i = 0; i < mol.getNumAtoms(); ++i) conf->setAtomPos(i, withH.getConformer().getAtomPos(i));
        mol.clearConformers();
        mol.addConformer(conf, true);
    }
    // RDKit packs separate fragments by their atoms alone, so the labels of ions and small
    // molecules run into each other (Na⁺Cl⁻, 3 H₂O): line them up left to right instead.
    std::vector<std::vector<int>> frags;
    if (RDKit::MolOps::getMolFrags(mol, frags) < 2) return;
    auto& conf = mol.getConformer();
    double x = 0;
    for (const auto& f : frags) {  // in SMILES order
        double lo = 1e9, hi = -1e9, top = 1e9, bottom = -1e9;
        for (int i : f) {
            const auto* a = mol.getAtomWithIdx(i);
            const auto& p = conf.getAtomPos(i);
            // ponytail: a label's extent guessed at 0.65 Å a character, its first letter centred on the atom.
            int chars = 0;
            if (a->getAtomicNum() != 6 || a->getFormalCharge() || !a->getDegree()) {
                const int h = int(a->getTotalNumHs()), q = std::abs(a->getFormalCharge());
                chars = int(a->getSymbol().size()) + (h ? 1 + (h > 1) : 0) + (q ? 1 + (q > 1) : 0);
            }
            lo = std::min(lo, p.x - (chars ? 0.35 : 0));
            hi = std::max(hi, p.x + (chars ? 0.65 * chars - 0.3 : 0));
            top = std::min(top, p.y), bottom = std::max(bottom, p.y);
        }
        for (int i : f) conf.getAtomPos(i) += RDGeom::Point3D(x - lo, -(top + bottom) / 2, 0);
        x += hi - lo + 0.6;
    }
}

// scale 0: normalise whatever bond length the source used to ours.
static Document fromRDKit(RWMol& mol, double scale = 0) {
    try {
        RDKit::MolOps::Kekulize(mol, true);  // draw explicit double bonds
    } catch (...) {
    }
    const auto& conf = mol.getConformer();
    if (!scale) {
        double sum = 0;
        for (const auto* b : mol.bonds())
            sum += (conf.getAtomPos(b->getBeginAtomIdx()) - conf.getAtomPos(b->getEndAtomIdx())).length();
        scale = mol.getNumBonds() && sum > 1e-6 ? kBondLength * mol.getNumBonds() / sum : kScale;
    }

    Document doc;
    for (const auto* a : mol.atoms()) {
        const auto& p = conf.getAtomPos(a->getIdx());
        QString label;  // a generic atom's: its alias (X, Ar), or R-group number
        unsigned r = 0;
        std::string alias;
        if (a->getAtomicNum() == 0 && a->getPropIfPresent(RDKit::common_properties::molFileAlias, alias) && !alias.empty())
            label = QString::fromStdString(alias);
        else if (a->getAtomicNum() == 0 && a->getPropIfPresent(RDKit::common_properties::_MolFileRLabel, r) && r)
            label = QString("R%1").arg(r);
        doc.atoms.push_back({QPointF(p.x * scale, -p.y * scale), int(a->getAtomicNum()),
                             a->getFormalCharge(), label, {}, int(a->getAtomMapNum()), 0,
                             int(std::min(2u, a->getNumRadicalElectrons()))});
        if (a->getAtomicNum() > 0) doc.atoms.back().isotope = int(a->getIsotope());  // on a dummy it's an R-group number
    }
    for (const auto* b : mol.bonds()) {
        Bond out{int(b->getBeginAtomIdx()), int(b->getEndAtomIdx())};
        out.order = b->getBondType() == RDKit::Bond::DOUBLE   ? 2
                    : b->getBondType() == RDKit::Bond::TRIPLE ? 3
                                                              : 1;
        out.stereo = b->getBondDir() == RDKit::Bond::BEGINWEDGE  ? BondStereo::Wedge
                     : b->getBondDir() == RDKit::Bond::BEGINDASH ? BondStereo::Hash
                     : b->getBondDir() == RDKit::Bond::UNKNOWN   ? BondStereo::Wavy
                                                                 : BondStereo::None;
        doc.bonds.push_back(out);
    }
    auto groups = mol.getStereoGroups();
    RDKit::assignStereoGroupIds(groups);  // numbers for groups the source left unnumbered
    for (const auto& g : groups) {
        const auto type = g.getGroupType();
        for (const auto* a : g.getAtoms()) {
            Atom& out = doc.atoms[a->getIdx()];
            out.stereoGroup = type == RDKit::StereoGroupType::STEREO_ABSOLUTE ? StereoGroup::Abs
                              : type == RDKit::StereoGroupType::STEREO_AND    ? StereoGroup::And
                                                                             : StereoGroup::Or;
            out.stereoGroupNumber = type == RDKit::StereoGroupType::STEREO_ABSOLUTE ? 0 : int(g.getWriteId());
        }
    }
    return doc;
}

std::optional<Document> fromSmiles(const std::string& smiles) {
    std::unique_ptr<RWMol> mol;
    try {
        mol.reset(RDKit::SmilesToMol(smiles));
    } catch (...) {
        return std::nullopt;
    }
    if (!mol) return std::nullopt;
    layout(*mol);
    RDKit::Chirality::wedgeMolBonds(*mol, &mol->getConformer());
    return fromRDKit(*mol);
}

std::optional<Document> fromSequence(const QString& sequence) {
    static const QHash<QString, QChar> codes{
        {"ala", 'A'}, {"arg", 'R'}, {"asn", 'N'}, {"asp", 'D'}, {"cys", 'C'}, {"gln", 'Q'}, {"glu", 'E'},
        {"gly", 'G'}, {"his", 'H'}, {"ile", 'I'}, {"leu", 'L'}, {"lys", 'K'}, {"met", 'M'}, {"phe", 'F'},
        {"pro", 'P'}, {"ser", 'S'}, {"thr", 'T'}, {"trp", 'W'}, {"tyr", 'Y'}, {"val", 'V'}};
    // Three-letter codes: Gly-Phe or GLY PHE, with H-…-OH free ends and D- or L- before a residue; a lone
    // unhyphenated one only as written (Gly), since GLY is also Gly-Leu-Tyr.
    QStringList words = sequence.split(QRegularExpression("[\\s-]+"), Qt::SkipEmptyParts);
    const bool hyphens = sequence.contains('-');
    if (hyphens && words.size() > 1 && words.first() == "H") words.removeFirst();
    if (hyphens && words.size() > 1 && words.last() == "OH") words.removeLast();
    QString letters;
    bool three = !words.isEmpty();
    for (int i = 0; i < words.size() && three; ++i) {
        const QString w = words[i].toLower();
        const bool prefix = (words[i] == "D" || words[i] == "L") && i + 1 < words.size() && codes.contains(words[i + 1].toLower());
        three = prefix || (codes.contains(w) && (hyphens || words.size() > 1 || words[i] == words[i].left(1).toUpper() + w.mid(1)));
        if (three && !prefix) letters += i > 0 && words[i - 1] == "D" ? codes[w].toLower() : codes[w];
    }
    if (!three && hyphens) return std::nullopt;  // Ala-X-Phe: refused, not read as one-letter codes
    if (!three) letters = sequence.split(QRegularExpression("\\s+"), Qt::SkipEmptyParts).join("");  // one-letter, lower case for D
    if (letters.isEmpty()) return std::nullopt;
    std::unique_ptr<RWMol> mol;
    try {
        mol.reset(RDKit::SequenceToMol(letters.toStdString(), true, true));  // lower case: D
    } catch (...) {
        return std::nullopt;
    }
    if (!mol || !mol->getNumAtoms()) return std::nullopt;
    layout(*mol);
    RDKit::Chirality::wedgeMolBonds(*mol, &mol->getConformer());
    return fromRDKit(*mol);
}

std::optional<Document> fromInchi(const std::string& inchi) {
    RDKit::ExtraInchiReturnValues rv;
    std::unique_ptr<RWMol> mol;
    try {
        mol.reset(RDKit::InchiToMol(inchi, rv));
    } catch (...) {
        return std::nullopt;
    }
    if (!mol) return std::nullopt;
    layout(*mol);
    RDKit::Chirality::wedgeMolBonds(*mol, &mol->getConformer());
    return fromRDKit(*mol);
}

// Sanitizing marks rings aromatic, and fromRDKit's Kekulize then picks its own
// double bonds. Put back the ones drawn in the file (#323): all of a molecule's
// aromatic bonds, or none when the file has no single/double form for one.
// `key` matches atoms between the sanitized and unsanitized readings.
template <class Key>
using DrawnBonds = std::map<std::pair<std::invoke_result_t<Key, const RWMol&, const RDKit::Atom*>,
                                      std::invoke_result_t<Key, const RWMol&, const RDKit::Atom*>>,
                            RDKit::Bond::BondType>;

template <class Key>
static DrawnBonds<Key> bondTypes(const RWMol& mol, Key key) {
    DrawnBonds<Key> out;
    for (const auto* b : mol.bonds()) {
        auto [it, fresh] = out.emplace(std::minmax(key(mol, b->getBeginAtom()), key(mol, b->getEndAtom())), b->getBondType());
        if (!fresh && it->second != b->getBondType()) it->second = RDKit::Bond::UNSPECIFIED;  // overlapping atoms (#365)
    }
    return out;
}

template <class Key>
static void keepDrawnKekule(RWMol& mol, const DrawnBonds<Key>& drawn, Key key) {
    std::vector<std::pair<RDKit::Bond*, RDKit::Bond::BondType>> set;
    for (auto* b : mol.bonds()) {
        if (!b->getIsAromatic()) continue;
        auto it = drawn.find(std::minmax(key(mol, b->getBeginAtom()), key(mol, b->getEndAtom())));
        if (it == drawn.end() || (it->second != RDKit::Bond::SINGLE && it->second != RDKit::Bond::DOUBLE)) return;
        set.push_back({b, it->second});
    }
    for (auto [b, type] : set) b->setBondType(type), b->setIsAromatic(false);
    for (auto* a : mol.atoms()) a->setIsAromatic(false);
}

std::optional<Document> fromMolBlock(const std::string& block) {
    std::unique_ptr<RWMol> mol;
    for (bool sanitize : {true, false}) {
        try {
            mol.reset(RDKit::MolBlockToMol(block, sanitize, false));
        } catch (...) {
        }
        if (mol) break;
    }
    if (!mol || !mol->getNumAtoms()) return std::nullopt;
    auto index = [](const RWMol&, const RDKit::Atom* a) { return int(a->getIdx()); };
    try {
        if (std::unique_ptr<RWMol> drawn(RDKit::MolBlockToMol(block, false, false)); drawn)
            keepDrawnKekule(*mol, bondTypes(*drawn, index), index);
    } catch (...) {
    }
    RDKit::Chirality::reapplyMolBlockWedging(*mol);
    return fromRDKit(*mol);
}

// A ChemDraw label node (nickname, generic R, or free text drawn as an atom).
struct LabelNode {
    int id = 0;
    QPointF pos, textPos;
    QString text;
    double textScale = 1;
    std::vector<QPointF> inner;  // positions of the atoms inside its fragment
    int headCharge = 0;  // the charge on the fragment's atom under the label (RDKit drops it)
};

// The colours the file gives atoms and bonds, by node position: RDKit reads none of them.
struct FileInks {
    std::vector<std::pair<QPointF, QColor>> atoms;
    std::vector<std::tuple<QPointF, QPointF, QColor>> bonds;  // from, to
};

// Arrows, free text and label nodes from the CDXML itself; RDKit only reads the
// molecules. Returns the label nodes and how many bonds each node id has.
// ponytail: plain lines, brackets, shapes and binary .cdx graphics are skipped.
static std::vector<LabelNode> chemDrawGraphics(const QByteArray& xml, Document& doc, QHash<int, int>& bondCount,
                                               std::vector<QPointF>& lonePairs, FileInks& inks) {
    QXmlStreamReader r(xml);
    double scale = kBondLength / 30;  // CDXML's default BondLength
    double lineWidth = 1;  // and LineWidth, in points
    struct Open { QString tag; int label = -1; };  // label: index into `labels` for label <n>s
    std::vector<Open> stack;
    std::vector<LabelNode> labels;
    auto point = [&](QStringView s) {
        auto v = s.split(' ');
        return v.size() >= 2 ? QPointF(v[0].toDouble(), v[1].toDouble()) * scale : QPointF();
    };
    // Where characters go: a free Text, or a label node's own text (not its inner fragment's).
    Text* text = nullptr;
    int labelText = -1;
    bool superscript = false;  // the current style run is superscript (face 64)
    auto inside = [&](const char* tag) {
        return std::any_of(stack.begin(), stack.end(), [&](const Open& o) { return o.tag == tag; });
    };
    // Colour numbers: 0 and 1 are the hard-wired black and white, the file's table starts at 2
    // (its first two entries are the page background and default ink).
    // Black is the ink (it follows the theme), so it comes back as no colour.
    std::vector<QColor> table;
    QHash<int, QPointF> nodePos;
    auto ink = [&](const QXmlStreamAttributes& at) {
        const int i = at.value("color").toInt() - 2;
        const QColor c = i >= 0 && i < int(table.size()) ? table[i] : QColor();
        return c == QColor(Qt::black) ? QColor() : c;
    };
    static const QStringList labelTypes{"Fragment",    "Nickname",  "GenericNickname", "Unspecified", "Anonymous",
                                        "AnonymousAlternativeGroup", "NamedAlternativeGroup", "Variable"};
    // Layers: a graphic whose Z is below every atom's is drawn behind the molecule.
    double moleculeZ = std::numeric_limits<double>::infinity();
    std::vector<double> arrowZ;
    while (!r.atEnd()) {
        auto tok = r.readNext();
        if (tok == QXmlStreamReader::EndElement) {
            if (r.name() == u"t") text = nullptr, labelText = -1;
            stack.pop_back();
            continue;
        }
        if (tok == QXmlStreamReader::Characters && !stack.empty() && stack.back().tag == "s") {
            // A superscript run comes back marked (SO4^2-) only where the plain text wouldn't
            // raise it already (Cu2+ does), so the ^ appears only where it's needed (#371).
            static const QRegularExpression charge("^[0-9+\\-\u2212]+$");
            if (text && superscript && !text->text.isEmpty() && charge.match(r.text()).hasMatch()) {
                const qsizetype from = std::max(text->text.lastIndexOf('\n'), text->text.lastIndexOf('\r')) + 1;
                const QString line = text->text.mid(from) + r.text();
                const auto sc = scripts(line);
                if (!std::all_of(sc.end() - r.text().size(), sc.end(), [](Script s) { return s == Script::Super; }))
                    text->text += '^';
            }
            if (text) text->text += r.text();
            if (labelText >= 0) labels[labelText].text += r.text();
            continue;
        }
        if (tok != QXmlStreamReader::StartElement) continue;
        const QString tag = r.name().toString();
        const auto at = r.attributes();
        auto pushArrow = [&](Arrow a) {
            a.color = ink(at);
            doc.arrows.push_back(a);
            arrowZ.push_back(at.hasAttribute("Z") ? at.value("Z").toDouble() : std::numeric_limits<double>::quiet_NaN());
        };
        const int parentLabel = stack.empty() ? -1 : stack.back().label;
        stack.push_back({tag});
        const bool inLabel = std::any_of(stack.begin(), stack.end(), [](const Open& o) { return o.label >= 0; });
        if (tag == "color")
            table.push_back(QColor(qRound(at.value("r").toDouble() * 255), qRound(at.value("g").toDouble() * 255),
                                   qRound(at.value("b").toDouble() * 255)));  // 8-bit, as the colours are chosen
        if (tag == "CDXML") {
            const double bond = at.hasAttribute("BondLength") ? std::max(1.0, at.value("BondLength").toDouble()) : 30;
            scale = kBondLength / bond;
            if (const double w = at.value("LineWidth").toDouble(); w > 0) lineWidth = w;
            // The file's labels relative to its bonds, which may be far from our styles' (#203).
            if (at.hasAttribute("LabelSize")) doc.labelRatio = at.value("LabelSize").toDouble() / bond;
        } else if (tag == "b") {
            const int from = at.value("B").toInt(), to = at.value("E").toInt();
            ++bondCount[from], ++bondCount[to];
            if (const QColor c = ink(at); c.isValid() && !inLabel && nodePos.contains(from) && nodePos.contains(to))
                inks.bonds.push_back({nodePos[from], nodePos[to], c});
        } else if (tag == "n") {
            if (at.hasAttribute("Z")) moleculeZ = std::min(moleculeZ, at.value("Z").toDouble());
            const int id = at.value("id").toInt();
            nodePos[id] = point(at.value("p"));
            if (const QColor c = ink(at); c.isValid() && !inLabel) inks.atoms.push_back({nodePos[id], c});
            for (const Open& o : stack)  // an atom inside a label's fragment (its link out isn't one)
                if (o.label >= 0 && at.value("NodeType") != u"ExternalConnectionPoint") {
                    LabelNode& l = labels[o.label];
                    l.inner.push_back(point(at.value("p")));
                    if (len(l.inner.back() - l.pos) < 0.6) l.headCharge = at.value("Charge").toInt();
                }
            if (labelTypes.contains(at.value("NodeType").toString())) {
                labels.push_back({id, point(at.value("p"))});
                stack.back().label = int(labels.size()) - 1;
            }
        } else if (tag == "t" && parentLabel >= 0) {
            labels[parentLabel].textPos = point(at.value("p"));
            auto box = at.value("BoundingBox").split(' ');
            if (box.size() == 4) labels[parentLabel].textPos.setX(std::min(box[0].toDouble(), box[2].toDouble()) * scale);
            labelText = parentLabel;
        }
        if (tag == "s") superscript = at.value("face").toInt() & 64;
        if (tag == "s" && text && !text->color.isValid()) text->color = ink(at);
        if (tag == "s" && (text || labelText >= 0) && at.hasAttribute("size")) {
            const double rel = at.value("size").toDouble() * scale / 10;  // 10 pt: the default (ACS) label size
            if (rel > 0) (text ? text->scale : labels[labelText].textScale) = rel;  // size="0" would save unopenable (#316)
        } else if (tag == "t" && !inside("n") && !inside("fragment")) {
            // p is the first baseline; the bounding box gives the left edge whatever the justification.
            QPointF p = point(at.value("p"));
            auto box = at.value("BoundingBox").split(' ');
            if (box.size() == 4) p.setX(std::min(box[0].toDouble(), box[2].toDouble()) * scale);
            doc.texts.push_back({p, {}});
            text = &doc.texts.back();
            text->color = ink(at);
        } else if (tag == "graphic") {
            // Plain lines, boxes and ellipses (not filled ones), and orbitals. A graphic ChemDraw
            // marks SupersededBy is the old copy of an <arrow> read above: skip it.
            const auto type = at.value("GraphicType");
            if (type == u"Orbital") {  // s, p, a single lobe or a hybrid; not ovals or d orbitals
                // Centre, the main lobe's tip (an s orbital's rim), and a half-width we don't keep.
                // ChemDraw's Filled is our solid Shaded; its Shaded is our Gradient.
                const auto kind = at.value("OrbitalType");
                Arrow a{point(at.value("Center3D")), point(at.value("MajorAxisEnd3D"))};
                a.kind = kind.startsWith(u"s") ? ArrowKind::SOrbital
                         : kind.startsWith(u"lobe") ? ArrowKind::Lobe
                         : kind.startsWith(u"p") ? ArrowKind::POrbital
                         : kind.startsWith(u"hybrid") ? ArrowKind::HybridOrbital
                                                      : ArrowKind::Line;
                a.look = kind.endsWith(u"Filled") ? OrbitalLook::Shaded
                         : kind.endsWith(u"Shaded") ? OrbitalLook::Gradient
                                                    : OrbitalLook::Outline;
                if (isOrbital(a.kind) && len(a.to - a.from) > 1e-6) pushArrow(a);
                continue;
            }
            if (at.hasAttribute("SupersededBy") || !at.value("ArrowType").isEmpty() || at.value("OvalType").contains(u"Filled"))
                continue;
            auto box = at.value("BoundingBox").split(' ');
            if (box.size() != 4) continue;
            if (type == u"Symbol") {  // charges and radical dots are the atoms' own; lone pairs are marks
                if (at.value("SymbolType") == u"LonePair") lonePairs.push_back(QPointF(box[0].toDouble(), box[1].toDouble()) * scale);
                continue;
            }
            Arrow a{QPointF(box[0].toDouble(), box[1].toDouble()) * scale, QPointF(box[2].toDouble(), box[3].toDouble()) * scale};
            a.dashed = at.value("LineType").contains(u"Dash");
            if (type == u"Line") {
                a.kind = ArrowKind::Line;
            } else if (type == u"Rectangle") {
                a.kind = at.value("RectangleType").contains(u"RoundEdge") ? ArrowKind::RoundedBox : ArrowKind::Box;
            } else if (type == u"Oval") {  // centre and two axis ends; the bounding box can be just a radius
                const QPointF c = point(at.value("Center3D")), major = point(at.value("MajorAxisEnd3D")),
                              minor = point(at.value("MinorAxisEnd3D"));
                const QPointF half(std::max(std::abs(major.x() - c.x()), std::abs(minor.x() - c.x())),
                                   std::max(std::abs(major.y() - c.y()), std::abs(minor.y() - c.y())));
                a = {c - half, c + half, ArrowKind::Ellipse, 0, {}, a.dashed};
            } else {
                continue;
            }
            pushArrow(a);
        } else if (tag == "arrow") {
            const auto head = at.value("ArrowheadHead"), tail = at.value("ArrowheadTail");
            Arrow a{point(at.value("Tail3D")), point(at.value("Head3D"))};
            a.crossed = at.value("NoGo") == u"Cross";
            // ChemDraw's head is HeadSize % of the line width long (1000 when not given).
            const double width = at.value("LineWidth").toDouble() > 0 ? at.value("LineWidth").toDouble() : lineWidth;
            const double size = (at.hasAttribute("HeadSize") ? at.value("HeadSize").toDouble() : 1000) / 100 * width * scale / kHeadLength;
            if (size > 0) a.head = std::clamp(size, 0.25, 4.0);
            if (head.isEmpty() && tail.isEmpty()) {  // a plain line
                a.kind = ArrowKind::Line;
                a.dashed = at.value("LineType").contains(u"Dash");
                pushArrow(a);
                continue;
            }
            if (!head.isEmpty() && !tail.isEmpty())
                a.kind = at.hasAttribute("ArrowShaftSpacing") ? ArrowKind::Equilibrium : ArrowKind::Resonance;
            else if (at.value("ArrowheadType") == u"Hollow")
                a.kind = ArrowKind::Retro;
            else if (head.startsWith(u"Half") || tail.startsWith(u"Half"))
                a.kind = ArrowKind::Fishhook;
            if (head.isEmpty()) std::swap(a.from, a.to);
            if (double deg = std::abs(at.value("AngularSize").toDouble()); deg > 1 && a.kind != ArrowKind::Equilibrium) {
                // Circular arc about Center3D: bend is how far its midpoint sits off the chord.
                QPointF c = point(at.value("Center3D")), mid = (a.from + a.to) / 2, d = a.to - a.from;
                double radius = std::hypot(point(at.value("MajorAxisEnd3D")).x() - c.x(),
                                           point(at.value("MajorAxisEnd3D")).y() - c.y());
                QPointF out = mid - c;
                double l = std::hypot(out.x(), out.y()), dl = std::hypot(d.x(), d.y());
                if (l > 1e-6 && dl > 1e-6) {
                    QPointF arcMid = c + out / l * radius * (deg <= 180 ? 1 : -1);
                    QPointF n(-d.y() / dl, d.x() / dl);
                    a.bend = -QPointF::dotProduct(arcMid - mid, n);
                }
            }
            pushArrow(a);
        }
    }
    for (size_t k = 0; k < arrowZ.size(); ++k) doc.arrows[k].behind = std::isfinite(moleculeZ) && arrowZ[k] < moleculeZ;
    std::erase_if(doc.texts, [](const Text& t) { return t.text.trimmed().isEmpty(); });
    for (auto& t : doc.texts) t.text = t.text.trimmed().replace('\r', '\n');
    for (auto& l : labels) l.text = l.text.trimmed().replace('\r', '\n');
    return labels;
}

// RDKit expands label nodes itself, at their fragment's own coordinates (often
// far off the page), and leaves unbonded ones as stray atoms. Put them back the
// way ChemDraw shows them (#86):
//  - an unbonded label is text (reagents such as "LiBr, acetone");
//  - a bonded label is one labelled atom, as ChemDraw draws it: a known
//    abbreviation keeps its chemistry, an unknown one (SCoA) becomes an
//    unknown group (*) rather than pretending to be its attachment atom;
//  - a label RDKit kept as one atom (R, X) keeps its text as the atom's label.
// RDKit keeps no node ids on the atoms it returns, but keeps every atom at its
// CDXML position (even inside a nickname's own fragment), so match by position (#372).
static void placeLabels(Document& doc, const std::vector<LabelNode>& labels,
                        const QHash<int, int>& bondCount) {
    std::vector<int> drop;
    auto at = [](QPointF a, QPointF b) { return std::abs(a.x() - b.x()) < 0.6 && std::abs(a.y() - b.y()) < 0.6; };
    for (const LabelNode& l : labels) {
        std::vector<int> atoms;
        for (int i = 0; i < int(doc.atoms.size()); ++i) {
            const QPointF p = doc.atoms[i].pos;
            if (at(p, l.pos) ||
                std::any_of(l.inner.begin(), l.inner.end(), [&](QPointF q) { return at(p, q); }))
                atoms.push_back(i);
        }
        const auto head = abbreviationHead(l.text);
        if (bondCount.value(l.id) == 0 && !l.inner.empty() && head) {  // a lone group (Boc, N3-): one labelled atom (#384)
            auto at0 = std::find_if(atoms.begin(), atoms.end(), [&](int i) { return at(doc.atoms[i].pos, l.pos); });
            if (at0 == atoms.end()) continue;
            Atom& a = doc.atoms[*at0];
            a.z = head->z, a.charge = l.headCharge, a.label = l.text, a.pos = l.pos;  // N3-: the file's charge
            for (int i : atoms)
                if (i != *at0) drop.push_back(i);
            continue;
        }
        if (bondCount.value(l.id) == 0) {
            if (!l.text.isEmpty()) doc.texts.push_back({l.textPos, l.text, l.textScale});
            drop.insert(drop.end(), atoms.begin(), atoms.end());
            continue;
        }
        if (atoms.empty()) continue;
        auto self = std::find_if(atoms.begin(), atoms.end(), [&](int i) {
            return at(doc.atoms[i].pos, l.pos) && l.inner.empty();
        });
        if (self != atoms.end()) {  // kept as a single atom
            doc.atoms[*self].pos = l.pos;
            if (!l.text.isEmpty()) doc.atoms[*self].label = l.text;
            continue;
        }
        const QSet<int> group(atoms.begin(), atoms.end());
        auto attach = std::find_if(atoms.begin(), atoms.end(), [&](int i) {
            auto n = doc.neighbors(i);
            return std::any_of(n.begin(), n.end(), [&](int j) { return !group.contains(j); });
        });
        if (attach == atoms.end()) continue;
        Atom& a = doc.atoms[*attach];
        a.z = head ? head->z : 0, a.charge = head ? head->charge : 0, a.label = l.text, a.pos = l.pos;
        for (int i : atoms)
            if (i != *attach) drop.push_back(i);
    }
    doc.removeAtoms(drop);
}

std::optional<Document> fromChemDraw(const QByteArray& data) {
    if (data.startsWith("\xEF\xBB\xBF")) return fromChemDraw(data.mid(3));  // UTF-8 byte order mark
    if (!data.trimmed().startsWith('<')) {  // binary CDX: read as the CDXML it converts to
        const QByteArray xml = cdxToCdxml(data);
        return xml.isEmpty() ? std::nullopt : fromChemDraw(xml);
    }
    std::vector<std::unique_ptr<RWMol>> mols;
    try {
        mols = RDKit::v2::CDXMLParser::MolsFromCDXML(data.toStdString());
    } catch (...) {
        return std::nullopt;
    }
    // The two readings share coordinates, not atom order (only one drops hydrogens).
    auto at = [](const RWMol& mol, const RDKit::Atom* a) {  // not getOwningMol: Windows links no RDKit logger
        const auto p = mol.getConformer().getAtomPos(a->getIdx());
        return std::pair{std::lround(p.x * 1000), std::lround(p.y * 1000)};
    };
    DrawnBonds<decltype(at)> drawn;  // over all fragments
    try {
        const RDKit::v2::CDXMLParser::CDXMLParserParams asDrawn(false, false, RDKit::v2::CDXMLParser::CDXMLFormat::CDXML);
        for (auto& m : RDKit::v2::CDXMLParser::MolsFromCDXML(data.toStdString(), asDrawn))
            if (m->getNumConformers()) drawn.merge(bondTypes(*m, at));
    } catch (...) {
    }
    Document doc;
    for (auto& mol : mols) {
        if (!mol->getNumConformers()) continue;
        if (!mol->getNumAtoms()) continue;
        keepDrawnKekule(*mol, drawn, at);
        RDKit::Chirality::wedgeMolBonds(*mol, &mol->getConformer());
        doc.append(fromRDKit(*mol, kScale));  // RDKit scales CDXML to 1.5 Å bonds, like MOL
    }
    if (data.trimmed().startsWith('<')) {
        Document graphics;
        QHash<int, int> bondCount;
        std::vector<QPointF> lonePairs;
        FileInks inks;
        placeLabels(doc, chemDrawGraphics(data, graphics, bondCount, lonePairs, inks), bondCount);
        auto atomAt = [&](QPointF p) {
            for (int i = 0; i < int(doc.atoms.size()); ++i)
                if (QLineF(doc.atoms[i].pos, p).length() < 0.6) return i;
            return -1;
        };
        for (const auto& [p, c] : inks.atoms)
            if (const int i = atomAt(p); i >= 0) doc.atoms[i].color = c;
        for (const auto& [from, to, c] : inks.bonds)
            for (Bond& b : doc.bonds)
                if ((atomAt(from) == b.a && atomAt(to) == b.b) || (atomAt(from) == b.b && atomAt(to) == b.a)) b.color = c;
        for (QPointF lp : lonePairs) {  // onto the nearest atom
            int best = -1;
            double bestDist = 1.4 * kBondLength;
            for (int i = 0; i < int(doc.atoms.size()); ++i)
                if (const double d = QLineF(doc.atoms[i].pos, lp).length(); d < bestDist) best = i, bestDist = d;
            if (best >= 0) doc.atoms[best].lonePairs = std::min(3, doc.atoms[best].lonePairs + 1);
        }
        doc.arrows = graphics.arrows;
        doc.labelRatio = graphics.labelRatio;
        doc.texts.insert(doc.texts.end(), graphics.texts.begin(), graphics.texts.end());
    }
    // Any dummy atom left unbonded and unlabelled is an artefact of the expansion.
    std::vector<int> strays;
    for (int i = 0; i < int(doc.atoms.size()); ++i)
        if (doc.atoms[i].z == 0 && doc.atoms[i].label.isEmpty() && doc.neighbors(i).empty()) strays.push_back(i);
    doc.removeAtoms(strays);
    if (doc.empty()) return std::nullopt;
    return doc;
}

std::vector<Record> readRecords(const QString& path) {
    const QFileInfo info(path);
    const QString base = info.completeBaseName(), ext = info.suffix().toLower();
    QFile f(path);
    if (!f.open(QIODevice::ReadOnly)) return {};
    std::vector<Record> out;
    const QString text = QString::fromUtf8(f.readAll()).remove('\r');
    const QStringList parts = ext == "sdf" ? text.split("$$$$") : text.split('\n');
    for (const QString& raw : parts) {
        const QString part = ext == "sdf" ? raw : raw.trimmed();
        if (part.trimmed().isEmpty() || part.startsWith('#')) continue;
        const QString n = QString("%1-%2").arg(base).arg(out.size() + 1);
        if (ext == "sdf") {
            // The newline after "$$$$"; the first record's leading newline is its empty title (#317).
            const QString block = &raw != &parts.front() && part.startsWith('\n') ? part.mid(1) : part;
            const QString title = block.section('\n', 0, 0).trimmed();  // the molfile's name line
            out.push_back({title.isEmpty() ? n : title, fromMolBlock(block.toStdString())});
        } else {
            const QStringList cols = part.split(QRegularExpression("\\s+"));
            if (out.empty() && cols[0].compare("SMILES", Qt::CaseInsensitive) == 0) continue;  // a header row (#329)
            const std::string first = cols[0].toStdString();
            out.push_back({cols.size() > 1 ? cols[1] : n, ext == "inchi"                        ? fromInchi(first)
                                                          : first.find('>') != std::string::npos ? fromReactionSmiles(first)
                                                                                                 : fromSmiles(first)});
        }
    }
    return out;
}

static QRectF atomBox(const Document& d) {
    QPolygonF pts;
    for (const auto& a : d.atoms) pts << a.pos;
    return pts.boundingRect();
}

// Records side by side, row by row, each centred in a cell as big as the largest.
static std::optional<Document> grid(const std::vector<Record>& records) {
    std::vector<Document> docs;
    for (const auto& r : records)
        if (r.doc && !r.doc->empty()) docs.push_back(*r.doc);
    if (docs.size() <= 1) return docs.empty() ? std::nullopt : std::optional(docs[0]);
    auto box = atomBox;
    QSizeF cell;
    for (const auto& d : docs) cell = cell.expandedTo(box(d).size());
    cell += QSizeF(2 * kBondLength, 2 * kBondLength);
    const int cols = int(std::ceil(std::sqrt(double(docs.size()))));
    Document out;
    for (size_t k = 0; k < docs.size(); ++k) {
        const QPointF centre((k % cols) * cell.width(), (k / cols) * cell.height());
        out.append(docs[k], centre - box(docs[k]).center());
    }
    return out;
}

// ---------------------------------------------------------------- reactions

std::vector<Document> molecules(const Document& doc) {
    std::vector<Document> out;
    std::vector<int> seen(doc.atoms.size(), 0);
    for (int start = 0; start < int(doc.atoms.size()); ++start) {
        if (seen[start]) continue;
        std::vector<int> stack{start}, keep;
        seen[start] = 1;
        while (!stack.empty()) {
            const int i = stack.back();
            stack.pop_back();
            keep.push_back(i);
            for (int nb : doc.neighbors(i))
                if (!seen[nb]) seen[nb] = 1, stack.push_back(nb);
        }
        Document m;
        m.atoms = doc.atoms;
        m.bonds = doc.bonds;
        std::vector<int> drop;
        std::sort(keep.begin(), keep.end());
        for (int i = 0; i < int(doc.atoms.size()); ++i)
            if (!std::binary_search(keep.begin(), keep.end(), i)) drop.push_back(i);
        m.removeAtoms(drop);
        out.push_back(std::move(m));
    }
    return out;
}

// One reaction per straight reaction or equilibrium arrow, in reading order (#334).
// A molecule nearer an arrow's line than any arrow's end is one of its agents; otherwise it's a reactant of
// the arrow whose tail is nearest ahead of it, and a product of the arrow whose
// head is nearest behind it, so B in A -> B -> C is both. A step with no reactants
// takes the previous step's products, so a scheme that wraps (B ends one row, the
// next arrow starts the next) keeps B as the link between rows (#388).
std::vector<Reaction> reactionsOf(const Document& doc) {
    std::vector<const Arrow*> arrows;
    for (const Arrow& a : doc.arrows)
        if (a.bend == 0 && (a.kind == ArrowKind::Reaction || a.kind == ArrowKind::Equilibrium) && a.from != a.to)
            arrows.push_back(&a);
    // Reading order: rows (a new one where arrows are over 2 bond lengths apart in y), then left to right.
    auto mid = [](const Arrow* a) { return (a->from + a->to) / 2; };
    std::sort(arrows.begin(), arrows.end(), [&](const Arrow* a, const Arrow* b) { return mid(a).y() < mid(b).y(); });
    std::map<const Arrow*, int> row;
    for (size_t k = 0; k < arrows.size(); ++k)
        row[arrows[k]] = k && mid(arrows[k]).y() - mid(arrows[k - 1]).y() > 2 * kBondLength ? row[arrows[k - 1]] + 1
                         : k ? row[arrows[k - 1]] : 0;
    std::stable_sort(arrows.begin(), arrows.end(), [&](const Arrow* a, const Arrow* b) {
        return std::pair{row[a], mid(a).x()} < std::pair{row[b], mid(b).x()};
    });
    std::vector<Reaction> out(arrows.size());
    for (auto& m : molecules(doc)) {
        const QPointF c = atomBox(m).center();
        int agentOf = -1, reactantOf = -1, productOf = -1;
        double agent = 1e300, reactant = 1e300, product = 1e300;
        for (int k = 0; k < int(arrows.size()); ++k) {
            const QPointF d = arrows[k]->to - arrows[k]->from;
            const double t = QPointF::dotProduct(c - arrows[k]->from, d) / QPointF::dotProduct(d, d);
            if (t >= 0 && t <= 1) {
                const QPointF off = c - (arrows[k]->from + d * t);
                if (double dist = std::hypot(off.x(), off.y()); dist < agent) agent = dist, agentOf = k;
            } else if (t < 0) {
                const QPointF off = c - arrows[k]->from;
                if (double dist = std::hypot(off.x(), off.y()); dist < reactant) reactant = dist, reactantOf = k;
            } else {
                const QPointF off = c - arrows[k]->to;
                if (double dist = std::hypot(off.x(), off.y()); dist < product) product = dist, productOf = k;
            }
        }
        if (agentOf >= 0 && agent < std::min(reactant, product)) {  // nearer that arrow than any end
            out[agentOf].agents.push_back(m);
            continue;
        }
        if (reactantOf >= 0) out[reactantOf].reactants.push_back(m);
        if (productOf >= 0) out[productOf].products.push_back(m);
    }
    for (size_t k = 1; k < out.size(); ++k)
        if (out[k].reactants.empty()) out[k].reactants = out[k - 1].products;
    return out;
}

std::optional<Reaction> reactionOf(const Document& doc) {
    auto all = reactionsOf(doc);
    return all.empty() ? std::nullopt : std::optional(all.front());
}

std::string toReactionSmiles(const std::vector<Reaction>& steps) {
    std::string out;
    for (const auto& r : steps) {
        const std::string line = toReactionSmiles(r);
        if (line.empty()) return "";
        out += (out.empty() ? "" : "\n") + line;
    }
    return out;
}

std::string toReactionSmiles(const Reaction& r) {
    bool valid = true;
    auto side = [&valid](const std::vector<Document>& ms) {
        std::string s;
        for (const auto& m : ms) {
            const std::string smiles = toSmiles(m);
            if (smiles.empty()) valid = false;
            s += (s.empty() ? "" : ".") + smiles.substr(0, smiles.find(' '));  // no CXSMILES extension mid-reaction
        }
        return s;
    };
    const std::string line = side(r.reactants) + ">" + side(r.agents) + ">" + side(r.products);
    return valid ? line : "";  // as toSmiles: an invalid molecule gives nothing, not a reaction without it
}

std::string toRxn(const Reaction& r) {
    std::string out = "$RXN\n\n  Penzene\n\n" + QString("%1%2").arg(r.reactants.size(), 3).arg(r.products.size(), 3).toStdString() + "\n";
    for (const auto* side : {&r.reactants, &r.products})
        for (const auto& m : *side) out += "$MOL\n" + toMolBlock(m);
    return out;
}

std::string toRdf(const std::vector<Reaction>& steps) {
    std::string out = "$RDFILE 1\n$DATM    " + QDateTime::currentDateTime().toString("MM/dd/yy HH:mm").toStdString() + "\n";
    for (const auto& r : steps) out += "$RFMT\n" + toRxn(r);
    return out;
}

// Reactants + … → (agents above the arrow) → products, left to right; a step whose
// reactants are the previous step's products carries on from them (#387).
Document layoutReaction(const std::vector<Reaction>& steps) {
    Document out;
    const double gap = kBondLength;
    double x = 0;
    auto place = [&](const std::vector<Document>& ms) {
        for (size_t k = 0; k < ms.size(); ++k) {
            if (k) {  // a "+" between molecules
                out.texts.push_back({{x + gap, kBondLength * 0.25}, "+"});
                x += 3 * gap;
            }
            const QRectF b = atomBox(ms[k]);
            out.append(ms[k], QPointF(x - b.left(), -b.center().y()));
            x += b.width();
        }
    };
    auto smiles = [](const std::vector<Document>& ms) {
        std::vector<std::string> s;
        for (const auto& m : ms) s.push_back(toSmiles(m));
        std::sort(s.begin(), s.end());
        return s;
    };
    for (size_t k = 0; k < steps.size(); ++k) {
        const Reaction& r = steps[k];
        if (!k || smiles(r.reactants) != smiles(steps[k - 1].products)) {
            if (k) x += 2 * gap;
            place(r.reactants);
        }
        double above = 0, width = 0;
        for (const auto& m : r.agents) width += atomBox(m).width() + gap, above = std::max(above, atomBox(m).height());
        const double len = std::max(3 * kBondLength, width + gap);
        out.arrows.push_back({{x + gap, 0}, {x + gap + len, 0}});
        double ax = x + gap + (len - width + gap) / 2;
        for (const auto& m : r.agents) {
            const QRectF b = atomBox(m);
            out.append(m, QPointF(ax - b.left(), -gap - above / 2 - b.center().y()));
            ax += b.width() + gap;
        }
        x += 2 * gap + len;
        place(r.products);
    }
    return out;
}

std::optional<Document> fromReactionSmiles(const std::string& smiles) {
    const QStringList sides = QString::fromStdString(smiles).trimmed().split('>');
    if (sides.size() != 3) return std::nullopt;
    Reaction r;
    std::vector<Document>* into[] = {&r.reactants, &r.agents, &r.products};
    for (int k = 0; k < 3; ++k)
        for (const QString& part : sides[k].split('.', Qt::SkipEmptyParts)) {
            auto m = fromSmiles(part.toStdString());
            if (!m) return std::nullopt;
            into[k]->push_back(std::move(*m));
        }
    if (r.reactants.empty() && r.products.empty()) return std::nullopt;
    return layoutReaction({r});
}

static std::optional<Reaction> readRxn(const QString& s) {
    if (!s.startsWith("$RXN")) return std::nullopt;
    const QString counts = s.section('\n', 4, 4);
    const int nr = counts.mid(0, 3).trimmed().toInt(), np = counts.mid(3, 3).trimmed().toInt();
    QStringList blocks = s.split("$MOL\n");
    blocks.removeFirst();
    if (nr + np == 0 || int(blocks.size()) < nr + np) return std::nullopt;
    Reaction r;
    for (int k = 0; k < nr + np; ++k) {
        auto m = fromMolBlock(blocks[k].toStdString());
        if (!m) return std::nullopt;
        (k < nr ? r.reactants : r.products).push_back(std::move(*m));
    }
    return r;
}

std::optional<Document> fromRxn(const std::string& text) {
    auto r = readRxn(QString::fromStdString(text).remove('\r'));
    return r ? std::optional(layoutReaction({*r})) : std::nullopt;
}

// Each $RXN record, up to the next $ line other than $MOL ($RFMT, $DTYPE, …).
std::optional<Document> fromRdf(const std::string& text) {
    QStringList records = QString::fromStdString(text).remove('\r').split("$RXN");
    records.removeFirst();
    std::vector<Reaction> steps;
    for (QString rec : records) {
        rec.prepend("$RXN");
        for (qsizetype at = rec.indexOf("\n$"); at >= 0; at = rec.indexOf("\n$", at + 1))
            if (!rec.mid(at + 1).startsWith("$MOL")) {
                rec.truncate(at + 1);
                break;
            }
        auto r = readRxn(rec);
        if (!r) return std::nullopt;
        steps.push_back(std::move(*r));
    }
    if (steps.empty()) return std::nullopt;
    return layoutReaction(steps);
}

std::optional<Document> readFile(const QString& path) {
    const QString ext = QFileInfo(path).suffix().toLower();
    if (ext == "sdf" || ext == "smi" || ext == "inchi") return grid(readRecords(path));
    QFile f(path);
    if (!f.open(QIODevice::ReadOnly)) return std::nullopt;
    const QByteArray data = f.readAll();
    if (ext == "penz") return Document::fromJson(data);
    if (ext == "rxn") return fromRxn(data.toStdString());
    if (ext == "rdf") return fromRdf(data.toStdString());
    if (ext == "png" || ext == "svg" || ext == "pdf") return Document::fromEmbedded(data);
    if (ext == "cdxml" || ext == "cdx") return fromChemDraw(data);
    return fromMolBlock(data.toStdString());
}

// CDXML in our own coordinates (BondLength = ours, y down, as ChemDraw), so
// ChemDraw and chemDrawGraphics read it back unscaled.
// ponytail: abbreviations are written expanded, free-text labels as generic
// nicknames; ChemDraw's own Fragment/Nickname nodes would keep "OMe" as a label.
QByteArray toCdxml(const Document& doc) {
    QByteArray out;
    QXmlStreamWriter w(&out);
    w.setAutoFormatting(true);
    w.writeStartDocument();
    w.writeDTD(R"(<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd">)");
    int id = 1;
    // Placed as ChemDraw places its own: across the middle of a US Letter page, an inch from the top (#443).
    const QRectF box = documentBounds(doc);
    const QPointF shift(std::round(std::max(306 - box.center().x(), 72 - box.left())), std::round(72 - box.top()));  // whole points: coordinates keep their decimals
    auto pt = [shift](QPointF p) { p += shift; return QString("%1 %2").arg(p.x(), 0, 'f', 2).arg(p.y(), 0, 'f', 2); };
    auto pt3 = [&](QPointF p) { return pt(p) + " 0"; };
    // The colour table starts at colour 2 (0 and 1 are black and white); its first two entries are
    // ChemDraw's page background and default ink, so white and black lead it and the colours used follow.
    std::vector<QColor> palette{Qt::white, Qt::black};
    auto colorNo = [&](const QColor& c) {
        auto it = std::find(palette.begin(), palette.end(), c);
        if (it == palette.end()) it = palette.insert(palette.end(), c);
        return int(it - palette.begin()) + 2;
    };
    auto paint = [&](const QColor& c) {
        if (c.isValid()) w.writeAttribute("color", QString::number(colorNo(c)));
    };
    for (const Atom& a : doc.atoms) if (a.color.isValid()) colorNo(a.color);
    for (const Bond& b : doc.bonds) if (b.color.isValid()) colorNo(b.color);
    for (const Text& tx : doc.texts) if (tx.color.isValid()) colorNo(tx.color);
    for (const Arrow& a : doc.arrows) if (a.color.isValid()) colorNo(a.color);
    w.writeStartElement("CDXML");
    w.writeAttribute("BondLength", QString::number(kBondLength));
    if (doc.labelRatio > 0) w.writeAttribute("LabelSize", QString::number(doc.labelRatio * kBondLength));
    w.writeAttribute("CreationProgram", "Penzene");
    if (palette.size() > 2) {
        w.writeStartElement("colortable");
        for (const QColor& c : palette) {
            w.writeEmptyElement("color");
            w.writeAttribute("r", QString::number(c.redF(), 'f', 4));
            w.writeAttribute("g", QString::number(c.greenF(), 'f', 4));
            w.writeAttribute("b", QString::number(c.blueF(), 'f', 4));
        }
        w.writeEndElement();
    }
    // ChemDraw draws a style run in its colour only when the run names a font.
    w.writeStartElement("fonttable");
    w.writeEmptyElement("font");
    w.writeAttribute("id", "3");
    w.writeAttribute("charset", "iso-8859-1");
    w.writeAttribute("name", "Arial");
    w.writeEndElement();
    w.writeStartElement("page");
    w.writeAttribute("id", QString::number(id++));
    // Stacking (Z, and document order for readers without it): arrows sent behind the molecule first.
    int z = 1;
    auto writeArrow = [&](const Arrow& a) {
        if (isShape(a.kind)) {  // ChemDraw graphics: Line, Rectangle, Oval
            w.writeStartElement("graphic");
            w.writeAttribute("id", QString::number(id++));
        w.writeAttribute("Z", QString::number(z++));
            w.writeAttribute("BoundingBox", pt(a.from) + " " + pt(a.to));
            if (isOrbital(a.kind)) {
                static const char* kType[] = {"s", "p", "lobe", "hybridPlus"};
                static const char* kLook[] = {"", "Filled", "Shaded"};
                const QPointF side = a.from + perp(a.to - a.from) * (a.kind == ArrowKind::SOrbital ? 1 : 0.5);
                w.writeAttribute("GraphicType", "Orbital");
                w.writeAttribute("OrbitalType", QString(kType[int(a.kind) - int(ArrowKind::SOrbital)]) +
                                                    kLook[a.kind == ArrowKind::HybridOrbital && a.look == OrbitalLook::Gradient
                                                              ? 1  // ChemDraw has no shaded hybrid
                                                              : int(a.look)]);
                w.writeAttribute("Center3D", pt3(a.from));
                w.writeAttribute("MajorAxisEnd3D", pt3(a.to));
                w.writeAttribute("MinorAxisEnd3D", pt3(side));
            } else if (a.kind == ArrowKind::Line) {
                w.writeAttribute("GraphicType", "Line");
            } else if (a.kind == ArrowKind::Ellipse) {
                const QRectF r = QRectF(a.from, a.to).normalized();
                w.writeAttribute("GraphicType", "Oval");
                w.writeAttribute("Center3D", pt3(r.center()));
                w.writeAttribute("MajorAxisEnd3D", pt3({r.right(), r.center().y()}));
                w.writeAttribute("MinorAxisEnd3D", pt3({r.center().x(), r.bottom()}));
            } else {
                w.writeAttribute("GraphicType", "Rectangle");
                if (a.kind == ArrowKind::RoundedBox) w.writeAttribute("RectangleType", "RoundEdge");
            }
            if (a.dashed) w.writeAttribute("LineType", "Dashed");
            paint(a.color);
            w.writeEndElement();
            return;
        }
        w.writeStartElement("arrow");
        w.writeAttribute("id", QString::number(id++));
        w.writeAttribute("Z", QString::number(z++));
        w.writeAttribute("Head3D", pt3(a.to));
        w.writeAttribute("Tail3D", pt3(a.from));
        switch (a.kind) {
        case ArrowKind::Reaction: w.writeAttribute("ArrowheadHead", "Full"); break;
        case ArrowKind::Retro:
            w.writeAttribute("ArrowheadHead", "Full");
            w.writeAttribute("ArrowheadType", "Hollow");
            break;
        case ArrowKind::Resonance:
            w.writeAttribute("ArrowheadHead", "Full");
            w.writeAttribute("ArrowheadTail", "Full");
            break;
        case ArrowKind::Equilibrium:
            w.writeAttribute("ArrowheadHead", "HalfLeft");
            w.writeAttribute("ArrowheadTail", "HalfLeft");
            w.writeAttribute("ArrowShaftSpacing", "4");
            break;
        case ArrowKind::Fishhook: w.writeAttribute("ArrowheadHead", "HalfLeft"); break;
        }
        if (a.kind != ArrowKind::Retro) w.writeAttribute("ArrowheadType", "Solid");
        if (a.crossed) w.writeAttribute("NoGo", "Cross");
        // At ChemDraw's 1 pt line (no LineWidth is written), a head kHeadLength long; notch and width in its proportions.
        const double head = a.head * kHeadLength * 100;
        for (auto [name, share] : {std::pair{"HeadSize", 1.0}, {"ArrowheadCenterSize", 0.875}, {"ArrowheadWidth", 0.25}})
            w.writeAttribute(name, QString::number(std::lround(head * share)));
        paint(a.color);
        if (std::abs(a.bend) > 1e-6 && a.kind != ArrowKind::Equilibrium) {
            // The circle through both ends and the arc's midpoint (bend off the chord, as read back).
            const QPointF d = a.to - a.from, mid = (a.from + a.to) / 2;
            const double c = std::hypot(d.x(), d.y()) / 2, s = std::abs(a.bend);
            const QPointF n(-d.y() / (2 * c), d.x() / (2 * c)), arcMid = mid - a.bend * n;
            const double r = (c * c + s * s) / (2 * s);
            const QPointF centre = arcMid + (mid - arcMid) / s * r;
            auto angle = [&](QPointF p) { return std::atan2(p.y() - centre.y(), p.x() - centre.x()) * 180 / std::numbers::pi; };
            auto wrap = [](double x) { return std::fmod(std::fmod(x, 360) + 360, 360); };
            double sweep = wrap(angle(a.to) - angle(a.from));
            if (wrap(angle(arcMid) - angle(a.from)) > sweep) sweep -= 360;
            w.writeAttribute("Center3D", pt3(centre));
            w.writeAttribute("MajorAxisEnd3D", pt3(centre + QPointF(r, 0)));
            w.writeAttribute("MinorAxisEnd3D", pt3(centre + QPointF(0, r)));
            w.writeAttribute("AngularSize", QString::number(-sweep, 'f', 2));  // ChemDraw: head = tail turned by -AngularSize
        }
        w.writeEndElement();
    };
    for (const Arrow& a : doc.arrows)
        if (a.behind) writeArrow(a);
    if (!doc.atoms.empty()) {
        w.writeStartElement("fragment");
        w.writeAttribute("id", QString::number(id++));
        std::vector<int> node(doc.atoms.size());  // each atom's id, for the bonds
        auto element = [&](const Atom& a) {
            if (a.z != 6) w.writeAttribute("Element", QString::number(a.z));
            if (a.charge) w.writeAttribute("Charge", QString::number(a.charge));
            if (a.radicals) w.writeAttribute("Radical", a.radicals == 1 ? "Doublet" : "Triplet");
            if (a.isotope && a.z > 0) w.writeAttribute("Isotope", QString::number(a.isotope));
            if (a.stereoGroup != StereoGroup::None) {
                w.writeAttribute("EnhancedStereoType", QStringList{"", "Absolute", "And", "Or"}[int(a.stereoGroup)]);
                if (a.stereoGroupNumber) w.writeAttribute("EnhancedStereoGroupNum", QString::number(a.stereoGroupNumber));
            }
        };
        for (size_t i = 0; i < doc.atoms.size(); ++i) {
            const Atom& a = doc.atoms[i];
            node[i] = id++;
            w.writeStartElement("n");
            w.writeAttribute("id", QString::number(node[i]));
            w.writeAttribute("Z", QString::number(z));
            w.writeAttribute("p", pt(a.pos));
            paint(a.color);
            if (a.label.isEmpty()) {
                element(a);
                w.writeEndElement();
                continue;
            }
            // An abbreviation is a ChemDraw nickname: the label, over the group it stands
            // for, joined to the rest through connection points (#333). Other text (R, X,
            // MgEt) is a generic nickname, a label with no chemistry.
            // Expanded on its own: the atom (0) and its neighbours, which set its direction (#375).
            Document group;
            group.atoms.push_back(a);
            for (int nb : doc.neighbors(int(i))) group.atoms.push_back(doc.atoms[nb]), group.bonds.push_back({0, int(group.atoms.size()) - 1});
            const int atomsBefore = int(group.atoms.size()), bondsBefore = int(group.bonds.size());
            if (const auto head = abbreviationHead(a.label); head && attach(group, 0, a.label.toStdString())) {
                group.atoms[0].charge += a.charge - head->charge;  // a charged group: N3- (#382)
                w.writeAttribute("NodeType", "Nickname");
                w.writeStartElement("fragment");
                w.writeAttribute("id", QString::number(id++));
                std::map<int, int> inner{{0, id++}};  // group atom -> id
                w.writeStartElement("n");
                w.writeAttribute("id", QString::number(inner[0]));
                w.writeAttribute("p", pt(a.pos));
                element(group.atoms[0]);
                w.writeEndElement();
                for (int k = atomsBefore; k < int(group.atoms.size()); ++k) {
                    inner[k] = id++;
                    w.writeStartElement("n");
                    w.writeAttribute("id", QString::number(inner[k]));
                    w.writeAttribute("p", pt(group.atoms[k].pos));
                    element(group.atoms[k]);
                    w.writeEndElement();
                }
                std::vector<std::pair<int, int>> links;  // connection point -> head
                for (int nb : doc.neighbors(int(i))) {
                    links.push_back({id++, inner[0]});
                    w.writeStartElement("n");
                    w.writeAttribute("id", QString::number(links.back().first));
                    w.writeAttribute("p", pt(doc.atoms[nb].pos));
                    w.writeAttribute("NodeType", "ExternalConnectionPoint");
                    w.writeEndElement();
                }
                for (int k = bondsBefore; k < int(group.bonds.size()); ++k) {
                    const Bond& b = group.bonds[k];
                    w.writeStartElement("b");
                    w.writeAttribute("id", QString::number(id++));
                    w.writeAttribute("B", QString::number(inner[b.a]));
                    w.writeAttribute("E", QString::number(inner[b.b]));
                    if (b.order > 1) w.writeAttribute("Order", QString::number(b.order));
                    w.writeEndElement();
                }
                for (auto [from, to] : links) {
                    w.writeStartElement("b");
                    w.writeAttribute("id", QString::number(id++));
                    w.writeAttribute("B", QString::number(from));
                    w.writeAttribute("E", QString::number(to));
                    w.writeEndElement();
                }
                w.writeEndElement();  // fragment
            } else {
                w.writeAttribute("NodeType", "GenericNickname");
                w.writeAttribute("GenericNickname", a.label);
            }
            w.writeStartElement("t");
            w.writeAttribute("p", pt(a.pos + QPointF(-3, 4)));
            w.writeStartElement("s");
            w.writeAttribute("font", "3");
            paint(a.color);
            w.writeCharacters(a.label);
            w.writeEndElement();
            w.writeEndElement();
            w.writeEndElement();
        }
        static const char* display[] = {nullptr, "WedgeBegin", "WedgedHashBegin", "Bold", "Dash", "Wavy", nullptr, "Dash"};
        static const char* side[] = {nullptr, "Left", "Center", "Right"};
        for (const Bond& b : doc.bonds) {
            w.writeStartElement("b");
            w.writeAttribute("id", QString::number(id++));
            w.writeAttribute("B", QString::number(node[b.a]));
            w.writeAttribute("E", QString::number(node[b.b]));
            // ChemDraw's own hydrogen-bond order, and half orders for partial bonds.
            if (b.stereo == BondStereo::Interaction) w.writeAttribute("Order", "hydrogen");
            else if (b.stereo == BondStereo::Partial) w.writeAttribute("Order", b.order == 2 ? "1.5" : "0.5");
            else if (b.order > 1) w.writeAttribute("Order", QString::number(b.order));
            if (display[int(b.stereo)]) w.writeAttribute("Display", display[int(b.stereo)]);
            if (side[int(b.position)]) w.writeAttribute("DoublePosition", side[int(b.position)]);
            paint(b.color);
            w.writeEndElement();
        }
        w.writeEndElement();
    }
    for (const Text& t : doc.texts) {
        w.writeStartElement("t");
        w.writeAttribute("id", QString::number(id++));
        w.writeAttribute("p", pt(t.pos));
        paint(t.color);  // on the text and its runs: the binary CDX keeps only the text's
        // One run per script, as ChemDraw styles them: subscript face 32, superscript 64.
        auto run = [&](Script s, const QString& chars) {
            w.writeStartElement("s");
            w.writeAttribute("font", "3");
            w.writeAttribute("size", QString::number(10 * t.scale));  // 10 pt: the ACS label size
            paint(t.color);
            if (s != Script::Base || t.compound) w.writeAttribute("face", QString::number((t.compound ? 1 : 0) | (s == Script::Sub ? 32 : s == Script::Super ? 64 : 0)));
            w.writeCharacters(chars);
            w.writeEndElement();
        };
        const QStringList lines = t.text.split('\n');
        QString chars;
        Script current = Script::Base;
        for (int li = 0; li < lines.size(); ++li) {
            const QString& line = lines[li];
            const auto sc = scripts(line);
            for (int i = 0; i < line.size(); ++i) {
                if (sc[i] == Script::Hidden) continue;
                if (sc[i] != current && !chars.isEmpty()) run(current, chars), chars.clear();
                current = sc[i], chars += line[i];
            }
            if (li + 1 < lines.size()) {  // the line break, as plain text
                if (current != Script::Base && !chars.isEmpty()) run(current, chars), chars.clear();
                current = Script::Base, chars += '\r';
                continue;  // the next line's plain text joins it
            }
            if (!chars.isEmpty()) run(current, chars);
        }
        w.writeEndElement();
    }
    ++z;  // the molecule's layer
    for (const Arrow& a : doc.arrows)
        if (!a.behind) writeArrow(a);
    w.writeEndDocument();
    return out;
}

// Binary CDX: our CDXML, converted by the ChemDraw file library.
QByteArray toCdx(const Document& doc) { return cdxmlToCdx(toCdxml(doc)); }

std::string toMolBlock(const Document& doc, bool v3000) {
    auto mol = toRDKit(doc);
    perceive(*mol, false);  // the drawn Kekulé bonds, not aromatic type 4, a query-only type (#321)
    RDKit::Chirality::reapplyMolBlockWedging(*mol);  // keep the user's wedges
    if (!v3000) mol->setStereoGroups({});  // V2000 has no place for them; RDKit would switch to V3000 (and RXN files are V2000)
    return RDKit::MolToMolBlock(*mol, true, -1, false, v3000);
}

std::string toSdf(const Document& doc, bool v3000) {
    std::string out;
    for (const Document& m : molecules(doc)) out += toMolBlock(m, v3000) + "$$$$\n";
    return out;
}

std::string toSmiles(const Document& doc) {
    auto mol = toRDKit(doc);
    if (!perceive(*mol)) return "";
    // Stereo groups (&1, or1) have no place in plain SMILES: write CXSMILES with only them then (#402).
    auto write = [](const RWMol& m) {
        return m.getStereoGroups().empty()
                   ? RDKit::MolToSmiles(m)
                   : RDKit::MolToCXSmiles(m, RDKit::SmilesWriteParams(), RDKit::SmilesWrite::CXSmilesFields::CX_ENHANCEDSTEREO);
    };
    std::string smiles = write(*mol);
    auto parses = [](const std::string& s) {
        try {
            return std::unique_ptr<RWMol>(RDKit::SmilesToMol(s)) != nullptr;
        } catch (...) {
            return false;
        }
    };
    // RDKit can call a ring with an odd charged atom ([b-2]1ccccc1) aromatic yet not
    // read that SMILES back; write it in Kekulé form then.
    if (parses(smiles)) return smiles;
    mol = toRDKit(doc);  // as drawn, without aromaticity
    if (perceive(*mol, false)) smiles = write(*mol);
    return parses(smiles) ? smiles : "";  // e.g. a hydrogen with four bonds (#266)
}

std::optional<Properties> properties(const Document& doc) {
    if (doc.atoms.empty()) return std::nullopt;
    auto mol = toRDKit(doc);
    if (!perceive(*mol)) return std::nullopt;
    return Properties{RDKit::Descriptors::calcMolFormula(*mol, true),  // isotopes apart: C[13C]H6O, CDH3
                       RDKit::Descriptors::calcAMW(*mol),
                      RDKit::Descriptors::calcExactMW(*mol)};
}

std::optional<Profile> profile(const Document& doc) {
    auto basic = properties(doc);
    if (!basic) return std::nullopt;
    auto mol = toRDKit(doc);
    if (!perceive(*mol)) return std::nullopt;
    Profile p{*basic};
    double mr = 0;
    RDKit::Descriptors::calcCrippenDescriptors(*mol, p.logP, mr);
    p.tpsa = RDKit::Descriptors::calcTPSA(*mol);
    p.hbd = int(RDKit::Descriptors::calcNumHBD(*mol));
    p.hba = int(RDKit::Descriptors::calcNumHBA(*mol));
    p.rotatable = int(RDKit::Descriptors::calcNumRotatableBonds(*mol));
    p.heavyAtoms = int(RDKit::Descriptors::calcNumHeavyAtoms(*mol));
    // Mass fractions per element, hydrogens included; masses via Atom (PeriodicTable's
    // inline methods don't link on Windows).
    RDKit::RWMol withH(*mol);
    RDKit::MolOps::addHs(withH);
    std::map<std::string, double> mass;
    for (const auto* a : withH.atoms()) mass[a->getSymbol()] += a->getMass();
    double total = 0;
    for (const auto& [s, m] : mass) total += m;
    std::vector<std::string> order;  // Hill: C, H, then alphabetical
    for (const char* s : {"C", "H"})
        if (mass.count(s)) order.push_back(s);
    for (const auto& [s, m] : mass)
        if (s != "C" && s != "H") order.push_back(s);
    for (const auto& s : order) p.elemental.push_back({s, 100 * mass[s] / total});
    p.lipinskiViolations = (p.basic.mw > 500) + (p.logP > 5) + (p.hbd > 5) + (p.hba > 10);
    p.veber = p.rotatable <= 10 && p.tpsa <= 140;
    return p;
}

// Natural isotopes per atomic number: (mass, abundance as a fraction), from RDKit's isotope table.
static const std::vector<std::vector<std::pair<double, double>>>& naturalIsotopes() {
    static const auto table = [] {
        std::vector<std::vector<std::pair<double, double>>> t(RDKit::elementNames.size() + 1);
        QByteArray all;
        for (int i = 0; RDKit::isotopesAtomData[i] != "EOS"; ++i) all += RDKit::isotopesAtomData[i];
        for (const QByteArray& line : all.split('\n')) {
            const auto f = line.simplified().split(' ');  // z, symbol, isotope, mass, abundance (%)
            bool ok = f.size() >= 5;
            const int z = ok ? f[0].toInt(&ok) : 0;
            const double abundance = ok ? f[4].toDouble(&ok) : 0;  // C locale, whatever the user's
            if (ok && abundance > 0 && z > 0 && z < int(t.size())) t[z].push_back({f[3].toDouble(), abundance});
        }
        for (auto& isotopes : t) {
            double sum = 0;
            for (const auto& [m, a] : isotopes) sum += a;
            for (auto& [m, a] : isotopes) a /= sum;
        }
        return t;
    }();
    return table;
}

static double isotopeMass(int z, unsigned isotope) {
    RDKit::Atom a(z);
    a.setIsotope(isotope);
    return a.getMass();
}

// Sticks closer than 0.001 Da become one, at their weighted mean; the faintest are dropped.
static std::vector<Peak> merged(std::vector<Peak> peaks, double floor) {
    std::sort(peaks.begin(), peaks.end(), [](const Peak& a, const Peak& b) { return a.mz < b.mz; });
    double top = 0;
    for (const Peak& p : peaks) top = std::max(top, p.intensity);
    std::vector<Peak> out;
    for (const Peak& p : peaks) {
        if (!out.empty() && p.mz - out.back().mz < 0.001) {
            Peak& q = out.back();
            q.mz = (q.mz * q.intensity + p.mz * p.intensity) / (q.intensity + p.intensity);
            q.intensity += p.intensity;
        } else {
            out.push_back(p);
        }
    }
    std::erase_if(out, [&](const Peak& p) { return p.intensity < floor * top; });
    return out;
}

static std::vector<Peak> convolve(const std::vector<Peak>& a, const std::vector<Peak>& b) {
    std::vector<Peak> out;
    for (const Peak& p : a)
        for (const Peak& q : b) out.push_back({p.mz + q.mz, p.intensity * q.intensity});
    return merged(std::move(out), 1e-7);  // pruned: the faint combinations never reach 0.1%
}

constexpr double electron = 0.000548579909;  // Da

std::vector<Peak> isotopePattern(const Document& doc, Ion ion) {
    if (doc.atoms.empty()) return {};
    auto mol = toRDKit(doc);
    if (!perceive(*mol)) return {};
    std::map<int, int> natural;  // atomic number → count
    std::vector<Peak> pattern{{0, 1}};
    int charge = 0;
    for (const auto* a : mol->atoms()) {
        if (a->getAtomicNum() == 0) return {};  // R groups have no mass
        charge += a->getFormalCharge();
        natural[1] += int(a->getTotalNumHs());
        if (a->getIsotope()) pattern[0].mz += isotopeMass(a->getAtomicNum(), a->getIsotope());
        else ++natural[a->getAtomicNum()];
    }
    if (ion == Ion::M && charge == 0) charge = 1;  // a neutral molecule is seen as M+•, one electron lighter (#556)
    if (ion == Ion::MplusH) ++natural[1], ++charge;
    if (ion == Ion::MplusNa) ++natural[11], ++charge;
    if (ion == Ion::MminusH) --natural[1], --charge;
    if (natural[1] < 0) return {};
    // ponytail: one atom at a time; square-and-multiply if thousand-atom polymers get slow.
    for (const auto& [z, n] : natural) {
        if (n == 0) continue;
        std::vector<Peak> one;
        for (const auto& [m, a] : naturalIsotopes()[z]) one.push_back({m, a});
        if (one.empty()) return {};  // no stable isotope (Tc, Pm…)
        for (int i = 0; i < n; ++i) pattern = convolve(pattern, one);
    }
    pattern = merged(std::move(pattern), 0.001);
    const double top = std::max_element(pattern.begin(), pattern.end(), [](const Peak& a, const Peak& b) {
                           return a.intensity < b.intensity;
                       })->intensity;
    for (Peak& p : pattern) {
        p.mz = (p.mz - charge * electron) / std::max(1, std::abs(charge));
        p.intensity *= 100 / top;
    }
    return pattern;
}

QString hrmsLine(const Document& doc, Ion ion) {
    if (doc.atoms.empty()) return {};
    auto mol = toRDKit(doc);
    if (!perceive(*mol)) return {};
    for (const auto* a : mol->atoms())
        if (a->getAtomicNum() == 0) return {};
    const int drawn = RDKit::MolOps::getFormalCharge(*mol);  // calcExactMW takes its electrons off already
    int charge = drawn;
    auto add = [&](int z) {
        auto* a = new RDKit::Atom(z);
        a->setNoImplicit(true);
        mol->addAtom(a, false, true);
        ++charge;
    };
    if (ion == Ion::MplusH) add(1);
    if (ion == Ion::MplusNa) add(11);
    if (ion == Ion::MminusH) {
        auto it = std::find_if(mol->atoms().begin(), mol->atoms().end(), [](const RDKit::Atom* a) { return a->getTotalNumHs() > 0; });
        if (it == mol->atoms().end()) return {};
        (*it)->setNumExplicitHs((*it)->getTotalNumHs() - 1);
        (*it)->setNoImplicit(true);
        --charge;
    }
    const bool ei = ion == Ion::M && charge == 0;  // the molecular ion M+ of electron impact
    if (ei) charge = 1;
    mol->updatePropertyCache(false);
    const double mz = (RDKit::Descriptors::calcExactMW(*mol) - (charge - drawn) * electron) / std::abs(charge);
    const char* label[] = {"[M]", "[M+H]+", "[M+Na]+", "[M-H]-"};
    QString what = label[int(ion)];
    if (ion == Ion::M) what += QString(std::abs(charge), charge > 0 ? '+' : '-');
    return QString("HRMS (%1) m/z: %2 calcd for %3 %4")
        .arg(ei ? "EI" : "ESI", what, QString::fromStdString(RDKit::Descriptors::calcMolFormula(*mol, true)))
        .arg(mz, 0, 'f', 4);
}

using Counts = std::map<int, int>;  // atomic number → how many

static QString hill(const Counts& c) {
    QString out;
    auto put = [&](int z) {
        if (c.count(z) && c.at(z) > 0) out += QString::fromStdString(symbol(z)) + (c.at(z) > 1 ? QString::number(c.at(z)) : "");
    };
    put(6), put(1);
    std::vector<std::pair<std::string, int>> rest;  // alphabetical, as Hill orders the others
    for (const auto& [z, n] : c)
        if (z != 1 && z != 6) rest.push_back({symbol(z), z});
    std::sort(rest.begin(), rest.end());
    for (const auto& r : rest) put(r.second);
    return out;
}

std::vector<EiIon> eiIons(const Document& doc) {
    if (doc.atoms.empty()) return {};
    auto mol = toRDKit(doc);
    if (!perceive(*mol) || RDKit::MolOps::getFormalCharge(*mol)) return {};  // EI ionises neutral molecules
    Counts m;
    auto add = [&](Counts& c, const RDKit::Atom* a) { ++c[a->getAtomicNum()], c[1] += int(a->getTotalNumHs()); };
    for (const auto* a : mol->atoms()) {
        if (a->getAtomicNum() == 0) return {};
        add(m, a);
    }
    auto mass = [](const Counts& c) {
        double sum = 0;
        for (const auto& [z, n] : c) {
            const auto& iso = naturalIsotopes()[z];
            if (iso.empty()) return 0.0;
            sum += n * std::max_element(iso.begin(), iso.end(), [](auto& x, auto& y) { return x.second < y.second; })->first;
        }
        return sum;
    };
    std::vector<EiIon> out;
    auto ion = [&](const Counts& c, const QString& from) {
        const QString formula = hill(c) + "+";
        if (mass(c) > 0 && std::none_of(out.begin(), out.end(), [&](const EiIon& i) { return i.formula == formula; }))
            out.push_back({mass(c) - electron, formula, from});
    };
    auto loss = [&](const Counts& lost, const QString& from) {
        Counts c = m;
        for (const auto& [z, n] : lost)
            if ((c[z] -= n) < 0) return;
        ion(c, "M − " + hill(lost) + ", " + from);
    };
    auto has = [&](const char* smarts) {
        const std::unique_ptr<RWMol> q(RDKit::SmartsToMol(smarts));
        return RDKit::SubstructMatch(*mol, *q);
    };
    ion(m, "M⁺•");
    if (!has("[CH3][CX4]").empty()) loss({{6, 1}, {1, 3}}, "methyl");
    if (!has("[CX4][OX2H]").empty()) loss({{1, 2}, {8, 1}}, "alcohol");
    if (!has("[CX3](=O)[OX2H]").empty()) loss({{1, 1}, {8, 1}}, "acid"), loss({{6, 1}, {1, 1}, {8, 2}}, "acid");
    if (!has("[CX3](=O)[OX2][CH3]").empty()) loss({{6, 1}, {1, 3}, {8, 1}}, "methyl ester");
    if (!has("[CX3](=O)[OX2][CH2][CH3]").empty()) loss({{6, 2}, {1, 5}, {8, 1}}, "ethyl ester");
    if (!has("[CX4]([CH3])([CH3])[CH3]").empty()) ion({{6, 4}, {1, 9}}, "tert-butyl"), loss({{6, 4}, {1, 9}}, "tert-butyl");
    if (!has("[Cl]").empty()) loss({{17, 1}}, "chlorine");
    if (!has("[Br]").empty()) loss({{35, 1}}, "bromine");
    if (!has("[CH3][CX3]=O").empty()) ion({{6, 2}, {1, 3}, {8, 1}}, "acetyl");
    if (!has("[cH]1[cH][cH][cH][cH]c1[CH2]").empty()) ion({{6, 7}, {1, 7}}, "benzyl (tropylium)");
    if (!has("[cH]1[cH][cH][cH][cH]c1[CX3]=O").empty()) ion({{6, 7}, {1, 5}, {8, 1}}, "benzoyl"), ion({{6, 6}, {1, 5}}, "phenyl");
    // McLafferty: a γ-H moves to the carbonyl O and the alkene beyond Cα leaves.
    for (const auto& match : has("[CX3](=O)[CX4]!@[CX4][#6;!H0]")) {
        const int alpha = match[2].second, beta = match[3].second;
        Counts alkene;
        std::set<int> seen{alpha, beta};
        std::vector<int> stack{beta};
        while (!stack.empty()) {
            const auto* a = mol->getAtomWithIdx(stack.back());
            stack.pop_back();
            add(alkene, a);
            for (const auto* n : mol->atomNeighbors(a))
                if (seen.insert(int(n->getIdx())).second) stack.push_back(int(n->getIdx()));
        }
        if (--alkene[1] >= 0) loss(alkene, "McLafferty");
    }
    return out;
}

QStringList descriptorColumns() {
    return {"id", "name", "smiles", "inchikey", "formula", "mw", "exact_mass", "clogp", "tpsa", "hbd", "hba",
            "rotatable_bonds", "heavy_atoms", "aromatic_rings", "lipinski_violations", "veber", "error"};
}

std::string descriptorsCsv(const std::vector<Record>& records, const QStringList& columns) {
    const QStringList cols = columns.isEmpty() ? descriptorColumns() : columns;
    auto quote = [](QString s) {
        static const QRegularExpression special("[\",\n]");
        return s.contains(special) ? '"' + s.replace('"', "\"\"") + '"' : s;
    };
    QString out = cols.join(',') + '\n';
    int id = 0;
    for (const auto& r : records) {
        QHash<QString, QString> v{{"id", QString::number(++id)}, {"name", r.name}};
        const bool reaction = r.doc && !r.doc->arrows.empty();  // its molecules' descriptors don't add up to one row
        const auto p = r.doc && !reaction ? profile(*r.doc) : std::nullopt;
        if (!p) {
            v["error"] = reaction ? "reaction" : r.doc ? "not valid chemistry" : "unreadable";
        } else {
            const Document& d = *r.doc;
            auto num = [](double x, int places) { return QString::number(x, 'f', places); };
            v["smiles"] = QString::fromStdString(toSmiles(d));
            v["inchikey"] = QString::fromStdString(toInchiKey(d));
            v["formula"] = QString::fromStdString(p->basic.formula);
            v["mw"] = num(p->basic.mw, 2), v["exact_mass"] = num(p->basic.exactMass, 4);
            v["clogp"] = num(p->logP, 2), v["tpsa"] = num(p->tpsa, 2);
            v["hbd"] = QString::number(p->hbd), v["hba"] = QString::number(p->hba);
            v["rotatable_bonds"] = QString::number(p->rotatable), v["heavy_atoms"] = QString::number(p->heavyAtoms);
            v["aromatic_rings"] = QString::number(aromaticRings(d).size());
            v["lipinski_violations"] = QString::number(p->lipinskiViolations);
            v["veber"] = p->veber ? "true" : "false";
        }
        QStringList row;
        for (const QString& c : cols) row << quote(v.value(c));
        out += row.join(',') + '\n';
    }
    return out.toStdString();
}

std::string toInchi(const Document& doc) {
    if (doc.atoms.empty()) return "";
    auto mol = toRDKit(doc);
    if (!perceive(*mol)) return "";
    RDKit::ExtraInchiReturnValues rv;
    return RDKit::MolToInchi(*mol, rv);
}

std::string toInchiKey(const Document& doc) {
    std::string inchi = toInchi(doc);
    return inchi.empty() ? "" : RDKit::InchiToInchiKey(inchi);
}

// One connected fragment, laid out around its old centroid.
static Document cleanFragment(const Document& doc) {
    // Abbreviations stay single nodes, as in ChemDraw: expanding a ring onto a
    // crowded atom squeezes the depictor's layout (#85).
    auto mol = toRDKit(doc, false);
    perceive(*mol, false);  // keeps the drawn Kekulé form (#322)
    layout(*mol);
    RDKit::Chirality::wedgeMolBonds(*mol, &mol->getConformer());
    Document out = fromRDKit(*mol);

    auto centroid = [](const Document& d) {
        QPointF c;
        for (const auto& a : d.atoms) c += a.pos;
        return c / double(d.atoms.size());
    };
    QPointF shift = centroid(doc) - centroid(out);
    for (auto& a : out.atoms) a.pos += shift;
    return out;
}

// Each fragment is cleaned in place, so a reaction scheme keeps its layout;
// arrows and text pass through untouched.
Document clean2D(const Document& doc, const std::vector<int>& only) {
    const int n = int(doc.atoms.size());
    std::vector<int> comp(n, -1);
    int count = 0;
    for (int s = 0; s < n; ++s) {
        if (comp[s] >= 0) continue;
        std::vector<int> stack{s};
        comp[s] = count;
        while (!stack.empty()) {
            int i = stack.back();
            stack.pop_back();
            for (int nb : doc.neighbors(i))
                if (comp[nb] < 0) comp[nb] = count, stack.push_back(nb);
        }
        ++count;
    }
    std::vector<bool> wanted(count, only.empty());
    for (int i : only)
        if (i >= 0 && i < n) wanted[comp[i]] = true;
    Document out = doc;
    out.bonds.clear();
    for (const auto& b : doc.bonds)  // molecules left alone keep their bonds as drawn
        if (!wanted[comp[b.a]]) out.bonds.push_back(b);
    for (int c = 0; c < count; ++c) {
        if (!wanted[c]) continue;
        std::vector<int> ids, drop;
        for (int i = 0; i < n; ++i) (comp[i] == c ? ids : drop).push_back(i);
        Document frag = doc;
        frag.arrows.clear(), frag.texts.clear();
        frag.removeAtoms(drop);  // keeps order: frag atom k is doc atom ids[k]
        Document clean = cleanFragment(frag);
        for (size_t k = 0; k < ids.size(); ++k) out.atoms[ids[k]].pos = clean.atoms[k].pos;
        const int m = int(ids.size());  // defensive: only bonds between this fragment's atoms
        for (auto b : clean.bonds) {
            if (b.a >= m || b.b >= m) continue;
            b.a = ids[b.a], b.b = ids[b.b];
            if (int o = doc.bondBetween(b.a, b.b); o >= 0 && chemicalOrder(doc.bonds[o]) != doc.bonds[o].order)
                continue;  // RDKit's view of a partial bond; the drawn one goes back below
            // RDKit only knows order and wedges: keep what the user chose for display.
            if (int o = doc.bondBetween(b.a, b.b); o >= 0) {
                const Bond& was = doc.bonds[o];
                const bool styled = was.stereo == BondStereo::Bold || was.stereo == BondStereo::Dashed ||
                                    was.stereo == BondStereo::Partial;
                if (styled && b.stereo == BondStereo::None) b.stereo = was.stereo;
                if (b.order == 2 && was.order == 2) b.position = was.position;
                b.color = was.color;
            }
            out.bonds.push_back(b);
        }
        // Interactions and partial bonds, which RDKit never saw, as drawn.
        for (const auto& b : doc.bonds)
            if (comp[b.a] == c && chemicalOrder(b) != b.order) out.bonds.push_back(b);
    }
    return out;
}

Document addHydrogens(const Document& doc) {
    Document out = doc;
    const int n = int(doc.atoms.size());
    if (!n) return out;
    auto mol = toRDKit(doc);
    perceive(*mol);
    const int before = int(mol->getNumAtoms());
    try {
        RDKit::MolOps::addHs(*mol, false, true);  // with 2D coordinates
    } catch (...) {
        return out;
    }
    const auto& conf = mol->getConformer();
    for (int i = before; i < int(mol->getNumAtoms()); ++i) {
        const auto* h = mol->getAtomWithIdx(i);
        const int heavy = int((*mol->getAtomNeighbors(h).first));
        if (heavy >= n || !doc.atoms[heavy].label.isEmpty()) continue;  // not onto abbreviations
        const auto& p = conf.getAtomPos(i);
        const int j = out.addAtom({p.x * kScale, -p.y * kScale}, 1);
        out.bonds.push_back({heavy, j});
    }
    return out;
}

Document removeHydrogens(const Document& doc) {
    Document out = doc;
    std::vector<int> drop;
    for (int i = 0; i < int(doc.atoms.size()); ++i) {
        const Atom& a = doc.atoms[i];
        if (a.z != 1 || a.isotope || a.charge || !a.label.isEmpty() || doc.neighbors(i).size() != 1) continue;  // D stays
        const Bond& b = doc.bonds[doc.bondBetween(i, doc.neighbors(i)[0])];
        if (b.stereo == BondStereo::Wedge || b.stereo == BondStereo::Hash) continue;  // stereo H stays
        if (doc.atoms[doc.neighbors(i)[0]].z == 1) continue;  // H2 stays
        drop.push_back(i);
    }
    out.removeAtoms(drop);
    return out;
}

std::vector<StereoLabel> stereoLabels(const Document& doc) {
    std::vector<StereoLabel> out;
    if (doc.atoms.empty()) return out;
    auto mol = toRDKit(doc);
    if (!perceive(*mol)) return out;
    try {
        RDKit::CIPLabeler::assignCIPLabels(*mol);
    } catch (...) {
        return out;  // e.g. the labeller's work limit on huge symmetric molecules
    }
    const int n = int(doc.atoms.size());
    std::string code;
    for (const auto* a : mol->atoms())
        if (int(a->getIdx()) < n && a->getPropIfPresent(RDKit::common_properties::_CIPCode, code))
            out.push_back({int(a->getIdx()), -1, QString::fromStdString(code)});
    for (const auto* b : mol->bonds()) {
        const int i = int(b->getBeginAtomIdx()), j = int(b->getEndAtomIdx());
        if (i < n && j < n && b->getPropIfPresent(RDKit::common_properties::_CIPCode, code))
            out.push_back({-1, doc.bondBetween(i, j), QString::fromStdString(code)});
    }
    return out;
}

std::vector<Problem> checkStructure(const Document& doc) {
    std::vector<Problem> out;
    const int n = int(doc.atoms.size());
    auto name = [&](int i) {
        const Atom& a = doc.atoms[i];
        return a.label.isEmpty() ? QString::fromStdString(symbol(a.z)) : a.label;
    };
    const auto info = atomInfo(doc);
    for (int i = 0; i < n; ++i)
        if (info[i].valenceError)
            out.push_back({QObject::tr("Valence error: %1 has too many bonds").arg(name(i)), {i}});
    for (int i = 0; i < n; ++i)
        if (!doc.atoms[i].label.isEmpty() && !abbreviationHead(doc.atoms[i].label))
            out.push_back({QObject::tr("Unknown label \"%1\": drawn, but treated as an unknown group")
                               .arg(doc.atoms[i].label),
                           {i}});
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            if (QLineF(doc.atoms[i].pos, doc.atoms[j].pos).length() < 0.3 * kBondLength)
                out.push_back({QObject::tr("Overlapping atoms: %1 and %2").arg(name(i), name(j)), {i, j}});

    auto mol = toRDKit(doc);
    if (n && perceive(*mol)) {
        QSet<int> centres;
        for (const auto& s : RDKit::Chirality::findPotentialStereo(*mol)) {
            if (s.type != RDKit::Chirality::StereoType::Atom_Tetrahedral || int(s.centeredOn) >= n) continue;
            centres.insert(int(s.centeredOn));
            if (s.specified == RDKit::Chirality::StereoSpecified::Unspecified)
                out.push_back({QObject::tr("Stereocentre %1 has no wedge or hash").arg(name(int(s.centeredOn))),
                               {int(s.centeredOn)}});
        }
        for (const Bond& b : doc.bonds)
            if ((b.stereo == BondStereo::Wedge || b.stereo == BondStereo::Hash) && !centres.contains(b.a))
                out.push_back({QObject::tr("Wedge or hash starts at %1, which is not a stereocentre").arg(name(b.a)),
                               {b.a, b.b}});
    }
    return out;
}

std::vector<std::vector<int>> aromaticRings(const Document& doc) {
    std::vector<std::vector<int>> out;
    const int n = int(doc.atoms.size());
    if (!n) return out;
    auto mol = toRDKit(doc);
    if (!perceive(*mol)) return out;
    for (const auto& r : mol->getRingInfo()->atomRings())
        if (std::all_of(r.begin(), r.end(), [&](int i) { return i < n && mol->getAtomWithIdx(i)->getIsAromatic(); }))
            out.push_back(r);
    return out;
}

std::vector<std::vector<int>> rings(const Document& doc) {
    auto mol = toRDKit(doc);
    std::vector<std::vector<int>> out;
    try {
        RDKit::MolOps::findSSSR(*mol);
    } catch (...) {
        return out;
    }
    const int n = int(doc.atoms.size());
    for (const auto& r : mol->getRingInfo()->atomRings())
        if (std::all_of(r.begin(), r.end(), [n](int i) { return i < n; })) out.push_back(r);
    return out;
}

std::vector<AtomInfo> atomInfo(const Document& doc) {
    auto mol = toRDKit(doc);
    std::vector<AtomInfo> info(doc.atoms.size());
    for (auto* a : mol->atoms()) {
        if (a->getIdx() >= info.size()) break;  // expanded abbreviation atoms
        auto& i = info[a->getIdx()];
        try {
            a->updatePropertyCache(true);
        } catch (...) {
            i.valenceError = true;
            continue;
        }
        i.hydrogens = int(a->getNumImplicitHs());
    }
    return info;
}

// Via Atom, not PeriodicTable: its inline methods reference a logger global
// that the Windows RDKit DLL doesn't export.
std::string symbol(int z) {
    return RDKit::Atom(z).getSymbol();
}

// A table rather than PeriodicTable::getElementName, for the same reason as symbol().
std::string elementName(int z) {
    static const char* const names[] = {
        "Hydrogen", "Helium", "Lithium", "Beryllium", "Boron", "Carbon",
        "Nitrogen", "Oxygen", "Fluorine", "Neon", "Sodium", "Magnesium",
        "Aluminium", "Silicon", "Phosphorus", "Sulfur", "Chlorine", "Argon",
        "Potassium", "Calcium", "Scandium", "Titanium", "Vanadium", "Chromium",
        "Manganese", "Iron", "Cobalt", "Nickel", "Copper", "Zinc", "Gallium",
        "Germanium", "Arsenic", "Selenium", "Bromine", "Krypton", "Rubidium",
        "Strontium", "Yttrium", "Zirconium", "Niobium", "Molybdenum",
        "Technetium", "Ruthenium", "Rhodium", "Palladium", "Silver", "Cadmium",
        "Indium", "Tin", "Antimony", "Tellurium", "Iodine", "Xenon", "Caesium",
        "Barium", "Lanthanum", "Cerium", "Praseodymium", "Neodymium",
        "Promethium", "Samarium", "Europium", "Gadolinium", "Terbium",
        "Dysprosium", "Holmium", "Erbium", "Thulium", "Ytterbium", "Lutetium",
        "Hafnium", "Tantalum", "Tungsten", "Rhenium", "Osmium", "Iridium",
        "Platinum", "Gold", "Mercury", "Thallium", "Lead", "Bismuth",
        "Polonium", "Astatine", "Radon", "Francium", "Radium", "Actinium",
        "Thorium", "Protactinium", "Uranium", "Neptunium", "Plutonium",
        "Americium", "Curium", "Berkelium", "Californium", "Einsteinium",
        "Fermium", "Mendelevium", "Nobelium", "Lawrencium", "Rutherfordium",
        "Dubnium", "Seaborgium", "Bohrium", "Hassium", "Meitnerium",
        "Darmstadtium", "Roentgenium", "Copernicium", "Nihonium", "Flerovium",
        "Moscovium", "Livermorium", "Tennessine", "Oganesson",
    };
    return z >= 1 && z <= int(std::size(names)) ? names[z - 1] : "";
}

int atomicNumber(const std::string& sym) {
    try {
        return int(RDKit::Atom(sym).getAtomicNum());
    } catch (...) {
        return 0;
    }
}

}  // namespace chem

namespace chem {

// ponytail: depth always comes from an embedding; a MOL file's own z isn't kept
// on import, so a drawn 3D structure gets a fresh conformer. Keep z to use it.
std::optional<Pose3D> pose3D(const Document& doc, const std::vector<int>& atoms) {
    if (atoms.empty()) return std::nullopt;
    auto mol = toRDKit(doc);
    perceive(*mol);  // chiral tags from the wedges, which the embedding keeps
    RWMol withH(*mol);
    try {
        RDKit::MolOps::addHs(withH);
        // ETKDGv3, set through JSON: the preset is a DLL data symbol Windows can't import.
        RDKit::DGeomHelpers::EmbedParameters params;
        RDKit::DGeomHelpers::updateEmbedParametersFromJSON(params, R"({"useExpTorsionAnglePrefs": true,
            "useBasicKnowledge": true, "ETversion": 2, "useMacrocycleTorsions": true, "useMacrocycle14config": true,
            "enforceChirality": true, "randomSeed": 12648430})");  // a fixed seed: the same pose every time
        if (RDKit::DGeomHelpers::EmbedMolecule(withH, params) < 0) {
            params.useRandomCoords = true;
            if (RDKit::DGeomHelpers::EmbedMolecule(withH, params) < 0) return std::nullopt;
        }
    } catch (...) {
        return std::nullopt;
    }
    // Kabsch: the rotation that best lays the conformer onto the drawing (z = 0).
    const auto& conf = withH.getConformer();
    const int n = int(atoms.size());
    Eigen::MatrixXd P(3, n), Q(3, n);
    for (int k = 0; k < n; ++k) {
        const auto& p = conf.getAtomPos(atoms[k]);
        P.col(k) << p.x, p.y, p.z;
        Q.col(k) << doc.atoms[atoms[k]].pos.x() / kScale, -doc.atoms[atoms[k]].pos.y() / kScale, 0;
    }
    const Eigen::Vector3d pc = P.rowwise().mean(), qc = Q.rowwise().mean();
    P.colwise() -= pc, Q.colwise() -= qc;
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(Q * P.transpose(), Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix3d fix = Eigen::Matrix3d::Identity();
    fix(2, 2) = (svd.matrixU() * svd.matrixV().transpose()).determinant() < 0 ? -1 : 1;  // a rotation, not a mirror
    const Eigen::Matrix3d R = svd.matrixU() * fix * svd.matrixV().transpose();
    Pose3D pose{atoms, {}, {}, {qc.x(), qc.y(), 0}};
    for (int k = 0; k < n; ++k) {
        const Eigen::Vector3d r = R * P.col(k);
        pose.drawn.push_back({Q(0, k), Q(1, k), r.z()});
        pose.conformer.push_back({r.x(), r.y(), r.z()});
    }
    return pose;
}

Document project3D(const Document& doc, const Pose3D& pose, double aboutX, double aboutY) {
    Document out = doc;
    const double ax = qDegreesToRadians(aboutX), ay = qDegreesToRadians(aboutY);
    const double w = std::min(1.0, std::max(std::abs(aboutX), std::abs(aboutY)) / 30);  // drawing → conformer
    for (size_t k = 0; k < pose.atoms.size(); ++k) {
        const double x = (1 - w) * pose.drawn[k][0] + w * pose.conformer[k][0];
        const double y = (1 - w) * pose.drawn[k][1] + w * pose.conformer[k][1];
        const double z = pose.conformer[k][2];
        const double y1 = y * std::cos(ax) - z * std::sin(ax), z1 = y * std::sin(ax) + z * std::cos(ax);  // about x
        const double x2 = x * std::cos(ay) + z1 * std::sin(ay);                                              // about y
        out.atoms[pose.atoms[k]].pos = QPointF((x2 + pose.centre[0]) * kScale, -(y1 + pose.centre[1]) * kScale);
    }
    // Wedges for the new view: the drawing's own stereo (chiral tags), placed by RDKit.
    auto mol = toRDKit(doc);
    if (!perceive(*mol)) return out;
    auto& conf = mol->getConformer();
    for (int i = 0; i < int(out.atoms.size()); ++i)
        conf.setAtomPos(i, {out.atoms[i].pos.x() / kScale, -out.atoms[i].pos.y() / kScale, 0});
    for (auto* b : mol->bonds()) b->setBondDir(RDKit::Bond::NONE);
    RDKit::Chirality::wedgeMolBonds(*mol, &conf);
    std::vector<bool> moved(out.atoms.size());
    for (int i : pose.atoms) moved[i] = true;
    for (auto& b : out.bonds)
        if (moved[b.a] && moved[b.b] && (b.stereo == BondStereo::Wedge || b.stereo == BondStereo::Hash)) b.stereo = BondStereo::None;
    for (const auto* b : mol->bonds()) {
        const int i = int(b->getBeginAtomIdx()), j = int(b->getEndAtomIdx());
        if (i >= int(out.atoms.size()) || j >= int(out.atoms.size()) || !moved[i]) continue;
        const int k = out.bondBetween(i, j);
        if (k < 0 || (b->getBondDir() != RDKit::Bond::BEGINWEDGE && b->getBondDir() != RDKit::Bond::BEGINDASH)) continue;
        out.bonds[k].a = i, out.bonds[k].b = j;
        out.bonds[k].stereo = b->getBondDir() == RDKit::Bond::BEGINWEDGE ? BondStereo::Wedge : BondStereo::Hash;
    }
    return out;
}

}  // namespace chem

namespace chem {

// HOSE codes (Bremser 1978), following the description of CDK's HOSECodeGenerator (default mode):
// spheres of neighbours out from the root, each sphere's branches ranked by element, bond and the
// degree of their children, with & for an atom met before (a ring closure) and "," between branches.
// cmake/nmr-table.py builds the shift table with this same function, through the Python module.
namespace {
struct HoseNode {
    int atom = -1;    // -1: the empty branch after a terminal atom
    int parent = -1;  // node index; -1: the root
    int bond = 0;     // 1-3, 4 aromatic, 0 for an empty branch
    int degree = 0;
    long score = 0;
    int ranking = 0, sortOrder = 1;
    bool stopper = false;
};

QString hoseCharge(int charge) {
    if (!charge) return {};
    if (std::abs(charge) == 1) return charge < 0 ? "-" : "+";
    return "'" + QString(charge > 0 ? "+" : "") + QString::number(charge) + "'";
}

long hoseRank(int z) {  // CDK's element ranks; other elements by their mass
    switch (z) {
    case 6: return 9000; case 8: return 8900; case 7: return 8800; case 16: return 8700; case 15: return 8600;
    case 14: return 8500; case 5: return 8400; case 9: return 8300; case 17: return 8200; case 35: return 8100;
    case 53: return 7900; case 0: return 800000;
    }
    // ponytail: rounded average mass, not the major isotope's mass number (PeriodicTable's inline
    // lookups aren't exported by the Windows DLL); it only orders rare elements, the same way everywhere.
    return 800000 - std::lround(RDKit::Atom(z).getMass());
}

std::string hoseSymbol(int z) {  // the one-letter stand-ins inside the spheres
    return z == 14 ? "Q" : z == 17 ? "X" : z == 35 ? "Y" : symbol(z);
}

// One atom's code to `spheres` spheres. `mol` has no H atoms; `ranks` are its canonical ranks.
std::string hoseCode(const RDKit::ROMol& mol, int root, int spheres, const std::vector<unsigned>& ranks) {
    static const char* delimiters[] = {"(", "/", "/", ")", "/", "/", "/", "/", "/", "/", "/", "/"};
    static const char* bondSymbols[] = {"", "", "=", "%", "*"};
    static const long bondRanks[] = {0, 0, 200000, 300000, 100000};
    auto neighbours = [&](int a) {
        std::vector<int> out;
        for (const auto* n : mol.atomNeighbors(mol.getAtomWithIdx(a))) out.push_back(int(n->getIdx()));
        return out;
    };
    auto bondType = [&](int a, int b) {
        const auto* bond = mol.getBondBetweenAtoms(a, b);
        if (bond->getIsAromatic()) return 4;
        const auto t = bond->getBondType();
        return t == RDKit::Bond::DOUBLE ? 2 : t == RDKit::Bond::TRIPLE ? 3 : 1;
    };
    std::vector<HoseNode> nodes;
    auto atomOf = [&](int parent) { return parent < 0 ? root : nodes[parent].atom; };
    auto byRank = [&](std::vector<int>& sphere) {  // the canonical order, which ties keep
        std::stable_sort(sphere.begin(), sphere.end(), [&](int x, int y) {
            const long rx = nodes[x].atom < 0 ? -1 : long(ranks[nodes[x].atom]);
            const long ry = nodes[y].atom < 0 ? -1 : long(ranks[nodes[y].atom]);
            return rx < ry;
        });
    };
    auto add = [&](int atom, int parent, int bond) {
        nodes.push_back({atom, parent, bond, atom < 0 ? 0 : int(mol.getAtomWithIdx(atom)->getDegree())});
        return int(nodes.size()) - 1;
    };
    // Breadth first, one sphere past the last written (its degrees rank the last one).
    std::vector<std::vector<int>> sphere(spheres + 1);
    for (int n : neighbours(root)) sphere[0].push_back(add(n, -1, bondType(root, n)));
    byRank(sphere[0]);
    for (int s = 1; s <= spheres; ++s) {
        for (int k : sphere[s - 1]) {
            const int atom = nodes[k].atom;
            if (atom < 0) continue;
            const auto next = neighbours(atom);
            if (next.size() == 1) {
                sphere[s].push_back(add(-1, k, 0));
                continue;
            }
            for (int n : next)
                if (n != atomOf(nodes[k].parent)) sphere[s].push_back(add(n, k, bondType(atom, n)));
        }
        byRank(sphere[s]);
    }
    for (int s = spheres; s >= 1; --s)
        for (int k : sphere[s]) nodes[nodes[k].parent].ranking += nodes[k].degree;
    std::vector<bool> visited(mol.getNumAtoms());
    auto parentOrder = [&](int k) { return nodes[k].parent < 0 ? 1 : nodes[nodes[k].parent].sortOrder; };
    for (int s = 0; s < spheres; ++s) {
        for (int k : sphere[s]) {
            HoseNode& n = nodes[k];
            n.score += n.atom >= 0 && visited[n.atom] ? 1100 : n.atom >= 0 ? hoseRank(mol.getAtomWithIdx(n.atom)->getAtomicNum()) : 1000;
            n.score += bondRanks[n.bond] + n.ranking;
        }
        for (int k : sphere[s])
            if (nodes[k].atom >= 0) visited[nodes[k].atom] = true;
        std::stable_sort(sphere[s].begin(), sphere[s].end(), [&](int x, int y) {
            if (parentOrder(x) != parentOrder(y)) return parentOrder(x) > parentOrder(y);
            return nodes[x].score > nodes[y].score;
        });
        for (size_t i = 0; i < sphere[s].size(); ++i) nodes[sphere[s][i]].sortOrder = int(sphere[s].size() - i);
    }
    const auto* r = mol.getAtomWithIdx(root);
    QString code = QString::fromStdString(symbol(r->getAtomicNum())) + "-" +
                   QString::number(r->getDegree() + r->getTotalNumHs()) + hoseCharge(r->getFormalCharge()) + ";";
    std::fill(visited.begin(), visited.end(), false);
    for (int s = 0; s < spheres; ++s) {
        if (!sphere[s].empty()) {
            int branch = atomOf(nodes[sphere[s][0]].parent);
            for (int k : sphere[s]) {
                HoseNode& n = nodes[k];
                const bool parentStops = n.parent >= 0 && nodes[n.parent].stopper;
                if (!parentStops) {
                    if (atomOf(n.parent) != branch) branch = atomOf(n.parent), code += ',';
                    code += bondSymbols[n.bond];
                    if (n.atom >= 0) {
                        const auto* a = mol.getAtomWithIdx(n.atom);
                        if (visited[n.atom]) code += '&', n.stopper = true;
                        else code += QString::fromStdString(hoseSymbol(a->getAtomicNum()));
                        code += hoseCharge(a->getFormalCharge());
                    }
                }
                if (n.atom >= 0) visited[n.atom] = true;
                if (parentStops) n.stopper = true;
            }
        }
        code += delimiters[s];
    }
    for (int s = spheres; s < 4; ++s) code += delimiters[s];
    return code.toStdString();
}

// Every atom's codes, on a copy without H atoms (explicit or not, the same codes); with `classes`,
// also each atom's symmetry class (equal for atoms the molecule can't tell apart; -1 for H atoms).
std::vector<std::vector<std::string>> hoseCodesOf(const RWMol& in, int maxSpheres, std::vector<int>* classes = nullptr) {
    std::vector<std::vector<std::string>> out(in.getNumAtoms());
    if (classes) classes->assign(in.getNumAtoms(), -1);
    RWMol mol(in);
    for (auto* a : mol.atoms()) a->setProp("penzeneIndex", int(a->getIdx()));
    try {
        RDKit::MolOps::removeAllHs(mol, false);
        mol.updatePropertyCache(false);
        std::vector<unsigned> ranks;
        RDKit::Canon::rankMolAtoms(mol, ranks, true, false, false, false, false, false);
        std::vector<unsigned> symmetry;  // ties kept: equivalent atoms share a rank
        if (classes) RDKit::Canon::rankMolAtoms(mol, symmetry, false, false, false, false, false, false);
        for (const auto* a : mol.atoms()) {
            if (classes) (*classes)[a->getProp<int>("penzeneIndex")] = int(symmetry[a->getIdx()]);
            auto& codes = out[a->getProp<int>("penzeneIndex")];
            for (int s = 1; s <= maxSpheres; ++s) codes.push_back(hoseCode(mol, int(a->getIdx()), s, ranks));
        }
    } catch (...) {
        return {};
    }
    return out;
}

// The table: sorted "<spheres> <code>\t13C\tcount\t1H\tcount" lines, found by binary search.
struct ShiftTable {
    QByteArray text;
    std::vector<std::string_view> lines;
    ShiftTable() {
        QFile f(":/nmr/hose.tsv");
        if (!f.open(QIODevice::ReadOnly)) return;
        text = f.readAll();
        std::string_view all(text.constData(), size_t(text.size()));
        for (size_t at = 0; at < all.size();) {
            const size_t end = std::min(all.find('\n', at), all.size());
            if (end > at && all[at] != '#') lines.push_back(all.substr(at, end - at));
            at = end + 1;
        }
    }
    // The 13C (column 0) or 1H (column 1) shift for a code, if the table has one.
    std::optional<double> find(const std::string& key, int column) const {
        auto keyOf = [](std::string_view line) { return line.substr(0, line.find('\t')); };
        auto it = std::lower_bound(lines.begin(), lines.end(), key,
                                   [&](std::string_view line, const std::string& k) { return keyOf(line) < k; });
        if (it == lines.end() || keyOf(*it) != key) return std::nullopt;
        const auto fields = QByteArray(it->data(), qsizetype(it->size())).split('\t');
        bool ok = false;
        const double v = fields.value(1 + 2 * column).toDouble(&ok);
        return ok ? std::optional(v) : std::nullopt;
    }
};
}  // namespace

std::vector<std::vector<std::string>> hoseCodes(const std::string& molBlock, int maxSpheres) {
    try {
        std::unique_ptr<RWMol> mol(RDKit::MolBlockToMol(molBlock, true, false));
        return mol ? hoseCodesOf(*mol, maxSpheres) : std::vector<std::vector<std::string>>{};
    } catch (...) {
        return {};
    }
}

// Molecule by molecule, so one RDKit can't read leaves the others their codes.
static std::vector<std::vector<std::string>> hoseCodes(const Document& doc, int maxSpheres, std::vector<int>* classes) {
    const int n = int(doc.atoms.size());
    std::vector<std::vector<std::string>> out(n);
    if (classes) classes->assign(n, -1);
    std::vector<int> piece(n, -1);
    int pieces = 0;
    for (int i = 0; i < n; ++i) {
        if (piece[i] >= 0) continue;
        std::vector<int> todo{i};
        piece[i] = pieces;
        while (!todo.empty()) {
            const int a = todo.back();
            todo.pop_back();
            for (int b : doc.neighbors(a))
                if (piece[b] < 0) piece[b] = pieces, todo.push_back(b);
        }
        ++pieces;
    }
    for (int p = 0; p < pieces; ++p) {
        std::vector<int> kept, others;  // the part's atom j is the drawing's kept[j]
        for (int i = 0; i < n; ++i) (piece[i] == p ? kept : others).push_back(i);
        Document part = doc;
        part.removeAtoms(others);
        auto mol = toRDKit(part);
        if (!perceive(*mol)) continue;  // half-perceived: codes the table can't match
        std::vector<int> ranks;
        auto codes = hoseCodesOf(*mol, maxSpheres, classes ? &ranks : nullptr);
        for (size_t j = 0; j < kept.size() && j < codes.size(); ++j) {  // not expanded abbreviations' atoms
            out[kept[j]] = std::move(codes[j]);
            if (classes) (*classes)[kept[j]] = p * n + ranks[j];  // ranks restart in each molecule
        }
    }
    return out;
}

std::vector<std::vector<std::string>> hoseCodes(const Document& doc, int maxSpheres) {
    return hoseCodes(doc, maxSpheres, nullptr);
}

std::vector<Shift> predictShifts(const Document& doc) {
    static const ShiftTable table;
    // Shifts follow the bonding, not where atoms are drawn, so dragging an atom reuses the last answer.
    std::string key;
    for (const Atom& a : doc.atoms)
        key += std::to_string(a.z) + ',' + std::to_string(a.charge) + ',' + std::to_string(a.radicals) + ',' +
               std::to_string(a.isotope) + ',' + a.label.toStdString() + ';';
    for (const Bond& b : doc.bonds)
        key += std::to_string(b.a) + '-' + std::to_string(b.b) + ',' + std::to_string(b.order) + ',' +
               std::to_string(int(b.stereo)) + ';';
    static thread_local std::string lastKey;
    static thread_local std::vector<Shift> last;
    if (!last.empty() && key == lastKey) return last;
    std::vector<Shift> out;
    std::vector<int> classes;
    const auto codes = hoseCodes(doc, 4, &classes);
    const auto info = atomInfo(doc);
    std::vector<int> hydrogens(codes.size());
    for (size_t i = 0; i < codes.size(); ++i) hydrogens[i] = info[i].hydrogens;
    for (const Bond& b : doc.bonds)  // drawn H atoms too
        for (auto [h, to] : {std::pair{b.a, b.b}, std::pair{b.b, b.a}})
            if (doc.atoms[h].z == 1 && size_t(to) < codes.size()) ++hydrogens[to];
    for (size_t i = 0; i < codes.size(); ++i) {
        Shift s{int(i)};
        s.hydrogens = hydrogens[i], s.symmetry = classes[i];
        if (doc.atoms[i].z == 6)  // H on O or N exchange, so split nothing and aren't split
            for (int nb : doc.neighbors(int(i)))
                if (doc.atoms[nb].z == 6 && classes[nb] != classes[i]) s.coupled += hydrogens[nb];
        auto look = [&](int column, double& value, int& spheres) {
            for (int n = int(codes[i].size()); n >= 1 && !spheres; --n)
                if (auto v = table.find(std::to_string(n) + " " + codes[i][n - 1], column)) value = *v, spheres = n;
        };
        if (codes[i].empty()) continue;
        if (doc.atoms[i].z == 6) look(0, s.carbon, s.carbonSpheres);
        if (hydrogens[i] && doc.atoms[i].z != 1) look(1, s.proton, s.protonSpheres);
        if (s.carbonSpheres || s.protonSpheres) out.push_back(s);
    }
    lastKey = key, last = out;
    return out;
}

std::vector<NmrStick> nmrSticks(const Document& doc, bool proton, const std::vector<int>& only) {
    const std::set<int> wanted(only.begin(), only.end());
    std::map<int, NmrStick> byClass;  // atoms the molecule can't tell apart make one stick
    for (const Shift& s : predictShifts(doc)) {
        if (!(proton ? s.protonSpheres : s.carbonSpheres) || (!wanted.empty() && !wanted.count(s.atom))) continue;
        NmrStick& k = byClass[s.symmetry];
        k.ppm = proton ? s.proton : s.carbon;
        k.count += proton ? s.hydrogens : 1;
        k.atoms.push_back(s.atom);
        k.weak = (proton ? s.protonSpheres : s.carbonSpheres) < 3;
        k.coupled = s.coupled;
    }
    std::vector<NmrStick> out;
    for (auto& [c, k] : byClass) out.push_back(std::move(k));
    std::sort(out.begin(), out.end(), [](const NmrStick& a, const NmrStick& b) { return a.ppm > b.ppm; });  // as a spectrum reads, left to right
    return out;
}

QString nmrLine(const Document& doc, bool proton, const std::vector<int>& only) {
    QStringList parts;
    for (const auto& k : nmrSticks(doc, proton, only))
        parts << (proton ? QString("%1 (%2, %3H)").arg(k.ppm, 0, 'f', 2).arg(k.multiplicity()).arg(k.count) : QString::number(k.ppm, 'f', 1));
    return parts.isEmpty() ? QString() : (proton ? "1H" : "13C") + QString(" NMR (predicted) δ ") + parts.join(", ") + ".";
}

}  // namespace chem
