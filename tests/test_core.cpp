#include "Chem.h"
#include "Geometry.h"
#include "Edit.h"
#include "Render.h"

#include <QDir>
#include <QFile>
#include <QImage>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QPainter>
#include <QTemporaryDir>
#include <set>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <algorithm>
#include <cmath>
#include <numeric>

TEST_CASE("SMILES gives 2D coordinates") {
    auto doc = chem::fromSmiles("c1ccccc1O");
    REQUIRE(doc);
    REQUIRE(doc->atoms.size() == 7);
    REQUIRE(doc->bonds.size() == 7);
    int doubles = 0;
    for (auto& b : doc->bonds) doubles += b.order == 2;
    CHECK(doubles == 3);  // kekulized
    auto& b = doc->bonds[0];
    auto d = doc->atoms[b.a].pos - doc->atoms[b.b].pos;
    CHECK(std::abs(std::hypot(d.x(), d.y()) - kBondLength) < 0.5);
}

TEST_CASE("bad SMILES is rejected") {
    CHECK_FALSE(chem::fromSmiles("C1CC"));
}

TEST_CASE("deleting a bond removes only endpoints left isolated") {
    Document chain;
    chain.atoms = {{{0, 0}}, {{kBondLength, 0}}, {{2 * kBondLength, 0}}};
    chain.bonds = {{0, 1}, {1, 2}};
    chain.removeBond(1);
    REQUIRE(chain.atoms.size() == 2);
    REQUIRE(chain.bonds.size() == 1);
    CHECK(chain.bonds[0].a == 0);
    CHECK(chain.bonds[0].b == 1);

    Document bridge;
    for (int i = 0; i < 4; ++i) bridge.addAtom({i * kBondLength, 0});
    bridge.bonds = {{0, 1}, {1, 2}, {2, 3}};
    bridge.removeBond(1);
    CHECK(bridge.atoms.size() == 4);
    CHECK(bridge.bonds.size() == 2);

    bridge.removeBond(0);
    CHECK(bridge.atoms.size() == 2);
    CHECK(bridge.bonds.size() == 1);
}

TEST_CASE(".penz round-trips") {
    auto doc = chem::fromSmiles("C[C@H](N)C(=O)[O-]");
    REQUIRE(doc);
    doc->bonds[0].stereo = BondStereo::Wedge;
    doc->arrows.push_back({{0, 0}, {40, 0}, ArrowKind::Equilibrium});
    doc->arrows.push_back({{0, 10}, {20, 10}, ArrowKind::Fishhook, -6});
    doc->texts.push_back({{5, -8}, "Pd(PPh3)4\n80 °C"});
    doc->style = "JDP";
    auto back = Document::fromJson(doc->toJson());
    REQUIRE(back);
    CHECK(*back == *doc);
}

TEST_CASE(".penz rejects junk and dangling bonds") {
    CHECK_FALSE(Document::fromJson("not json"));
    CHECK_FALSE(Document::fromJson(R"({"format":"penzene","version":1,
        "atoms":[{"x":0,"y":0,"z":6}],"bonds":[{"a":0,"b":5}]})"));
}

TEST_CASE("MOL round-trips atoms, bonds, coordinates and wedges") {
    auto doc = chem::fromSmiles("C[C@H](N)C(=O)O");
    REQUIRE(doc);
    int wedged = 0;
    for (auto& b : doc->bonds) wedged += b.stereo != BondStereo::None;
    CHECK(wedged == 1);
    auto back = chem::fromMolBlock(chem::toMolBlock(*doc));
    REQUIRE(back);
    REQUIRE(back->atoms.size() == doc->atoms.size());
    REQUIRE(back->bonds.size() == doc->bonds.size());
    for (size_t i = 0; i < doc->atoms.size(); ++i) {
        CHECK(back->atoms[i].z == doc->atoms[i].z);
        auto d = back->atoms[i].pos - doc->atoms[i].pos;
        CHECK(std::hypot(d.x(), d.y()) < 0.01);
    }
    for (size_t i = 0; i < doc->bonds.size(); ++i) CHECK(back->bonds[i] == doc->bonds[i]);
    CHECK(chem::toSmiles(*back) == chem::toSmiles(*doc));
    CHECK(chem::toSmiles(*doc) == "C[C@H](N)C(=O)O");
}

TEST_CASE("implicit hydrogens and valence errors") {
    Document d;
    d.atoms = {{{0, 0}, 8}, {{kBondLength, 0}, 6}};
    d.bonds = {{0, 1}};
    auto info = chem::atomInfo(d);
    CHECK(info[0].hydrogens == 1);  // OH
    CHECK(info[1].hydrogens == 3);  // CH3
    d.bonds[0].order = 3;           // O#C: oxygen over valence
    CHECK(chem::atomInfo(d)[0].valenceError);
}

TEST_CASE("clean keeps atom order and centroid") {
    Document d;
    for (int i = 0; i < 6; ++i) d.atoms.push_back({QPointF(100 + i * 3, 50 + i * i)});
    for (int i = 0; i < 6; ++i) d.bonds.push_back({i, (i + 1) % 6, i % 2 ? 2 : 1});
    auto c = chem::clean2D(d);
    REQUIRE(c.atoms.size() == 6);
    auto e = c.atoms[0].pos - c.atoms[1].pos;
    CHECK(std::abs(std::hypot(e.x(), e.y()) - kBondLength) < 0.5);
    QPointF c0, c1;
    for (int i = 0; i < 6; ++i) c0 += d.atoms[i].pos, c1 += c.atoms[i].pos;
    CHECK(std::hypot((c0 - c1).x(), (c0 - c1).y()) < 0.01);
}

TEST_CASE("element symbols") {
    CHECK(chem::symbol(17) == "Cl");
    CHECK(chem::atomicNumber("Br") == 35);
    CHECK(chem::atomicNumber("Xx") == 0);
}

TEST_CASE(".penz rejects bad arrows; v0.1 files still load") {
    CHECK_FALSE(Document::fromJson(R"({"format":"penzene","version":1,"arrows":[{"kind":"wiggly"}]})"));
    auto old = Document::fromJson(R"({"format":"penzene","version":1,"atoms":[{"x":0,"y":0,"z":6}],"bonds":[]})");
    REQUIRE(old);
    CHECK(old->arrows.empty());
}

TEST_CASE("clean lays out each fragment in place and keeps arrows and text") {
    auto left = chem::fromSmiles("CCO"), right = chem::fromSmiles("CC=O");
    REQUIRE(left);
    REQUIRE(right);
    Document scheme;
    scheme.append(*left, {-100, 0});
    scheme.append(*right, {100, 0});
    scheme.arrows.push_back({{-30, 0}, {30, 0}});
    scheme.texts.push_back({{-20, -5}, "PCC"});
    for (auto& a : scheme.atoms) a.pos += QPointF(0, a.pos.x() * 0.1);  // skew it
    Document clean = chem::clean2D(scheme);
    CHECK(clean.arrows == scheme.arrows);
    CHECK(clean.texts == scheme.texts);
    CHECK(clean.bonds.size() == scheme.bonds.size());
    for (int i : {0, 1, 2}) CHECK(clean.atoms[i].pos.x() < -50);  // reactant stays left
    for (int i : {3, 4, 5}) CHECK(clean.atoms[i].pos.x() > 50);

    // Only the molecule with a selected atom moves; the other is untouched.
    Document partial = chem::clean2D(scheme, {4});
    for (int i : {0, 1, 2}) CHECK(partial.atoms[i] == scheme.atoms[i]);
    CHECK(partial.atoms[3] == clean.atoms[3]);  // the product is laid out as in a full clean
    CHECK(partial.atoms[5] == clean.atoms[5]);
    CHECK(partial.bonds.size() == scheme.bonds.size());
}

TEST_CASE("formula, weights and InChI") {
    auto aspirin = chem::fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
    REQUIRE(aspirin);
    auto p = chem::properties(*aspirin);
    REQUIRE(p);
    CHECK(p->formula == "C9H8O4");
    CHECK(std::abs(p->mw - 180.159) < 0.01);
    CHECK(std::abs(p->exactMass - 180.0423) < 0.001);
    CHECK(chem::toInchiKey(*aspirin) == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N");
    CHECK(chem::toInchi(*aspirin).rfind("InChI=1S/C9H8O4/", 0) == 0);
    CHECK_FALSE(chem::properties(Document{}));
}

TEST_CASE("CDXML import: molecules, arrows and text in place") {
    auto doc = chem::readFile(QString(PENZENE_TEST_DATA) + "/scheme.cdxml");
    REQUIRE(doc);
    CHECK(chem::toSmiles(*doc) == "CCO");
    REQUIRE(doc->atoms.size() == 3);
    CHECK(std::abs(doc->atoms[0].pos.x() - 100) < 0.5);  // 14.4 pt bonds: points map 1:1
    CHECK(std::abs(doc->atoms[0].pos.y() - 100) < 0.5);
    REQUIRE(doc->arrows.size() == 3);
    CHECK(doc->arrows[0].kind == ArrowKind::Reaction);
    CHECK(doc->arrows[0].to == QPointF(190, 104));
    CHECK(doc->arrows[1].kind == ArrowKind::Equilibrium);
    // Quarter circle of radius 25 about (165,175): its top is 7.3 pt above the chord.
    CHECK(std::abs(doc->arrows[2].bend - 7.32) < 0.1);
    REQUIRE(doc->texts.size() == 1);
    CHECK(doc->texts[0].text == "PCC");
}

TEST_CASE("hotkeys without a canvas: ChemDraw's dipeptide example") {
    Document doc;
    int n = doc.addAtom({0, 0}, 7);
    edit::link(doc, n, doc.addAtom({kBondLength, 0}));
    edit::Hotspot h{1, -1};
    for (QChar k : QString("42n152o")) {
        h = edit::hotkey(doc, h, k);
        REQUIRE(h.valid());
    }
    std::string smi = chem::toSmiles(doc);
    std::erase(smi, '@');
    CHECK(chem::toSmiles(*chem::fromSmiles(smi)) == chem::toSmiles(*chem::fromSmiles("CC(N)C(=O)NC(C)C(=O)O")));
    CHECK_FALSE(edit::hotkey(doc, {0, -1}, "~").valid());  // not a hotkey
}

TEST_CASE("descriptor table: a row per record, invalid ones kept with the reason (#152)") {
    const std::vector<chem::Record> records{
        {"aspirin", chem::fromSmiles("CC(=O)Oc1ccccc1C(=O)O")}, {"broken", std::nullopt}, {"a, \"quoted\" name", chem::fromSmiles("c1ccccc1")}};
    const QStringList lines = QString::fromStdString(chem::descriptorsCsv(records)).split('\n', Qt::SkipEmptyParts);
    REQUIRE(lines.size() == 4);
    CHECK(lines[0] == chem::descriptorColumns().join(','));
    CHECK(lines[1].startsWith("1,aspirin,CC(=O)Oc1ccccc1C(=O)O,BSYNRYMUTXBXSQ-UHFFFAOYSA-N,C9H8O4,180.16,180.0423,1.31,63.60,1,3,2,13,1,0,true,"));
    CHECK(lines[2] == "2,broken" + QString(",").repeated(chem::descriptorColumns().size() - 3) + ",unreadable");
    CHECK(lines[3].startsWith("3,\"a, \"\"quoted\"\" name\",c1ccccc1,"));
    // Chosen columns, in the order given.
    CHECK(QString::fromStdString(chem::descriptorsCsv(records, {"name", "tpsa"})).startsWith("name,tpsa\naspirin,63.60\nbroken,\n"));
    // The molecules of a drawing, each its own row.
    Document two = *chem::fromSmiles("CCO.c1ccccc1");
    CHECK(chem::molecules(two).size() == 2);
}

TEST_CASE("SMILES export never writes text that doesn't parse (#266)") {
    // From the fuzzer: an H label with four bonds, a charge and a stereocentre.
    Document doc;
    edit::Hotspot h{doc.addAtom({0, 0}), -1};
    for (QChar k : QString("i0JKl4dy:+d")) h = edit::hotkey(doc, h, k);
    const std::string smi = chem::toSmiles(doc);
    INFO(smi);
    CHECK((smi.empty() || chem::fromSmiles(smi)));
}

TEST_CASE("a ChemDraw file's label size relative to its bonds is kept (#203)") {
    auto cdxml = [](const QString& attrs) {
        return QString(R"(<?xml version="1.0"?><CDXML %1><page><fragment>
            <n id="1" p="100 100"/><n id="2" p="126 115" NodeType="Element" Element="8"><t><s>OH</s></t></n>
            <b B="1" E="2"/></fragment></page></CDXML>)").arg(attrs).toUtf8();
    };
    auto small = chem::fromChemDraw(cdxml(R"(BondLength="30" LabelSize="6")"));
    REQUIRE(small);
    CHECK(small->labelRatio == Catch::Approx(0.2));
    CHECK(documentStyle(*small).fontSize == Catch::Approx(0.2 * kBondLength));  // 6 pt labels on 30 pt bonds
    CHECK(documentStyle(*small).labelRadius < drawingStyle("").labelRadius);
    CHECK(chem::fromChemDraw(cdxml(R"(BondLength="30")"))->labelRatio == 0);  // no LabelSize: the style's own
    // Kept through .penz and CDXML.
    CHECK(Document::fromJson(small->toJson())->labelRatio == Catch::Approx(0.2));
    CHECK(chem::fromChemDraw(chem::toCdxml(*small))->labelRatio == Catch::Approx(0.2));
    CHECK_FALSE(Document().toJson().contains("labelRatio"));  // absent unless set
}

TEST_CASE("j onto a ring already there makes that ring Cp⁻ (#254)") {
    Document doc;
    doc.addAtom({0, 0});
    // A plain five-ring exactly where j puts its Cp: j merges into it, adding no ring atoms.
    edit::ringAt(doc, doc.atoms[0].pos + doc.awayDirection(0) * (1.6 * kBondLength), 5, false);
    REQUIRE(edit::hotkey(doc, {0, -1}, "j").valid());
    CHECK(doc.atoms.size() == 7);  // the atom, the ring and its centroid
    CHECK(std::count_if(doc.atoms.begin(), doc.atoms.end(), [](const Atom& a) { return a.charge == -1; }) == 1);
    CHECK(std::count_if(doc.bonds.begin(), doc.bonds.end(), [](const Bond& b) { return b.order == 2; }) == 2);
}

TEST_CASE("clean keeps bond display styles and double-bond positions (#84)") {
    auto d = chem::fromSmiles("CC=CC(C)C");
    REQUIRE(d);
    int dbl = -1, single = -1;
    for (int i = 0; i < int(d->bonds.size()); ++i)
        (d->bonds[i].order == 2 ? dbl : single) = i;
    d->bonds[dbl].position = BondPosition::Right;
    d->bonds[single].stereo = BondStereo::Bold;
    Document clean = chem::clean2D(*d);
    const Bond& b1 = clean.bonds[clean.bondBetween(d->bonds[dbl].a, d->bonds[dbl].b)];
    const Bond& b2 = clean.bonds[clean.bondBetween(d->bonds[single].a, d->bonds[single].b)];
    CHECK(b1.position == BondPosition::Right);
    CHECK(b2.stereo == BondStereo::Bold);
}

TEST_CASE("clean lays abbreviations out as single nodes, so bonds stay even (#85)") {
    // Found by fuzzing: "79P7" puts a Ph label on a crowded atom. Expanding the
    // ring before layout squeezed some bonds to about half length.
    auto d = std::make_optional<Document>();
    d->addAtom({0, 0});
    edit::Hotspot h{0, -1};
    for (QChar k : QString("79P7")) h = edit::hotkey(*d, h, k);
    Document clean = chem::clean2D(*d);
    REQUIRE(clean.atoms.size() == d->atoms.size());
    for (const auto& b : clean.bonds) {
        QPointF v = clean.atoms[b.a].pos - clean.atoms[b.b].pos;
        CHECK(std::abs(std::hypot(v.x(), v.y()) - kBondLength) < 0.2 * kBondLength);
    }
}

TEST_CASE("CDXML label nodes: reagent labels become text, R groups stay labelled (#86)") {
    auto doc = chem::readFile(QString(PENZENE_TEST_DATA) + "/labels.cdxml");
    REQUIRE(doc);
    REQUIRE(doc->atoms.size() == 2);  // the carbon and R; no stray atoms from the reagent label
    CHECK(doc->atoms[1].label == "R");
    CHECK(doc->bonds.size() == 1);
    REQUIRE(doc->texts.size() == 1);
    CHECK(doc->texts[0].text == "LiBr, acetone");
    CHECK(std::abs(doc->texts[0].scale - 0.7) < 0.01);  // 7 pt against the 10 pt default
    CHECK(std::abs(doc->texts[0].pos.x() - 150) < 0.5);
    auto back = Document::fromJson(doc->toJson());  // text scale survives .penz
    REQUIRE(back);
    CHECK(*back == *doc);
}

TEST_CASE("explicit hydrogens and carbon/H display options (#97)") {
    auto eth = chem::fromSmiles("CCO");
    REQUIRE(eth);
    Document withH = chem::addHydrogens(*eth);
    CHECK(withH.atoms.size() == 9);  // C2H6O: 3 heavy + 6 H
    CHECK(chem::properties(withH)->formula == "C2H6O");  // chemistry unchanged
    for (const auto& b : withH.bonds) {  // placed at a bond's length, not piled up
        QPointF v = withH.atoms[b.a].pos - withH.atoms[b.b].pos;
        CHECK(std::hypot(v.x(), v.y()) > 0.4 * kBondLength);
    }
    Document without = chem::removeHydrogens(withH);
    CHECK(without.atoms.size() == 3);

    Document wedgedH = withH;  // a wedged H carries stereo, so it stays
    for (auto& b : wedgedH.bonds)
        if (wedgedH.atoms[b.b].z == 1) { b.stereo = BondStereo::Wedge; break; }
    CHECK(chem::removeHydrogens(wedgedH).atoms.size() == 4);

    eth->carbonLabels = Document::CarbonLabels::Terminal;
    eth->hideImplicitH = true;
    auto back = Document::fromJson(eth->toJson());
    REQUIRE(back);
    CHECK(back->carbonLabels == Document::CarbonLabels::Terminal);
    CHECK(back->hideImplicitH);
}

TEST_CASE("CIP stereo labels, and E/Z read from the drawing (#94)") {
    auto ala = chem::fromSmiles("C[C@H](N)C(=O)O");  // L-alanine
    REQUIRE(ala);
    auto labels = chem::stereoLabels(*ala);
    REQUIRE(labels.size() == 1);
    CHECK(labels[0].atom == 1);
    CHECK(labels[0].text == "S");

    auto ene = chem::fromSmiles("C/C=C/C");
    REQUIRE(ene);
    labels = chem::stereoLabels(*ene);
    REQUIRE(labels.size() == 1);
    CHECK(labels[0].bond >= 0);
    CHECK(labels[0].text == "E");

    // A plain zig-zag drawing of 2-butene is trans: SMILES now says so.
    Document drawn;
    drawn.atoms = {{{0, 0}}, {{12.47, -7.2}}, {{24.94, 0}}, {{37.41, -7.2}}};
    drawn.bonds = {{0, 1}, {1, 2, 2}, {2, 3}};
    CHECK(chem::toSmiles(drawn) == "C/C=C/C");

    ala->showStereo = true;
    auto back = Document::fromJson(ala->toJson());
    REQUIRE(back);
    CHECK(back->showStereo);
}

TEST_CASE("check structure finds valence, stereo, label and overlap problems (#95)") {
    auto has = [](const std::vector<chem::Problem>& ps, const QString& text) {
        return std::any_of(ps.begin(), ps.end(), [&](const auto& p) { return p.message.contains(text); });
    };
    auto butanol = chem::fromSmiles("CCC(C)O");  // a stereocentre drawn without a wedge
    REQUIRE(butanol);
    auto ps = chem::checkStructure(*butanol);
    CHECK(has(ps, "no wedge"));

    Document d;  // ethane with a wedge (not a stereocentre), a bad label, overlapping atoms
    d.atoms = {{{0, 0}}, {{kBondLength, 0}}, {{40, 0}}, {{40.5, 0}}};
    d.bonds = {{0, 1, 1, BondStereo::Wedge}};
    d.atoms[2].label = "Xyz";
    ps = chem::checkStructure(d);
    CHECK(has(ps, "not a stereocentre"));
    CHECK(has(ps, "Unknown label"));
    CHECK(has(ps, "Overlapping"));

    auto pentavalent = chem::fromSmiles("C");
    Document c5 = *pentavalent;
    for (int k = 0; k < 5; ++k) edit::link(c5, 0, c5.addAtom({10.0 * k, 10}));
    CHECK(has(chem::checkStructure(c5), "Valence error"));

    CHECK(chem::checkStructure(*chem::fromSmiles("CC(=O)Oc1ccccc1C(=O)O")).empty());  // aspirin is fine
}

TEST_CASE("aromatic circles preserve chemistry and survive save/load (#98)") {
    auto benzene = chem::fromSmiles("c1ccccc1");
    auto naphthalene = chem::fromSmiles("c1ccc2ccccc2c1");
    auto pyridine = chem::fromSmiles("c1ccncc1");
    auto cyclohexane = chem::fromSmiles("C1CCCCC1");
    REQUIRE(benzene);
    REQUIRE(naphthalene);
    REQUIRE(pyridine);
    REQUIRE(cyclohexane);
    CHECK(chem::aromaticRings(*benzene).size() == 1);
    CHECK(chem::aromaticRings(*naphthalene).size() == 2);
    CHECK(chem::aromaticRings(*pyridine).size() == 1);
    CHECK(chem::aromaticRings(*cyclohexane).empty());

    const auto smiles = chem::toSmiles(*naphthalene);
    naphthalene->aromaticCircles = true;
    auto back = Document::fromJson(naphthalene->toJson());
    REQUIRE(back);
    CHECK(back->aromaticCircles);
    CHECK(chem::toSmiles(*back) == smiles);

    auto ring = chem::aromaticRings(*benzene).front();
    std::sort(ring.begin(), ring.end());
    benzene->aromaticCircleOverrides.push_back(ring);
    Document pair = *benzene;
    pair.append(*benzene, {100, 0});
    REQUIRE(pair.aromaticCircleOverrides.size() == 2);
    pair.removeAtoms({0, 1, 2, 3, 4, 5});
    REQUIRE(pair.aromaticCircleOverrides.size() == 1);
    CHECK(pair.aromaticCircleOverrides.front() == ring);

    if (auto prefix = qgetenv("PENZENE_CIRCLE_SHOTS"); !prefix.isEmpty()) {
        pyridine->aromaticCircles = true;
        CHECK(exportDocument(*naphthalene, QString::fromUtf8(prefix) + "-naphthalene.png", {150, Qt::white}));
        CHECK(exportDocument(*pyridine, QString::fromUtf8(prefix) + "-pyridine.png", {150, Qt::white}));
    }
}

TEST_CASE("properties panel profile for aspirin (#96)") {
    auto p = chem::profile(*chem::fromSmiles("CC(=O)Oc1ccccc1C(=O)O"));
    REQUIRE(p);
    INFO("logP " << p->logP << " tpsa " << p->tpsa << " hbd " << p->hbd << " hba " << p->hba << " rot " << p->rotatable);
    CHECK(p->basic.formula == "C9H8O4");
    CHECK(std::abs(p->logP - 1.31) < 0.01);   // Crippen
    CHECK(std::abs(p->tpsa - 63.6) < 0.1);
    CHECK(p->hbd == 1);
    CHECK(p->hba == 3);
    CHECK(p->rotatable == 2);
    CHECK(p->heavyAtoms == 13);
    REQUIRE(p->elemental.size() == 3);   // C, H, O in Hill order
    CHECK(p->elemental[0].first == "C");
    CHECK(std::abs(p->elemental[0].second - 60.00) < 0.01);
    CHECK(std::abs(p->elemental[1].second - 4.48) < 0.01);
    CHECK(p->elemental[2].first == "O");
    CHECK(p->lipinskiViolations == 0);
    CHECK(p->veber);
    CHECK_FALSE(chem::profile(Document{}));
}

TEST_CASE("atom-map numbers survive SMILES, .penz and the ' hotkey (#99)") {
    auto doc = chem::fromSmiles("[CH3:1][OH:2]");
    REQUIRE(doc);
    CHECK(doc->atoms[0].map + doc->atoms[1].map == 3);
    CHECK(chem::toSmiles(*doc) == "[CH3:1][OH:2]");
    doc->showAtomNumbers = true;
    auto back = Document::fromJson(doc->toJson());
    REQUIRE(back);
    CHECK(*back == *doc);

    auto ethanol = *chem::fromSmiles("CCO");
    CHECK(edit::hotkey(ethanol, {2, -1}, "'").atom == 2);
    CHECK(edit::hotkey(ethanol, {0, -1}, "'").valid());
    CHECK(ethanol.atoms[2].map == 1);
    CHECK(ethanol.atoms[0].map == 2);
    edit::hotkey(ethanol, {2, -1}, "'");
    CHECK(ethanol.atoms[2].map == 0);
    CHECK(chem::toSmiles(ethanol).find("[CH3:2]") != std::string::npos);
}

TEST_CASE("multi-record SDF, .smi and .inchi open as a grid; MOL V3000; InChI (#102)") {
    const std::string aspirin = "CC(=O)Oc1ccccc1C(=O)O";
    const auto inchi = chem::toInchi(*chem::fromSmiles(aspirin));
    auto fromInchi = chem::fromInchi(inchi);
    REQUIRE(fromInchi);
    CHECK(chem::toSmiles(*fromInchi) == chem::toSmiles(*chem::fromSmiles(aspirin)));
    CHECK_FALSE(chem::fromInchi("InChI=nonsense"));

    const auto v3000 = chem::toMolBlock(*fromInchi, true);
    CHECK(v3000.find("V3000") != std::string::npos);
    CHECK(chem::toSmiles(*chem::fromMolBlock(v3000)) == chem::toSmiles(*fromInchi));

    QTemporaryDir dir;
    auto write = [&](const QString& name, const std::string& text) {
        QFile f(dir.filePath(name));
        REQUIRE(f.open(QIODevice::WriteOnly));
        f.write(text.c_str());
        return dir.filePath(name);
    };
    const std::string mol = chem::toMolBlock(*fromInchi), ethanol = chem::toMolBlock(*chem::fromSmiles("CCO"));
    const QString sdf = write("two.sdf", "aspirin" + mol.substr(mol.find('\n')) + "$$$$\n" + ethanol + "$$$$\n");
    const QString smi = write("two.smi", aspirin + " aspirin\nCCO ethanol\n");
    const QString inchis = write("two.inchi", inchi + "\n" + chem::toInchi(*chem::fromSmiles("CCO")) + "\n");
    for (const QString& path : {sdf, smi, inchis}) {
        INFO(path.toStdString());
        auto doc = chem::readFile(path);
        REQUIRE(doc);
        CHECK(doc->atoms.size() == 16);  // 13 + 3, side by side
        CHECK(chem::toSmiles(*doc).find('.') != std::string::npos);
        CHECK(chem::readRecords(path).size() == 2);
    }
    CHECK(chem::readRecords(sdf)[0].name == "aspirin");
    CHECK(chem::readRecords(smi)[1].name == "ethanol");
}

TEST_CASE("reactions: reaction SMILES and RXN, both ways (#101)") {
    // Aspirin synthesis: salicylic acid + acetic anhydride, with pyridine over the arrow.
    const std::string rsmi = "OC(=O)c1ccccc1O.CC(=O)OC(C)=O>c1ccncc1>CC(=O)Oc1ccccc1C(=O)O.CC(=O)O";
    auto doc = chem::fromReactionSmiles(rsmi);
    REQUIRE(doc);
    REQUIRE(doc->arrows.size() == 1);
    CHECK(doc->texts.size() == 2);  // the "+" signs
    auto r = chem::reactionOf(*doc);
    REQUIRE(r);
    CHECK(r->reactants.size() == 2);
    CHECK(r->agents.size() == 1);
    CHECK(r->products.size() == 2);
    auto canon = [](const std::string& s) { return chem::toSmiles(*chem::fromSmiles(s)); };
    CHECK(chem::toReactionSmiles(*r) == canon("OC(=O)c1ccccc1O") + "." + canon("CC(=O)OC(C)=O") + ">" +
                                           canon("c1ccncc1") + ">" + canon("CC(=O)Oc1ccccc1C(=O)O") + "." +
                                           canon("CC(=O)O"));

    const std::string rxn = chem::toRxn(*r);
    CHECK(rxn.starts_with("$RXN"));
    auto back = chem::fromRxn(rxn);
    REQUIRE(back);
    auto rb = chem::reactionOf(*back);
    REQUIRE(rb);
    CHECK(rb->reactants.size() == 2);
    CHECK(rb->products.size() == 2);  // RXN V2000 carries no agents
    CHECK(chem::toSmiles(rb->products[0]) == canon("CC(=O)Oc1ccccc1C(=O)O"));

    CHECK_FALSE(chem::reactionOf(*chem::fromSmiles("CCO")));
    CHECK_FALSE(chem::fromReactionSmiles("CCO"));
    CHECK_FALSE(chem::fromRxn("not an rxn"));
}

TEST_CASE("CDXML export reads back: molecules, wedges, arrows and text (#29)") {
    Document doc = *chem::fromSmiles("C[C@H](N)C(=O)O");  // L-alanine, wedged
    const std::string smiles = chem::toSmiles(doc);
    const QPointF right(60, 0);
    doc.arrows.push_back({right, right + QPointF(40, 0)});
    doc.arrows.push_back({right + QPointF(0, 30), right + QPointF(40, 30), ArrowKind::Reaction, 10});  // curved
    doc.texts.push_back({{0, 40}, "L-alanine"});
    const QByteArray cdxml = chem::toCdxml(doc);
    CHECK(cdxml.contains("<CDXML"));
    auto back = chem::fromChemDraw(cdxml);
    REQUIRE(back);
    CHECK(chem::toSmiles(*back) == smiles);  // stereo survives
    REQUIRE(back->arrows.size() == 2);
    CHECK(QLineF(back->arrows[0].from, doc.arrows[0].from).length() < 0.1);
    CHECK(std::abs(back->arrows[1].bend - 10) < 0.1);  // the arc comes back on the same side
    REQUIRE(back->texts.size() == 1);
    CHECK(back->texts[0].text == "L-alanine");
    for (size_t i = 0; i < doc.atoms.size(); ++i) CHECK(QLineF(back->atoms[i].pos, doc.atoms[i].pos).length() < 0.1);

    Document r = *chem::fromSmiles("CC");
    edit::applyLabel(r, 1, "R", true);
    auto rb = chem::fromChemDraw(chem::toCdxml(r));
    REQUIRE(rb);
    CHECK(rb->atoms.size() == 2);
    CHECK(rb->atoms[1].label == "R");

    // Binary CDX on every platform (#185), carrying arrows and text too (#181).
    const QByteArray cdx = chem::toCdx(doc);
    REQUIRE(cdx.startsWith("VjCD0100"));
    auto fromCdx = chem::fromChemDraw(cdx);
    REQUIRE(fromCdx);
    CHECK(chem::toSmiles(*fromCdx) == smiles);
    CHECK(fromCdx->arrows.size() == 2);
    REQUIRE(fromCdx->texts.size() == 1);
    INFO(fromCdx->texts[0].text.toStdString());
    CHECK(fromCdx->texts[0].text == "L-alanine");
    CHECK(chem::cdxToCdxml("not a CDX file").isEmpty());
    CHECK(chem::cdxmlToCdx("<not closed").isEmpty());
}

TEST_CASE("a real ChemDraw file survives CDXML → CDX → CDXML (#185)") {
    QFile f(QString(PENZENE_TEST_DATA) + "/scheme.cdxml");
    REQUIRE(f.open(QIODevice::ReadOnly));
    const QByteArray cdxml = f.readAll();
    auto direct = chem::fromChemDraw(cdxml);
    auto viaCdx = chem::fromChemDraw(chem::cdxmlToCdx(cdxml));
    REQUIRE(direct);
    REQUIRE(viaCdx);
    CHECK(chem::toSmiles(*viaCdx) == chem::toSmiles(*direct));
    CHECK(viaCdx->arrows.size() == direct->arrows.size());
    CHECK(viaCdx->texts.size() == direct->texts.size());
}

TEST_CASE("SMILES of a ring with a charged boron or phosphorus reads back (fuzz)") {
    for (int z : {5, 15})
        for (int charge : {-2, 2}) {
            Document ring = *chem::fromSmiles("C1=CC=CC=C1");
            ring.atoms[0].z = z;
            ring.atoms[0].charge = charge;
            const std::string smiles = chem::toSmiles(ring);
            INFO(z << " " << charge << " " << smiles);
            if (!smiles.empty()) CHECK(chem::fromSmiles(smiles));
        }
}

TEST_CASE("interaction and partial bonds are drawn, not chemistry (#169)") {
    // Water dimer: an H-bond from one water's H to the other's O.
    Document dimer = chem::addHydrogens(*chem::fromSmiles("O.O"));
    int h = -1, o2 = -1;  // an H on one water, and the other water's O
    for (int i = 0; i < int(dimer.atoms.size()); ++i)
        if (dimer.atoms[i].z == 1 && h < 0) h = i;
    REQUIRE(h >= 0);
    for (int i = 0; i < int(dimer.atoms.size()); ++i)
        if (dimer.atoms[i].z == 8 && dimer.bondBetween(i, h) < 0) o2 = i;
    INFO(chem::toSmiles(dimer));
    REQUIRE(h >= 0);
    REQUIRE(o2 >= 0);
    const std::string water = chem::toSmiles(dimer);
    dimer.bonds.push_back({h, o2});
    CHECK(chem::toSmiles(dimer) != water);  // as a covalent bond: not water any more
    REQUIRE(edit::hotkey(dimer, {-1, int(dimer.bonds.size()) - 1}, "i").valid());
    CHECK(dimer.bonds.back().stereo == BondStereo::Interaction);
    CHECK(chem::toSmiles(dimer) == water);
    CHECK(chem::properties(dimer)->formula == "H4O2");
    auto cleaned = chem::clean2D(dimer);
    CHECK(cleaned.bondBetween(h, o2) >= 0);  // Clean keeps the interaction
    auto back = Document::fromJson(dimer.toJson());
    REQUIRE(back);
    CHECK(*back == dimer);

    // Diels–Alder transition state: butadiene + ethylene, two forming bonds (partial singles),
    // the diene's and dienophile's π bonds partial doubles.
    Document ts = *chem::fromSmiles("C=CC=C.C=C");
    const int c1 = 0, c4 = 3, c5 = 4, c6 = 5;
    for (int bi = 0; bi < int(ts.bonds.size()); ++bi)
        if (ts.bonds[bi].order == 2) REQUIRE(edit::hotkey(ts, {-1, bi}, "P").valid());
    ts.bonds.push_back({c1, c6});
    ts.bonds.push_back({c4, c5});
    for (int bi : {int(ts.bonds.size()) - 2, int(ts.bonds.size()) - 1}) REQUIRE(edit::hotkey(ts, {-1, bi}, "p").valid());
    CHECK(chem::properties(ts)->formula == "C6H10");  // the reactants' atoms, nothing invented
    CHECK(chem::toSmiles(ts) == chem::toSmiles(*chem::fromSmiles("C=CC=C.C=C")));  // the forming bonds don't count
    CHECK(chem::toCdxml(ts).contains("Order=\"0.5\""));
    CHECK(chem::toCdxml(dimer).contains("Order=\"hydrogen\""));
}
TEST_CASE("Clean keeps partial doubles partial (#169)") {
    Document ts = *chem::fromSmiles("C=CC=C");
    for (int bi = 0; bi < int(ts.bonds.size()); ++bi)
        if (ts.bonds[bi].order == 2) edit::hotkey(ts, {-1, bi}, "P");
    int partial = 0;
    for (const auto& b : chem::clean2D(ts).bonds) partial += b.stereo == BondStereo::Partial && b.order == 2;
    CHECK(partial == 2);
}

TEST_CASE("R-groups and generic atoms export and read back (#34)") {
    Document doc = *chem::fromSmiles("CC(=O)CC");  // R1 and R2 on either side of a ketone
    REQUIRE(edit::applyLabel(doc, 0, "R1", true));
    REQUIRE(edit::applyLabel(doc, 4, "R2", true));
    doc.append(*chem::fromSmiles("c1ccccc1C"), {60, 0});  // and a halobenzene: X on the ring
    REQUIRE(edit::applyLabel(doc, 11, "X", true));
    const std::string smiles = chem::toSmiles(doc);
    INFO(smiles);
    CHECK(smiles.find("[1*]") != std::string::npos);
    CHECK(smiles.find("[2*]") != std::string::npos);
    for (bool v3000 : {false, true}) {
        const std::string mol = chem::toMolBlock(doc, v3000);
        INFO(mol);
        CHECK(mol.find("R#") != std::string::npos);
        if (v3000) CHECK(mol.find("RGROUPS") != std::string::npos);
        else CHECK(mol.find("M  RGP") != std::string::npos);
        auto back = chem::fromMolBlock(mol);
        REQUIRE(back);
        QStringList labels;
        for (const auto& a : back->atoms)
            if (!a.label.isEmpty()) labels << a.label;
        labels.sort();
        CHECK(labels.join(",").toStdString() == (v3000 ? "R1,R2" : "R1,R2,X"));  // V3000 has no alias block
    }
    Document ar = *chem::fromSmiles("CC");
    REQUIRE(edit::applyLabel(ar, 1, "Ar", true));
    CHECK(ar.atoms[1].z == 0);  // aryl, not argon
    CHECK(ar.atoms[1].label == "Ar");
    REQUIRE(edit::applyLabel(ar, 1, "Ar"));
    CHECK(ar.atoms[1].z == 18);  // the strict (Python) path: the element
}

TEST_CASE("radicals are chemistry; lone pairs and δ are drawn (#106)") {
    Document methyl = *chem::fromSmiles("C");
    methyl.atoms[0].radicals = 1;
    CHECK(chem::properties(methyl)->formula == "CH3");
    CHECK(chem::toSmiles(methyl) == "[CH3]");
    CHECK(chem::fromSmiles("[CH3]")->atoms[0].radicals == 1);  // read back from SMILES

    Document water = *chem::fromSmiles("O");
    water.atoms[0].lonePairs = 2;
    water.atoms[0].partial = -1;
    CHECK(chem::toSmiles(water) == "O");  // display only
    auto back = Document::fromJson(water.toJson());
    REQUIRE(back);
    CHECK(*back == water);
}

TEST_CASE("shapes and lines: saved, exported and read from CDXML (#107)") {
    Document doc = *chem::fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
    doc.arrows.push_back({{-40, -40}, {40, 40}, ArrowKind::RoundedBox, 0, {}, true});  // a dashed box around it
    doc.arrows.push_back({{-40, 50}, {40, 50}, ArrowKind::Line});
    doc.arrows.push_back({{50, -20}, {90, 10}, ArrowKind::Ellipse});
    auto back = Document::fromJson(doc.toJson());
    REQUIRE(back);
    CHECK(*back == doc);
    CHECK_FALSE(chem::reactionOf(doc));  // shapes aren't reaction arrows
    CHECK(documentBounds(doc).contains(QRectF(-40, -40, 80, 80)));

    auto cdx = chem::fromChemDraw(chem::toCdxml(doc));
    REQUIRE(cdx);
    REQUIRE(cdx->arrows.size() == 3);
    for (size_t k = 0; k < 3; ++k) {
        INFO(k);
        CHECK(cdx->arrows[k].kind == doc.arrows[k].kind);
        CHECK(cdx->arrows[k].dashed == doc.arrows[k].dashed);
        CHECK(QRectF(cdx->arrows[k].from, cdx->arrows[k].to).normalized() ==
              QRectF(doc.arrows[k].from, doc.arrows[k].to).normalized());
    }
}

TEST_CASE("3D rotation keeps stereo: (R)-alanine turned over (#173)") {
    Document ala = *chem::fromSmiles("C[C@@H](C(=O)O)N");  // L-alanine, (S)
    const std::string smiles = chem::toSmiles(ala);
    auto cip = [](const Document& d) {
        for (const auto& l : chem::stereoLabels(d))
            if (l.atom >= 0) return l.text.toStdString();
        return std::string();
    };
    REQUIRE(cip(ala) == "S");
    std::vector<int> all(ala.atoms.size());
    std::iota(all.begin(), all.end(), 0);
    auto pose = chem::pose3D(ala, all);
    REQUIRE(pose);
    const Document same = chem::project3D(ala, *pose, 0, 0);  // no turn: the drawing
    for (size_t i = 0; i < ala.atoms.size(); ++i) CHECK(QLineF(same.atoms[i].pos, ala.atoms[i].pos).length() < 1e-6);
    const Document over = chem::project3D(ala, *pose, 0, 180);  // turned over left to right
    CHECK(chem::toSmiles(over) == smiles);
    CHECK(cip(over) == "S");
    auto cx = [](const Document& d) {
        double x = 0;
        for (const auto& a : d.atoms) x += a.pos.x() / double(d.atoms.size());
        return x;
    };
    int flipped = 0;  // seen from behind: left and right swap
    for (size_t i = 0; i < ala.atoms.size(); ++i)
        flipped += (ala.atoms[i].pos.x() - cx(ala)) * (over.atoms[i].pos.x() - cx(over)) < 0;
    CHECK(flipped >= int(ala.atoms.size()) - 2);
    const Document tilted = chem::project3D(ala, *pose, 35, 20);
    CHECK(chem::toSmiles(tilted) == smiles);
}

TEST_CASE("CDXML: superseded graphics aren't doubled; lone-pair symbols land on atoms (#193)") {
    const QByteArray cdxml = R"(<?xml version="1.0" encoding="UTF-8" ?>
<CDXML BondLength="14.4"><page id="1">
<fragment id="2"><n id="3" p="0 0"/><n id="4" p="14.4 0" Element="8"/><b id="5" B="3" E="4"/></fragment>
<arrow id="6" Head3D="40 0 0" Tail3D="80 0 0"/>
<graphic id="7" SupersededBy="6" BoundingBox="40 0 80 0" GraphicType="Line"/>
<graphic id="8" BoundingBox="18 -6 21 -6" GraphicType="Symbol" SymbolType="LonePair"/>
</page></CDXML>)";
    auto doc = chem::fromChemDraw(cdxml);
    REQUIRE(doc);
    CHECK(doc->arrows.size() == 1);  // the line once, not twice
    CHECK(doc->arrows[0].kind == ArrowKind::Line);
    int pairs = 0;
    for (const auto& a : doc->atoms) pairs += a.z == 8 ? a.lonePairs : 0;
    CHECK(pairs == 1);
}

TEST_CASE("element names for accessibility, without RDKit's logger (#112)") {
    CHECK(chem::elementName(1) == "Hydrogen");
    CHECK(chem::elementName(6) == "Carbon");
    CHECK(chem::elementName(118) == "Oganesson");
    CHECK(chem::elementName(0).empty());
    CHECK(chem::elementName(119).empty());
}

TEST_CASE("a CDXML file that starts with a byte order mark opens (#244)") {
    const QByteArray cdxml = chem::toCdxml(*chem::fromSmiles("CC(=O)Oc1ccccc1C(=O)O"));
    auto doc = chem::fromChemDraw("\xEF\xBB\xBF" + cdxml);
    REQUIRE(doc);
    CHECK(chem::toSmiles(*doc) == "CC(=O)Oc1ccccc1C(=O)O");
}

TEST_CASE("arrow heads: even on a tight curve, and a half head has no sliver (#221)") {
    // Ink either side of the head's axis, which runs from where the shaft is 4.8 pt (the head's
    // notch) back from the tip. `side` +1 or -1; samples the head, not the shaft behind it.
    auto inkBeside = [](const Arrow& a, int side) {
        Document d;
        d.arrows.push_back(a);
        const double k = 20;
        QImage img(800, 800, QImage::Format_ARGB32);
        img.fill(Qt::white);
        QPainter p(&img);
        p.translate(400, 400);
        p.scale(k, k);
        p.translate(-a.to);
        paintDocument(p, d, {Qt::black, Qt::black, 0.6});
        p.end();
        const QPainterPath path = arrowPath(a);
        const QPointF back = path.pointAtPercent(path.percentAtLength(path.length() - 4.8));
        const QPointF dir = (a.to - back) / QLineF(back, a.to).length(), n(-dir.y(), dir.x());
        int ink = 0;
        for (double t = 0.3; t < 4.8; t += 0.1)
            for (double s = 0.4; s < 2.4; s += 0.1) {
                const QPointF q = (a.to - dir * t + n * (side * s) - a.to) * k + QPointF(400, 400);
                ink += qGray(img.pixel(q.toPoint())) < 128;
            }
        return ink;
    };
    const Arrow curved{{0, 0}, {16, 0}, ArrowKind::Reaction, 10};  // the tool icon's
    const int left = inkBeside(curved, 1), right = inkBeside(curved, -1);
    CHECK(std::min(left, right) > 0.7 * std::max(left, right));
    // A fishhook's single barb is on one side; the other side stays clear of the head.
    const Arrow hook{{0, 0}, {40, 0}, ArrowKind::Fishhook, 12};
    CHECK(std::min(inkBeside(hook, 1), inkBeside(hook, -1)) == 0);
}

TEST_CASE(".penz files from every release still open, and save back the same (#117)") {
    const QStringList files = QDir(QString(PENZENE_TEST_DATA) + "/penz").entryList({"v*.penz"});
    CHECK(files.size() >= 8);  // v0.1.0 to v0.8.0, plus one per later release
    for (const QString& name : files) {
        INFO(name.toStdString());
        QFile f(QString(PENZENE_TEST_DATA) + "/penz/" + name);
        REQUIRE(f.open(QIODevice::ReadOnly));
        auto doc = Document::fromJson(f.readAll());
        REQUIRE(doc);
        CHECK(doc->atoms.size() == 13);  // aspirin
        auto again = Document::fromJson(doc->toJson());
        REQUIRE(again);
        CHECK(again->toJson() == doc->toJson());
    }
}

TEST_CASE(".penz rejects what would abort the app later: unknown elements, duplicate bonds, non-finite atoms") {
    auto doc = [](const char* atoms, const char* bonds = "[]") {
        return Document::fromJson(QString(R"({"format":"penzene","version":1,"atoms":%1,"bonds":%2})").arg(atoms, bonds).toUtf8());
    };
    CHECK(doc(R"([{"x":0,"y":0,"z":118}])"));
    CHECK(doc(R"([{"x":0,"y":0,"z":0,"label":"R"}])"));
    CHECK_FALSE(doc(R"([{"x":0,"y":0,"z":119}])"));
    CHECK_FALSE(doc(R"([{"x":0,"y":0,"z":-3}])"));
    CHECK_FALSE(doc(R"([{"x":1e999,"y":0}])"));
    const char* two = R"([{"x":0,"y":0},{"x":14,"y":0,"z":8}])";
    CHECK(doc(two, R"([{"a":0,"b":1,"order":2}])"));
    CHECK_FALSE(doc(two, R"([{"a":0,"b":1},{"a":1,"b":0,"order":2}])"));
    for (int z : {1, 118}) CHECK_FALSE(chem::symbol(z).empty());  // the whole accepted range is looked up safely
}

TEST_CASE("formula charges come after the counts") {
    CHECK(chem::properties(*chem::fromSmiles("[O-]S(=O)(=O)[O-]"))->formula == "O4S-2");
}

TEST_CASE("text charges and counts reach ChemDraw as super- and subscripts, and come back (#284)") {
    Document doc;
    doc.addAtom({0, 0});
    doc.texts.push_back({{0, 30}, "SO4^2- and NH4+\nH2O"});
    const QByteArray cdxml = chem::toCdxml(doc);
    CHECK_FALSE(cdxml.contains('^'));
    CHECK(cdxml.contains("face=\"32\">4<"));   // the count
    CHECK(cdxml.contains("face=\"64\">2-<"));  // the charge
    CHECK(cdxml.contains("face=\"64\">+<"));
    auto back = chem::fromChemDraw(cdxml);
    REQUIRE(back);
    REQUIRE(back->texts.size() == 1);
    CHECK(back->texts[0].text == "SO4^2- and NH4^+\nH2O");  // the same drawing: + was a charge already
}

TEST_CASE("merging atoms keeps the dropped atom's brackets and ring overrides (#303)") {
    Document doc;
    for (int i = 0; i < 4; ++i) doc.addAtom({i * kBondLength, 0});
    doc.bonds = {{0, 1}, {2, 3}};
    doc.brackets = {{{2, 3}}};
    doc.aromaticCircleOverrides = {{1, 3}};
    edit::mergeAtoms(doc, {{0, 3}});  // atom 3 fuses onto atom 0
    REQUIRE(doc.atoms.size() == 3);
    CHECK(doc.brackets.at(0).atoms == std::vector<int>{0, 2});
    CHECK(doc.aromaticCircleOverrides == std::vector<std::vector<int>>{{0, 1}});
}

TEST_CASE("a .penz file keeps its pages, and older versions still open the first (#219)") {
    Document one, two;
    one.addAtom({0, 0});
    two.addAtom({0, 0}, 8), two.addAtom({kBondLength, 0}, 7);
    two.bonds = {{0, 1}};
    const std::vector<Sheet> sheets{{"Scheme", one}, {"Mechanism", two}};
    const QByteArray json = sheetsToJson(sheets);
    CHECK(sheetsFromJson(json) == sheets);
    CHECK(Document::fromJson(json) == one);  // what a version without pages reads
    const auto single = sheetsFromJson(one.toJson());
    REQUIRE(single.size() == 1);
    CHECK(single[0].name == "Page 1");
    QJsonObject root = QJsonDocument::fromJson(json).object();
    QJsonArray pages = root["pages"].toArray();
    QJsonObject second = pages[0].toObject();
    second["version"] = 9;  // a page this version can't read
    pages[0] = second;
    root["pages"] = pages;
    CHECK(sheetsFromJson(QJsonDocument(root).toJson()).empty());
}

TEST_CASE("a ChemDraw text of size 0 saves a .penz that opens again (#316)") {
    const QByteArray cdxml = R"(<?xml version="1.0"?><CDXML><page>
        <t p="100 100"><s size="0">heat</s></t></page></CDXML>)";
    auto doc = chem::fromChemDraw(cdxml);
    REQUIRE(doc);
    REQUIRE(doc->texts.size() == 1);
    CHECK(doc->texts[0].scale > 0);
    CHECK(Document::fromJson(doc->toJson()));
    // A file saved before the fix opens, at the normal size.
    Document old;
    old.texts = {{{0, 0}, "heat"}};
    QByteArray json = old.toJson();
    json.replace("\"text\": \"heat\"", "\"scale\": 0, \"text\": \"heat\"");
    REQUIRE(json.contains("\"scale\": 0"));
    auto reopened = Document::fromJson(json);
    REQUIRE(reopened);
    CHECK(reopened->texts[0].scale == 1);
}

TEST_CASE("an SDF whose first record has no title keeps that record (#317)") {
    QTemporaryDir dir;
    const QString path = dir.filePath("untitled.sdf");
    QFile f(path);
    REQUIRE(f.open(QIODevice::WriteOnly));
    const std::string ethanol = chem::toMolBlock(*chem::fromSmiles("CCO"));
    REQUIRE(ethanol.front() == '\n');  // an empty name line
    f.write(QByteArray::fromStdString(ethanol + "$$$$\n" + chem::toMolBlock(*chem::fromSmiles("N")) + "$$$$\n"));
    f.close();
    const auto records = chem::readRecords(path);
    REQUIRE(records.size() == 2);
    REQUIRE(records[0].doc);
    CHECK(chem::toSmiles(*records[0].doc) == "CCO");
    REQUIRE(records[1].doc);
    CHECK(chem::toSmiles(*records[1].doc) == "N");
}

TEST_CASE("a .smi header row is not read as a structure (#329)") {
    QTemporaryDir dir;
    const QString path = dir.filePath("lib.smi");
    QFile f(path);
    REQUIRE(f.open(QIODevice::WriteOnly));
    f.write("SMILES Name\nCCO ethanol\nN ammonia\n");
    f.close();
    const auto records = chem::readRecords(path);
    REQUIRE(records.size() == 2);
    CHECK(records[0].name == "ethanol");
    CHECK(records[0].doc);
}

TEST_CASE("ChemDraw export keeps radicals (#319)") {
    for (const char* smiles : {"[CH2]C", "[CH]C"}) {  // doublet and triplet
        auto doc = chem::fromSmiles(smiles);
        REQUIRE(doc);
        auto back = chem::fromChemDraw(chem::toCdxml(*doc));
        REQUIRE(back);
        CHECK(chem::toSmiles(*back) == chem::toSmiles(*doc));
    }
}

TEST_CASE("SDF export writes one record per molecule (#330)") {
    auto doc = chem::fromSmiles("CCO.N");
    REQUIRE(doc);
    QTemporaryDir dir;
    const QString path = dir.filePath("two.sdf");
    QFile f(path);
    REQUIRE(f.open(QIODevice::WriteOnly));
    f.write(QByteArray::fromStdString(chem::toSdf(*doc)));
    f.close();
    const auto records = chem::readRecords(path);
    REQUIRE(records.size() == 2);
    std::set<std::string> smiles;
    for (const auto& r : records) {
        REQUIRE(r.doc);
        smiles.insert(chem::toSmiles(*r.doc));
    }
    CHECK(smiles == std::set<std::string>{"CCO", "N"});
    CHECK(chem::toSdf(Document{}).empty());
}

TEST_CASE("MOL export writes Kekulé bonds, not query bond type 4 (#321)") {
    for (bool v3000 : {false, true}) {
        const std::string mol = chem::toMolBlock(*chem::fromSmiles("c1ccccc1"), v3000);
        INFO(mol);
        CHECK(mol.find(v3000 ? "M  V30 1 4 " : "  1  2  4  0") == std::string::npos);
        auto back = chem::fromMolBlock(mol);
        REQUIRE(back);
        CHECK(chem::toSmiles(*back) == "c1ccccc1");
}
}

TEST_CASE("Clean keeps bond colours and the drawn Kekulé structure (#322)") {
}

TEST_CASE("MOL and CDXML imports keep the drawn Kekulé form (#323)") {
    auto orders = [](const Document& d) {
        std::vector<std::tuple<int, int, int>> out;
        for (const auto& b : d.bonds) out.push_back({std::min(b.a, b.b), std::max(b.a, b.b), b.order});
        std::sort(out.begin(), out.end());
        return out;
    };
    for (bool flip : {false, true}) {  // either Kekulé form stays as drawn
        Document d = *chem::fromSmiles("c1ccccc1");
        for (auto& b : d.bonds) b.order = flip ? 3 - b.order : b.order, b.color = QColor("#ff0000");
        const Document c = chem::clean2D(d);
        CHECK(orders(c) == orders(d));
        for (const auto& b : c.bonds) CHECK(b.color == QColor("#ff0000"));
}

    for (int form : {0, 1}) {  // both Kekulé forms of benzene
        Document d;
        for (int k = 0; k < 6; ++k)
            d.addAtom({kBondLength * std::cos(M_PI / 3 * k), kBondLength * std::sin(M_PI / 3 * k)});
        for (int k = 0; k < 6; ++k) d.bonds.push_back({k, (k + 1) % 6, (k + form) % 2 ? 1 : 2});
        auto cdxml = chem::fromChemDraw(chem::toCdxml(d));
        REQUIRE(cdxml);
        CHECK(orders(*cdxml) == orders(d));
        std::string mol = "\n  hand\n\n  6  6  0  0  0  0  0  0  0  0999 V2000\n";
        for (int k = 0; k < 6; ++k)
            mol += QString("%1%2    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n")
                       .arg(1.5 * std::cos(M_PI / 3 * k), 10, 'f', 4)
                       .arg(1.5 * std::sin(M_PI / 3 * k), 10, 'f', 4)
                       .toStdString();
        for (const auto& b : d.bonds)
            mol += QString("%1%2%3  0\n").arg(b.a + 1, 3).arg(b.b + 1, 3).arg(b.order, 3).toStdString();
        auto fromMol = chem::fromMolBlock(mol + "M  END\n");
        REQUIRE(fromMol);
        CHECK(orders(*fromMol) == orders(d));
}
}

TEST_CASE("a typed charged label sets the element and charge (#324)") {
    for (auto [label, z, charge] : {std::tuple{"NH3+", 7, 1}, {"O-", 8, -1}, {"Na+", 11, 1}, {"Fe3+", 26, 3},
                                    {"O2-", 8, -2}, {"NH4+", 7, 1}, {"O−", 8, -1}, {"S+2", 16, 2}, {"15NH4+", 7, 1}}) {  // the last with an isotope too
        Document d;
        d.addAtom({0, 0});
        REQUIRE(edit::applyLabel(d, 0, QString::fromUtf8(label), true));
        INFO(label << " -> z " << d.atoms[0].z << " label '" << d.atoms[0].label.toStdString() << "'");
        CHECK(d.atoms[0].z == z);
        CHECK(d.atoms[0].charge == charge);
        CHECK(d.atoms[0].label.isEmpty());
    }
}

TEST_CASE("CDXML import's Kekulé restore isn't fooled by overlapping atoms (#365)") {
    Document d;
    for (int k = 0; k < 6; ++k) d.addAtom({kBondLength * std::cos(M_PI / 3 * k), kBondLength * std::sin(M_PI / 3 * k)});
    for (int k = 0; k < 6; ++k) d.bonds.push_back({k, (k + 1) % 6, k % 2 ? 1 : 2});  // 1-2 is single
    const int a = d.addAtom(d.atoms[1].pos), b = d.addAtom(d.atoms[2].pos);  // an ethylene on top of it
    d.bonds.push_back({a, b, 2});
    auto back = chem::fromChemDraw(chem::toCdxml(d));
    REQUIRE(back);
    for (int i = 0; i < int(back->atoms.size()); ++i) {
        int doubles = 0;
        for (const auto& bond : back->bonds) doubles += (bond.a == i || bond.b == i) && bond.order == 2;
        CHECK(doubles <= 1);
    }
}

TEST_CASE("orbitals: kind, look and colour survive .penz and CDXML (#204)") {
    const QByteArray cdxml = R"(<?xml version="1.0" encoding="UTF-8" ?>
<CDXML BondLength="14.4"><page id="1">
<fragment id="2"><n id="3" p="0 0" Z="10"/><n id="4" p="14.4 0" Z="11"/><b id="5" B="3" E="4"/></fragment>
<graphic id="6" Z="5" BoundingBox="6 0 0 0" GraphicType="Orbital" OvalType="Circle Shaded" OrbitalType="sShaded" Center3D="0 0 0" MajorAxisEnd3D="6 0 0" MinorAxisEnd3D="0 6 0"/>
<graphic id="7" BoundingBox="14.4 -12 14.4 0" GraphicType="Orbital" OrbitalType="pFilled" Center3D="14.4 0 0" MajorAxisEnd3D="14.4 -12 0" MinorAxisEnd3D="20.4 0 0"/>
<graphic id="8" BoundingBox="0 12 0 0" GraphicType="Orbital" OrbitalType="hybridPlus" Center3D="0 0 0" MajorAxisEnd3D="0 12 0" MinorAxisEnd3D="6 0 0"/>
<graphic id="9" BoundingBox="0 12 0 0" GraphicType="Orbital" OrbitalType="dz2Plus" Center3D="0 0 0" MajorAxisEnd3D="0 12 0" MinorAxisEnd3D="6 0 0"/>
</page></CDXML>)";
    auto doc = chem::fromChemDraw(cdxml);
    REQUIRE(doc);
    REQUIRE(doc->arrows.size() == 3);  // no d orbitals
    // ChemDraw's Shaded is our gradient, its Filled our solid shading.
    CHECK(doc->arrows[0].kind == ArrowKind::SOrbital);
    CHECK(doc->arrows[0].look == OrbitalLook::Gradient);
    CHECK(doc->arrows[0].behind);  // below the atoms' Z
    CHECK(!doc->arrows[1].behind);
    CHECK(doc->arrows[1].kind == ArrowKind::POrbital);
    CHECK(doc->arrows[1].look == OrbitalLook::Shaded);
    CHECK(doc->arrows[1].from == doc->atoms[1].pos);
    CHECK(doc->arrows[1].to.y() < doc->arrows[1].from.y());
    CHECK(doc->arrows[2].kind == ArrowKind::HybridOrbital);
    CHECK(doc->arrows[2].look == OrbitalLook::Outline);

    doc->arrows[1].color = QColor(200, 30, 30);
    auto again = Document::fromJson(doc->toJson());
    REQUIRE(again);
    CHECK(again->arrows == doc->arrows);

    auto back = chem::fromChemDraw(chem::toCdxml(*doc));
    REQUIRE(back);
    REQUIRE(back->arrows.size() == 3);
    for (int i = 0; i < 3; ++i) {
        CHECK(back->arrows[i].kind == doc->arrows[i].kind);
        CHECK(back->arrows[i].look == doc->arrows[i].look);
        CHECK(back->arrows[i].behind == doc->arrows[i].behind);
        CHECK(len(back->arrows[i].to - doc->arrows[i].to) < 0.05);
    }
}

TEST_CASE("arrange: arrows and orbitals stack around the molecule's layer (#204)") {
    Document doc;
    for (int i = 0; i < 3; ++i) doc.arrows.push_back({{double(i), 0}, {double(i), 5}, ArrowKind::Line});
    doc.arrows[0].behind = true;  // stack: 0 | molecule | 1 2
    auto xs = [&] {
        std::vector<std::pair<double, bool>> v;
        for (const auto& a : doc.arrows) v.push_back({a.from.x(), a.behind});
        return v;
    };
    CHECK(edit::restack(doc, {2}, edit::Restack::Back) == std::vector<int>{0});
    CHECK(xs() == std::vector<std::pair<double, bool>>{{2, true}, {0, true}, {1, false}});
    CHECK(edit::restack(doc, {1}, edit::Restack::Forward) == std::vector<int>{1});  // past the molecule
    CHECK(xs() == std::vector<std::pair<double, bool>>{{2, true}, {0, false}, {1, false}});
    CHECK(edit::restack(doc, {1}, edit::Restack::Backward) == std::vector<int>{1});  // back under it
    CHECK(xs() == std::vector<std::pair<double, bool>>{{2, true}, {0, true}, {1, false}});
    CHECK(edit::restack(doc, {0}, edit::Restack::Front) == std::vector<int>{2});
    CHECK(xs() == std::vector<std::pair<double, bool>>{{0, true}, {1, false}, {2, false}});
}

TEST_CASE("more abbreviations: acids, alkyls, amines and protecting groups (#331)") {
    // Each alone, as its parent compound (the group plus H).
    for (auto [label, formula] : {std::pair{"COOH", "CH2O2"}, {"Bu", "C4H10"}, {"iBu", "C4H10"}, {"OEt", "C2H6O"},
                                  {"NHBoc", "C5H11NO2"}, {"NMe2", "C2H7N"}, {"SO3H", "H2O3S"},
                                  {"TIPS", "C9H22Si"}, {"MOM", "C2H6O"}, {"THP", "C5H10O"}}) {
        Document d;
        d.addAtom({0, 0});
        REQUIRE(edit::applyLabel(d, 0, label, true));
        INFO(label);
        CHECK(d.atoms[0].label == label);
        auto p = chem::properties(d);
        REQUIRE(p);
        CHECK(p->formula == formula);
    }
}
