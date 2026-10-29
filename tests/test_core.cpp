#include "Chem.h"
#include "Geometry.h"
#include "Edit.h"
#include "Render.h"
#include "Styles.h"
#include "Jcamp.h"

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/MolOps.h>

#include <QDir>
#include <QFile>
#include <QImage>
#include <QFontMetricsF>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QPainter>
#include <QPicture>
#include <QTemporaryDir>
#include <QXmlStreamReader>
#include <set>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <algorithm>
#include <cmath>
#include <numeric>

TEST_CASE("JCAMP decodes numeric and compressed samples with checked difference checkpoints (#553)") {
    const QByteArray header = "##TITLE=synthetic\n##JCAMP-DX=5.01\n##DATA TYPE=NMR SPECTRUM\n"
        "##.OBSERVE NUCLEUS=^1H\n##XUNITS=PPM\n##FIRSTX=4\n##LASTX=0\n##NPOINTS=5\n"
        "##XFACTOR=0.1\n##YFACTOR=2\n##XYDATA=(X++(Y..Y))\n";
    for (const QByteArray& data : {QByteArray("40 1 2 3 6 10\n"), QByteArray("40 A B C F A0\n"),
                                  QByteArray("40AJT\n20CLM\n0A0\n"), QByteArray("40+1+2+3+6+10\n")}) {
        QString error;
        const auto spectra = readJcamp(header + data + "##END=\n", &error);
        INFO(error.toStdString());
        REQUIRE(spectra.size() == 1);
        CHECK(spectra[0].proton);
        CHECK(spectra[0].points == std::vector<QPointF>{{4, 2}, {3, 4}, {2, 6}, {1, 12}, {0, 20}});
    }
    auto gaps = readJcamp(header + "40 A?BC@\n##END=\n");
    REQUIRE(gaps.size() == 1);
    CHECK(std::isnan(gaps[0].points[1].y()));
    const auto scientific = readJcamp(header + "40 1e0 2E+0 3.0e0 6 10\n##END=\n");
    REQUIRE(scientific.size() == 1);
    CHECK(scientific[0].points[1].y() == 4);
    const auto packedE = readJcamp(header + "40E3BCFA0\n##END=\n");
    REQUIRE(packedE.size() == 1);
    CHECK(packedE[0].points[0].y() == 106);
    const auto repeated = readJcamp(header + "40 A U b @\n##END=\n");
    REQUIRE(repeated.size() == 1);
    CHECK(repeated[0].points == std::vector<QPointF>{{4, 2}, {3, 2}, {2, 2}, {1, -4}, {0, 0}});
    for (const QByteArray& bad : {QByteArray("40AJT\n20DLM\n0A0\n"), QByteArray("40AJT\n21CLM\n0A0\n"),
                                 QByteArray("40AJT\n20CLM\n"), QByteArray("40AS999999999\n"), QByteArray("40 A B\n"),
                                 QByteArray("40 A T T @\n"), QByteArray("40 A ? J B @\n")}) {
        QString error;
        CHECK(readJcamp(header + bad + "##END=\n", &error).empty());
        CHECK_FALSE(error.isEmpty());
    }
    CHECK(readJcamp(header + "40 A B C F A0\n").empty());
}

TEST_CASE("JCAMP converts Hz to referenced ppm and reads NTUPLES real pages and linked spectra (#553)") {
    const QByteArray header = "##TITLE=carbon\n##JCAMP-DX=5.01\n##DATA TYPE=NMR SPECTRUM\n"
        "##.OBSERVE NUCLEUS=13C\n##.OBSERVE FREQUENCY=100\n##XUNITS=HZ\n"
        "##NPOINTS=3\n##.SHIFT REFERENCE=INTERNAL, CDCl3, 2, 77\n##XYPOINTS=(XY..XY)\n14000,1;12000,2;10000,1\n##END=\n";
    const auto spectra = readJcamp(header);
    REQUIRE(spectra.size() == 1);
    CHECK_FALSE(spectra[0].proton);
    CHECK(spectra[0].points == std::vector<QPointF>{{97, 1}, {77, 2}, {57, 1}});
    QByteArray offset = header;
    offset.replace("##.SHIFT REFERENCE=INTERNAL, CDCl3, 2, 77", "##FIRSTX=14000\n##$OFFSET=97");
    CHECK(readJcamp(offset)[0].points == spectra[0].points);
    offset.replace("14000,1;12000,2;10000,1", "10000,1;12000,2;14000,1");
    CHECK(readJcamp(offset)[0].points == std::vector<QPointF>{{57, 1}, {77, 2}, {97, 1}});
    const auto peaks = readJcamp("##TITLE=negative shifts\n##DATA TYPE=NMR PEAK TABLE\n##.OBSERVE NUCLEUS=1H\n"
                               "##XUNITS=PPM\n##$OFFSET=1\n##PEAK TABLE=(XY..XY)\n-2,-1;-1,2\n##END=\n");
    REQUIRE(peaks.size() == 1);
    CHECK(peaks[0].peaks);
    CHECK(peaks[0].points == std::vector<QPointF>{{0, -1}, {1, 2}});
    const QByteArray ntuples = "##TITLE=real and imaginary\n##DATA TYPE=NMR SPECTRUM\n##.OBSERVE NUCLEUS=1H\n"
        "##.OBSERVE FREQUENCY=400\n##NUM DIM=1\n##NTUPLES=NMR SPECTRUM\n##SYMBOL=X,R,I\n"
        "##VAR_TYPE=INDEPENDENT,DEPENDENT,DEPENDENT\n##VAR_DIM=3,3,3\n##UNITS=HZ,ARBITRARY UNITS,ARBITRARY UNITS\n"
        "##FACTOR=1,2,3\n##FIRST=800,0,0\n##LAST=0,0,0\n##PAGE=N=1\n"
        "##DATA TABLE=(X++(R..R)),XYDATA\n800 A B C\n##PAGE=N=2\n##DATA TABLE=(X++(I..I)),XYDATA\n800 ZZZ\n##END NTUPLES=\n##END=\n";
    const auto real = readJcamp(ntuples);
    REQUIRE(real.size() == 1);
    CHECK(real[0].points == std::vector<QPointF>{{2, 2}, {1, 4}, {0, 6}});
    QByteArray generic = ntuples;
    generic.replace("SYMBOL=X,R,I", "SYMBOL=F1,R,I");
    generic.replace("(X++", "(F1++");
    CHECK(readJcamp(generic)[0].points == real[0].points);
    const auto linked = readJcamp("##TITLE=two spectra\n##DATA TYPE=LINK\n" + header + ntuples + "##END=\n");
    CHECK(linked.size() == 2);
    QByteArray bad = header;
    bad.replace("##.OBSERVE FREQUENCY=100\n", "");
    CHECK(readJcamp(bad).empty());
    bad = header;
    bad.replace("CDCl3, 2, 77", "CDCl3, 0, 77");
    CHECK(readJcamp(bad).empty());
    bad = ntuples;
    bad.replace("##NUM DIM=1", "##NUM DIM=2");
    CHECK(readJcamp(bad).empty());
    bad = header;
    bad.replace("NMR SPECTRUM", "NMR FID");
    CHECK(readJcamp(bad).empty());
    CHECK(readJcamp("not a spectrum").empty());
}

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

TEST_CASE("themes have a validated versioned JSON round trip (#475)") {
    for (const Theme& t : themes()) {
        INFO(t.name.toStdString());
        const auto back = Theme::fromJson(t.toJson());
        REQUIRE(back);
        CHECK(back->toJson() == t.toJson());
        CHECK(chrome(*back).border.rgba() == chrome(t).border.rgba());
        CHECK(chrome(*back).secondary.rgba() == chrome(t).secondary.rgba());
        CHECK(chrome(*back).accentBg.rgba() == chrome(t).accentBg.rgba());
    }
    const auto valid = QJsonDocument::fromJson(theme("Light").toJson()).object();
    for (auto [key, value] : {std::pair<QString, QJsonValue>{"version", 2}, {"version", "1"}, {"format", "penzene"},
                              {"name", " "}, {"dark", "false"}, {"colors", "#ffffff"}}) {
        auto bad = valid;
        bad[key] = value;
        CHECK_FALSE(Theme::fromJson(QJsonDocument(bad).toJson()));
    }
    for (QJsonValue badColor : {QJsonValue("invalid"), QJsonValue(12), QJsonValue(QJsonValue::Null)}) {
        auto bad = valid;
        auto colors = bad["colors"].toObject();
        colors["ink"] = badColor;
        bad["colors"] = colors;
        CHECK_FALSE(Theme::fromJson(QJsonDocument(bad).toJson()));
    }
    auto missing = valid;
    auto colors = missing["colors"].toObject();
    colors.remove("accent");
    missing["colors"] = colors;
    CHECK_FALSE(Theme::fromJson(QJsonDocument(missing).toJson()));
    CHECK_FALSE(Theme::fromJson("{}"));
    CHECK_FALSE(Theme::fromJson("not JSON"));
    QFile example(QString(PENZENE_TEST_DATA) + "/../../docs/_static/theme-example.json");
    REQUIRE(example.open(QIODevice::ReadOnly));
    const auto shared = Theme::fromJson(example.readAll());
    REQUIRE(shared);
    CHECK(shared->name == "Shared Light");
    CHECK(chrome(*shared).border == chrome(theme("Light")).border);
}

TEST_CASE("substructure search finds unique chemistry matches and maps abbreviations (#485)") {
    const Document doc = *chem::fromSmiles("CC(=O)OC.CC=O");
    const auto matches = chem::substructureMatches(doc, "C=O");
    REQUIRE(matches);
    CHECK(*matches == std::vector<std::vector<int>>{{1, 2}, {6, 7}});
    CHECK(chem::substructureMatches(doc, "[#8]")->size() == 3);
    CHECK(chem::substructureMatches(doc, "[$(C=O)]")->size() == 2);
    CHECK(chem::substructureMatches(doc, "N")->empty());
    CHECK_FALSE(chem::substructureMatches(doc, "C1"));
    CHECK_FALSE(chem::substructureMatches(doc, ""));
    const auto overlap = chem::substructureMatches(*chem::fromSmiles("CCC"), "CC");
    REQUIRE(overlap);
    CHECK(*overlap == std::vector<std::vector<int>>{{0, 1}, {1, 2}});
    Document abbrev = *chem::fromSmiles("CO");
    REQUIRE(edit::applyLabel(abbrev, 0, "tBu", true));
    const auto methyls = chem::substructureMatches(abbrev, "[CH3]");
    REQUIRE(methyls);
    CHECK(*methyls == std::vector<std::vector<int>>{{0}, {0}, {0}});
    const Document stereo = *chem::fromSmiles("F[C@H](Cl)Br");
    CHECK(chem::substructureMatches(stereo, "F[C@H](Cl)Br")->size() == 1);
    CHECK(chem::substructureMatches(stereo, "F[C@@H](Cl)Br")->empty());
    Document many;
    for (int i = 0; i < 1100; ++i) many.addAtom({double(i), 0}, 8);
    CHECK(chem::substructureMatches(many, "O")->size() == 1100);
    CHECK(chem::substructureMatches(Document{}, "O")->empty());
}

TEST_CASE("bad SMILES is rejected") {
    CHECK_FALSE(chem::fromSmiles("C1CC"));
}

TEST_CASE("drawing style definitions preserve presets and reject invalid data (#625)") {
    for (const auto& preset : builtInDrawingStyles()) {
        const auto s = styleSettings(preset);
        const auto back = styleFromJson(styleToJson(s));
        REQUIRE(back);
        CHECK(*back == s);
        Document original = *chem::fromSmiles("C[C@H](N)C(=O)O");
        original.style = preset.name;
        Document portable = original;
        portable.styleSettings = s;
        CHECK(renderImage(original) == renderImage(portable));
        CHECK(documentStyle(portable).lineWidth == Catch::Approx(preset.lineWidth));
        CHECK(exportScale(portable) == Catch::Approx(preset.bondLength / kBondLength));
        CHECK(Document::fromJson(portable.toJson()) == portable);
    }
    const auto good = styleToJson(styleSettings(drawingStyle("RSC")));
    for (const auto& key : {"format", "version", "name", "font", "weight", "fontSize", "bondLength", "lineWidth",
                           "boldWidth", "wedgeWidth", "hashSpacing", "bondSpacing", "labelRadius", "labelMargin"}) {
        auto missing = good;
        missing.remove(key);
        CHECK_FALSE(styleFromJson(missing));
    }
    for (const auto& pair : {std::pair{"bondLength", QJsonValue(0)}, {"lineWidth", QJsonValue(-1)},
          {"fontSize", QJsonValue(0)}, {"weight", QJsonValue(400.5)}, {"name", QJsonValue("  ")},
          {"bondSpacing", QJsonValue(2)}, {"version", QJsonValue(2)}, {"hashSpacing", QJsonValue("2")}}) {
        auto bad = good;
        bad[pair.first] = pair.second;
        CHECK_FALSE(styleFromJson(bad));
    }
}

TEST_CASE("custom styles are snapshots shared by files, pages and exports (#625)") {
    auto s = styleSettings(drawingStyle("RSC"));
    s.name = "Penzene style regression";
    deleteDrawingStyle(s.name);
    struct Cleanup { QString name; ~Cleanup() { deleteDrawingStyle(name); } } cleanup{s.name};
    REQUIRE(saveDrawingStyle(s));
    CHECK_FALSE(saveDrawingStyle(s));
    auto duplicate = s;
    duplicate.name = s.name.toUpper();
    CHECK_FALSE(saveDrawingStyle(duplicate));
    duplicate.name = "ACS 1996";
    CHECK_FALSE(saveDrawingStyle(duplicate));
    CHECK_FALSE(deleteDrawingStyle("RSC"));
    Document d = *chem::fromSmiles("[NH3]->[Pt+2](<-[NH3])(Cl)Cl");
    d.labelRatio = 2;
    REQUIRE(applyDrawingStyle(d, s.name));
    CHECK(d.labelRatio == 0);
    CHECK(d.styleSettings == s);
    const auto png = renderImage(d);
    const auto sheets = sheetsToJson({{"A", d}, {"B", d}});
    const auto before = d.toJson();
    auto changed = s;
    changed.bondLength = 24;
    REQUIRE(saveDrawingStyle(changed, s.name));
    CHECK(renderImage(d) == png);
    Document newlyApplied = d;
    REQUIRE(applyDrawingStyle(newlyApplied, s.name));
    CHECK(exportScale(newlyApplied) == Catch::Approx(24 / kBondLength));
    REQUIRE(deleteDrawingStyle(s.name));
    CHECK_FALSE(applyDrawingStyle(newlyApplied, s.name));
    REQUIRE(Document::fromJson(before));
    CHECK(renderImage(*Document::fromJson(before)) == png);
    const std::vector<Sheet> expected{{"A", d}, {"B", d}};
    CHECK(sheetsFromJson(sheets) == expected);
    CHECK(sheetsFromJson(sheetsToJsonV1({{"A", d}, {"B", d}}))[1].doc == d);
    CHECK(Document::fromEmbedded(renderSvg(d)) == d);
    CHECK(Document::fromEmbedded(renderPng(d)) == d);
    CHECK(Document::fromEmbedded(renderPdf(d)) == d);
    CHECK(chem::toCdxml(d).contains("BondLength=\"12.2\""));
    auto broken = QJsonDocument::fromJson(before).object();
    auto settings = broken["styleSettings"].toObject();
    settings["fontSize"] = -1;
    broken["styleSettings"] = settings;
    CHECK_FALSE(Document::fromJson(QJsonDocument(broken).toJson()));
}

TEST_CASE("circled charges cover signs, magnitudes, bare carbons and portable styles (#623)") {
    for (const std::string smiles : {"[Na+]", "[Cl-]", "[Pt+2]", "[O-2]", "C[CH2+]", "[N+](C)(C)(C)C", "[O-][N+](=O)c1ccccc1"}) {
        INFO(smiles);
        Document d = *chem::fromSmiles(smiles);
        REQUIRE(chem::properties(d));
        const auto formula = chem::properties(d)->formula;
        const auto chemistry = chem::toSmiles(d);
        const QRectF plain = outputBounds(d);
        const auto image = renderImage(d);
        d.circledCharges = true;
        CHECK(chem::properties(d)->formula == formula);
        CHECK(chem::toSmiles(d) == chemistry);
        const auto singles = std::count_if(d.atoms.begin(), d.atoms.end(), [](const Atom& a) { return std::abs(a.charge) == 1; });
        CHECK((renderImage(d) != image) == (singles > 0));
        CHECK(outputBounds(d).width() >= plain.width());
        if (d.atoms.size() == 1 && singles) CHECK(outputBounds(d).width() > plain.width());
        CHECK(renderSvg(d).count("<circle") == singles);
        CHECK(Document::fromJson(d.toJson()) == d);
        CHECK(Document::fromEmbedded(renderPng(d)) == d);
        const auto pages = sheetsFromJson(sheetsToJson({{"A", d}, {"B", d}}));
        REQUIRE(pages.size() == 2);
        CHECK(pages[1].doc.circledCharges);
        QPicture picture;
        QPainter p(&picture);
        std::vector<QRectF> labels;
        RenderStyle style;
        style.labels = &labels;
        paintDocument(p, d, style);
        p.end();
        REQUIRE(labels.size() == d.atoms.size());
        for (size_t i = 0; i < d.atoms.size(); ++i)
            if (d.atoms[i].charge) CHECK_FALSE(labels[i].isEmpty());
        CHECK(outputBounds(d).contains(QRectF(picture.boundingRect())));
        auto settings = styleSettings(drawingStyle("RSC"));
        settings.name = "Circled custom style";
        settings.fontSize = 12;
        d.style = settings.name;
        d.styleSettings = settings;
        CHECK(renderSvg(d).count("<circle") == singles);
        CHECK(Document::fromEmbedded(renderSvg(d)) == d);
    }
    Document group;
    const int atom = group.addAtom({}, 8);
    REQUIRE(edit::applyLabel(group, atom, "OMe"));
    group.atoms[atom].charge = 1;
    group.circledCharges = true;
    CHECK(renderSvg(group).contains("<circle"));
    for (int charge : {-12, -2, 2, 12}) {
        group.atoms[atom].charge = charge;
        CHECK_FALSE(renderSvg(group).contains("<circle"));
        Document plain = group;
        plain.circledCharges = false;
        CHECK(renderImage(group) == renderImage(plain));
    }
    CHECK_FALSE(Document::fromJson(R"({"format":"penzene","version":1,"atoms":[],"bonds":[]})")->circledCharges);
}

TEST_CASE("circled charges preserve visible bonds around a crowded metal centre (#623)") {
    for (const auto& preset : builtInDrawingStyles()) {
        INFO(preset.name.toStdString());
        Document d = *chem::fromSmiles("[NH3]->[Pt+](<-[NH3])(Cl)Cl");
        REQUIRE(applyDrawingStyle(d, preset.name));
        // Font metrics can trim a crowded bond; adding a circle must preserve the normal rendering.
        const auto visible = renderSvg(d).count("<polyline");
        REQUIRE(visible > 0);
        d.circledCharges = true;
        const auto circled = renderSvg(d);
        REQUIRE(circled.contains("<circle"));
        CHECK(circled.count("<polyline") == visible);
    }
}

TEST_CASE("a bond stops at an abbreviation's text, not inside it (#677)") {
    for (const char* text : {"CF3", "Me", "Fmoc", "MgBr"}) {
        for (double side : {-1.0, 1.0}) {
            INFO(text << " bond on the " << (side > 0 ? "right" : "left"));
            Document d;
            const int label = d.addAtom({0, 0}, 6);
            REQUIRE(edit::applyLabel(d, label, text));
            const int carbon = d.addAtom({side * kBondLength, 0});
            d.bonds.push_back({label, carbon});
            auto matches = QRegularExpression(R"re(<polyline[^>]*points="([^"]+)")re").globalMatch(renderSvg(d));
            REQUIRE(matches.hasNext());
            const QStringList points = matches.next().captured(1).trimmed().split(' ');
            REQUIRE(points.size() == 2);
            const double begin = points.front().section(',', 0, 0).toDouble();  // at the label atom
            CHECK(side * begin > 1);  // clear of the attaching letter, not through it
        }
    }
}

TEST_CASE("journal presets match publisher measurements and ChemDraw physical sizes (#473)") {
    for (const auto& [name, bond, font, line, bold, margin, hash] : {
        std::tuple{"Nature", 10.8, 6.0, 0.5952755906, 1.5590551181, 1.1905511811, 1.7007874016},
        {"Angewandte", 17.0 * 2 / 3, 8.0, 0.5, 1.734375, 4.0 / 3, 1.734375},
        {"SYNLETT/SYNTHESIS", 10.1196850394, 6.0, 0.5669291339, 1.3889763780, 0.9070866142, 1.7574803150}}) {
        INFO(name);
        const auto s = styleSettings(drawingStyle(name));
        CHECK(s.name == name);
        CHECK(s.bondLength == Catch::Approx(bond).margin(1e-6));
        CHECK(s.fontSize == Catch::Approx(font));
        CHECK(s.lineWidth == Catch::Approx(line).margin(1e-6));
        CHECK(s.boldWidth == Catch::Approx(bold).margin(1e-6));
        CHECK(s.labelMargin == Catch::Approx(margin).margin(1e-6));
        CHECK(s.hashSpacing == Catch::Approx(hash).margin(1e-6));
        CHECK(s.bondSpacing == Catch::Approx(0.18));
        CHECK(s.wedgeWidth == Catch::Approx(1.5 * bold));
        CHECK(s.labelRadius == 0);
        CHECK(s.chainAngle == 120);
        CHECK_FALSE(saveDrawingStyle(s));
        Document d = *chem::fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
        d.circledCharges = true;
        REQUIRE(applyDrawingStyle(d, name));
        CHECK(d.circledCharges == (QString(name) != "Nature"));  // Nature applies its explicit publisher convention
        CHECK_FALSE(d.styleSettings);
        CHECK(exportScale(d) == Catch::Approx(bond / kBondLength));
        const auto xml = chem::toCdxml(d);
        QXmlStreamReader reader(xml);
        REQUIRE(reader.readNextStartElement());
        CHECK(reader.attributes().value("BondLength").toDouble() == Catch::Approx(bond));
        CHECK(reader.attributes().value("LabelSize").toDouble() == Catch::Approx(font));
        CHECK(reader.attributes().value("LineWidth").toDouble() == Catch::Approx(line));
        const auto back = chem::fromChemDraw(xml);
        REQUIRE(back);
        CHECK(back->style == name);
        CHECK(chem::toSmiles(*back) == chem::toSmiles(d));
        CHECK(Document::fromJson(d.toJson()) == d);
        CHECK_FALSE(renderPng(d).isEmpty());
        CHECK_FALSE(renderPdf(d).isEmpty());
    }
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

TEST_CASE("a peptide from its one- or three-letter sequence (#503)") {
    const auto gfls = chem::fromSequence("GFLS");
    REQUIRE(gfls);
    CHECK(chem::properties(*gfls)->formula == "C20H30N4O6");
    for (const char* same : {"Gly-Phe-Leu-Ser", "GLY PHE LEU SER", " G F L S ", "H-Gly-Phe-Leu-Ser-OH", "Gly-L-Phe-Leu-Ser"}) {
        INFO(same);
        const auto doc = chem::fromSequence(same);
        REQUIRE(doc);
        CHECK(chem::toSmiles(*doc) == chem::toSmiles(*gfls));
    }
    CHECK(chem::properties(*chem::fromSequence("Gly"))->formula == "C2H5NO2");  // glycine, not Gly-Leu-Tyr
    CHECK(chem::properties(*chem::fromSequence("GLY"))->formula == "C17H25N3O5");  // Gly-Leu-Tyr
    CHECK(chem::toSmiles(*chem::fromSequence("a")) != chem::toSmiles(*chem::fromSequence("A")));  // D- and L-alanine
    CHECK_FALSE(chem::fromSequence(""));
    CHECK_FALSE(chem::fromSequence("G1S"));
    CHECK(chem::toSmiles(*chem::fromSequence("Ala-D-Phe")) == chem::toSmiles(*chem::fromSequence("Af")));  // D-phenylalanine (#570)
    CHECK(chem::properties(*chem::fromSequence("H-GLY-OH"))->formula == "C2H5NO2");
    for (const char* bad : {"Gly-Xyz", "Ac-NH2", "Gly-Ac", "NH2-Gly", "Ac-Ac-Gly", "Ac-H-Gly", "Gly-OH-NH2",
                            "Gly-NH2-NH2", "Ala-D-D-Phe", "G-F-L"})
        CHECK_FALSE(chem::fromSequence(bad));
}

TEST_CASE("peptide caps change only the termini and keep stereochemistry (#572)") {
    for (auto [sequence, smiles] : {std::pair{"Ac-Gly", "CC(=O)NCC(=O)O"}, {"Gly-NH2", "NCC(N)=O"},
                                    {"Ac-Gly-Phe-NH2", "CC(=O)NCC(=O)N[C@@H](Cc1ccccc1)C(N)=O"},
                                    {"Ac-Pro-NH2", "CC(=O)N1CCC[C@H]1C(N)=O"},
                                    {"Ac-Lys-Asp-NH2", "CC(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CC(=O)O)C(N)=O"}}) {
        INFO(sequence);
        const auto doc = chem::fromSequence(sequence);
        REQUIRE(doc);
        CHECK(chem::toSmiles(*doc) == chem::toSmiles(*chem::fromSmiles(smiles)));
        CHECK(chem::properties(*doc)->formula == chem::properties(*chem::fromSmiles(smiles))->formula);
        CHECK(chem::checkStructure(*doc).empty());
    }
    CHECK(chem::toSmiles(*chem::fromSequence("Ac-Gly-D-Phe-NH2")) !=
          chem::toSmiles(*chem::fromSequence("Ac-Gly-Phe-NH2")));
    CHECK(chem::toSmiles(*chem::fromSequence("H-Gly-NH2")) == chem::toSmiles(*chem::fromSequence("Gly-NH2")));
    CHECK(chem::toSmiles(*chem::fromSequence("Ac-Gly-OH")) == chem::toSmiles(*chem::fromSequence("Ac-Gly")));
    CHECK(chem::toSmiles(*chem::fromSequence("Ac-GF-NH2")) == chem::toSmiles(*chem::fromSequence("Ac-Gly-Phe-NH2")));
    if (const QString preview = qEnvironmentVariable("PENZENE_PEPTIDE_PREVIEW"); !preview.isEmpty())
        REQUIRE(exportDocument(*chem::fromSequence("Ac-Gly-Phe-NH2"), preview, {.background = Qt::white, .margin = 10}));
}

TEST_CASE("a peptide's side chains keep clear of the backbone's N-H (#559)") {
    for (const char* seq : {"YGGFL", "ACDEFGHIKLMNPQRSTVWY"}) {
        INFO(seq);
        const Document d = *chem::fromSequence(seq);
        const auto info = chem::atomInfo(d);
        double closest = 1e9;  // from where an N-H's H goes (away from its two bonds) to the nearest atom
        for (int i = 0; i < int(d.atoms.size()); ++i) {
            const auto nbs = d.neighbors(i);
            if (d.atoms[i].z != 7 || nbs.size() != 2 || info[i].hydrogens != 1) continue;
            QPointF out;
            for (int n : nbs) out -= unit(d.atoms[n].pos - d.atoms[i].pos);
            const QPointF h = d.atoms[i].pos + unit(out) * kBondLength;
            for (int j = 0; j < int(d.atoms.size()); ++j)
                if (j != i) closest = std::min(closest, len(d.atoms[j].pos - h) / kBondLength);
        }
        CHECK(closest > 0.6);
    }
}

TEST_CASE("macrocycles and peptides are laid out cleanly (#502)") {
    // A cyclophane's bridges leave their rings well clear of the ring bonds, not squeezed against one.
    const Document phane = *chem::fromSmiles("C1Cc2ccc(cc2)CCc2ccc1cc2");  // [2.2]paracyclophane
    double tightest = 180;  // degrees between a bridge bond and a ring bond at the same atom
    for (const Bond& b : phane.bonds) {
        for (int end : {b.a, b.b}) {
            const int other = end == b.a ? b.b : b.a;
            if (phane.neighbors(end).size() != 3 || phane.neighbors(other).size() != 2) continue;  // ring atom, bridge CH2
            for (int n : phane.neighbors(end))
                if (n != other)
                    tightest = std::min(tightest, qRadiansToDegrees(std::acos(QPointF::dotProduct(
                        unit(phane.atoms[other].pos - phane.atoms[end].pos), unit(phane.atoms[n].pos - phane.atoms[end].pos)))));
        }
    }
    CHECK(tightest > 80);

    // A peptide's backbone runs out end to end instead of folding back, where each C=O met the next NH.
    const Document gfls = *chem::fromSmiles("NCC(=O)NC(Cc1ccccc1)C(=O)NC(CC(C)C)C(=O)NC(CO)C(=O)O");  // Gly-Phe-Leu-Ser
    CHECK(len(gfls.atoms.front().pos - gfls.atoms.back().pos) / kBondLength > 9);  // N-terminus to the C-terminal OH
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

TEST_CASE("the 1 hotkey grows one straight chain (#676)") {
    Document d;
    d.addAtom({0, 0});
    edit::Hotspot h{0, -1};
    for (int i = 0; i < 8; ++i) h = edit::hotkey(d, h, "1");
    REQUIRE(d.atoms.size() == 9);
    for (int i = 1; i < 9; ++i) {
        CHECK(d.atoms[i].pos.x() > d.atoms[i - 1].pos.x());
        CHECK(len(d.atoms[i].pos - d.atoms[i - 1].pos) == Catch::Approx(kBondLength));
    }
    const double rise = std::abs(d.atoms[1].pos.y() - d.atoms[0].pos.y());
    for (int i = 2; i < 9; ++i)
        CHECK(std::abs(std::abs(d.atoms[i].pos.y() - d.atoms[i - 1].pos.y()) - rise) < 1e-9);
}

TEST_CASE("0 closes a hotkey chain into a ring (#676)") {
    Document d;
    d.addAtom({0, 0});
    edit::Hotspot h{0, -1};
    for (int i = 0; i < 6; ++i) h = edit::hotkey(d, h, "1");
    REQUIRE(d.atoms.size() == 7);
    const QPointF hotspot = d.atoms[h.atom].pos;
    h = edit::hotkey(d, h, "0");
    CHECK(h.valid());
    CHECK(d.atoms.size() == 7);
    CHECK(d.bonds.size() == 7);
    CHECK(d.bondBetween(0, 6) >= 0);
    CHECK(len(d.atoms[h.atom].pos - hotspot) < 1e-9);  // the hotspot atom stays put
    for (const Bond& b : d.bonds)
        CHECK(len(d.atoms[b.a].pos - d.atoms[b.b].pos) == Catch::Approx(kBondLength));
    for (int i = 0; i < 7; ++i) CHECK(d.neighbors(i).size() == 2);
}

TEST_CASE("chair hotkey on two bonds of one ring fuses cleanly (#417)") {
    auto clean = [](const Document& doc, bool unitBonds = true) {
        for (int i = 0; i < int(doc.atoms.size()); ++i) {
            CHECK(doc.neighbors(i).size() >= 2);  // no dangling bond
            for (int j = 0; j < i; ++j)           // nothing drawn over or next to another atom
                if (doc.bondBetween(i, j) < 0) CHECK(len(doc.atoms[i].pos - doc.atoms[j].pos) > 0.45 * kBondLength);
        }
        for (const Bond& x : doc.bonds) {
            if (unitBonds) CHECK(std::abs(len(doc.atoms[x.a].pos - doc.atoms[x.b].pos) - kBondLength) < 0.15 * kBondLength);
            for (const Bond& y : doc.bonds)
                if (x.a != y.a && x.a != y.b && x.b != y.a && x.b != y.b)
                    CHECK(QLineF(doc.atoms[x.a].pos, doc.atoms[x.b].pos)
                              .intersects(QLineF(doc.atoms[y.a].pos, doc.atoms[y.b].pos)) != QLineF::BoundedIntersection);
        }
    };
    std::vector<QPointF> firstChair[2];
    for (int k : {0, 1})
        for (int second = 1; second < 6; ++second) {
            const QString key = k ? "0" : "9";
            Document doc;
            edit::ringAt(doc, {0, 0}, 6, false);
            edit::hotkey(doc, {-1, 0}, key);
            INFO(key.toStdString() << " then bond " << second);
            REQUIRE(doc.atoms.size() == 10);
            clean(doc);
            if (second == 1)
                for (const Atom& a : doc.atoms) firstChair[k].push_back(a.pos);
            edit::hotkey(doc, {-1, second}, key);
            CHECK(doc.atoms.size() == 14);
            clean(doc);
        }
    CHECK(firstChair[0] != firstChair[1]);  // 9 and 0 are mirror images

    // A chair fused onto any bond of another chair is just as clean (#449).
    for (const QString first : {"9", "0"})
        for (const QString key : {"9", "0"})
            for (int bond = 1; bond < 6; ++bond) {
                Document doc;
                const int a = doc.addAtom({0, 0}), b = doc.addAtom({kBondLength, 0});
                edit::link(doc, a, b);
                edit::hotkey(doc, {-1, 0}, first);  // a chair on the single bond
                INFO(first.toStdString() << " then " << key.toStdString() << " on chair bond " << bond);
                REQUIRE(doc.atoms.size() == 6);
                edit::hotkey(doc, {-1, bond}, key);
                CHECK(doc.atoms.size() == 10);
                clean(doc, false);  // the template's bonds aren't all one length, and fusing on one carries that on
            }
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
    // A reaction is a row that says so, not one molecule's numbers (#501).
    Document scheme = *chem::fromReactionSmiles("CC(=O)O>>CC(=O)OC");
    CHECK(QString::fromStdString(chem::descriptorsCsv({{"ester", scheme}}, {"name", "error"})) == "name,error\nester,reaction\n");
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

TEST_CASE("ChemDraw keeps drawing styles at their physical size (#624)") {
    for (const auto& style : drawingStyles())
        for (double ratio : {0.0, 0.9}) {
            INFO(style.name.toStdString() << " " << ratio);
            Document d = *chem::fromSmiles("CCO");
            REQUIRE(edit::applyLabel(d, 2, "OMe", true));
            d.style = style.name;
            d.labelRatio = ratio;
            d.texts.push_back({{0, 50}, "heat", 1.2});
            d.arrows.push_back({{50, 0}, {90, 0}, ArrowKind::Reaction, 6});
            d.arrows.back().head = 1.3;
            d.arrows.push_back({{50, 30}, {80, 60}, ArrowKind::Ellipse});
            d.arrows.push_back({{50, 80}, {90, 80}, ArrowKind::Equilibrium});
            const QByteArray xml = chem::toCdxml(d);
            QXmlStreamReader r(xml);
            std::vector<QPointF> nodes;
            while (!r.atEnd()) {
                r.readNext();
                if (!r.isStartElement()) continue;
                const auto attrs = r.attributes();
                if (r.name() == u"CDXML") {
                    CHECK(attrs.value("BondLength").toDouble() == Catch::Approx(style.bondLength));
                    CHECK(attrs.value("LabelSize").toDouble() == Catch::Approx(documentStyle(d).fontSize * exportScale(d)));
                    CHECK(attrs.value("LineWidth").toDouble() == Catch::Approx(style.lineWidth * exportScale(d)));
                }
                if (r.name() == u"n") {
                    const auto xy = attrs.value("p").toString().split(' ');
                    REQUIRE(xy.size() == 2);
                    nodes.push_back({xy[0].toDouble(), xy[1].toDouble()});
                }
            }
            REQUIRE(nodes.size() > 2);
            CHECK(len(nodes[1] - nodes[0]) == Catch::Approx(style.bondLength).margin(0.02));
            for (const QByteArray& file : {xml, chem::toCdx(d)}) {
                INFO((file.startsWith("VjCD") ? "CDX" : "CDXML"));
                auto back = chem::fromChemDraw(file);
                REQUIRE(back);
                CHECK(drawingStyle(back->style).bondLength == Catch::Approx(style.bondLength));
                CHECK(documentStyle(*back).fontSize == Catch::Approx(documentStyle(d).fontSize).margin(0.06));
                CHECK(chem::properties(*back)->formula == chem::properties(d)->formula);
                REQUIRE(back->atoms.size() == d.atoms.size());
                const QPointF moved = back->atoms[0].pos - d.atoms[0].pos;
                for (size_t i = 0; i < d.atoms.size(); ++i) CHECK(len(back->atoms[i].pos - moved - d.atoms[i].pos) < 0.04);
                REQUIRE(back->texts.size() == 1);
                CHECK(back->texts[0].scale == Catch::Approx(1.2).margin(0.01));
                CHECK(len(back->texts[0].pos - moved - d.texts[0].pos) < 0.04);
                REQUIRE(back->arrows.size() == d.arrows.size());
                for (const auto& arrow : d.arrows) {
                    // Binary ChemDraw groups graphics and arrows separately.
                    const auto found = std::find_if(back->arrows.begin(), back->arrows.end(), [&](const Arrow& a) { return a.kind == arrow.kind; });
                    REQUIRE(found != back->arrows.end());
                    CHECK(len(found->from - moved - arrow.from) < 0.04);
                    CHECK(len(found->to - moved - arrow.to) < 0.04);
                    CHECK(found->bend == Catch::Approx(arrow.bend).margin(0.04));
                    // CDX rounds stroke widths to hundredths of a point; compare physical head lengths.
                    if (!isShape(arrow.kind)) CHECK(found->head * kHeadLength * exportScale(d) ==
                        Catch::Approx(arrow.head * kHeadLength * exportScale(d)).margin(0.05));
                }
            }
        }
    auto drawing = [](double bond) {
        return QString(R"(<CDXML BondLength="%1" LabelSize="9"><page><fragment><n id="1" p="0 0"/><n id="2" p="%1 0"/><b B="1" E="2"/></fragment></page></CDXML>)").arg(bond).toUtf8();
    };
    for (double bond : {12.15, 12.2, 12.25}) {
        auto matched = chem::fromChemDraw(drawing(bond));
        REQUIRE(matched);
        CHECK(matched->style == "RSC");
    }
    auto unmatched = chem::fromChemDraw(drawing(12.26));
    REQUIRE(unmatched);
    CHECK(unmatched->style.isEmpty());
    CHECK(unmatched->labelRatio == Catch::Approx(9 / 12.26));
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
    CHECK(QFontMetricsF(labelFont(documentStyle(*small))).capHeight() < QFontMetricsF(labelFont(drawingStyle(""))).capHeight());
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

TEST_CASE("a ChemDraw label node without an id leaves the molecules alone (#372)") {
    // tests/data/labels.cdxml, with the reagent label's fragment and node left without ids.
    const QByteArray cdxml = R"(<?xml version="1.0"?><CDXML BondLength="14.40"><page id="1">
        <fragment id="10"><n id="11" p="100 100"/>
        <n id="12" p="112.47 107.20" NodeType="GenericNickname" GenericNickname="R"><t p="109 111"><s size="10">R</s></t></n>
        <b id="13" B="11" E="12"/></fragment>
        <fragment><n p="160 80" NodeType="Unspecified"><t p="150 83" BoundingBox="150 76 190 83"><s size="7">LiBr, acetone</s></t></n></fragment>
        </page></CDXML>)";
    auto doc = chem::fromChemDraw(cdxml);
    REQUIRE(doc);
    CHECK(doc->atoms.size() == 2);  // CH3-R
    CHECK(doc->bonds.size() == 1);
    REQUIRE(doc->texts.size() == 1);
    CHECK(doc->texts[0].text == "LiBr, acetone");
}

TEST_CASE("CDXML export keeps abbreviations as ChemDraw nicknames (#333)") {
    Document d = *chem::fromSmiles("c1ccccc1");
    const int boc = d.addAtom({d.atoms[0].pos.x() + kBondLength, d.atoms[0].pos.y()});
    d.bonds.push_back({0, boc});
    REQUIRE(edit::applyLabel(d, boc, "Boc", true));
    const int ome = d.addAtom({d.atoms[3].pos.x() - kBondLength, d.atoms[3].pos.y()});
    d.bonds.push_back({3, ome});
    REQUIRE(edit::applyLabel(d, ome, "OMe", true));
    const QByteArray cdxml = chem::toCdxml(d);
    CHECK(cdxml.count("NodeType=\"Nickname\"") == 2);  // not drawn out, and not a label without chemistry
    CHECK(cdxml.count("ExternalConnectionPoint") == 2);
    for (const QByteArray& file : {cdxml, chem::toCdx(d)}) {
        if (file.isEmpty()) continue;  // a build without binary CDX
        auto back = chem::fromChemDraw(file);
        REQUIRE(back);
        CHECK(back->atoms.size() == d.atoms.size());
        CHECK(chem::properties(*back)->formula == chem::properties(d)->formula);
        QStringList labels;
        for (const auto& a : back->atoms)
            if (!a.label.isEmpty()) labels << a.label;
        labels.sort();
        CHECK(labels == QStringList{"Boc", "OMe"});
    }
}

TEST_CASE("ChemDraw export keeps a charged abbreviation's charge (#382)") {
    Document d;
    d.addAtom({0, 0});
    REQUIRE(edit::applyLabel(d, 0, "N3-", true));
    const QByteArray cdxml = chem::toCdxml(d);
    CHECK(cdxml.count("Charge=\"-1\"") == 2);  // N-=N+=N-: the attaching N carries the typed charge too
    CHECK(cdxml.count("Charge=\"1\"") == 1);
}

TEST_CASE("a ChemDraw-made file keeps its Boc nickname and OMe fragment labels (#391)") {
    auto doc = chem::readFile(QString(PENZENE_TEST_DATA) + "/chemdraw-nicknames.cdxml");  // ChemDraw 22
    REQUIRE(doc);
    REQUIRE(doc->atoms.size() == 8);  // the ring plus one atom per label
    CHECK(doc->atoms[6].label == "Boc");
    CHECK(doc->atoms[7].label == "OMe");
    CHECK(chem::properties(*doc)->formula == "C12H16O3");
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
    QJsonObject old = QJsonDocument::fromJson(eth->toJson()).object();
    old["hideImplicitH"] = true;  // saved by 1.4: the labels keep their hydrogens now
    auto back = Document::fromJson(QJsonDocument(old).toJson());
    REQUIRE(back);
    CHECK(back->carbonLabels == Document::CarbonLabels::Terminal);
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

TEST_CASE("isotope pattern: aspirin's ions and a chlorine M+2 (#397)") {
    const Document aspirin = *chem::fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
    // Intensity summed per nominal mass above the monoisotopic peak (13C, 2H and 17O sticks apart).
    auto bins = [](const std::vector<chem::Peak>& peaks) {
        std::map<int, double> sum;
        for (const auto& p : peaks) sum[int(std::lround(p.mz - peaks[0].mz))] += p.intensity;
        return sum;
    };
    const auto m = chem::isotopePattern(aspirin, chem::Ion::M);
    REQUIRE_FALSE(m.empty());
    CHECK(std::abs(m[0].mz - 180.0417) < 0.0001);  // M+•: the neutral 180.0423 less an electron (#556)

    const auto mh = chem::isotopePattern(aspirin, chem::Ion::MplusH);
    REQUIRE(mh.size() > 2);
    CHECK(std::abs(mh[0].mz - 181.0495) < 0.0001);
    CHECK(mh[0].intensity == 100);
    // C9H9O4: M+1/M = Σ n·a(M+1)/a(M) over C, H and O.
    const double theory = 9 * 1.07 / 98.93 + 9 * 0.0115 / 99.9885 + 4 * 0.038 / 99.757;
    const double ratio = bins(mh)[1] / bins(mh)[0];
    INFO("M+1/M " << ratio << " theory " << theory);
    CHECK(std::abs(ratio / theory - 1) < 0.005);
    for (const auto& p : mh) CHECK(p.intensity >= 0.1);

    const auto mna = chem::isotopePattern(aspirin, chem::Ion::MplusNa);
    REQUIRE_FALSE(mna.empty());
    CHECK(std::abs(mna[0].mz - 203.0315) < 0.0001);
    const auto deprotonated = chem::isotopePattern(aspirin, chem::Ion::MminusH);
    REQUIRE_FALSE(deprotonated.empty());
    CHECK(std::abs(deprotonated[0].mz - 179.0350) < 0.0001);

    const auto cl = bins(chem::isotopePattern(*chem::fromSmiles("Clc1ccccc1"), chem::Ion::M));
    INFO("M+2/M " << cl.at(2) / cl.at(0));
    CHECK(std::abs(100 * cl.at(2) / cl.at(0) - 32) < 1);
    // A drawn 13C is that isotope only: one mass unit up, and no more 13C to spread.
    const auto labelled = chem::isotopePattern(*chem::fromSmiles("[13CH4]"), chem::Ion::M);
    REQUIRE_FALSE(labelled.empty());
    CHECK(std::abs(labelled[0].mz - 17.0341) < 0.0001);
    CHECK(labelled.size() == 1);
    CHECK(chem::isotopePattern(Document{}, chem::Ion::M).empty());
    CHECK(chem::isotopePattern(*chem::fromSmiles("[Na+].[Cl-]"), chem::Ion::MminusH).empty());
}

TEST_CASE("EI ions to look for, from the groups present (#554)") {
    auto find = [](const std::vector<chem::EiIon>& ions, const QString& formula) {
        auto it = std::find_if(ions.begin(), ions.end(), [&](const chem::EiIon& i) { return i.formula == formula; });
        return it == ions.end() ? 0.0 : it->mz;
    };
    const auto ester = chem::eiIons(*chem::fromSmiles("COC(=O)Cc1ccccc1"));  // methyl phenylacetate
    REQUIRE_FALSE(ester.empty());
    CHECK(ester[0].formula == "C9H10O2+");
    CHECK(std::abs(find(ester, "C7H7+") - 91.0542) < 1e-4);   // tropylium
    CHECK(std::abs(find(ester, "C8H7O+") - 119.0491) < 1e-4);  // M − OMe
    CHECK(find(ester, "C8H7O2+") == 0);                        // its methyl is on O, so no M − 15
    const auto ketone = chem::eiIons(*chem::fromSmiles("CC(=O)CCCC"));  // 2-hexanone
    CHECK(std::abs(find(ketone, "C3H6O+") - 58.0413) < 1e-4);  // McLafferty: propene lost
    CHECK(std::abs(find(ketone, "C2H3O+") - 43.0178) < 1e-4);  // acetyl
    const auto ethanol = chem::eiIons(*chem::fromSmiles("CCO"));
    CHECK(find(ethanol, "C7H7+") == 0);
    CHECK(find(ethanol, "C2H4+") > 0);  // M − H2O
    CHECK(chem::eiIons(Document{}).empty());
    CHECK(chem::eiIons(*chem::fromSmiles("CC(=O)[O-]")).empty());
}

TEST_CASE("the HRMS line for a supporting-information entry (#398)") {
    const Document aspirin = *chem::fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
    CHECK(chem::hrmsLine(aspirin, chem::Ion::MplusH).toStdString() == "HRMS (ESI) m/z: [M+H]+ calcd for C9H9O4 181.0495");
    CHECK(chem::hrmsLine(aspirin, chem::Ion::MplusNa).toStdString() == "HRMS (ESI) m/z: [M+Na]+ calcd for C9H8NaO4 203.0315");
    CHECK(chem::hrmsLine(aspirin, chem::Ion::MminusH).toStdString() == "HRMS (ESI) m/z: [M-H]- calcd for C9H7O4 179.0350");
    CHECK(chem::hrmsLine(aspirin, chem::Ion::M).toStdString() == "HRMS (EI) m/z: [M]+ calcd for C9H8O4 180.0417");
    CHECK(chem::hrmsLine(*chem::fromSmiles("C[N+](C)(C)C"), chem::Ion::M).toStdString() == "HRMS (ESI) m/z: [M]+ calcd for C4H12N+ 74.0964");
    CHECK(chem::hrmsLine(Document{}, chem::Ion::MplusH).isEmpty());
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
    // A reaction SMILES line is a scheme, not an unreadable molecule (#501).
    const auto rxn = chem::readRecords(write("rxn.smi", "CC(=O)O.OC>>CC(=O)OC ester\n"));
    REQUIRE(rxn.size() == 1);
    REQUIRE(rxn[0].doc);
    CHECK(rxn[0].doc->arrows.size() == 1);
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

TEST_CASE("multi-step reaction SMILES, one step a line, read back as one scheme (#591)") {
    const std::string text = "CCO>>CC=O\r\nCC=O>>CC(=O)O\n";  // as copied, with Windows line ends too
    auto doc = chem::fromReactionSmiles(text);
    REQUIRE(doc);
    CHECK(doc->arrows.size() == 2);
    const auto steps = chem::reactionsOf(*doc);
    REQUIRE(steps.size() == 2);
    CHECK(doc->atoms.size() == 3 + 3 + 4);  // the shared intermediate drawn once
    CHECK(chem::toReactionSmiles(steps) == chem::toReactionSmiles(chem::reactionsOf(*chem::fromReactionSmiles(chem::toReactionSmiles(steps)))));
    CHECK_FALSE(chem::fromReactionSmiles("CCO>>CC=O\nnot a reaction"));  // any bad line: nothing
}

TEST_CASE("an RD file opens as one scheme (#387)") {
    // Ethanol -> acetaldehyde -> acetic acid, written as two Rxnfiles.
    auto mol = [](const char* s) { return *chem::fromSmiles(s); };
    const std::vector<chem::Reaction> steps{{{mol("CCO")}, {}, {mol("CC=O")}}, {{mol("CC=O")}, {}, {mol("CC(=O)O")}}};
    QTemporaryDir dir;
    const QString path = dir.filePath("scheme.rdf");
    QFile f(path);
    REQUIRE(f.open(QIODevice::WriteOnly));
    f.write(QByteArray::fromStdString(chem::toRdf(steps)));
    f.close();
    auto doc = chem::readFile(path);
    REQUIRE(doc);
    CHECK(doc->arrows.size() == 2);
    CHECK(chem::molecules(*doc).size() == 3);  // the aldehyde once, shared by both steps
    CHECK(chem::toReactionSmiles(chem::reactionsOf(*doc)) == chem::toReactionSmiles(steps));
    CHECK_FALSE(chem::fromRdf("$RDFILE 1\n"));
}

TEST_CASE("a multi-step scheme exports every step (#334)") {
    // Ethanol -> acetaldehyde -> acetic acid, in a row.
    Document doc = *chem::fromReactionSmiles("CCO>>CC=O");
    double right = -1e9, mid = 0;
    for (const auto& a : doc.atoms) right = std::max(right, a.pos.x());
    for (const auto& a : doc.arrows) mid = a.from.y();
    doc.arrows.push_back({{right + kBondLength, mid}, {right + 4 * kBondLength, mid}});
    Document acid = *chem::fromSmiles("CC(=O)O");
    double left = 1e9, y = 0;
    for (const auto& a : acid.atoms) left = std::min(left, a.pos.x()), y += a.pos.y() / acid.atoms.size();
    doc.append(acid, QPointF(right + 5 * kBondLength - left, mid - y));

    const auto steps = chem::reactionsOf(doc);
    REQUIRE(steps.size() == 2);
    auto canon = [](const std::string& s) { return chem::toSmiles(*chem::fromSmiles(s)); };
    CHECK(chem::toReactionSmiles(steps) ==
          canon("CCO") + ">>" + canon("CC=O") + "\n" + canon("CC=O") + ">>" + canon("CC(=O)O"));  // the aldehyde is both
    const std::string rdf = chem::toRdf(steps);
    CHECK(rdf.starts_with("$RDFILE 1"));
    size_t rxns = 0;
    for (size_t at = rdf.find("$RFMT\n$RXN"); at != std::string::npos; at = rdf.find("$RFMT\n$RXN", at + 1)) ++rxns;
    CHECK(rxns == 2);
    CHECK(chem::reactionOf(doc)->products.size() == 1);  // the first step, as before
    // Arrows a little apart in height are still one row, read left to right.
    Document straddle = doc;
    straddle.arrows[0].from.ry() += 1, straddle.arrows[0].to.ry() += 1;
    for (double shift : {0.0, 4 * kBondLength - 0.5}) {  // wherever the row sits
        Document moved = straddle;
        for (auto& a : moved.atoms) a.pos.ry() += shift;
        for (auto& a : moved.arrows) a.from.ry() += shift, a.to.ry() += shift;
        CHECK(chem::toReactionSmiles(chem::reactionsOf(moved)) == chem::toReactionSmiles(steps));
    }
}

TEST_CASE("a scheme wrapped onto two rows keeps the link between them (#388)") {
    // Ethanol -> acetaldehyde ending row 1; -> acetic acid starting row 2.
    Document doc = *chem::fromReactionSmiles("CCO>>CC=O");
    const double below = doc.arrows[0].from.y() + 6 * kBondLength;
    doc.arrows.push_back({{0, below}, {3 * kBondLength, below}});
    Document acid = *chem::fromSmiles("CC(=O)O");
    double left = 1e9, y = 0;
    for (const auto& a : acid.atoms) left = std::min(left, a.pos.x()), y += a.pos.y() / acid.atoms.size();
    doc.append(acid, QPointF(4 * kBondLength - left, below - y));

    const auto steps = chem::reactionsOf(doc);
    REQUIRE(steps.size() == 2);
    auto canon = [](const std::string& s) { return chem::toSmiles(*chem::fromSmiles(s)); };
    CHECK(chem::toReactionSmiles(steps) ==
          canon("CCO") + ">>" + canon("CC=O") + "\n" + canon("CC=O") + ">>" + canon("CC(=O)O"));  // the aldehyde is both
}

TEST_CASE("double dashed bonds keep their order through ChemDraw round trips (#668)") {
    Document single;
    single.atoms = {{{0, 0}}, {{kBondLength, 0}}};
    single.bonds = {{0, 1, 1, BondStereo::DoubleDashed}};
    CHECK(chem::toSmiles(single) == "CC");
    auto json = Document::fromJson(single.toJson());
    REQUIRE(json);
    CHECK(json->bonds[0].stereo == BondStereo::DoubleDashed);
    CHECK(json->bonds[0].order == 1);
    auto back = chem::fromChemDraw(chem::toCdxml(single));
    REQUIRE(back);
    REQUIRE(back->bonds.size() == 1);
    CHECK(back->bonds[0].stereo == BondStereo::DoubleDashed);
    CHECK(back->bonds[0].order == 1);
    CHECK(chem::toSmiles(*back) == "CC");
    CHECK(renderSvg(single).contains("stroke-dasharray"));

    Document doubleBond;
    doubleBond.atoms = {{{0, 0}}, {{kBondLength, 0}}};
    doubleBond.bonds = {{0, 1, 2, BondStereo::DoubleDashed}};
    CHECK(chem::toSmiles(doubleBond) == "C=C");
    CHECK(renderSvg(doubleBond).contains("stroke-dasharray"));
    back = chem::fromChemDraw(chem::toCdxml(doubleBond));
    REQUIRE(back);
    REQUIRE(back->bonds.size() == 1);
    CHECK(back->bonds[0].stereo == BondStereo::DoubleDashed);
    CHECK(back->bonds[0].order == 2);
    CHECK(chem::toSmiles(*back) == "C=C");

    Document dashed = doubleBond;
    dashed.bonds[0].stereo = BondStereo::Dashed;  // the existing dashed double is not the new style
    back = chem::fromChemDraw(chem::toCdxml(dashed));
    REQUIRE(back);
    REQUIRE(back->bonds.size() == 1);
    CHECK(back->bonds[0].stereo == BondStereo::Dashed);
    CHECK(back->bonds[0].order == 2);
    CHECK(chem::toSmiles(*back) == "C=C");
}

TEST_CASE("a quadruple bond is saved, drawn and survives ChemDraw") {
    Document d;
    d.atoms = {{{0, 0}}, {{kBondLength, 0}}};
    d.bonds = {{0, 1, 4}};
    CHECK(chem::toSmiles(d) == "C$C");
    auto json = Document::fromJson(d.toJson());
    REQUIRE(json);
    CHECK(json->bonds[0].order == 4);
    auto back = chem::fromChemDraw(chem::toCdxml(d));
    REQUIRE(back);
    REQUIRE(back->bonds.size() == 1);
    CHECK(back->bonds[0].order == 4);
    CHECK(chem::toSmiles(*back) == "C$C");
    CHECK(renderSvg(d).count("<polyline") >= 4);
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
    const QPointF moved = back->atoms[0].pos - doc.atoms[0].pos;  // onto the page (#443)
    CHECK(QLineF(back->arrows[0].from - moved, doc.arrows[0].from).length() < 0.1);
    CHECK(std::abs(back->arrows[1].bend - 10) < 0.1);  // the arc comes back on the same side
    REQUIRE(back->texts.size() == 1);
    CHECK(back->texts[0].text == "L-alanine");
    for (size_t i = 0; i < doc.atoms.size(); ++i) CHECK(QLineF(back->atoms[i].pos - moved, doc.atoms[i].pos).length() < 0.1);

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

TEST_CASE("ChemDraw save and reopen keeps custom colours (#426)") {
    Document d = *chem::fromSmiles("CCO");
    const QColor red(255, 0, 0), blue(0, 0, 255), green(0, 128, 0), orange(255, 128, 0);
    for (auto& a : d.atoms)
        if (a.z == 8) a.color = red;
    d.bonds[0].color = blue;  // C-C
    d.texts.push_back({{0, 60}, "note"});
    d.texts.back().color = green;
    d.arrows.push_back({{60, 0}, {100, 0}});
    d.arrows.back().color = orange;
    d.arrows.push_back({{60, 20}, {100, 20}});  // left as the ink
    for (const QByteArray& file : {chem::toCdxml(d), chem::toCdx(d)}) {
        if (file.isEmpty()) continue;  // a build without binary CDX
        auto back = chem::fromChemDraw(file);
        REQUIRE(back);
        int coloured = 0;
        for (const auto& a : back->atoms) {
            CHECK(a.color == (a.z == 8 ? red : QColor()));
            coloured += a.color.isValid();
        }
        CHECK(coloured == 1);
        int blueBonds = 0;
        for (const auto& b : back->bonds) blueBonds += b.color == blue;
        CHECK(blueBonds == 1);
        REQUIRE(back->texts.size() == 1);
        CHECK(back->texts[0].color == green);
        REQUIRE(back->arrows.size() == 2);
        CHECK(back->arrows[0].color == orange);
        CHECK(!back->arrows[1].color.isValid());
    }
    // ChemDraw takes the table's first two entries as the page and the ink: white, then black, then the colours used.
    CHECK(chem::toCdxml(d).simplified().replace("> <", "><").contains(R"(<colortable><color r="1.0000" g="1.0000" b="1.0000"/><color r="0.0000" g="0.0000" b="0.0000"/><color r="1.0000" g="0.0000" b="0.0000"/>)"));
    // and colours a text run only when the run names a font from the font table.
    CHECK(chem::toCdxml(d).contains(R"(<font id="3")"));
    CHECK(chem::toCdxml(d).contains(R"(<s font="3")"));
    CHECK(!chem::toCdxml(*chem::fromSmiles("CCO")).contains("colortable"));  // nothing coloured, nothing written
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
TEST_CASE("atoms on partial bonds keep the reactants' hydrogens (#655)") {
    // Isocyanide to nitrile through a three-centre TS: Ca...N breaking, Ca...C forming, N#C kept.
    Document ts = *chem::fromSmiles("CC(C)(C)[N+]#[C-]");  // tBu-NC, C5H9N
    const int ca = 1, n = 4, c = 5;
    ts.atoms[n].charge = ts.atoms[c].charge = 0;
    ts.bonds[ts.bondBetween(ca, n)].stereo = BondStereo::Partial;
    ts.bonds.push_back({ca, c});
    ts.bonds.back().stereo = BondStereo::Partial;
    CHECK(chem::properties(ts)->formula == "C5H9N");  // no H invented on Ca or the cyanide C
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

TEST_CASE("labelled charges and radicals share the free upper corner (#721)") {
    for (const QString font : {QString("Arial"), QString("Courier New")})
      for (const QString label : {QString("R"), QString("Ar")})
        for (const auto [charge, radicals] : {std::pair{-1, 0}, std::pair{0, 1}, std::pair{-1, 1}, std::pair{1, 1}})
            for (bool circled : {false, true})
                for (double direction : {-1.0, 1.0}) {
                    CAPTURE(font, label, charge, radicals, circled, direction);
                    Document doc;
                    doc.styleSettings = StyleSettings{};
                    doc.styleSettings->font = font;
                    doc.atoms = {{{0, 0}, 0, charge, label}, {{direction * kBondLength, 0}, 7}};
                    doc.atoms[0].radicals = radicals;
                    doc.bonds = {{0, 1}};
                    doc.circledCharges = circled;
                    auto render = [](const Document& d) {
                        QImage image(800, 800, QImage::Format_ARGB32);
                        image.fill(Qt::white);
                        QPainter painter(&image);
                        painter.translate(400, 400);
                        painter.scale(20, 20);
                        paintDocument(painter, d, {Qt::black, Qt::black, 0.6});
                        return image;
                    };
                    const QImage with = render(doc);
                    doc.atoms[0].charge = doc.atoms[0].radicals = 0;
                    const QImage without = render(doc);
                    QRect marks;
                    std::vector<QPoint> ink;
                    for (int y = 0; y < with.height(); ++y)
                        for (int x = 0; x < with.width(); ++x)
                            if (qGray(with.pixel(x, y)) < 128 && qGray(without.pixel(x, y)) >= 128) {
                                marks |= QRect(x, y, 1, 1);
                                ink.push_back({x, y});
                            }
                    REQUIRE(marks.isValid());
                    CHECK(marks.bottom() < 400);
                    CHECK((marks.center().x() - 400) * direction < 0);
                    int clash = 0;
                    for (QPoint pixel : ink) {
                        const QRect near = QRect(pixel, QSize(1, 1)).adjusted(-8, -8, 8, 8) & without.rect();
                        for (int y = near.top(); y <= near.bottom(); ++y)
                            for (int x = near.left(); x <= near.right(); ++x) clash += qGray(without.pixel(x, y)) < 128;
                    }
                    CHECK(clash == 0);
                }
}

TEST_CASE("a lone labelled radical sits at the label's upper right (#666)") {
    Document radical;
    Atom atom;
    atom.pos = {0, 0};
    atom.z = 0;
    atom.label = "R";
    atom.radicals = 1;
    radical.atoms.push_back(atom);
    Document plain = radical;
    plain.atoms[0].radicals = 0;
    auto render = [](const Document& doc) {
        QImage image(800, 800, QImage::Format_ARGB32);
        image.fill(Qt::white);
        QPainter painter(&image);
        painter.translate(400, 400);
        painter.scale(20, 20);
        paintDocument(painter, doc, {Qt::black, Qt::black, 0.6});
        return image;
    };
    const QImage with = render(radical), without = render(plain);
    QPointF sum;
    int n = 0;
    for (int y = 0; y < with.height(); ++y)
        for (int x = 0; x < with.width(); ++x)
            if (qGray(with.pixel(x, y)) < 128 && qGray(without.pixel(x, y)) >= 128)
                sum += QPointF(x, y), ++n;
    REQUIRE(n > 0);
    const QPointF centre = sum / n - QPointF(with.width(), with.height()) / 2;
    CHECK(centre.x() > std::abs(centre.y()) * 0.4);
    CHECK(centre.y() < 0);
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
    auto box = [](const Arrow& a) { return QRectF(a.from, a.to).normalized(); };
    const QPointF moved = box(cdx->arrows[0]).topLeft() - box(doc.arrows[0]).topLeft();  // onto the page (#443)
    for (size_t k = 0; k < 3; ++k) {
        INFO(k);
        CHECK(cdx->arrows[k].kind == doc.arrows[k].kind);
        CHECK(cdx->arrows[k].dashed == doc.arrows[k].dashed);
        CHECK(box(cdx->arrows[k]) == box(doc.arrows[k]).translated(moved));
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

TEST_CASE("a CDXML export sits inside the page, not at its top-left corner (#443)") {
    const QRectF box = documentBounds(*chem::fromChemDraw(chem::toCdxml(*chem::fromSmiles("CC(=O)Oc1ccccc1C(=O)O"))));
    CHECK(box.center().x() == Catch::Approx(306).margin(1));  // across a US Letter page
    CHECK(box.top() == Catch::Approx(72).margin(1));          // an inch down
    Document wide = *chem::fromSmiles("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC");
    CHECK(documentBounds(*chem::fromChemDraw(chem::toCdxml(wide))).left() == Catch::Approx(72).margin(1));  // wider than the page
}

TEST_CASE("separate ions and molecules from SMILES don't overlap") {
    // Each single-atom fragment's label, drawn alone, keeps clear of the others'.
    for (const char* smiles : {"[Na+].[Cl-]", "O.O.O", "[Li+].[Al+3].[H-].[H-].[H-].[H-]"}) {
        INFO(smiles);
        const auto doc = chem::fromSmiles(smiles);
        REQUIRE(doc);
        std::vector<QRect> boxes;
        for (const Atom& a : doc->atoms) {
            Document one;
            one.atoms = {a};
            QImage img(1600, 400, QImage::Format_ARGB32);
            img.fill(Qt::white);
            QPainter p(&img);
            p.translate(100, 200);
            p.scale(4, 4);
            p.translate(-doc->atoms[0].pos);
            paintDocument(p, one, {Qt::black, Qt::black, 0.6});
            p.end();
            QRect ink;
            for (int y = 0; y < img.height(); ++y)
                for (int x = 0; x < img.width(); ++x)
                    if (qGray(img.pixel(x, y)) < 128) ink |= QRect(x, y, 1, 1);
            for (const QRect& b : boxes) CHECK_FALSE(b.intersects(ink));
            boxes.push_back(ink);
        }
    }
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

TEST_CASE("a charge sits clear of the bonds around its atom (#494)") {
    // The box of the charge's own pixels (ink only with it), grown by 0.4 pt, holds no other ink.
    for (const char* smiles : {"C[N+](C)(C)C", "C[N+](=O)[O-]", "C=[N+]=[N-]", "CC[N+](C)(C)CC", "CC(C)(C)[C+](C)C", "C[n+]1ccccc1"}) {
        INFO(smiles);
        auto doc = chem::fromSmiles(smiles);
        REQUIRE(doc);
        const auto at = std::find_if(doc->atoms.begin(), doc->atoms.end(), [](const Atom& a) { return a.charge > 0; });
        REQUIRE(at != doc->atoms.end());
        const double k = 20;
        auto render = [&](const Document& d) {
            QImage img(800, 800, QImage::Format_ARGB32);
            img.fill(Qt::white);
            QPainter p(&img);
            p.translate(400, 400);
            p.scale(k, k);
            p.translate(-at->pos);
            paintDocument(p, d, {Qt::black, Qt::black, 0.6});
            return img;
        };
        const QImage with = render(*doc);
        Document bare = *doc;
        bare.atoms[at - doc->atoms.begin()].charge = 0;
        const QImage without = render(bare);
        QRect glyph;
        for (int y = 0; y < 800; ++y)
            for (int x = 0; x < 800; ++x)
                if (qGray(with.pixel(x, y)) < 128 && qGray(without.pixel(x, y)) >= 128) glyph |= QRect(x, y, 1, 1);
        REQUIRE(glyph.isValid());
        int clash = 0;
        const QRect near = glyph.adjusted(-8, -8, 8, 8) & with.rect();
        for (int y = near.top(); y <= near.bottom(); ++y)
            for (int x = near.left(); x <= near.right(); ++x) clash += qGray(without.pixel(x, y)) < 128;
        CHECK(clash == 0);
    }
}

TEST_CASE("a shown carbon's charge keeps off its bond (#656)") {
    // Cyanide with its carbon shown and the triple bond to its right: the charge is placed by the
    // same rule as an element's, more than 50 degrees from the bond, not raised after the label.
    auto doc = Document::fromJson(R"({"format":"penzene","version":1,"atoms":[
        {"x":0,"y":0,"z":6,"label":"C","charge":-1},{"x":14.4,"y":0,"z":7}],"bonds":[{"a":0,"b":1,"order":3}]})");
    REQUIRE(doc);
    const double k = 20;
    auto render = [&](const Document& d) {
        QImage img(800, 800, QImage::Format_ARGB32);
        img.fill(Qt::white);
        QPainter p(&img);
        p.translate(400, 400);
        p.scale(k, k);
        paintDocument(p, d, {Qt::black, Qt::black, 0.6});
        return img;
    };
    Document bare = *doc;
    bare.atoms[0].charge = 0, bare.atoms[0].hydrogens = 0;  // the same C, without its charge or an H
    const QImage with = render(*doc), without = render(bare);
    QPointF sum;
    int n = 0;
    for (int y = 0; y < 800; ++y)
        for (int x = 0; x < 800; ++x)
            if (qGray(with.pixel(x, y)) < 128 && qGray(without.pixel(x, y)) >= 128) sum += QPointF(x, y), ++n;
    REQUIRE(n > 0);
    const QPointF d = sum / n - QPointF(400, 400);  // from the carbon to the charge's centre
    CHECK(std::abs(std::atan2(d.y(), d.x())) > qDegreesToRadians(50.0));
}

TEST_CASE("bonds leave a margin around wide element symbols in every drawing style (#617)") {
    for (const auto& st : drawingStyles())
        for (int z : {78, 26, 44, 17, 35, 42})
            for (QPointF toward : {QPointF(1, 0), QPointF(-1, 0), QPointF(0, 1), QPointF(0, -1)}) {
                INFO(st.name.toStdString() << " " << chem::symbol(z) << " toward " << toward.x() << "," << toward.y());
                Document d;
                d.style = st.name;
                d.atoms = {{{0, 0}, z}, {{toward * kBondLength}, 6}};
                d.bonds = {{0, 1}};
                d.bonds[0].color = Qt::blue;
                QFontMetricsF fm(labelFont(st));
                QPainterPath symbol;
                const QString text = QString::fromStdString(chem::symbol(z));
                symbol.addText({-fm.horizontalAdvance(text) / 2, fm.capHeight() / 2}, labelFont(st), text);
                const QRectF box = symbol.boundingRect();
                const double edge = toward.x() > 0 ? box.right() : toward.x() < 0 ? -box.left() :
                                    toward.y() > 0 ? box.bottom() : -box.top();
                QImage img(800, 800, QImage::Format_ARGB32);
                img.fill(Qt::white);
                QPainter p(&img);
                p.translate(400, 400);
                p.scale(20, 20);
                paintDocument(p, d, {Qt::black, Qt::black});
                p.end();
                double first = kBondLength;
                int pixels = 0;
                for (int y = 0; y < img.height(); ++y)
                    for (int x = 0; x < img.width(); ++x) {
                        const QColor c = img.pixelColor(x, y);
                        if (c.blue() > 200 && c.red() < 100 && c.green() < 100) {
                            first = std::min(first, QPointF::dotProduct((QPointF(x, y) - QPointF(400, 400)) / 20, toward));
                            ++pixels;
                        }
                    }
                REQUIRE(pixels > 0);
                CHECK(first >= edge + 1 - st.lineWidth / 2 - 0.1);  // a point of margin, less the round cap and raster rounding
                CHECK(first >= st.labelRadius - st.lineWidth / 2 - 0.1);
            }
}

TEST_CASE("a bond between two labels closer than their clearances isn't drawn over them (#496)") {
    auto render = [](const Document& d) {
        QImage img(400, 400, QImage::Format_ARGB32);
        img.fill(Qt::white);
        QPainter p(&img);
        p.translate(200, 200);
        p.scale(20, 20);
        paintDocument(p, d, {Qt::black, Qt::black, 0.6});
        return img;
    };
    Document d;
    d.atoms = {{{-1.5, 0}, 10}, {{1.5, 0}, 10}};  // neon: no H either way
    const QImage alone = render(d);
    d.bonds.push_back({0, 1});
    CHECK(render(d) == alone);
}

TEST_CASE("an aldehyde's or chain-end alkene's double bond sits toward its one neighbour; a ketone's is centred (#542)") {
    // A zigzag chain from the origin; the last atom is `end`, doubly bonded to the one before it.
    auto chain = [](int atoms, int end, std::vector<QPointF> extra = {}) {
        Document d;
        for (int i = 0; i < atoms; ++i) {
            Atom a;
            a.pos = QPointF(i * 12.47, i % 2 ? -7.2 : 0);
            if (i == atoms - 1) a.z = end;
            d.atoms.push_back(a);
            if (i) d.bonds.push_back({i - 1, i, i == atoms - 1 ? 2 : 1});
        }
        for (QPointF p : extra) {  // more atoms on the double bond's inner end
            Atom a;
            a.pos = p;
            d.atoms.push_back(a);
            d.bonds.push_back({atoms - 2, int(d.atoms.size()) - 1});
        }
        return d;
    };
    auto towardNeighbour = [](const Document& d) {  // the second line on the side of the atom before
        const Bond& b = d.bonds.back();
        const QPointF n = perp(unit(d.atoms[b.b].pos - d.atoms[b.a].pos)) * doubleBondSide(d, b);
        return QPointF::dotProduct(n, d.atoms[0].pos - d.atoms[b.a].pos) > 0;
    };
    CHECK(doubleBondSide(chain(3, 8), chain(3, 8).bonds.back()) != 0);  // CC=O
    CHECK(towardNeighbour(chain(3, 8)));
    CHECK(towardNeighbour(chain(3, 6)));  // CC=C
    CHECK(towardNeighbour(chain(4, 8)));  // CCC=O
    CHECK(doubleBondSide(chain(2, 8), chain(2, 8).bonds.back()) == 0);  // C=O
    const Document ketone = chain(3, 8, {{12.47, -21.6}});  // CC(C)=O
    CHECK(doubleBondSide(ketone, ketone.bonds[1]) == 0);
}

TEST_CASE("a wedge's wide end lies along the bonds beside it, not past them or into a double bond (#509)") {
    const QPointF b(14.4, 0);
    const double k = 20;
    auto render = [&](const Document& d) {
        QImage img(800, 800, QImage::Format_ARGB32);
        img.fill(Qt::white);
        QPainter p(&img);
        p.translate(400, 400);
        p.scale(k, k);
        p.translate(-b);
        paintDocument(p, d, {Qt::black, Qt::black, 0.6});
        return img;
    };
    auto inked = [&](const QImage& img, QPointF at) { return qGray(img.pixel(((at - b) * k + QPointF(400, 400)).toPoint())) < 128; };
    // The wedge ends at b, whose one other bond turns back 60° (a cyclopropane) or 120° (a chain or
    // six-membered ring): no ink across that bond's line, either way along it, on both of the wedge's sides.
    for (double turn : {60.0, 120.0}) {
        INFO(turn);
        const QPointF u(-std::cos(qDegreesToRadians(turn)), std::sin(qDegreesToRadians(turn)));
        Document d;
        d.atoms = {{QPointF(0, 0)}, {b}, {b + u * 14.4}};
        d.bonds = {{0, 1, 1, BondStereo::Wedge}, {1, 2}};
        const QImage img = render(d);
        const QPointF out(u.y(), -u.x());  // across the b–c bond, away from the wedge's narrow end
        REQUIRE(QPointF::dotProduct(out, QPointF(0, 0) - b) < 0);
        int ink = 0;
        for (double t = -4; t < 4; t += 0.1)
            for (double s = 0.5; s < 1.5; s += 0.1) ink += inked(img, b + u * t + out * s);
        CHECK(ink == 0);
    }
    // A carboxyl: the wedge meets the C=O's nearer line, and the space between its two lines stays clear.
    const QPointF u(0.5, std::sqrt(3.0) / 2);
    Document d;
    d.atoms = {{QPointF(0, 0)}, {b}, {b + u * 14.4, 8}, {b + QPointF(u.x(), -u.y()) * 14.4, 8}};
    d.bonds = {{0, 1, 1, BondStereo::Wedge}, {1, 2, 2}, {1, 3}};
    const QImage img = render(d);
    int ink = 0;
    for (double t = 1; t < 4; t += 0.1) ink += inked(img, b + u * t);
    CHECK(ink == 0);
}

TEST_CASE("a charge on an atom with one bond sits above it, not across from the bond") {
    // Straight across from its only bond, a minus reads as another bond (−F–C).
    auto render = [](const Document& d) {
        QImage img(400, 400, QImage::Format_ARGB32);
        img.fill(Qt::white);
        QPainter p(&img);
        p.translate(200, 200);
        p.scale(20, 20);
        paintDocument(p, d, {Qt::black, Qt::black, 0.6});
        return img;
    };
    Document d;
    d.atoms = {{{0, 0}, 9, -1}, {{kBondLength, 0}}};
    d.bonds = {{0, 1}};
    const QImage with = render(d);
    d.atoms[0].charge = 0;
    const QImage without = render(d);
    QRect glyph;
    for (int y = 0; y < 400; ++y)
        for (int x = 0; x < 400; ++x)
            if (qGray(with.pixel(x, y)) < 128 && qGray(without.pixel(x, y)) >= 128) glyph |= QRect(x, y, 1, 1);
    REQUIRE(glyph.isValid());
    CHECK(glyph.bottom() < 200);
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
        f.seek(0);
        const auto sheets = sheetsFromJson(f.readAll());  // and every page, through this version's format
        REQUIRE_FALSE(sheets.empty());
        CHECK(sheetsFromJson(sheetsToJson(sheets)) == sheets);
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
    CHECK(back->texts[0].text == "SO4^2- and NH4+\nH2O");  // as typed: ^ only where it's needed (#371)
}

TEST_CASE("ChemDraw superscripts come in with a ^ only where it's needed (#371)") {
    auto read = [](const char* runs) {
        const QByteArray xml = QByteArray(R"(<?xml version="1.0"?><CDXML><page><t p="100 100">)") + runs +
                               "</t></page></CDXML>";
        auto doc = chem::fromChemDraw(xml);
        return doc && !doc->texts.empty() ? doc->texts[0].text : QString();
    };
    CHECK(read(R"(<s>Cu</s><s face="64">2+</s>)") == "Cu2+");      // raised as typed
    CHECK(read(R"(<s>SO</s><s face="32">4</s><s face="64">2-</s>)") == "SO4^2-");  // 42 would read as a count
    CHECK(read(R"(<s>10</s><s face="64">5</s><s> M</s>)") == "10^5 M");  // plain, it wouldn't be raised
}

TEST_CASE("brackets enclose their atoms' δ marks (#657)") {
    // A three-centre TS in square brackets: R···N, R···C partial, δ+ on R, δ− on C.
    auto doc = Document::fromJson(R"({"format":"penzene","version":1,"atoms":[
        {"x":0,"y":0,"z":0,"label":"R","partial":1},{"x":16,"y":-7.6,"z":7},{"x":16,"y":7.6,"z":6,"label":"C","partial":-1}],
        "bonds":[{"a":0,"b":1,"stereo":"partial"},{"a":1,"b":2,"order":3},{"a":0,"b":2,"stereo":"partial"}],
        "brackets":[{"atoms":[0,1,2]}]})");
    REQUIRE(doc);
    const double k = 10;
    auto render = [&](const Document& d) {
        QImage img(600, 600, QImage::Format_ARGB32);
        img.fill(Qt::white);
        QPainter p(&img);
        p.translate(300, 300);
        p.scale(k, k);
        p.translate(-8, 0);
        paintDocument(p, d, {Qt::black, Qt::black, 0.6});
        return img;
    };
    // Ink only in `a`, not in `b`.
    auto only = [](const QImage& a, const QImage& b) {
        QRect r;
        for (int y = 0; y < a.height(); ++y)
            for (int x = 0; x < a.width(); ++x)
                if (qGray(a.pixel(x, y)) < 128 && qGray(b.pixel(x, y)) >= 128) r |= QRect(x, y, 1, 1);
        return r;
    };
    Document unbracketed = *doc, unmarked = *doc;
    unbracketed.brackets.clear();
    for (Atom& a : unmarked.atoms) a.partial = 0;
    const QImage drawn = render(*doc);
    const QRect brackets = only(drawn, render(unbracketed)), marks = only(render(unbracketed), render([&] {
        Document d = unmarked;
        d.brackets.clear();
        return d;
    }()));
    REQUIRE(brackets.isValid());
    REQUIRE(marks.isValid());
    CHECK(brackets.left() < marks.left());
    CHECK(brackets.right() > marks.right());
    CHECK(brackets.top() < marks.top());
    CHECK(brackets.bottom() > marks.bottom());
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

TEST_CASE(".penz version 2: pages as documents, shared settings at the top (#404)") {
    Document a, b, c;
    a.addAtom({0, 0}), b.addAtom({0, 0}, 8), c.addAtom({0, 0}, 7);
    a.style = b.style = c.style = "RSC";
    a.showStereo = b.showStereo = true;  // not c's: stays on the pages
    c.page = "A4";
    const std::vector<Sheet> sheets{{"One", a}, {"Two", b}, {"Three", c}};
    const QByteArray json = sheetsToJson(sheets);
    const QJsonObject root = QJsonDocument::fromJson(json).object();
    CHECK(root["format"].toString() == "penzene");
    CHECK(root["version"].toInt() == 2);
    CHECK(root["style"].toString() == "RSC");
    CHECK_FALSE(root.contains("showStereo"));
    const QJsonArray pages = root["pages"].toArray();
    REQUIRE(pages.size() == 3);
    CHECK_FALSE(pages[0].toObject().contains("style"));
    CHECK(pages[0].toObject()["showStereo"].toBool());
    CHECK(pages[2].toObject()["name"].toString() == "Three");
    CHECK(sheetsFromJson(json) == sheets);
    CHECK(Document::fromJson(json) == a);  // a single document (pz.read, the command line) is page 1
    CHECK(sheetsFromJson(R"({"format":"penzene","version":2,"pages":[]})").empty());
    CHECK(sheetsFromJson(R"({"format":"penzene","version":3,"pages":[{}]})").empty());
}

TEST_CASE(".penz version 1 files open with every page, and save for Penzene 1 (#404)") {
    QFile f(QString(PENZENE_TEST_DATA) + "/penz/v1.4.0.penz");  // two pages, the #219 layout
    REQUIRE(f.open(QIODevice::ReadOnly));
    const auto sheets = sheetsFromJson(f.readAll());
    REQUIRE(sheets.size() == 2);
    CHECK(sheets[0].name == "Aspirin");
    CHECK(sheets[1].name == "Salicylic acid");
    CHECK(sheets[0].doc.style == "ACS 1996");
    CHECK(sheets[1].doc.atoms.size() == 11);

    const QByteArray v1 = sheetsToJsonV1(sheets);
    const QJsonObject root = QJsonDocument::fromJson(v1).object();
    CHECK(root["version"].toInt() == 1);
    for (const auto& p : root["pages"].toArray()) CHECK(p.toObject()["version"].toInt() == 1);
    CHECK(Document::fromJson(v1) == sheets[0].doc);  // what Penzene before pages reads
    CHECK(sheetsFromJson(v1) == sheets);             // and 1.2 to 1.4, every page
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

// Bond (lower atom, higher atom, order) triples, sorted, to compare Kekulé forms.
static std::vector<std::tuple<int, int, int>> bondOrders(const Document& d) {
    std::vector<std::tuple<int, int, int>> out;
    for (const auto& b : d.bonds) out.push_back({std::min(b.a, b.b), std::max(b.a, b.b), b.order});
    std::sort(out.begin(), out.end());
    return out;
}

TEST_CASE("Clean keeps bond colours and the drawn Kekulé structure (#322)") {
    for (bool flip : {false, true}) {  // either Kekulé form stays as drawn
        Document d = *chem::fromSmiles("c1ccccc1");
        for (auto& b : d.bonds) b.order = flip ? 3 - b.order : b.order, b.color = QColor("#ff0000");
        const Document c = chem::clean2D(d);
        CHECK(bondOrders(c) == bondOrders(d));
        for (const auto& b : c.bonds) CHECK(b.color == QColor("#ff0000"));
    }
}

TEST_CASE("MOL and CDXML imports keep the drawn Kekulé form (#323)") {
    auto orders = bondOrders;
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

TEST_CASE("a racemic (&1) centre keeps its stereo group through .penz, MOL V3000, CDXML and CDX (#389)") {
    auto doc = chem::fromSmiles("C[C@H](N)C(=O)O |&1:1|");  // alanine, racemic, read from CXSMILES
    REQUIRE(doc);
    auto tagged = [](const Document& d) {
        std::vector<Atom> out;
        for (const Atom& a : d.atoms)
            if (a.stereoGroup != StereoGroup::None) out.push_back(a);
        return out;
    };
    auto isAnd1 = [&](const std::optional<Document>& d) {
        REQUIRE(d);
        const auto t = tagged(*d);
        REQUIRE(t.size() == 1);
        CHECK(t[0].stereoGroup == StereoGroup::And);
        CHECK(t[0].stereoGroupNumber == 1);
        CHECK(stereoGroupTag(t[0]) == "&1");
    };
    isAnd1(doc);
    isAnd1(Document::fromJson(doc->toJson()));
    isAnd1(chem::clean2D(*doc));
    const std::string mol = chem::toMolBlock(*doc, true);
    CHECK(mol.find("MDLV30/STERAC1") != std::string::npos);
    isAnd1(chem::fromMolBlock(mol));
    CHECK(chem::toMolBlock(*doc).find("V3000") == std::string::npos);  // V2000 stays V2000, without the group
    const QByteArray cdxml = chem::toCdxml(*doc);
    CHECK(cdxml.contains(R"(EnhancedStereoType="And")"));
    isAnd1(chem::fromChemDraw(cdxml));
    isAnd1(chem::fromChemDraw(chem::toCdx(*doc)));  // binary CDX
}

TEST_CASE("SMILES keeps stereo groups as CXSMILES, and only when there are some (#402)") {
    auto doc = chem::fromSmiles("C[C@H](O)[C@@H](C)N |&1:1,3|");
    REQUIRE(doc);
    const std::string smiles = chem::toSmiles(*doc);
    INFO(smiles);
    CHECK(smiles.find("|&1:") != std::string::npos);
    auto back = chem::fromSmiles(smiles);
    REQUIRE(back);
    int and1 = 0;
    for (const Atom& a : back->atoms) and1 += stereoGroupTag(a) == "&1";
    CHECK(and1 == 2);
    CHECK(chem::toSmiles(*chem::fromSmiles("C[C@H](O)[C@@H](C)N")) == "C[C@H](O)[C@@H](C)N");  // no groups: plain, as before
    chem::Reaction r;
    r.reactants.push_back(*doc);
    r.products.push_back(*doc);
    CHECK(chem::toReactionSmiles(r).find('|') == std::string::npos);  // no extension mid-reaction
}

TEST_CASE("a reaction with an invalid molecule gives no reaction SMILES (#427)") {
    Document pentavalent;  // a carbon with five single bonds
    const int c = pentavalent.addAtom({0, 0});
    for (int i = 0; i < 5; ++i) edit::link(pentavalent, c, pentavalent.addAtom({kBondLength * (i + 1), 0}));
    chem::Reaction r;
    r.reactants.push_back(pentavalent);
    r.reactants.push_back(*chem::fromSmiles("C"));
    r.products.push_back(*chem::fromSmiles("O"));
    CHECK(chem::toReactionSmiles(r).empty());  // not "C>>O"
    CHECK(chem::toReactionSmiles(std::vector{r}).empty());
    r.reactants.erase(r.reactants.begin());
    CHECK(chem::toReactionSmiles(r) == "C>>O");
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

TEST_CASE("typing an element clears a charge (#669)") {
    Document d;
    d.atoms = {{{0, 0}, 8, -1}};
    REQUIRE(edit::applyLabel(d, 0, "O", true));
    CHECK(d.atoms[0].z == 8);
    CHECK(d.atoms[0].charge == 0);
    d.atoms[0].charge = 1;
    REQUIRE(edit::applyLabel(d, 0, "N"));
    CHECK(d.atoms[0].z == 7);
    CHECK(d.atoms[0].charge == 0);
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
    const QPointF moved = back->atoms[0].pos - doc->atoms[0].pos;  // onto the page (#443)
    for (int i = 0; i < 3; ++i) {
        CHECK(back->arrows[i].kind == doc->arrows[i].kind);
        CHECK(back->arrows[i].look == doc->arrows[i].look);
        CHECK(back->arrows[i].behind == doc->arrows[i].behind);
        CHECK(len(back->arrows[i].to - moved - doc->arrows[i].to) < 0.05);
    }
}

TEST_CASE("ChemDraw groups keep whole fragments, arrows and text together (#579)") {
    Document doc = *chem::fromSmiles("CO.CN");
    REQUIRE(edit::applyLabel(doc, 1, "OMe", true));
    for (int i = 0; i < 4; ++i) doc.atoms[i].group = i < 2 ? 7 : 8;
    doc.texts = {{{0, 80}, "first"}, {{80, 80}, "second"}, {{160, 80}, "ungrouped"}};
    doc.texts[0].group = 7;
    doc.texts[1].group = 8;
    doc.arrows = {{{0, -40}, {30, -40}}, {{80, -40}, {110, -40}}};
    doc.arrows[0].group = 7, doc.arrows[0].behind = true;
    doc.arrows[1].group = 8;
    const QByteArray xml = chem::toCdxml(doc);
    CHECK(xml.count("<group ") == 2);
    for (const QByteArray& data : {xml, chem::toCdx(doc)}) {
        const auto back = chem::fromChemDraw(data);
        REQUIRE(back);
        CHECK(chem::toSmiles(*back) == chem::toSmiles(doc));
        REQUIRE(back->atoms.size() == 4);
        REQUIRE(back->arrows.size() == 2);
        REQUIRE(back->texts.size() == 3);
        auto groupOf = [&](const char* text) {
            const auto it = std::find_if(back->texts.begin(), back->texts.end(), [&](const Text& t) { return t.text == text; });
            REQUIRE(it != back->texts.end());
            return it->group;
        };
        const int first = groupOf("first"), second = groupOf("second");
        REQUIRE(first >= 0);
        REQUIRE(second >= 0);
        CHECK(first != second);
        CHECK(groupOf("ungrouped") == -1);
        CHECK(std::count_if(back->atoms.begin(), back->atoms.end(), [&](const Atom& a) { return a.group == first; }) == 2);
        CHECK(std::count_if(back->atoms.begin(), back->atoms.end(), [&](const Atom& a) { return a.group == second; }) == 2);
        CHECK(back->arrows[0].group == first);
        CHECK(back->arrows[0].behind);
        CHECK(back->arrows[1].group == second);
        CHECK_FALSE(back->arrows[1].behind);
    }
    doc.atoms[1].group = -1;  // grouping half a molecule cannot split its ChemDraw fragment
    const auto partial = chem::fromChemDraw(chem::toCdxml(doc));
    REQUIRE(partial);
    CHECK(chem::toSmiles(*partial) == chem::toSmiles(doc));
    CHECK(std::count_if(partial->atoms.begin(), partial->atoms.end(), [](const Atom& a) { return a.group < 0; }) == 2);
}

TEST_CASE("ChemDraw nested groups flatten while inline fragment groups are ignored (#579)") {
    const auto doc = chem::fromChemDraw(R"(<CDXML BondLength="30"><page id="1">
<group id="2"><fragment id="3"><n id="4" p="0 0"/><n id="5" p="30 0"/><b id="6" B="4" E="5"/></fragment>
<group id="7"><t id="8" p="0 30"><s>nested</s></t></group></group>
<fragment id="9"><group id="10"><n id="11" p="90 0"/></group><n id="12" p="120 0"/><b id="13" B="11" E="12"/></fragment>
</page></CDXML>)");
    REQUIRE(doc);
    REQUIRE(doc->atoms.size() == 4);
    REQUIRE(doc->texts.size() == 1);
    CHECK(doc->atoms[0].group == doc->texts[0].group);
    CHECK(doc->atoms[1].group == doc->texts[0].group);
    CHECK(doc->texts[0].group >= 0);
    CHECK(doc->atoms[2].group == -1);
    CHECK(doc->atoms[3].group == -1);
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

TEST_CASE(".penz reads an isotope lighter than its element as none (#368)") {
    auto doc = Document::fromJson(R"({"format":"penzene","version":1,"atoms":[{"x":0,"y":0,"z":6,"isotope":1},{"x":0,"y":0,"z":6,"isotope":13}]})");
    REQUIRE(doc);
    CHECK(doc->atoms[0].isotope == 0);
    CHECK(doc->atoms[1].isotope == 13);
}

TEST_CASE("a charged abbreviation keeps the group: N3- is azide, not nitride (#370)") {
    Document d;
    d.addAtom({0, 0});
    REQUIRE(edit::applyLabel(d, 0, "N3-", true));
    CHECK(d.atoms[0].label == "N3");
    CHECK(d.atoms[0].charge == -1);
    auto p = chem::properties(d);
    REQUIRE(p);
    CHECK(p->formula == "N3-");  // the azide anion
    Document fe;
    fe.addAtom({0, 0});
    REQUIRE(edit::applyLabel(fe, 0, "Fe3+", true));  // an element's digits are still its charge
    CHECK(fe.atoms[0].charge == 3);
}

TEST_CASE("an atom labelled with its own element keeps its charge and isotope (#654)") {
    auto d = Document::fromJson(R"({"format":"penzene","version":1,"atoms":[
        {"x":0,"y":0,"z":6,"label":"C","charge":-1,"isotope":13},{"x":14.4,"y":0,"z":7}],
        "bonds":[{"a":0,"b":1,"order":3}]})");
    REQUIRE(d);
    CHECK(chem::toSmiles(*d) == "[13C-]#N");
    CHECK(chem::properties(*d)->formula == "[13C]N-");
}

TEST_CASE("a lone abbreviation comes back from ChemDraw as a labelled atom (#384)") {
    for (const char* label : {"Boc", "N3-"}) {
        Document d;
        d.addAtom({0, 0});
        REQUIRE(edit::applyLabel(d, 0, label, true));
        auto back = chem::fromChemDraw(chem::toCdxml(d));
        REQUIRE(back);
        INFO(label);
        REQUIRE(back->atoms.size() == 1);
        CHECK(back->texts.empty());
        CHECK(back->atoms[0].label == d.atoms[0].label);
        CHECK(chem::properties(*back)->formula == chem::properties(d)->formula);
    }
}

TEST_CASE("HOSE codes match CDK's reference codes (#403)") {
    auto codesOf = [](const char* smiles) {
        std::vector<std::string> out;
        for (const auto& c : chem::hoseCodes(*chem::fromSmiles(smiles), 4)) out.push_back(c.empty() ? "" : c.back());
        return out;
    };
    // CDK's HOSECodeGeneratorTest, 4 spheres: a chain (Br drawn as Y), and indole's rings (& closes one).
    CHECK(codesOf("CC=CBr") == std::vector<std::string>{"C-4;C(=C/Y/)", "C-3;=CC(Y,//)", "C-3;=CY(C,//)", "Br-1;C(=C/C/)"});
    CHECK(codesOf("C1(C=CN2)=C2C=CC=C1") ==
          std::vector<std::string>{"C-3;*C*C*C(*C*N,*C,*C/*C,*&,*&,*&/*&)", "C-3;*C*C(*C*C,*N/*C*&,*C,*&/*C,*&)",
                                   "C-3;*C*N(*C,*C/*C*&,*C*&/*C,*C)", "N-3;*C*C(*C*C,*C/*C*&,*C,*&/*C,*&)",
                                   "C-3;*C*C*N(*C*C,*C,*C/*C,*&,*&,*&/*&)", "C-3;*C*C(*C*N,*C/*C*C,*C,*&/*&,*&,*&)",
                                   "C-3;*C*C(*C,*C/*C*N,*&/*C*&,*C)", "C-3;*C*C(*C,*C/*C*C,*&/*N*&,*C)",
                                   "C-3;*C*C(*C*C,*C/*C*N,*C,*&/*&,*&,*&)"});
    CHECK(codesOf("CCO")[0] == "C-4;C(O//)");
    CHECK(chem::hoseCodes(*chem::fromSmiles("CCO"), 1)[0][0] == "C-4;C(//)");
}

TEST_CASE("HOSE codes are the same from a drawing and from its MOL block, H drawn or not (#403)") {
    for (const char* smiles : {"c1ccccc1O", "C[N+](C)(C)CC(=O)[O-]", "CC(=O)Oc1ccccc1C(=O)O", "[H]OC([H])([H])C", "c1cc[nH]c1", "c1ccc2[nH]ccc2c1"}) {
        auto doc = chem::fromSmiles(smiles);
        REQUIRE(doc);
        const auto drawn = chem::hoseCodes(*doc);
        CHECK(drawn == chem::hoseCodes(chem::toMolBlock(*doc)));
        const auto withH = chem::hoseCodes(chem::addHydrogens(*doc));
        REQUIRE(withH.size() >= drawn.size());
        for (size_t i = 0; i < drawn.size(); ++i)
            if (doc->atoms[i].z != 1) CHECK(withH[i] == drawn[i]);
    }
    Document boc = *chem::fromSmiles("CN");  // an abbreviation: codes for the drawn atoms, of the expanded structure
    boc.atoms[1].label = "NHBoc";
    const auto codes = chem::hoseCodes(boc);
    REQUIRE(codes.size() == 2);
    CHECK(codes[0][0] == "C-4;N(//)");
    CHECK(codes[0][1].starts_with("C-4;N(C/"));
}

TEST_CASE("predicted 13C and 1H shifts: benzene, ethanol (#403)") {
    auto benzene = chem::predictShifts(*chem::fromSmiles("c1ccccc1"));
    REQUIRE(benzene.size() == 6);
    for (const auto& s : benzene) {
        CHECK(s.carbon == Catch::Approx(128.5).margin(2));
        CHECK(s.carbonSpheres == 5);
        CHECK(s.proton == Catch::Approx(7.3).margin(0.3));
    }
    auto ethanol = chem::predictShifts(*chem::fromSmiles("CCO"));
    REQUIRE(ethanol.size() == 3);
    CHECK(ethanol[0].carbon == Catch::Approx(18).margin(3));  // CH3
    CHECK(ethanol[1].carbon == Catch::Approx(58).margin(3));  // CH2
    CHECK(ethanol[0].proton == Catch::Approx(1.2).margin(0.3));
    CHECK(ethanol[1].proton == Catch::Approx(3.7).margin(0.3));
    CHECK(ethanol[2].carbonSpheres == 0);  // O: 1H only
    CHECK(ethanol[2].protonSpheres > 0);
}

TEST_CASE("a fifth HOSE sphere tells a primary aniline's ring from N in a ring (#682)") {
    // The 4-sphere code of aniline's para carbon is shared with ring N, which pulls its median to
    // 125 ppm; the fifth sphere (the terminal N's empty branch) gives 117.4, near the 118.5 measured.
    const auto shifts = chem::predictShifts(*chem::fromSmiles("Nc1ccccc1"));
    const auto para = std::find_if(shifts.begin(), shifts.end(), [](const chem::Shift& s) { return s.atom == 4; });
    REQUIRE(para != shifts.end());
    CHECK(para->carbonSpheres == 5);
    CHECK(para->carbon == Catch::Approx(117.4).margin(1.5));
}

TEST_CASE("predicted spectra: one stick per set of equivalent atoms, for the chosen molecule (#444)") {
    const Document d = *chem::fromSmiles("CCO.Cc1ccccc1");  // ethanol, then toluene (atoms 3-9)
    auto counts = [](const std::vector<chem::NmrStick>& sticks) {
        std::vector<int> out;
        for (const auto& k : sticks) out.push_back(k.count);
        return out;
    };
    const std::vector<int> toluene{3, 4, 5, 6, 7, 8, 9};
    const auto carbon = chem::nmrSticks(d, false, toluene);
    CHECK(counts(carbon) == std::vector<int>{1, 2, 2, 1, 1});  // ipso, ortho, meta, para, CH3 (highest ppm first)
    for (size_t k = 1; k < carbon.size(); ++k) CHECK(carbon[k - 1].ppm >= carbon[k].ppm);
    for (const auto& k : carbon)
        for (int a : k.atoms) CHECK(a >= 3);  // the drawing's own indices, so a stick can light its atoms
    CHECK(chem::nmrSticks(d, false).size() == carbon.size() + 2);  // the whole page: no stick shared across molecules
    const auto proton = chem::nmrSticks(d, true, toluene);
    int hydrogens = 0;
    for (const auto& k : proton) hydrogens += k.count;
    CHECK(hydrogens == 8);
    CHECK(proton.back().count == 3);  // the CH3, furthest upfield
    auto ethanol = counts(chem::nmrSticks(d, true, {0, 1, 2}));
    std::sort(ethanol.begin(), ethanol.end());
    CHECK(ethanol == std::vector<int>{1, 2, 3});  // OH, CH2, CH3
    QStringList split;  // first order: CH3 by CH2 a triplet, CH2 by CH3 a quartet (OH exchanges), OH a singlet
    for (const auto& k : chem::nmrSticks(d, true, {0, 1, 2})) split << QString::number(k.count) + k.multiplicity();
    split.sort();
    CHECK(split == QStringList{"1s", "2q", "3t"});
    CHECK(chem::nmrSticks(d, true, {3}).at(0).multiplicity() == "s");  // toluene's CH3: no H next door
    CHECK(chem::nmrSticks(*chem::fromSmiles("C1CCCCC1"), true).at(0).multiplicity() == "s");  // equivalent H don't split each other
    const QString h = chem::nmrLine(d, true, {0, 1, 2});  // for the SI (#566)
    CHECK(h.startsWith("1H NMR (predicted) δ "));
    CHECK(h.contains(QRegularExpression(R"(\d\.\d\d \(q, 2H\), .*\d\.\d\d \(t, 3H\)\.$)")));
    CHECK(QRegularExpression(R"(^13C NMR \(predicted\) δ (\d+\.\d, ){4}\d+\.\d\.$)").match(chem::nmrLine(d, false, {3, 4, 5, 6, 7, 8, 9})).hasMatch());
    CHECK(chem::nmrLine(*chem::fromSmiles("[Na+].[Cl-]"), true).isEmpty());
    CHECK(chem::nmrSticks(d, true).size() == proton.size() + 3);  // everything
}

TEST_CASE("higher NMR multiplicities require one equivalent coupling set (#574)") {
    for (auto [smiles, atom, expected] : {std::tuple{"CCCCC", 2, "quint"}, {"CC(O)C", 1, "sept"},
                                        {"CCCO", 1, "m"}}) {
        INFO(smiles);
        const Document d = *chem::fromSmiles(smiles);
        const auto sticks = chem::nmrSticks(d, true, {atom});
        REQUIRE(sticks.size() == 1);
        CHECK(sticks[0].multiplicity() == expected);
        CHECK(chem::nmrLine(d, true, {atom}).contains(QString("(%1,").arg(expected)));
        const auto withH = chem::addHydrogens(d);
        CHECK(chem::nmrSticks(withH, true, {atom})[0].multiplicity() == expected);
    }
    chem::NmrStick stick;
    stick.singleCouplingSet = true;
    stick.coupled = 5;
    CHECK(stick.multiplicity() == "sext");
    stick.singleCouplingSet = false;
    CHECK(stick.multiplicity() == "m");
    stick.singleCouplingSet = true;
    stick.coupled = 7;
    CHECK(stick.multiplicity() == "m");
}

TEST_CASE("a spectrum counts the atoms inside an abbreviation, the same as drawn out (#611)") {
    Document d = *chem::fromSmiles("CO");
    REQUIRE(edit::applyLabel(d, 0, "tBu", true));
    const Document full = chem::expanded(d);
    for (bool proton : {true, false}) {
        INFO(chem::nmrLine(d, proton).toStdString());
        CHECK(chem::nmrLine(d, proton) == chem::nmrLine(full, proton));
    }
    CHECK(chem::nmrLine(d, true).contains("9H"));
    const auto methyls = chem::nmrSticks(d, false, {0});  // the label's atom stands for its hidden ones
    REQUIRE(methyls.size() == 2);
    for (const auto& k : methyls) CHECK(k.atoms == std::vector<int>{0});
    CHECK(chem::predictShifts(d).size() == chem::predictShifts(d).back().atom + 1u);  // per drawn atom only
}

TEST_CASE("predicted shifts follow the bonds, not where the atoms are drawn") {
    auto doc = *chem::fromSmiles("CCO");
    const auto before = chem::predictShifts(doc);
    doc.atoms[1].pos += QPointF(30, -12);  // dragged
    const auto moved = chem::predictShifts(doc);
    REQUIRE(moved.size() == before.size());
    for (size_t i = 0; i < moved.size(); ++i) {
        CHECK(moved[i].carbon == before[i].carbon);
        CHECK(moved[i].proton == before[i].proton);
    }
    doc.atoms[2].z = 7;  // O -> N: a different molecule, a different answer
    CHECK(chem::predictShifts(doc)[1].carbon != before[1].carbon);
}

TEST_CASE("compound numbers keep scheme order, follow their molecules and keep series (#504)") {
    const double L = kBondLength;
    Document d;
    auto add = [&](const char* smiles, QPointF at) {  // a molecule and its number, anchored to its first atom
        const int first = int(d.atoms.size());
        d.append(*chem::fromSmiles(smiles), at);
        d.texts.push_back({at + QPointF(0, 2 * L), "", 1, {}, true, first});
    };
    add("CCO", {20 * L, 0});  // added out of order
    add("c1ccccc1", {0, 0});
    add("CC", {10 * L, 0});
    add("CN", {0, 10 * L});  // the next row
    auto numbers = [&] {
        QStringList out;
        for (const Text& t : d.texts) out << t.text;
        return out;
    };
    edit::renumberCompounds(d);
    CHECK(numbers() == QStringList{"3", "1", "2", "4"});
    Document again = d;
    edit::renumberCompounds(again);
    CHECK(again == d);  // no change, so no phantom undo step
    add("C", {5 * L, 0});  // a step inserted: the later numbers move up
    edit::renumberCompounds(d);
    CHECK(numbers() == QStringList{"4", "1", "3", "5", "2"});

    SECTION("a tall molecule's number sits lower but is in its row") {
        const int first = int(d.atoms.size());
        d.append(*chem::fromSmiles("CCCCCCC"), {-10 * L, 0});
        for (int i = first; i < int(d.atoms.size()); ++i) d.atoms[i].pos = {-10 * L, (i - first - 3) * L};  // upright, centred on the row
        d.texts.push_back({{-10 * L, 5 * L}, "", 1, {}, true, first});
        edit::renumberCompounds(d);
        CHECK(numbers() == QStringList{"5", "2", "4", "6", "3", "1"});
    }
    SECTION("2a and 2b keep one number between them") {
        d.texts[1].text = "7a", d.texts[4].text = "7b";
        edit::renumberCompounds(d);
        CHECK(numbers() == QStringList{"3", "1a", "2", "4", "1b"});
    }
    SECTION("a number follows its molecule, unless it was moved itself") {
        Document moved = d;
        for (int i : edit::moleculeOf(d, d.texts[0].anchor)) moved.atoms[i].pos += QPointF(-30 * L, L);
        for (int i : edit::moleculeOf(d, d.texts[2].anchor)) moved.atoms[i].pos += QPointF(0, L);  // dragged with its number
        moved.texts[2].pos += QPointF(0, L);
        Document after = moved;
        edit::followNumbers(d, after);
        CHECK(after.texts[0].pos == d.texts[0].pos + QPointF(-30 * L, L));
        CHECK(after.texts[2].pos == moved.texts[2].pos);
        CHECK(after.texts[1].pos == d.texts[1].pos);
        edit::renumberCompounds(after);
        CHECK(after.texts[0].text == "1");  // now first in the scheme
    }
    SECTION("anchors are renumbered with the atoms, and saved") {
        const int benzene = d.texts[1].anchor;
        d.removeAtoms({0, 1, 2});  // ethanol: its number stays, free
        CHECK(d.texts[0].anchor == -1);
        CHECK(d.texts[1].anchor == benzene - 3);
        Document cut = d;
        cut.removeAtoms({d.texts[1].anchor});  // just its atom: it takes a neighbour of the same molecule
        REQUIRE(cut.texts[1].anchor >= 0);
        CHECK(edit::moleculeOf(cut, cut.texts[1].anchor).size() == 5);
        Document two;
        two.append(d, {});
        two.append(d, {});
        CHECK(two.texts[6].anchor == d.texts[1].anchor + int(d.atoms.size()));
        const auto back = Document::fromJson(d.toJson());
        REQUIRE(back);
        CHECK(back->texts == d.texts);
        CHECK(chem::toCdxml(d).contains(R"(face="1")"));  // bold in ChemDraw too
    }
}

TEST_CASE("a spectrum's legend takes the emptier top corner and the sticks keep clear of it (#444)") {
    const QRectF plot(0, 0, 400, 200);
    const QSizeF size(100, 80);
    auto clear = [&](const Legend& l, const std::vector<QPointF>& sticks) {
        for (QPointF s : sticks)
            if (s.x() > l.rect.left() && s.x() < l.rect.right()) CHECK(plot.bottom() - s.y() * l.scale * plot.height() >= l.rect.bottom() + 10 - 1e-9);
    };
    const std::vector<QPointF> right{{350, 1}, {380, 0.5}};  // tall sticks on the right: the legend goes left, nothing shrinks
    Legend l = placeLegend(plot, size, right, 10);
    CHECK(l.rect.left() == 0);
    CHECK(l.scale == 1);
    const std::vector<QPointF> both{{20, 1}, {350, 0.8}};  // under either corner: the one that shrinks them less
    l = placeLegend(plot, size, both, 10);
    CHECK(l.rect.right() == 400);
    CHECK(l.scale == Catch::Approx((200 - 80 - 10) / (0.8 * 200)));
    clear(l, both);
}

TEST_CASE(".penz schema describes every field the app writes (#573)") {
    QFile f(QString(PENZENE_TEST_DATA) + "/../../docs/_static/penz.schema.json");
    REQUIRE(f.open(QIODevice::ReadOnly));
    const QJsonObject defs = QJsonDocument::fromJson(f.readAll()).object()["$defs"].toObject();
    auto properties = [&](const QString& def) { return defs[def].toObject()["properties"].toObject(); };
    auto described = [](const QJsonObject& written, const QJsonObject& schema) {
        for (const QString& key : written.keys()) {
            INFO(key.toStdString());
            CHECK(schema.contains(key));
        }
    };
    Document d = *chem::fromSmiles("C[C@H](N)C(=O)O");  // every feature, so every field is written
    Atom& a = d.atoms[1];
    a.charge = 1, a.label = "Me", a.color = Qt::red, a.map = 1, a.lonePairs = 1, a.radicals = 1, a.partial = 1;
    a.isotope = 13, a.stereoGroup = StereoGroup::And, a.stereoGroupNumber = 1, a.attachments = {0, 2}, a.group = 0, a.standsFor = "N, O";
    d.bonds[0].stereo = BondStereo::Wedge, d.bonds[0].position = BondPosition::Left, d.bonds[0].color = Qt::red;
    Arrow arrow;
    arrow.to = {50, 0}, arrow.bend = 5, arrow.color = Qt::red, arrow.dashed = arrow.behind = arrow.crossed = true;
    arrow.look = OrbitalLook::Shaded, arrow.head = 2, arrow.fromAt = {0, -1}, arrow.toAt = {1, 2}, arrow.group = 0;
    arrow.filled = true, arrow.fill = Qt::blue;
    d.arrows.push_back(arrow);
    Text text{{0, 40}, "1"};
    text.scale = 2, text.color = Qt::red, text.compound = text.legend = true, text.anchor = 0, text.group = 0;
    text.font = "Times New Roman", text.bold = text.italic = true;
    d.texts.push_back(text);
    d.fills.push_back({{0, 1, 2}, Qt::red});
    d.brackets.push_back({{0, 1}, false, "n"});
    d.style = "ACS 1996", d.carbonLabels = Document::CarbonLabels::All, d.labelRatio = 0.5;
    d.showStereo = d.showAtomNumbers = d.showShifts = d.aromaticCircles = true;
    d.page = "A4", d.aromaticCircleOverrides = {{0, 1, 2}};
    const QJsonObject root = QJsonDocument::fromJson(d.toJson()).object();
    QJsonObject page = properties("page");
    page["format"] = page["version"] = true;
    described(root, page);
    for (const auto& [key, def] : {std::pair{"atoms", "atom"}, {"bonds", "bond"}, {"arrows", "arrow"}, {"texts", "text"},
                                   {"fills", "fill"}, {"brackets", "bracket"}}) {
        REQUIRE(!root[key].toArray().isEmpty());
        for (const auto& o : root[key].toArray()) described(o.toObject(), properties(def));
    }
    described(root["page"].toObject(), page["page"].toObject()["properties"].toObject());
    QJsonObject v2 = properties("v2");
    v2["format"] = true;
    described(QJsonDocument::fromJson(sheetsToJson({{"A", d}, {"B", d}})).object(), v2);
}

TEST_CASE("a molecule's predicted shifts don't depend on the others on the page (#575)") {
    const Document ethanol = *chem::fromSmiles("CCO");
    Document d;
    const int c = d.addAtom({100, 0});  // first, a five-bonded carbon: no valid molecule
    for (int k = 0; k < 5; ++k) d.bonds.push_back({c, d.addAtom({100.0 + k, 10})});
    d.append(ethanol);  // atoms 6-8
    const auto alone = chem::nmrSticks(ethanol, false);
    const auto sticks = chem::nmrSticks(d, false, {6, 7, 8});
    REQUIRE(sticks.size() == 2);
    REQUIRE(sticks.size() == alone.size());
    for (size_t k = 0; k < sticks.size(); ++k) CHECK(sticks[k].ppm == alone[k].ppm);
    CHECK(chem::predictShifts(d).size() == chem::predictShifts(ethanol).size());  // the canvas's shift labels too
}

TEST_CASE("alternatives typed as one label make one variable atom, not a chain (#587)") {
    Document d = *chem::fromSmiles("c1ccccc1CCC");
    const size_t atoms = d.atoms.size(), bonds = d.bonds.size();
    REQUIRE(edit::applyLabel(d, 6, "N,O,S", true));
    REQUIRE(edit::applyLabel(d, 7, "N/O/S", true));
    REQUIRE(edit::applyLabel(d, 8, "[Cl, Br]", true));
    CHECK(d.atoms.size() == atoms);
    CHECK(d.bonds.size() == bonds);
    CHECK(d.atoms[6].label == "X");
    CHECK(d.atoms[7].label == "Y");  // the next free name
    CHECK(d.atoms[8].label == "Z");
    CHECK(d.atoms[6].standsFor == "N, O, S");
    CHECK(d.atoms[7].standsFor == "N, O, S");
    CHECK(d.atoms[8].standsFor == "Cl, Br");
    CHECK(d.atoms[6].z == 0);
    REQUIRE(edit::applyLabel(d, 6, "O,S", true));
    CHECK(d.atoms[6].label == "X");  // keeps its own name
    CHECK(d.atoms[6].standsFor == "O, S");
    REQUIRE(edit::applyLabel(d, 6, "OMe", true));
    CHECK(d.atoms[6].standsFor.isEmpty());  // a new label drops the old definition
    REQUIRE(edit::applyLabel(d, 7, "Y", true));
    CHECK(d.atoms[7].standsFor == "N, O, S");  // the same label keeps it
    Document strict = *chem::fromSmiles("CC");
    CHECK_FALSE(strict.atoms[1].label == "X");
    edit::applyLabel(strict, 1, "N,Q", true);  // Q isn't an element: not alternatives
    CHECK(strict.atoms[1].standsFor.isEmpty());
}

TEST_CASE("a variable label stands for elements or groups: a legend, and an atom list in MOL and SMARTS (#505)") {
    Document d = *chem::fromSmiles("Clc1ccc(cc1)C(=O)NCCO");
    REQUIRE(edit::applyLabel(d, 0, "R", true));
    REQUIRE(edit::applyLabel(d, 9, "X", true));
    const std::string plain = chem::toSmarts(d);
    CHECK(plain.find("#9") == std::string::npos);  // a bare X is a variable, not "any halogen"
    d.atoms[0].standsFor = "H, Me, OMe";  // groups: a note only
    d.atoms[9].standsFor = "N, O, S";     // elements: a query atom list
    edit::syncLegends(d);
    auto legends = [](const Document& doc) {
        QStringList out;
        for (const Text& t : doc.texts)
            if (t.legend) out << t.text;
        return out;
    };
    CHECK(legends(d) == QStringList{"R = H, Me, OMe\nX = N, O, S"});
    CHECK(d.texts[0].anchor >= 0);
    CHECK(*Document::fromJson(d.toJson()) == d);

    const std::string smarts = chem::toSmarts(d);
    INFO(smarts);
    CHECK(smarts.find("#7,#8,#16") != std::string::npos);
    for (bool v3000 : {false, true}) {
        auto back = chem::fromMolBlock(chem::toMolBlock(d, v3000));
        REQUIRE(back);
        CHECK(back->atoms[9].label == "X");  // the alias, or the first free of X, Y... (V3000)
        CHECK(back->atoms[9].standsFor == "N, O, S");
        CHECK(back->atoms[9].z == 0);
        edit::syncLegends(*back);
        CHECK(legends(*back).join("").contains("X = N, O, S"));
    }
    std::vector<QString> problems;
    for (const auto& p : chem::checkStructure(d)) problems.push_back(p.message);
    CHECK(std::find(problems.begin(), problems.end(), "Query atom X: one of N, O, S") != problems.end());
    CHECK(std::none_of(problems.begin(), problems.end(), [](const QString& p) { return p.startsWith("Unknown label"); }));

    d.atoms[0].standsFor.clear(), d.atoms[9].standsFor.clear();
    edit::syncLegends(d);
    CHECK(legends(d).isEmpty());  // gone with the last definition
    REQUIRE(edit::applyLabel(d, 12, "Q", true));  // MDL's heteroatom code, typed
    CHECK(chem::toSmarts(d).find("!#6") != std::string::npos);
    CHECK(chem::fromMolBlock(chem::toMolBlock(d))->atoms[12].label == "Q");
}

TEST_CASE("SMILES wildcards [1*] and [*:1] come in as R sites, keeping other map numbers (#432)") {
    auto site = [](const std::string& smiles) {
        auto doc = chem::fromSmiles(smiles);
        REQUIRE(doc);
        for (const auto& a : doc->atoms)
            if (a.z == 0) return std::pair{a.label.toStdString(), a.map};
        return std::pair{std::string(), -1};
    };
    CHECK(site("[1*]c1ccccc1") == std::pair{std::string("R1"), 0});
    CHECK(site("[*:2]c1ccccc1") == std::pair{std::string("R2"), 0});   // the map number was the site's
    CHECK(site("[1*:5]c1ccccc1") == std::pair{std::string("R1"), 5});  // here it's a real atom map
    CHECK(site("*c1ccccc1") == std::pair{std::string(), 0});           // a bare * is no site
    Document d = *chem::fromSmiles("Cc1ccccc1");
    REQUIRE(chem::attach(d, 0, "*OC"));  // the atom bonded to * attaches: methoxy, not a C-O-
    CHECK(chem::toSmiles(d) == chem::toSmiles(*chem::fromSmiles("COc1ccccc1")));
    CHECK_FALSE(chem::attach(d, 0, "*C*"));
}

TEST_CASE("a negated atom list stays negated through MOL and SMARTS (#610)") {
    const auto d = chem::fromMolBlock("Negated atom list\n     RDKit          2D\n\n"
                                      "  0  0  0  0  0  0  0  0  0  0999 V3000\n"
                                      "M  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 0 0 0\nM  V30 BEGIN ATOM\n"
                                      "M  V30 1 \"NOT [N,O]\" 0.000000 0.000000 0.000000 0\n"
                                      "M  V30 2 C 1.299038 0.750000 0.000000 0\nM  V30 END ATOM\n"
                                      "M  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n");
    REQUIRE(d);
    CHECK(d->atoms[0].standsFor == "not N, O");
    const std::string smarts = chem::toSmarts(*d);
    INFO(smarts);
    CHECK(smarts.find('!') != std::string::npos);
    for (bool v3000 : {false, true}) {  // V2000 writes it as a negated M  ALS
        auto back = chem::fromMolBlock(chem::toMolBlock(*d, v3000));
        REQUIRE(back);
        CHECK(back->atoms[0].standsFor == "not N, O");
    }
    std::vector<QString> problems;
    for (const auto& p : chem::checkStructure(*d)) problems.push_back(p.message);
    CHECK(std::find(problems.begin(), problems.end(), "Query atom X: none of N, O") != problems.end());
}

TEST_CASE("eta complexes have chemical formula and mass without drawing centroids (#619)") {
    Document doc;
    const int fe = doc.addAtom({0, 0}, 26);
    edit::hotkey(doc, {fe, -1}, "j");
    edit::hotkey(doc, {fe, -1}, "j");
    const QByteArray drawn = doc.toJson();
    const auto p = chem::properties(doc);
    REQUIRE(p);
    CHECK(p->formula == "C10H10Fe");
    CHECK(p->exactMass == Catch::Approx(186.0132).margin(0.0001));
    CHECK(chem::hrmsLine(doc, chem::Ion::M).contains("C10H10Fe 186.0126"));
    CHECK_FALSE(chem::isotopePattern(doc, chem::Ion::M).empty());
    CHECK(doc.toJson() == drawn);

    const auto ringIds = chem::rings(doc);
    REQUIRE(ringIds.size() == 2);
    Document clean = chem::clean2D(doc);
    CHECK(clean.atoms.size() == doc.atoms.size());
    CHECK(clean.bonds.size() == doc.bonds.size());
    for (const auto& ring : ringIds) {
        QPointF before, after;
        for (int i : ring) before += doc.atoms[i].pos, after += clean.atoms[i].pos;
        before /= double(ring.size()), after /= double(ring.size());
        auto dummy = std::find_if(doc.atoms.begin(), doc.atoms.end(), [&](const Atom& a) {
            return a.z == 0 && len(a.pos - before) < 0.5;
        });
        REQUIRE(dummy != doc.atoms.end());
        const int i = int(dummy - doc.atoms.begin());
        CHECK(len(clean.atoms[i].pos - after) < 0.5);
        CHECK(clean.atoms[i].z == 0);
        CHECK(clean.atoms[i].charge == 0);
    }
    REQUIRE(chem::properties(clean));
    CHECK(chem::properties(clean)->formula == "C10H10Fe");
    CHECK_FALSE(chem::hrmsLine(clean, chem::Ion::M).isEmpty());
}

TEST_CASE("eta recognition keeps real wildcard atoms and original drawing indices (#619)") {
    SECTION("mixed Cp and benzene on a metal: only Cp changes the oxidation state") {
        Document doc;
        doc.addAtom({0, 0}, 26);
        edit::hotkey(doc, {0, -1}, "j");
        edit::hotkey(doc, {0, -1}, "J");
        doc = *Document::fromJson(doc.toJson());  // existing .penz drawing convention, no new metadata
        REQUIRE(chem::properties(doc));
        CHECK(chem::properties(doc)->formula == "C11H11Fe");
        CHECK_FALSE(chem::hrmsLine(doc, chem::Ion::M).isEmpty());
        auto exported = chem::fromSmiles(chem::toSmiles(doc));
        REQUIRE(exported);
        CHECK(chem::properties(*exported)->formula == "C11H11Fe");
        CHECK(doc.atoms[0].charge == 0);
        CHECK(chem::rings(doc).size() == 2);
        CHECK(chem::atomInfo(doc).size() == doc.atoms.size());
        const int base = int(doc.atoms.size());
        doc.append(*chem::fromSmiles("F[C@H](Cl)Br"), {100, 0});
        const auto stereo = chem::stereoLabels(doc);
        REQUIRE(stereo.size() == 1);
        CHECK(stereo[0].atom == base + 1);
        const Document withH = chem::addHydrogens(doc);
        REQUIRE(chem::properties(withH));
        CHECK(chem::properties(withH)->formula == chem::properties(doc)->formula);
        CHECK(chem::checkStructure(doc).empty());
    }
    SECTION("a plain attachment point away from a ring is not an eta centroid") {
        Document doc;
        doc.addAtom({0, 0}, 26);
        edit::hotkey(doc, {0, -1}, "j");
        edit::hotkey(doc, {0, -1}, ".");
        REQUIRE(chem::properties(doc));
        CHECK(chem::properties(doc)->formula.find('*') != std::string::npos);
        CHECK(chem::hrmsLine(doc, chem::Ion::M).isEmpty());
    }
    SECTION("a named or variable attachment at a ring centre keeps its meaning") {
        for (bool variable : {false, true}) {
            Document doc;
            doc.addAtom({0, 0}, 26);
            edit::hotkey(doc, {0, -1}, "j");
            auto& dummy = doc.atoms.back();
            if (variable) dummy.attachments = {1, 2};
            else dummy.label = "R1";
            CHECK(chem::toSmiles(doc).find('*') != std::string::npos);
        }
    }
}

TEST_CASE("metal complexes keep their ligands' hydrogens through every format (#616)") {
    QFile f(QString(PENZENE_TEST_DATA) + "/metals.smi");
    REQUIRE(f.open(QIODevice::ReadOnly | QIODevice::Text));
    int read = 0;
    for (const QString& line : QString::fromUtf8(f.readAll()).split('\n', Qt::SkipEmptyParts)) {
        if (line.startsWith('#')) continue;
        const QStringList cols = line.split('\t');
        REQUIRE(cols.size() >= 3);
        INFO(cols[1].toStdString());
        auto doc = chem::fromSmiles(cols[0].toStdString());
        REQUIRE(doc);
        ++read;
        CHECK(chem::properties(*doc)->formula == cols[2]);
        const auto info = chem::atomInfo(*doc);
        QStringList donors;
        for (int i = 0; i < int(doc->atoms.size()); ++i) {
            CHECK_FALSE(info[i].valenceError);
            const int z = doc->atoms[i].z;
            const auto nbs = doc->neighbors(i);
            if ((z == 7 || z == 8 || z == 15) && std::any_of(nbs.begin(), nbs.end(), [&](int n) { return chem::isMetal(doc->atoms[n].z); }))
                donors << QString("%1=%2").arg(QString::fromStdString(chem::symbol(z))).arg(info[i].hydrogens);
        }
        donors.sort();
        CHECK(donors.join(',') == (cols.size() > 3 ? cols[3] : QString()));
        // Coordination to a metal is drawn as a plain line, and comes back the same from every format.
        for (const Bond& b : doc->bonds) CHECK(b.stereo != BondStereo::Dative);
        auto verify = [&](const Document& d) {
            REQUIRE(chem::properties(d));
            CHECK(chem::properties(d)->formula == cols[2]);
            const auto atoms = chem::atomInfo(d);
            QStringList hs;
            for (int i = 0; i < int(d.atoms.size()); ++i) {
                CHECK_FALSE(atoms[i].valenceError);
                const int z = d.atoms[i].z;
                const auto nbs = d.neighbors(i);
                if ((z == 7 || z == 8 || z == 15) && std::any_of(nbs.begin(), nbs.end(), [&](int n) { return chem::isMetal(d.atoms[n].z); }))
                    hs << QString("%1=%2").arg(QString::fromStdString(chem::symbol(z))).arg(atoms[i].hydrogens);
            }
            hs.sort();
            CHECK(hs == donors);
        };
        auto again = chem::fromSmiles(chem::toSmiles(*doc));
        REQUIRE(again);
        verify(*again);
        const auto block = chem::toMolBlock(*doc);
        auto mol = chem::fromMolBlock(block);
        REQUIRE(mol);
        verify(*mol);
        auto cdxml = chem::fromChemDraw(chem::toCdxml(*doc));
        REQUIRE(cdxml);
        verify(*cdxml);
        chem::Reaction reaction;
        reaction.reactants = {*doc}, reaction.products = {*doc};
        const auto rxn = chem::toRxn(reaction);
        CHECK(rxn.starts_with(block.find("V3000") != std::string::npos ? "$RXN V3000" : "$RXN\n"));
        for (const auto& scheme : {chem::fromRxn(rxn), chem::fromRdf(chem::toRdf({reaction}))}) {
            REQUIRE(scheme);
            const auto steps = chem::reactionsOf(*scheme);
            REQUIRE(steps.size() == 1);
            for (const auto* side : {&steps[0].reactants, &steps[0].products}) {
                REQUIRE_FALSE(side->empty());
                Document combined;
                for (const auto& part : *side) combined.append(part);
                verify(combined);
            }
        }
        auto penz = Document::fromJson(doc->toJson());
        REQUIRE(penz);
        CHECK(*penz == *doc);
    }
    CHECK(read >= 15);
}

TEST_CASE("a ligand typed with its hydrogens keeps them; a phosphine on a metal gains none (#616)") {
    Document doc;
    const int pt = doc.addAtom({0, 0}, 78), n = doc.addAtom({kBondLength, 0}, 7), o = doc.addAtom({-kBondLength, 0}, 8);
    doc.bonds.push_back({pt, n});
    doc.bonds.push_back({pt, o});
    CHECK(chem::atomInfo(doc)[n].hydrogens == 2);  // a plain line: an amide, as drawn
    REQUIRE(edit::applyLabel(doc, n, "NH3"));
    REQUIRE(edit::applyLabel(doc, o, "H2O"));
    CHECK(doc.atoms[o].z == 8);
    CHECK(edit::atomText(doc.atoms[n]) == "NH3");
    auto info = chem::atomInfo(doc);
    CHECK(info[n].hydrogens == 3);
    CHECK(info[o].hydrogens == 2);
    CHECK_FALSE(info[n].valenceError);
    CHECK(chem::properties(doc)->formula == "H5NOPt");
    REQUIRE(edit::applyLabel(doc, n, "N"));  // typed again without them: implicit again
    CHECK(chem::atomInfo(doc)[n].hydrogens == 2);
    // Off a metal, NH3 is still just N.
    Document chain;
    chain.addAtom({0, 0});
    const int end = chain.addAtom({kBondLength, 0});
    chain.bonds.push_back({0, end});
    REQUIRE(edit::applyLabel(chain, end, "NH3"));
    CHECK(chain.atoms[end].hydrogens == -1);

    REQUIRE(edit::applyLabel(doc, o, "PPh3"));
    const auto expanded = chem::expanded(doc);
    for (int i = 0; i < int(expanded.atoms.size()); ++i)
        if (expanded.atoms[i].z == 15) CHECK(chem::atomInfo(expanded)[i].hydrogens == 0);
    CHECK(chem::properties(doc)->formula == "C18H17NPPt");  // the N typed back to an amide
}

TEST_CASE("a square-planar metal drawn as a + isn't a Fischer projection (#616)") {
    Document doc;
    const int pt = doc.addAtom({0, 0}, 78);
    const QPointF at[] = {{kBondLength, 0}, {-kBondLength, 0}, {0, kBondLength}, {0, -kBondLength}};
    const int z[] = {17, 35, 53, 7};
    for (int k = 0; k < 4; ++k) doc.bonds.push_back({pt, doc.addAtom(at[k], z[k])});
    CHECK(chem::toSmiles(doc).find('@') == std::string::npos);
    const auto cleaned = chem::clean2D(doc);
    for (const Bond& b : cleaned.bonds) CHECK(b.stereo == BondStereo::None);
}

TEST_CASE("dative bonds: drawn as arrows, read from SMILES, kept in .penz and CDXML (#616)") {
    CHECK_FALSE(chem::isReactionSmiles("[NH3]->[Pt](<-[NH3])(Cl)Cl"));
    CHECK(chem::isReactionSmiles("CC=O>>CCO"));
    CHECK(chem::isReactionSmiles("[NH3]->[BH3]>>N.B"));
    Document doc;
    const int n = doc.addAtom({0, 0}, 7), b = doc.addAtom({kBondLength, 0}, 5);
    doc.bonds.push_back({b, n});
    REQUIRE(edit::hotkey(doc, {-1, 0}, "k").valid());
    CHECK(doc.bonds[0].stereo == BondStereo::Dative);
    REQUIRE(edit::hotkey(doc, {-1, 0}, "k").valid());  // again: the other way
    CHECK(doc.bonds[0].a == n);
    CHECK(chem::atomInfo(doc)[n].hydrogens == 3);  // the donor gives no hydrogen to the bond
    auto penz = Document::fromJson(doc.toJson());
    REQUIRE(penz);
    CHECK(penz->bonds[0].stereo == BondStereo::Dative);
    auto cdxml = chem::fromChemDraw(chem::toCdxml(doc));
    REQUIRE(cdxml);
    REQUIRE(cdxml->bonds.size() == 1);
    CHECK(cdxml->bonds[0].stereo == BondStereo::Dative);
    CHECK(chem::toSmiles(doc) == "B<-N");
    QImage image(200, 100, QImage::Format_ARGB32);
    image.fill(Qt::white);
    QPainter p(&image);
    paintDocument(p, doc);  // the arrowhead path
}

TEST_CASE("donor coordination respects charge, bridging chloride and covalent ligands (#616)") {
    for (const auto& [smiles, arrows] : std::vector<std::pair<std::string, int>>{
             {"[Pt]Cl[Pt]", 1}, {"C[N+](C)(C)[Pt]", 0}, {"C[O-][Pt]", 1},
             {"CN(C)[Ti]", 0}, {"N[Pt]", 0}, {"C[Mg]Br", 0},
             {"CP(C)(C)[Pd]", 1}, {"CN(C)(C)[Pt]", 1}, {"C=N(C)[Pt]", 1}}) {
        INFO(smiles);
        auto d = chem::fromSmiles(smiles);
        REQUIRE(d);
        const auto exported = chem::toSmiles(*d);
        CHECK_FALSE(exported.empty());
        int count = 0;
        for (size_t i = 0; i + 1 < exported.size(); ++i) count += exported.substr(i, 2) == "->" || exported.substr(i, 2) == "<-";
        CHECK(count == arrows);
        for (const auto& a : chem::atomInfo(*d)) CHECK_FALSE(a.valenceError);
    }
}

TEST_CASE("all metal ligand labels survive chemistry and CDXML (#616)") {
    for (const auto& [label, formula] : std::vector<std::pair<QString, std::string>>{
             {"CO", "COPt"}, {"PPh3", "C18H15PPt"}, {"PCy3", "C18H33PPt"}, {"PMe3", "C3H9PPt"},
             {"py", "C5H5NPt"}, {"THF", "C4H8OPt"}, {"MeCN", "C2H3NPt"},
             {"NH3", "H3NPt"}, {"H3N", "H3NPt"}, {"OH2", "H2OPt"}, {"H2O", "H2OPt"}, {"H3N+", "H3NPt+"}}) {
        INFO(label.toStdString());
        Document d;
        d.addAtom({0, 0}, 78), d.addAtom({kBondLength, 0}, 7);
        d.bonds.push_back({0, 1});
        REQUIRE(edit::applyLabel(d, 1, label));
        REQUIRE(chem::properties(d));
        CHECK(chem::properties(d)->formula == formula);
        auto roundtrip = chem::fromChemDraw(chem::toCdxml(d));
        REQUIRE(roundtrip);
        REQUIRE(chem::properties(*roundtrip));
        CHECK(chem::properties(*roundtrip)->formula == formula);
        for (const auto& a : chem::atomInfo(*roundtrip)) CHECK_FALSE(a.valenceError);
    }
}

TEST_CASE("CDXML keeps drawing atom indices with eta points and dative donors (#616)") {
    Document d;
    d.addAtom({0, 0}, 26);
    edit::hotkey(d, {0, -1}, "j");
    edit::hotkey(d, {0, -1}, "j");
    d.append(*chem::fromSmiles("[NH3]->[Pt](<-[NH3])(Cl)Cl"), {100, 0});
    auto roundtrip = chem::fromChemDraw(chem::toCdxml(d));
    REQUIRE(roundtrip);
    CHECK(roundtrip->atoms.size() == d.atoms.size());
    REQUIRE(chem::properties(*roundtrip));
    CHECK(chem::properties(*roundtrip)->formula == "C10H16Cl2FeN2Pt");
}

TEST_CASE("CDXML interaction bonds keep their drawing and hydrogen semantics (#616)") {
    Document d;
    d.addAtom({0, 0}, 78), d.addAtom({kBondLength, 0}, 7);
    d.bonds.push_back({0, 1, 1, BondStereo::Interaction});
    auto roundtrip = chem::fromChemDraw(chem::toCdxml(d));
    REQUIRE(roundtrip);
    REQUIRE(roundtrip->bonds.size() == 1);
    CHECK(roundtrip->bonds[0].stereo == BondStereo::Interaction);
    CHECK(chem::atomInfo(*roundtrip)[1].hydrogens == 3);
    CHECK(chem::properties(*roundtrip)->formula == "H3NPt");
}


static std::vector<QRectF> drawnLabels(const Document& doc) {
    QImage image(1, 1, QImage::Format_ARGB32);
    QPainter p(&image);
    std::vector<QRectF> labels(doc.atoms.size());
    paintDocument(p, doc, {Qt::black, Qt::black, 0, &labels});
    return labels;
}

TEST_CASE("metal coordination sites and Clean keep the drawing's topology (#618)") {
    auto angle = [](const Document& d, int m, int a, int b) {
        const auto u = unit(d.atoms[a].pos - d.atoms[m].pos), v = unit(d.atoms[b].pos - d.atoms[m].pos);
        return qRadiansToDegrees(std::acos(std::clamp(QPointF::dotProduct(u, v), -1.0, 1.0)));
    };
    SECTION("cisplatin's ammines occupy adjacent square-planar sites") {
        auto d = chem::fromSmiles("[NH3]->[Pt](<-[NH3])(Cl)Cl");
        REQUIRE(d);
        CHECK(angle(*d, 1, 0, 2) == Catch::Approx(90).margin(2));
        CHECK(angle(chem::clean2D(*d), 1, 0, 2) == Catch::Approx(90).margin(2));
    }
    SECTION("interaction bonds are layout connections and keep their drawing and chemistry") {
        Document d;
        d.addAtom({20, 30}, 78);
        for (QPointF v : {QPointF(1, 0), QPointF(0, 1), QPointF(-1, 0), QPointF(0, -1)}) {
            d.addAtom(QPointF(20, 30) + v * kBondLength, 7);
            d.atoms.back().hydrogens = 3;
            d.bonds.push_back({0, int(d.atoms.size()) - 1, 1, BondStereo::Interaction});
        }
        const auto smiles = chem::toSmiles(d);
        const Document clean = chem::clean2D(d);
        REQUIRE(clean.bonds.size() == 4);
        for (const auto& b : clean.bonds) {
            CHECK(b.stereo == BondStereo::Interaction);
            CHECK(len(clean.atoms[b.a].pos - clean.atoms[b.b].pos) == Catch::Approx(kBondLength).margin(0.5));
        }
        CHECK(chem::toSmiles(clean) == smiles);
        CHECK(angle(clean, 0, 1, 3) == Catch::Approx(180).margin(2));
    }
    SECTION("tagged square planar cis and trans survive Clean") {
        for (const auto& [smi, expected] : std::vector<std::pair<std::string, double>>{
                 {"Cl[Pt@SP1](Cl)([NH3])[NH3]", 90}, {"Cl[Pt@SP2](Cl)([NH3])[NH3]", 180},
                 {"Cl[Pt@SP1]([NH3])([OH2])Br", 90}}) {
            auto d = chem::fromSmiles(smi);
            REQUIRE(d);
            CHECK(angle(*d, 1, 0, 2) == Catch::Approx(expected).margin(2));
            CHECK(angle(chem::clean2D(*d), 1, 0, 2) == Catch::Approx(expected).margin(2));
        }
    }
    SECTION("tagged trigonal-bipyramidal and octahedral axial pairs survive Clean") {
        for (const char* smi : {"[Cl-]->[Fe@TB1](<-[NH3])(<-[OH2])(<-P(C)(C)C)<-[Br-]",
                               "[Cl-]->[Co@OH1](<-[NH3])(<-[OH2])(<-P(C)(C)C)(F)Br"}) {
            auto d = chem::fromSmiles(smi);
            REQUIRE(d);
            const int bromine = int(std::find_if(d->atoms.begin(), d->atoms.end(), [](const Atom& a) { return a.z == 35; }) - d->atoms.begin());
            CHECK(angle(*d, 1, 0, bromine) == Catch::Approx(180).margin(2));
            CHECK(angle(chem::clean2D(*d), 1, 0, bromine) == Catch::Approx(180).margin(2));
        }
    }
    SECTION("CX coordinates keep their supplied scale and orientation") {
        auto d = chem::fromSmiles("[NH3]->[Pt]<-[NH3] |(0,0,;2,0,;2,3,)|");
        REQUIRE(d);
        CHECK(d->atoms[0].pos == QPointF(0, 0));
        CHECK(d->atoms[1].pos.x() == Catch::Approx(2 * kBondLength / 1.5));
        CHECK(d->atoms[2].pos.y() == Catch::Approx(-3 * kBondLength / 1.5));
    }
    SECTION("octahedral drawings include an axial wedge and hash") {
        auto d = chem::fromSmiles("[NH3]->[Co+3](<-[NH3])(<-[NH3])(<-[NH3])(<-[NH3])<-[NH3]");
        REQUIRE(d);
        int wedges = 0, hashes = 0;
        for (const auto& b : d->bonds) wedges += b.stereo == BondStereo::Wedge, hashes += b.stereo == BondStereo::Hash;
        CHECK(wedges == 1);
        CHECK(hashes == 1);
        d->style = "RSC";
        d->labelRatio = 1.5;
        const auto labels = drawnLabels(chem::clean2D(*d));
        for (size_t i = 0; i < labels.size(); ++i)
            for (size_t j = i + 1; j < labels.size(); ++j) CHECK_FALSE(labels[i].intersects(labels[j]));
    }
}

TEST_CASE("metal corpus labels keep clear of one another (#618)") {
    QFile f(QString(PENZENE_TEST_DATA) + "/metals.smi");
    REQUIRE(f.open(QIODevice::ReadOnly | QIODevice::Text));
    for (const auto& line : QString::fromUtf8(f.readAll()).split('\n', Qt::SkipEmptyParts)) {
        if (line.startsWith('#')) continue;
        const auto cols = line.split('\t');
        INFO(cols[1].toStdString());
        auto doc = chem::fromSmiles(cols[0].toStdString());
        REQUIRE(doc);
        for (const Document& d : {*doc, chem::clean2D(*doc)}) {
            auto boxes = drawnLabels(d);
            for (size_t i = 0; i < boxes.size(); ++i)
                for (size_t j = i + 1; j < boxes.size(); ++j) {
                    INFO(i << " " << j);
                    CHECK_FALSE(boxes[i].intersects(boxes[j]));
                }
        }
    }
}

TEST_CASE("bulky metal ligands and coordination lines keep clear (#618)") {
    auto doc = chem::fromSmiles("Cl[Ru](Cl)(=Cc1ccccc1)([P](C2CCCCC2)(C2CCCCC2)C2CCCCC2)[P](C2CCCCC2)(C2CCCCC2)C2CCCCC2");
    REQUIRE(doc);
    for (const Document& d : {*doc, chem::clean2D(*doc)}) {
        RDKit::RWMol cut;
        for (const auto& a : d.atoms) cut.addAtom(new RDKit::Atom(a.z));
        for (const auto& b : d.bonds)
            if (b.a != 1 && b.b != 1) cut.addBond(b.a, b.b, RDKit::Bond::SINGLE);
        std::vector<int> membership;
        RDKit::MolOps::getMolFrags(cut, membership);
        const auto labels = drawnLabels(d);
        for (size_t i = 0; i < d.atoms.size(); ++i)
            for (size_t j = i + 1; j < d.atoms.size(); ++j)
                if (membership[i] != membership[j]) CHECK(len(d.atoms[i].pos - d.atoms[j].pos) >= 0.65 * kBondLength / 1.5);
        for (size_t i = 0; i < d.bonds.size(); ++i) {
            const auto& a = d.bonds[i];
            const int af = membership[a.a == 1 ? a.b : a.a];
            const QLineF line(d.atoms[a.a].pos, d.atoms[a.b].pos);
            for (size_t j = i + 1; j < d.bonds.size(); ++j) {
                const auto& b = d.bonds[j];
                const int bf = membership[b.a == 1 ? b.b : b.a];
                if (af == bf || a.a == b.a || a.a == b.b || a.b == b.a || a.b == b.b) continue;
                CHECK(line.intersects(QLineF(d.atoms[b.a].pos, d.atoms[b.b].pos)) != QLineF::BoundedIntersection);
            }
            for (size_t j = 0; j < labels.size(); ++j) {
                if (membership[j] == af || j == 1 || labels[j].isEmpty()) continue;
                const QRectF& box = labels[j];
                CHECK_FALSE(box.contains(line.p1()));
                CHECK_FALSE(box.contains(line.p2()));
                for (const auto& edge : {QLineF(box.topLeft(), box.topRight()), QLineF(box.topRight(), box.bottomRight()),
                                         QLineF(box.bottomRight(), box.bottomLeft()), QLineF(box.bottomLeft(), box.topLeft())})
                    CHECK(line.intersects(edge) != QLineF::BoundedIntersection);
            }
        }
    }
    QFile f(QString(PENZENE_TEST_DATA) + "/metal2d-readme.smi");
    REQUIRE(f.open(QIODevice::ReadOnly | QIODevice::Text));
    auto crowded = chem::fromSmiles(f.readLine().split('\t')[0].toStdString());
    REQUIRE(crowded);
    crowded->labelRatio = 1.4;  // also cover fonts wider than the local platform's
    const auto labels = drawnLabels(chem::clean2D(*crowded));
    for (size_t i = 0; i < labels.size(); ++i)
        for (size_t j = i + 1; j < labels.size(); ++j) CHECK_FALSE(labels[i].intersects(labels[j]));
}

TEST_CASE("metal2d README examples and the polynuclear fallback (#618)") {
    QFile f(QString(PENZENE_TEST_DATA) + "/metal2d-readme.smi");
    REQUIRE(f.open(QIODevice::ReadOnly | QIODevice::Text));
    for (const auto& line : QString::fromUtf8(f.readAll()).split('\n', Qt::SkipEmptyParts)) {
        const auto cols = line.split('\t');
        INFO(cols[1].toStdString());
        auto doc = chem::fromSmiles(cols[0].toStdString());
        REQUIRE(doc);
        const int metals = int(std::count_if(doc->atoms.begin(), doc->atoms.end(), [](const Atom& a) { return chem::isMetal(a.z); }));
        if (metals > 1) {
            // The fallback remains the platform's native depictor, whose CoordGen layouts vary across platforms.
            std::unique_ptr<RDKit::RWMol> native(RDKit::SmilesToMol(cols[0].toStdString()));
            REQUIRE(native);
            const auto& rings = native->getRingInfo()->atomRings();
            RDDepict::Compute2DCoordParameters params;
            params.canonOrient = true;
            params.useRingTemplates = true;
            params.forceRDKit = std::none_of(rings.begin(), rings.end(), [](const auto& r) { return r.size() >= 9; });
            RDDepict::compute2DCoords(*native, params);
            const auto& conf = native->getConformer();
            double sum = 0;
            for (const auto* bond : native->bonds())
                sum += (conf.getAtomPos(bond->getBeginAtomIdx()) - conf.getAtomPos(bond->getEndAtomIdx())).length();
            const double scale = kBondLength * native->getNumBonds() / sum;
            REQUIRE(native->getNumAtoms() == doc->atoms.size());
            for (size_t i = 0; i < doc->atoms.size(); ++i) {
                const auto& p = conf.getAtomPos(i);
                CHECK(len(QPointF(p.x * scale, -p.y * scale) - doc->atoms[i].pos) < 1e-5);
            }
        } else {
            for (const Document& d : {*doc, chem::clean2D(*doc)}) {
                const auto boxes = drawnLabels(d);
                for (size_t i = 0; i < boxes.size(); ++i)
                    for (size_t j = i + 1; j < boxes.size(); ++j) {
                        INFO(i << " " << j);
                        CHECK_FALSE(boxes[i].intersects(boxes[j]));
                    }
            }
        }
    }
}


TEST_CASE("format toolbar model: text and separate fills round trip (#621)") {
    Document d;
    d.style = "RSC";
    Arrow box{{0, 0}, {40, 20}, ArrowKind::RoundedBox};
    box.color = Qt::red, box.filled = true, box.fill = QColor(207, 227, 255);
    d.arrows.push_back(box);
    Text text{{0, 40}, "H2O and SO4^2-"};
    text.font = "Times New Roman", text.bold = text.italic = true, text.scale = 1.5, text.color = Qt::blue;
    d.texts.push_back(text);
    REQUIRE(Document::fromJson(d.toJson()));
    CHECK(*Document::fromJson(d.toJson()) == d);
    CHECK(sheetsFromJson(sheetsToJson({{"Sample", d}}))[0].doc == d);
    for (const QByteArray& bytes : {chem::toCdxml(d), chem::toCdx(d)}) {
        if (bytes.isEmpty()) continue;  // a build without binary CDX
        const auto back = chem::fromChemDraw(bytes);
        REQUIRE(back);
        REQUIRE(back->arrows.size() == 1);
        CHECK(back->arrows[0].filled);
        CHECK(back->arrows[0].fill == box.fill);
        CHECK(back->arrows[0].color == box.color);
        REQUIRE(back->texts.size() == 1);
        CHECK(back->texts[0].font == text.font);
        CHECK(back->texts[0].bold);
        CHECK(back->texts[0].italic);
        CHECK(back->texts[0].color == text.color);
        CHECK(std::abs(back->texts[0].scale - text.scale) < 1e-4);
        CHECK(back->texts[0].text == text.text);
    }
    const auto xml = chem::toCdxml(d);
    CHECK(xml.contains("bgcolor="));
    CHECK(xml.contains("face=\"3\""));
    Text styled = text;
    styled.text = "Format";
    QFont expected = labelFont(documentStyle(d), styled.scale);
    expected.setFamily(styled.font);
    expected.setBold(styled.bold);
    expected.setItalic(styled.italic);
    QPainterPath expectedPath;
    expectedPath.addText(styled.pos, expected, styled.text);
    CHECK(textPath(styled, documentStyle(d)) == expectedPath);
    CHECK(renderSvg(d).contains("#cfe3ff"));
}

TEST_CASE("old filled shapes retain their line colour (#621)") {
    const auto old = Document::fromJson(R"({"format":"penzene","version":1,"arrows":[{"x1":0,"y1":0,"x2":30,"y2":20,"kind":"box","filled":true,"color":"#ff0000"}],"texts":[{"x":0,"y":40,"text":"heat"}]})");
    REQUIRE(old);
    REQUIRE(old->arrows.size() == 1);
    CHECK(old->arrows[0].filled);
    CHECK(!old->arrows[0].fill.isValid());
    CHECK(old->texts[0].font.isEmpty());
    CHECK(!old->texts[0].bold);
    CHECK(!old->texts[0].italic);
    CHECK(old->texts[0].scale == 1);
    CHECK(renderSvg(*old).contains("#ff0000"));
    CHECK(Document::fromJson(old->toJson()) == old);
}

TEST_CASE("portable styles retain advanced typography, layout and legacy defaults") {
    auto settings = styleSettings(drawingStyles().front());
    settings.name = "Advanced style regression";
    settings.chainAngle = 100;
    settings.snapAngle = 15;
    settings.bondSpacingAbsolute = 1.25;
    settings.italic = true;
    settings.captionFont = "Helvetica";
    settings.captionSize = 6.5;
    settings.captionWeight = 600;
    settings.captionItalic = true;
    settings.scriptScale = 0.8;
    settings.subscriptOffset = 0.25;
    settings.superscriptOffset = 0.6;
    settings.chargesOnAtom = true;
    settings.rGroupSuperscript = true;
    settings.hideImplicitHydrogens = true;
    settings.carbonLabels = "terminal";
    CHECK(styleFromJson(styleToJson(settings)) == settings);
    const auto physical = styleSettings(modelStyle(settings));
    CHECK(physical == settings);
    Document doc = *chem::fromSmiles("CCO");
    doc.styleSettings = settings;
    REQUIRE(Document::fromJson(doc.toJson()));
    CHECK(Document::fromJson(doc.toJson())->styleSettings == settings);
    auto legacy = styleToJson(styleSettings(drawingStyles().front()));
    for (const auto& key : {"chainAngle", "snapAngle", "bondSpacingAbsolute", "italic", "captionFont", "captionSize", "captionWeight",
                           "captionItalic", "scriptScale", "subscriptOffset", "superscriptOffset", "chargesOnAtom",
                           "rGroupSuperscript", "hideImplicitHydrogens", "carbonLabels"}) legacy.remove(key);
    REQUIRE(styleFromJson(legacy));
    CHECK(styleFromJson(legacy)->chainAngle == 120);
    CHECK(styleFromJson(legacy)->captionSize == 0);
    CHECK_FALSE(styleFromJson(legacy)->hideImplicitHydrogens);
    for (const auto& key : {"chainAngle", "bondSpacingAbsolute", "captionSize", "scriptScale", "italic", "captionWeight", "captionFont", "carbonLabels"}) {
        auto bad = styleToJson(settings);
        bad[key] = QJsonArray{};
        CHECK_FALSE(styleFromJson(bad));
    }
    const QPointF back = unit(doc.atoms[0].pos - doc.atoms[1].pos);
    Document ethane = *chem::fromSmiles("CC");
    ethane.styleSettings = settings;
    const QPointF direction = edit::freeDirection(ethane, 1);
    CHECK(QPointF::dotProduct(direction, unit(ethane.atoms[0].pos - ethane.atoms[1].pos)) == Catch::Approx(std::cos(qDegreesToRadians(100.0))));
    CHECK(scripts("R1", true)[1] == Script::Super);
    CHECK(scripts("H2O", true)[1] == Script::Sub);
}

TEST_CASE("fractional label sizes and independent caption sizes survive physical scaling") {
    auto style = drawingStyle("RSC");
    style.fontSize = 8.1;
    const auto small = labelFont(style);
    style.fontSize = 8.9;
    const auto large = labelFont(style);
    CHECK(large.pointSizeF() * QFontMetricsF(large).fontDpi() / 72 == Catch::Approx(8.9));
    QPainterPath a, b;
    a.addText(0, 0, small, "NH");
    b.addText(0, 0, large, "NH");
    CHECK(b.boundingRect().height() > a.boundingRect().height() * 1.05);
    style.captionSize = 6.5;
    style.captionFont = "Arial";
    style.captionWeight = 700;
    style.captionItalic = true;
    const auto caption = captionStyle(style);
    CHECK(caption.fontSize == 6.5);
    CHECK(caption.font == "Arial");
    CHECK(labelFont(caption).bold());
    CHECK(labelFont(caption).italic());
}

TEST_CASE("publisher point dimensions are used without guessed wedge or hash multipliers") {
    const auto acs = styleSettings(drawingStyle("ACS 1996"));
    CHECK(acs.wedgeWidth == 3);
    CHECK(acs.hashSpacing == 2.5);
    CHECK(acs.labelMargin == 1.6);
    CHECK(acs.labelRadius == 0);
    const auto rsc = styleSettings(drawingStyle("RSC"));
    CHECK(rsc.lineWidth == Catch::Approx(0.5));
    CHECK(rsc.boldWidth == Catch::Approx(1.6));
    CHECK(rsc.wedgeWidth == Catch::Approx(1.6));
    CHECK(rsc.hashSpacing == Catch::Approx(1.8));
    CHECK(rsc.fontSize == Catch::Approx(7));
    CHECK(rsc.snapAngle == 15);
    const QPointF snapped = edit::snapped({0, 0}, dirAt(17), rsc.snapAngle);
    CHECK(snapped.x() == Catch::Approx(dirAt(15).x()));
    CHECK(snapped.y() == Catch::Approx(dirAt(15).y()));
}

TEST_CASE("ChemDraw carries physical style dimensions and separate caption typography") {
    Document doc = *chem::fromSmiles("CCO");
    auto style = styleSettings(drawingStyle("RSC"));
    style.captionFont = "Arial";
    style.captionSize = 6.5;
    style.captionWeight = 700;
    style.captionItalic = true;
    style.bondSpacingAbsolute = 1.25;
    style.chainAngle = 110;
    doc.styleSettings = style;
    doc.texts.push_back({{0, 40}, "H2O"});
    std::vector<QByteArray> files{chem::toCdxml(doc)};
if (const auto cdx = chem::toCdx(doc); !cdx.isEmpty()) files.push_back(chem::cdxToCdxml(cdx));
    for (const auto& xml : files) {
        QXmlStreamReader reader(xml);
        REQUIRE(reader.readNextStartElement());
        const auto attrs = reader.attributes();
        CHECK(attrs.value("MarginWidth").toDouble() == Catch::Approx(1.25));
        CHECK(attrs.value("ChainAngle").toDouble() == Catch::Approx(110));
        CHECK(attrs.value("CaptionSize").toDouble() == Catch::Approx(6.5));
        CHECK(attrs.value("CaptionFace").toInt() == 3);
        CHECK(attrs.value("BondSpacingAbs").toDouble() == Catch::Approx(1.25));
        CHECK(attrs.value("BondSpacingType").toInt() == 1);
        const auto back = chem::fromChemDraw(xml);
        REQUIRE(back);
        CHECK(styleSettings(documentStyle(*back)).labelMargin == Catch::Approx(1.25));
        CHECK(styleSettings(documentStyle(*back)).bondSpacingAbsolute == Catch::Approx(1.25));
        CHECK(styleSettings(documentStyle(*back)).chainAngle == Catch::Approx(110));
        CHECK(captionStyle(documentStyle(*back)).font == "Arial");
        CHECK(captionStyle(documentStyle(*back)).fontSize * exportScale(*back) == Catch::Approx(6.5));
        REQUIRE(back->texts.size() == 1);
        CHECK(back->texts[0].scale == Catch::Approx(1));
    }
}

TEST_CASE("absolute double spacing and hash spacing reach the actual rendered strokes") {
    Document doc;
    doc.atoms = {{{0, 0}}, {{2 * kBondLength, 0}}};
    doc.bonds = {{0, 1}};
    auto style = styleSettings(drawingStyle("ACS 1996"));
    style.bondSpacingAbsolute = 2;
    doc.styleSettings = style;
    auto ink = [&](QPointF point) {
        const double margin = 2;  // model units of room to probe beyond the bond, now that exports crop tight
        const QImage image = renderImage(doc, {720, Qt::white, 1, margin * exportScale(doc)});
        const QPoint pixel = ((point - outputBounds(doc).topLeft() + QPointF(margin, margin)) * exportScale(doc) * 10).toPoint();
        return image.pixelColor(pixel).lightness() < 160;
    };
    doc.bonds[0].order = 2;
    CHECK(ink({kBondLength, 1}));
    CHECK_FALSE(ink({kBondLength, 2.6}));
    doc.styleSettings->bondSpacingAbsolute = 0;
    CHECK_FALSE(ink({kBondLength, 1}));
    CHECK(ink({kBondLength, 0.18 * kBondLength}));
    doc.bonds[0].order = 1;
    doc.bonds[0].stereo = BondStereo::Hash;
    const auto image = renderImage(doc, {720, Qt::white});
    const auto bounds = outputBounds(doc);
    const int row = qRound(-bounds.top() * 10);
    std::vector<double> centres;
    int first = -1;
    for (int x = 0; x < image.width(); ++x) {
        if (image.pixelColor(x, row).lightness() < 160) { if (first < 0) first = x; }
        else if (first >= 0) { centres.push_back((first + x - 1) / 20.0); first = -1; }
    }
    REQUIRE(centres.size() > 5);
    for (size_t i = 2; i < centres.size(); ++i)
        CHECK(centres[i] - centres[i - 1] == Catch::Approx(style.hashSpacing).margin(0.1));
}

TEST_CASE("styles carry an optional charge-circle choice without locking document overrides") {
    auto style = styleSettings(drawingStyle("RSC"));
    style.name = "Charge style regression";
    style.circledCharges = false;
    deleteDrawingStyle(style.name);
    struct Cleanup { QString name; ~Cleanup() { deleteDrawingStyle(name); } } cleanup{style.name};
    REQUIRE(saveDrawingStyle(style));
    Document doc = *chem::fromSmiles("C[NH3+]");
    doc.circledCharges = true;
    REQUIRE(applyDrawingStyle(doc, style.name));
    CHECK_FALSE(doc.circledCharges);
    CHECK(doc.styleSettings->circledCharges == false);
    CHECK(Document::fromJson(doc.toJson()) == doc);
    doc.circledCharges = true;
    CHECK(Document::fromJson(doc.toJson())->circledCharges);
    auto json = styleToJson(style);
    json.remove("circledCharges");
    REQUIRE(styleFromJson(json));
    CHECK_FALSE(styleFromJson(json)->circledCharges.has_value());
    json["circledCharges"] = "false";
    CHECK_FALSE(styleFromJson(json));
}

TEST_CASE("caption typography does not change abbreviated atom labels") {
    Document doc = *chem::fromSmiles("CO");
    REQUIRE(edit::applyLabel(doc, 1, "OMe", true));
    auto settings = styleSettings(drawingStyle("RSC"));
    doc.styleSettings = settings;
    const auto before = renderImage(doc);
    doc.styleSettings->captionFont = "Times New Roman";
    doc.styleSettings->captionSize = 25;
    doc.styleSettings->captionWeight = 700;
    CHECK(renderImage(doc) == before);
}

TEST_CASE("Nature applies its complete supported label conventions and keeps them portable") {
    Document doc = *chem::fromSmiles("c1ccccc1[NH3+]");
    doc.circledCharges = doc.aromaticCircles = true;
    doc.aromaticCircleOverrides = {{0, 1, 2, 3, 4, 5}};
    doc.carbonLabels = Document::CarbonLabels::All;
    REQUIRE(applyDrawingStyle(doc, "Nature"));
    CHECK_FALSE(doc.circledCharges);
    CHECK_FALSE(doc.aromaticCircles);
    CHECK(doc.aromaticCircleOverrides.empty());
    CHECK(doc.carbonLabels == Document::CarbonLabels::None);
    const auto settings = styleSettings(documentStyle(doc));
    CHECK(settings.chargesOnAtom);
    CHECK(settings.rGroupSuperscript);
    CHECK_FALSE(settings.hideImplicitHydrogens);
    CHECK(settings.weight == 400);
    CHECK(settings.captionSize == 6);
    CHECK(styleFromJson(styleToJson(settings)) == settings);
    CHECK(Document::fromJson(doc.toJson()) == doc);
}

TEST_CASE("a bracket's superscript is saved and drawn at its top right (#658)") {
    auto doc = Document::fromJson(R"({"format":"penzene","version":1,"atoms":[{"x":0,"y":0,"z":7},{"x":14.4,"y":0,"z":6}],
        "bonds":[{"a":0,"b":1,"order":3}],"brackets":[{"atoms":[0,1],"label":"n","superscript":"‡"}]})");
    REQUIRE(doc);
    REQUIRE(doc->brackets.size() == 1);
    CHECK(doc->brackets[0].superscript == "‡");
    CHECK(*Document::fromJson(doc->toJson()) == *doc);
    const double k = 10;
    auto render = [&](const Document& d) {
        QImage img(600, 600, QImage::Format_ARGB32);
        img.fill(Qt::white);
        QPainter p(&img);
        p.translate(300, 300);
        p.scale(k, k);
        p.translate(-7.2, 0);
        paintDocument(p, d, {Qt::black, Qt::black, 0.6});
        return img;
    };
    Document plain = *doc;
    plain.brackets[0].superscript.clear();
    const QImage with = render(*doc), without = render(plain);
    QRect sup;
    for (int y = 0; y < 600; ++y)
        for (int x = 0; x < 600; ++x)
            if (qGray(with.pixel(x, y)) < 128 && qGray(without.pixel(x, y)) >= 128) sup |= QRect(x, y, 1, 1);
    REQUIRE(sup.isValid());
    CHECK(sup.left() > 300 + 0.55 * kBondLength * k);  // right of the closing bracket
    CHECK(sup.bottom() < 300);                          // in the upper half, where the subscript isn't
}
