#include "StyleEditor.h"
#include "Styles.h"
#include "Chem.h"
#include "NumberBox.h"
#include <QComboBox>
#include <QCheckBox>
#include <QSignalBlocker>
#include <QDialog>
#include <QDialogButtonBox>
#include <QFile>
#include <QFileDialog>
#include <QFontDatabase>
#include <QFormLayout>
#include <QHBoxLayout>
#include <QJsonDocument>
#include <QLabel>
#include <QLineEdit>
#include <QListWidget>
#include <QMessageBox>
#include <QPainter>
#include <QPushButton>
#include <QTabWidget>
#include <QVBoxLayout>
#include <algorithm>

namespace {
class StylePreview : public QWidget {
public:
    Document document;
    QString message;
    double scale = 2;
    bool fit = false;

    using QWidget::QWidget;

protected:
    void paintEvent(QPaintEvent*) override {
        QPainter p(this);
        p.fillRect(rect(), Qt::white);
        p.setRenderHint(QPainter::Antialiasing);
        if (!message.isEmpty()) {
            p.setPen(Qt::black);
            p.drawText(rect().adjusted(16, 16, -16, -16), Qt::AlignCenter | Qt::TextWordWrap, message);
            return;
        }
        const QRectF bounds = outputBounds(document);
        const double zoom = fit ? std::min((width() - 32) / std::max(1.0, bounds.width()),
                                          (height() - 32) / std::max(1.0, bounds.height()))
                                : scale * exportScale(document);
        p.translate(rect().center());
        p.scale(zoom, zoom);
        p.translate(-bounds.center());
        paintDocument(p, document);
    }
};
}

void editDrawingStyles(QWidget* parent, const Document& document) {
    QDialog dialog(parent);
    dialog.setObjectName("styleEditor");
    dialog.setWindowTitle(QObject::tr("Drawing Styles"));
    dialog.resize(980, 610);
    auto* outer = new QVBoxLayout(&dialog);
    outer->setContentsMargins(24, 24, 24, 20);
    outer->setSpacing(20);
    auto* columns = new QHBoxLayout;
    columns->setSpacing(24);
    outer->addLayout(columns, 1);
    auto* library = new QVBoxLayout;
    library->setSpacing(8);
    columns->addLayout(library);
    auto heading = [](const QString& title) {
        auto* label = new QLabel(title);
        label->setObjectName("sectionTitle");
        return label;
    };
    library->addWidget(heading(QObject::tr("Styles")));
    auto* list = new QListWidget;
    list->setObjectName("styleList");
    list->setAccessibleName(QObject::tr("Drawing styles"));
    list->setMinimumWidth(170);
    list->setMaximumWidth(210);
    list->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    library->addWidget(list, 1);
    auto* libraryActions = new QHBoxLayout;
    library->addLayout(libraryActions);
    auto* right = new QVBoxLayout;
    right->setSpacing(12);
    columns->addLayout(right, 2);
    auto* name = new QLineEdit;
    name->setObjectName("styleName");
    name->setAccessibleName(QObject::tr("Style name"));
    auto* identity = new QFormLayout;
    identity->setContentsMargins(0, 0, 0, 0);
    identity->addRow(QObject::tr("Name:"), name);
    right->addLayout(identity);
    auto* state = new QLabel;
    state->setObjectName("secondary");
    right->addWidget(state);
    auto* tabs = new QTabWidget;
    tabs->setObjectName("styleSections");
    tabs->setAccessibleName(QObject::tr("Style settings"));
    right->addWidget(tabs, 1);
    auto form = [&](const QString& title) {
        auto* page = new QWidget;
        auto* layout = new QFormLayout(page);
        layout->setContentsMargins(0, 16, 0, 0);
        layout->setVerticalSpacing(10);
        layout->setFieldGrowthPolicy(QFormLayout::AllNonFixedFieldsGrow);
        tabs->addTab(page, title);
        return layout;
    };
    auto* bonds = form(QObject::tr("Bonds"));
    auto* labels = form(QObject::tr("Labels"));
    auto* captions = form(QObject::tr("Captions"));
    auto* layout = form(QObject::tr("Layout"));
    auto* gapMode = new QComboBox;
    gapMode->setObjectName("bondSpacingMode");
    gapMode->addItems({QObject::tr("Relative (%)"), QObject::tr("Absolute (pt)")});
    bonds->addRow(QObject::tr("Double-bond spacing:"), gapMode);
    struct Measure { QFormLayout* section; const char* key; QString label; double StyleSettings::*member; double min, max; };
    const std::vector<Measure> measures{
        {bonds, "bondLength", QObject::tr("Bond length:"), &StyleSettings::bondLength, 0.1, 1000},
        {bonds, "lineWidth", QObject::tr("Line width:"), &StyleSettings::lineWidth, 0.001, 100},
        {bonds, "boldWidth", QObject::tr("Bold width:"), &StyleSettings::boldWidth, 0.001, 100},
        {bonds, "wedgeWidth", QObject::tr("Wedge width:"), &StyleSettings::wedgeWidth, 0.001, 100},
        {bonds, "hashSpacing", QObject::tr("Hash spacing:"), &StyleSettings::hashSpacing, 0.001, 100},
        {bonds, "doubleBondGap", QObject::tr("Double-bond gap:"), &StyleSettings::bondSpacing, 0, 1000},
        {labels, "labelMargin", QObject::tr("Label margin:"), &StyleSettings::labelMargin, 0, 100},
        {labels, "fontSize", QObject::tr("Label size:"), &StyleSettings::fontSize, 1, 200},
        {labels, "labelRadius", QObject::tr("Minimum label clearance:"), &StyleSettings::labelRadius, 0, 1000},
        {layout, "chainAngle", QObject::tr("Chain angle:"), &StyleSettings::chainAngle, 1, 179},
        {captions, "captionSize", QObject::tr("Caption size:"), &StyleSettings::captionSize, 0, 200},
        {layout, "scriptScale", QObject::tr("Script size / label size:"), &StyleSettings::scriptScale, 0.05, 2},
        {layout, "subscriptOffset", QObject::tr("Subscript drop / cap height:"), &StyleSettings::subscriptOffset, -2, 2},
        {layout, "superscriptOffset", QObject::tr("Superscript rise / cap height:"), &StyleSettings::superscriptOffset, -2, 2},
        {layout, "snapAngle", QObject::tr("Fixed angle increment:"), &StyleSettings::snapAngle, 1, 90}};
    std::vector<QDoubleSpinBox*> numbers;
    for (size_t i = 0; i < measures.size(); ++i) {
        const auto& m = measures[i];
        auto* box = new NumberBox;
        box->setObjectName(m.key);
        box->setAccessibleName(m.label);
        box->setDecimals(6);
        box->setRange(m.min, m.max);
        box->setSuffix(QObject::tr(" pt"));
        box->setSingleStep(0.1);
        if (m.section == layout) box->setSuffix({});
        if (m.member == &StyleSettings::chainAngle || m.member == &StyleSettings::snapAngle) box->setSuffix(QObject::tr("°"));
        if (m.member == &StyleSettings::captionSize) box->setSpecialValueText(QObject::tr("Same as labels"));
        m.section->addRow(m.label, box);
        numbers.push_back(box);
    }
    numbers[8]->setToolTip(QObject::tr("Minimum distance from an atom centre to a bond end. Set to zero to use only the label margin."));
    auto* font = new QComboBox;
    font->setEditable(true);
    font->addItems(QFontDatabase::families());
    font->setSizeAdjustPolicy(QComboBox::AdjustToMinimumContentsLengthWithIcon);
    font->setMinimumContentsLength(12);
    font->setObjectName("styleFont");
    labels->addRow(QObject::tr("Label font:"), font);
    auto* weight = new QComboBox;
    weight->setObjectName("styleWeight");
    for (const auto& pair : {std::pair{QObject::tr("Thin"), 100}, {QObject::tr("Extra light"), 200},
         {QObject::tr("Light"), 300}, {QObject::tr("Normal"), 400}, {QObject::tr("Medium"), 500},
         {QObject::tr("Demi bold"), 600}, {QObject::tr("Bold"), 700}, {QObject::tr("Extra bold"), 800}, {QObject::tr("Black"), 900}})
        weight->addItem(pair.first, pair.second);
    labels->addRow(QObject::tr("Label weight:"), weight);
    auto* italic = new QCheckBox(QObject::tr("Italic atom labels"));
    italic->setObjectName("styleItalic");
    labels->addRow(italic);
    auto* captionFont = new QComboBox;
    captionFont->setObjectName("captionFont");
    captionFont->setEditable(true);
    captionFont->setSizeAdjustPolicy(QComboBox::AdjustToMinimumContentsLengthWithIcon);
    captionFont->setMinimumContentsLength(12);
    captionFont->addItem(QObject::tr("Same as labels"));
    captionFont->addItems(QFontDatabase::families());
    captions->addRow(QObject::tr("Caption font:"), captionFont);
    auto* captionWeight = new QComboBox;
    captionWeight->setObjectName("captionWeight");
    for (int i = 0; i < weight->count(); ++i) captionWeight->addItem(weight->itemText(i), weight->itemData(i));
    captions->addRow(QObject::tr("Caption weight:"), captionWeight);
    auto* captionItalic = new QCheckBox(QObject::tr("Italic captions"));
    captionItalic->setObjectName("captionItalic");
    captions->addRow(captionItalic);
    auto* chargesOnAtom = new QCheckBox(QObject::tr("Place charges on the atom, before H"));
    chargesOnAtom->setObjectName("chargesOnAtom");
    auto* rGroupSuperscript = new QCheckBox(QObject::tr("Superscript R-group numbers (R¹, R²)"));
    rGroupSuperscript->setObjectName("rGroupSuperscript");
    auto* hideHydrogens = new QCheckBox(QObject::tr("Hide implicit hydrogens on labels"));
    hideHydrogens->setObjectName("hideImplicitHydrogens");
    layout->addRow(chargesOnAtom);
    layout->addRow(rGroupSuperscript);
    layout->addRow(hideHydrogens);
    auto* carbons = new QComboBox;
    carbons->setObjectName("styleCarbonLabels");
    for (const auto& pair : {std::pair{QObject::tr("Keep document choice"), ""}, {QObject::tr("Skeletal"), "none"},
                            {QObject::tr("Terminal"), "terminal"}, {QObject::tr("All carbons"), "all"}})
        carbons->addItem(pair.first, pair.second);
    layout->addRow(QObject::tr("Carbon labels:"), carbons);

    auto* fontNotice = new QLabel;
    fontNotice->setWordWrap(true);
    fontNotice->setObjectName("secondary");
    right->addWidget(fontNotice);
    auto* paper = new QVBoxLayout;
    paper->setSpacing(12);
    columns->addLayout(paper, 1);
    auto* previewHeader = new QHBoxLayout;
    previewHeader->addWidget(heading(QObject::tr("Preview")));
    auto* previewScale = new QComboBox;
    previewScale->setAccessibleName(QObject::tr("Preview scale"));
    previewScale->addItems({"200%", "100%", QObject::tr("Fit")});
    previewHeader->addWidget(previewScale);
    paper->addLayout(previewHeader);
    auto* preview = new StylePreview;
    preview->setObjectName("stylePreview");
    preview->setMinimumSize(230, 230);
    preview->setAccessibleName(QObject::tr("Live drawing style preview"));
    paper->addWidget(preview, 1);
    auto* previewNote = new QLabel(QObject::tr("100% uses one screen pixel per print point. Fit scales the preview only."));
    previewNote->setObjectName("secondary");
    previewNote->setWordWrap(true);
    paper->addWidget(previewNote);
    auto action = [&](const char* key, const QString& title) {
        auto* b = new QPushButton(title);
        b->setObjectName(key);
        b->setAutoDefault(false);
        return b;
    };
    auto* duplicate = action("duplicateStyle", QObject::tr("Duplicate"));
    auto* remove = action("deleteStyle", QObject::tr("Delete"));
    libraryActions->addWidget(duplicate);
    libraryActions->addWidget(remove);
    auto* actions = new QHBoxLayout;
    outer->addLayout(actions);
    auto* import = action("importStyle", QObject::tr("Import…"));
    auto* exportButton = action("exportStyle", QObject::tr("Export…"));
    actions->addWidget(import);
    actions->addWidget(exportButton);
    actions->addStretch();
    auto* save = action("saveStyle", QObject::tr("Save"));
    actions->addWidget(save);
    auto* buttons = new QDialogButtonBox(QDialogButtonBox::Close);
    buttons->button(QDialogButtonBox::Close)->setAutoDefault(false);
    actions->addWidget(buttons);
    QObject::connect(buttons, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);
    StyleSettings source;
    QString replacing;
    bool loading = false, editing = false;
    Document sample = chem::fromSmiles("[NH3]->[Pt+2](<-[NH3])(Cl)Cl").value_or(Document());
    sample.append(chem::fromSmiles("C[C@H](O)[C@@H](N)C(=O)O").value_or(Document()), {0, 55});
    sample.texts.push_back({{0, 100}, "R1 · conditions"});
    sample.circledCharges = document.circledCharges;
    auto draft = [&] {
        if (!editing) return source;
        StyleSettings s = source;
        s.name = name->text().trimmed();
        for (size_t i = 0; i < measures.size(); ++i)
            if (measures[i].member != &StyleSettings::bondSpacing) s.*measures[i].member = numbers[i]->value();
        s.bondSpacingAbsolute = gapMode->currentIndex() == 1 ? numbers[5]->value() : 0;
        if (gapMode->currentIndex() == 0) s.bondSpacing = numbers[5]->value() / 100;
        s.font = font->currentText().trimmed();
        s.weight = weight->currentData().toInt();
        s.italic = italic->isChecked();
        s.captionFont = captionFont->currentText() == QObject::tr("Same as labels") ? QString() : captionFont->currentText().trimmed();
        s.captionWeight = captionWeight->currentData().toInt();
        s.captionItalic = captionItalic->isChecked();
        s.chargesOnAtom = chargesOnAtom->isChecked();
        s.rGroupSuperscript = rGroupSuperscript->isChecked();
        s.hideImplicitHydrogens = hideHydrogens->isChecked();
        s.carbonLabels = carbons->currentData().toString();
        return s;
    };
    auto redraw = [&] {
        if (loading) return;
        const auto s = draft();
        if (!styleFromJson(styleToJson(s))) {
            preview->message = QObject::tr("Enter a name and valid measurements.");
            preview->update();
            return;
        }
        sample.style = s.name;
        sample.styleSettings = s;
        sample.carbonLabels = s.carbonLabels == "all" ? Document::CarbonLabels::All :
            s.carbonLabels == "terminal" ? Document::CarbonLabels::Terminal : Document::CarbonLabels::None;
        preview->document = sample;
        preview->message.clear();
        preview->update();
        fontNotice->setText(QFontDatabase::families().contains(s.font) && (s.captionFont.isEmpty() || QFontDatabase::families().contains(s.captionFont))
                            ? QString() : QObject::tr("This font is unavailable; a substitute is used."));
    };
    auto load = [&](const StyleSettings& s, bool editable, const QString& original) {
        loading = true;
        source = s;
        replacing = original;
        name->setText(s.name);
        gapMode->setCurrentIndex(s.bondSpacingAbsolute > 0 ? 1 : 0);
        numbers[5]->setSuffix(s.bondSpacingAbsolute > 0 ? QObject::tr(" pt") : QObject::tr("%"));
        for (size_t i = 0; i < measures.size(); ++i)
            numbers[i]->setValue(i == 5 ? (s.bondSpacingAbsolute > 0 ? s.bondSpacingAbsolute : s.bondSpacing * 100) : s.*measures[i].member);
        italic->setChecked(s.italic);
        captionFont->setCurrentText(s.captionFont.isEmpty() ? QObject::tr("Same as labels") : s.captionFont);
        if (captionWeight->findData(s.captionWeight) < 0) captionWeight->addItem(QString::number(s.captionWeight), s.captionWeight);
        captionWeight->setCurrentIndex(captionWeight->findData(s.captionWeight));
        captionItalic->setChecked(s.captionItalic);
        chargesOnAtom->setChecked(s.chargesOnAtom);
        rGroupSuperscript->setChecked(s.rGroupSuperscript);
        hideHydrogens->setChecked(s.hideImplicitHydrogens);
        carbons->setCurrentIndex(std::max(0, carbons->findData(s.carbonLabels)));
        font->setCurrentText(s.font);
        if (weight->findData(s.weight) < 0) weight->addItem(QString::number(s.weight), s.weight);
        weight->setCurrentIndex(weight->findData(s.weight));
        editing = editable;
        for (auto* box : numbers) box->setEnabled(editable);
        font->setEnabled(editable);
        weight->setEnabled(editable);
        for (QWidget* widget : std::initializer_list<QWidget*>{gapMode, italic, captionFont, captionWeight, captionItalic,
                                chargesOnAtom, rGroupSuperscript, hideHydrogens, carbons}) widget->setEnabled(editable);
        name->setReadOnly(!editable || !original.isEmpty());  // Duplicate supplies a new name; saved names stay stable.
        save->setEnabled(editable);
        save->setDefault(editable);
        duplicate->setDefault(!editable);
        state->setText(editable ? QObject::tr("Custom style") : QObject::tr("Built-in style · duplicate to edit"));
        remove->setEnabled(!original.isEmpty());
        loading = false;
        redraw();
    };
    auto refresh = [&](const QString& selected, bool includeDocument = false) {
        list->clear();
        for (const auto& s : drawingStyles()) {
            auto* item = new QListWidgetItem(s.name, list);
            item->setToolTip(s.name);
            item->setData(Qt::UserRole, styleToJson(styleSettings(s)));
        }
        if (includeDocument && document.styleSettings) {
            const auto& embedded = *document.styleSettings;
            const auto all = drawingStyles();
            const auto known = std::find_if(all.begin(), all.end(), [&](const auto& s) { return s.name == embedded.name; });
            if (known == all.end() || styleSettings(*known) != embedded) {
                auto* item = new QListWidgetItem(QObject::tr("%1 (document)").arg(embedded.name), list);
                auto draft = embedded;
                if (known != all.end()) draft.name += QObject::tr(" copy");
                item->setData(Qt::UserRole, styleToJson(draft));
                item->setData(Qt::UserRole + 1, true);
                list->setCurrentItem(item);
                return;
            }
        }
        const auto items = list->findItems(selected, Qt::MatchExactly);
        if (!items.isEmpty()) list->setCurrentItem(items.first());
        else list->setCurrentRow(0);
    };
    QObject::connect(list, &QListWidget::currentItemChanged, &dialog, [&](QListWidgetItem* item) {
        if (!item) return;
        const auto s = styleFromJson(item->data(Qt::UserRole).toJsonObject());
        if (!s) return;
        const bool embedded = item->data(Qt::UserRole + 1).toBool();
        const auto& builtins = builtInDrawingStyles();
        const bool builtin = !embedded && std::any_of(builtins.begin(), builtins.end(), [&](const auto& b) { return b.name == s->name; });
        load(*s, !builtin, builtin || embedded ? QString() : s->name);
    });
    QObject::connect(duplicate, &QPushButton::clicked, &dialog, [&] {
        auto s = draft();
        const QString base = s.name + QObject::tr(" copy");
        s.name = base;
        const auto all = drawingStyles();
        for (int n = 2; std::any_of(all.begin(), all.end(), [&](const auto& known) { return !known.name.compare(s.name, Qt::CaseInsensitive); }); ++n)
            s.name = base + " " + QString::number(n);
        load(s, true, {});
        name->setFocus();
        name->selectAll();
    });
    auto failed = [&](const QString& message) { QMessageBox::warning(&dialog, QObject::tr("Drawing Styles"), message); };
    QObject::connect(save, &QPushButton::clicked, &dialog, [&] {
        const auto s = draft();
        if (!saveDrawingStyle(s, replacing)) { failed(QObject::tr("Could not save the style. Check its name, measurements and folder permissions; names must be unique.")); return; }
        refresh(s.name);
    });
    QObject::connect(remove, &QPushButton::clicked, &dialog, [&] {
        if (QMessageBox::question(&dialog, QObject::tr("Delete style"), QObject::tr("Delete %1? Existing drawings keep their style.").arg(replacing)) != QMessageBox::Yes) return;
        if (!deleteDrawingStyle(replacing)) { failed(QObject::tr("Could not delete the style.")); return; }
        refresh({});
    });
    QObject::connect(import, &QPushButton::clicked, &dialog, [&] {
        const QString path = QFileDialog::getOpenFileName(&dialog, QObject::tr("Import drawing style"), {}, QObject::tr("Drawing styles (*.json)"));
        if (path.isEmpty()) return;
        QFile f(path);
        if (!f.open(QIODevice::ReadOnly) || f.size() > 64 * 1024) { failed(QObject::tr("Could not read the style file.")); return; }
        const auto s = styleFromJson(QJsonDocument::fromJson(f.readAll()).object());
        if (!s || !saveDrawingStyle(*s)) { failed(QObject::tr("Invalid style or duplicate name. Import requires a unique user style name.")); return; }
        refresh(s->name);
    });
    QObject::connect(exportButton, &QPushButton::clicked, &dialog, [&] {
        const auto s = draft();
        if (!styleFromJson(styleToJson(s))) { failed(QObject::tr("Enter a name and valid measurements.")); return; }
        const QString path = QFileDialog::getSaveFileName(&dialog, QObject::tr("Export drawing style"), {}, QObject::tr("Drawing styles (*.json)"));
        if (!path.isEmpty() && !writeWhole(path, QJsonDocument(styleToJson(s)).toJson())) failed(QObject::tr("Could not write the style file."));
    });
    QObject::connect(name, &QLineEdit::textChanged, &dialog, redraw);
    for (auto* box : numbers) QObject::connect(box, &QDoubleSpinBox::valueChanged, &dialog, redraw);
    QObject::connect(font, &QComboBox::editTextChanged, &dialog, redraw);
    QObject::connect(weight, &QComboBox::currentIndexChanged, &dialog, redraw);
    for (auto* box : {italic, captionItalic, chargesOnAtom, rGroupSuperscript, hideHydrogens})
        QObject::connect(box, &QCheckBox::toggled, &dialog, redraw);
    QObject::connect(captionFont, &QComboBox::editTextChanged, &dialog, redraw);
    QObject::connect(captionWeight, &QComboBox::currentIndexChanged, &dialog, redraw);
    QObject::connect(carbons, &QComboBox::currentIndexChanged, &dialog, redraw);
    QObject::connect(gapMode, &QComboBox::currentIndexChanged, &dialog, [&](int index) {
        if (loading) return;
        QSignalBlocker quiet(numbers[5]);
        const double length = numbers[0]->value();
        numbers[5]->setValue(index == 1 ? numbers[5]->value() * length / 100 : numbers[5]->value() / length * 100);
        numbers[5]->setSuffix(index == 1 ? QObject::tr(" pt") : QObject::tr("%"));
        redraw();
    });
    QObject::connect(previewScale, &QComboBox::currentIndexChanged, &dialog, [=](int index) {
        preview->fit = index == 2;
        preview->scale = index == 1 ? 1 : 2;
        preview->update();
    });
    refresh(documentStyle(document).name, true);
    dialog.exec();
}
