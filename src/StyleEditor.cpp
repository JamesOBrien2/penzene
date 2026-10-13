#include "StyleEditor.h"
#include "Styles.h"
#include "Chem.h"
#include <QComboBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QDoubleSpinBox>
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
#include <QPushButton>
#include <QVBoxLayout>
#include <algorithm>

void editDrawingStyles(QWidget* parent, const Document& document) {
    QDialog dialog(parent);
    dialog.setObjectName("styleEditor");
    dialog.setWindowTitle(QObject::tr("Drawing Styles"));
    dialog.resize(720, 650);
    auto* outer = new QVBoxLayout(&dialog);
    auto* columns = new QHBoxLayout;
    outer->addLayout(columns);
    auto* list = new QListWidget;
    list->setObjectName("styleList");
    list->setAccessibleName(QObject::tr("Drawing styles"));
    columns->addWidget(list, 1);
    auto* right = new QVBoxLayout;
    columns->addLayout(right, 2);
    auto* fields = new QWidget;
    auto* form = new QFormLayout(fields);
    right->addWidget(fields);
    auto* name = new QLineEdit;
    name->setObjectName("styleName");
    form->addRow(QObject::tr("Name:"), name);
    struct Measure { const char* key; QString label; double StyleSettings::*member; double min, max; };
    const std::vector<Measure> measures{
        {"bondLength", QObject::tr("Bond length:"), &StyleSettings::bondLength, 0.1, 1000},
        {"lineWidth", QObject::tr("Line width:"), &StyleSettings::lineWidth, 0.001, 100},
        {"boldWidth", QObject::tr("Bold width:"), &StyleSettings::boldWidth, 0.001, 100},
        {"wedgeWidth", QObject::tr("Wedge width:"), &StyleSettings::wedgeWidth, 0.001, 100},
        {"hashSpacing", QObject::tr("Hash spacing:"), &StyleSettings::hashSpacing, 0.001, 100},
        {"doubleBondGap", QObject::tr("Double-bond gap:"), &StyleSettings::bondSpacing, 0, 1000},
        {"labelMargin", QObject::tr("Label margin:"), &StyleSettings::labelMargin, 0, 100},
        {"fontSize", QObject::tr("Label size:"), &StyleSettings::fontSize, 1, 200}};
    std::vector<QDoubleSpinBox*> numbers;
    for (const auto& m : measures) {
        auto* box = new QDoubleSpinBox;
        box->setObjectName(m.key);
        box->setDecimals(6);
        box->setRange(m.min, m.max);
        box->setSuffix(QObject::tr(" pt"));
        form->addRow(m.label, box);
        numbers.push_back(box);
    }
    auto* font = new QComboBox;
    font->setEditable(true);
    font->addItems(QFontDatabase::families());
    font->setObjectName("styleFont");
    form->addRow(QObject::tr("Label font:"), font);
    auto* weight = new QComboBox;
    weight->setObjectName("styleWeight");
    for (const auto& pair : {std::pair{QObject::tr("Thin"), 100}, {QObject::tr("Extra light"), 200},
         {QObject::tr("Light"), 300}, {QObject::tr("Normal"), 400}, {QObject::tr("Medium"), 500},
         {QObject::tr("Demi bold"), 600}, {QObject::tr("Bold"), 700}, {QObject::tr("Extra bold"), 800}, {QObject::tr("Black"), 900}})
        weight->addItem(pair.first, pair.second);
    form->addRow(QObject::tr("Label weight:"), weight);
    auto* preview = new QLabel;
    preview->setObjectName("stylePreview");
    preview->setMinimumSize(300, 135);
    preview->setAlignment(Qt::AlignCenter);
    preview->setAccessibleName(QObject::tr("Live drawing style preview"));
    right->addWidget(preview);
    auto* fontNotice = new QLabel;
    fontNotice->setWordWrap(true);
    right->addWidget(fontNotice);
    auto* actions = new QHBoxLayout;
    outer->addLayout(actions);
    auto action = [&](const char* key, const QString& title) {
        auto* b = new QPushButton(title);
        b->setObjectName(key);
        actions->addWidget(b);
        return b;
    };
    auto* duplicate = action("duplicateStyle", QObject::tr("Duplicate"));
    auto* save = action("saveStyle", QObject::tr("Save"));
    auto* remove = action("deleteStyle", QObject::tr("Delete"));
    auto* import = action("importStyle", QObject::tr("Import…"));
    auto* exportButton = action("exportStyle", QObject::tr("Export…"));
    auto* buttons = new QDialogButtonBox(QDialogButtonBox::Close);
    outer->addWidget(buttons);
    QObject::connect(buttons, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);
    StyleSettings source;
    QString replacing;
    bool loading = false;
    Document sample = chem::fromSmiles("[NH3]->[Pt+2](<-[NH3])(Cl)Cl").value_or(Document());
    auto draft = [&] {
        if (!fields->isEnabled()) return source;
        StyleSettings s = source;
        s.name = name->text().trimmed();
        for (size_t i = 0; i < measures.size(); ++i)
            if (measures[i].member != &StyleSettings::bondSpacing) s.*measures[i].member = numbers[i]->value();
        s.bondSpacing = numbers[5]->value() / s.bondLength;
        s.labelRadius = std::max(0.0, source.labelRadius + 0.39 * (s.fontSize - source.fontSize) + s.labelMargin - source.labelMargin);
        s.font = font->currentText().trimmed();
        s.weight = weight->currentData().toInt();
        return s;
    };
    auto redraw = [&] {
        if (loading) return;
        const auto s = draft();
        if (!styleFromJson(styleToJson(s))) { preview->setText(QObject::tr("Enter a name and valid measurements.")); return; }
        sample.style = s.name;
        sample.styleSettings = s;
        const QImage image = renderImage(sample, {144, Qt::white});
        preview->setPixmap(QPixmap::fromImage(image).scaled(preview->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));
        fontNotice->setText(QFontDatabase::families().contains(s.font) ? QString() : QObject::tr("This font is unavailable; a substitute is used."));
    };
    auto load = [&](const StyleSettings& s, bool editable, const QString& original) {
        loading = true;
        source = s;
        replacing = original;
        name->setText(s.name);
        for (size_t i = 0; i < measures.size(); ++i)
            numbers[i]->setValue(measures[i].member == &StyleSettings::bondSpacing ? s.bondSpacing * s.bondLength : s.*measures[i].member);
        font->setCurrentText(s.font);
        if (weight->findData(s.weight) < 0) weight->addItem(QString::number(s.weight), s.weight);
        weight->setCurrentIndex(weight->findData(s.weight));
        fields->setEnabled(editable);
        name->setReadOnly(!original.isEmpty());  // Duplicate supplies a new name; saved names stay stable.
        save->setEnabled(editable);
        remove->setEnabled(!original.isEmpty());
        loading = false;
        redraw();
    };
    auto refresh = [&](const QString& selected) {
        list->clear();
        for (const auto& s : drawingStyles()) list->addItem(s.name);
        const auto items = list->findItems(selected, Qt::MatchExactly);
        if (!items.isEmpty()) list->setCurrentItem(items.first());
        else list->setCurrentRow(0);
    };
    QObject::connect(list, &QListWidget::currentTextChanged, &dialog, [&](const QString& text) {
        if (text.isEmpty()) return;
        const auto s = styleSettings(drawingStyle(text));
        const auto& builtins = builtInDrawingStyles();
        const bool builtin = std::any_of(builtins.begin(), builtins.end(), [&](const auto& b) { return b.name == text; });
        load(s, !builtin, builtin ? QString() : text);
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
    refresh(documentStyle(document).name);
    dialog.exec();
}
