#include "Stoichiometry.h"
#include "NumberBox.h"

#include <QApplication>
#include <QClipboard>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QFormLayout>
#include <QHeaderView>
#include <QLabel>
#include <QPushButton>
#include <QSignalBlocker>
#include <QTableWidget>
#include <QVBoxLayout>
#include <algorithm>
#include <cmath>

StoichiometryPanel::StoichiometryPanel(QWidget* parent) : QWidget(parent) {
    auto* layout = new QVBoxLayout(this);
    layout->setContentsMargins(12, 12, 12, 12);
    layout->setSpacing(10);
    auto* form = new QFormLayout;
    reaction_ = new QComboBox;
    reaction_->setObjectName("stoichReaction");
    limiting_ = new QComboBox;
    limiting_->setObjectName("stoichLimiting");
    limiting_->setPlaceholderText(tr("Choose a limiting reagent"));
    scale_ = new NumberBox;
    scale_->setObjectName("stoichScale");
    scale_->setDecimals(6);
    scale_->setRange(0, 1e9);
    scale_->setValue(1);
    form->addRow(tr("Reaction:"), reaction_);
    form->addRow(tr("Limiting reagent:"), limiting_);
    form->addRow(tr("Scale (mmol):"), scale_);
    layout->addLayout(form);
    table_ = new QTableWidget(0, 8);
    table_->setObjectName("stoichTable");
    table_->setAccessibleName(tr("Reaction quantities"));
    table_->setHorizontalHeaderLabels({tr("Component"), tr("Formula"), tr("MW (g/mol)"), tr("Equiv."),
                                      tr("mmol"), tr("Mass (mg)"), tr("Density (g/mL)"), tr("Volume (µL)")});
    table_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    table_->setShowGrid(false);
    table_->verticalHeader()->hide();
    table_->verticalHeader()->setDefaultSectionSize(38);
    table_->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    table_->horizontalHeader()->setStretchLastSection(true);
    layout->addWidget(table_);
    message_ = new QLabel;
    message_->setWordWrap(true);
    layout->addWidget(message_);
    auto* note = new QLabel(tr("Quantities are kept for this session. Copy as Text to keep a table in your lab notebook."));
    note->setWordWrap(true);
    note->setObjectName("secondary");
    layout->addWidget(note);
    copy_ = new QPushButton(tr("Copy as Text"));
    copy_->setObjectName("stoichCopy");
    copy_->setEnabled(false);
    layout->addWidget(copy_, 0, Qt::AlignRight);
    connect(copy_, &QPushButton::clicked, this, &StoichiometryPanel::copyText);
    connect(reaction_, &QComboBox::currentIndexChanged, this, [this](int i) {
        if (!page_ || i < 0) return;
        pages_[page_].reaction = reactionKeys_[i];
        rebuild();
    });
    connect(limiting_, &QComboBox::currentIndexChanged, this, [this](int i) {
        if (!page_ || i < 0) return;
        const int row = limiting_->itemData(i).toInt();
        currentStep().limiting = components_[row].key;
        calculate();
    });
    connect(scale_, &QDoubleSpinBox::valueChanged, this, [this](double value) {
        if (!page_ || reaction_->currentIndex() < 0) return;
        currentStep().scale = value;
        calculate();
    });
}

void StoichiometryPanel::reset() {
    for (auto& page : pages_) page = Page{};
}

void StoichiometryPanel::setDocument(const Document& doc, QUndoStack* page) {
    if (!pages_.contains(page)) {
        pages_[page] = Page{};
        connect(page, &QObject::destroyed, this, [this, page] {
            pages_.remove(page);
            if (!page_) {
                reaction_->clear();
                reaction_->setEnabled(false);
                rebuild();
            }
        });
    }
    page_ = page;
    reactions_ = chem::reactionsOf(doc);
    reactionKeys_.clear();
    QHash<QString, int> occurrences;
    for (auto r : reactions_) {
        r.products.clear();  // quantities follow the starting materials, not edits to a product
        const QString identity = QString::fromStdString(chem::toReactionSmiles(r));
        reactionKeys_.push_back(identity + ":" + QString::number(occurrences[identity]++));
    }
    QSignalBlocker quiet(reaction_);
    reaction_->clear();
    for (int i = 0; i < int(reactions_.size()); ++i) reaction_->addItem(tr("Reaction %1").arg(i + 1));
    const auto found = std::find(reactionKeys_.begin(), reactionKeys_.end(), pages_[page].reaction);
    const int step = found == reactionKeys_.end() ? 0 : int(found - reactionKeys_.begin());
    reaction_->setCurrentIndex(step);
    if (!reactionKeys_.empty()) pages_[page].reaction = reactionKeys_[step];
    reaction_->setEnabled(reactions_.size() > 1);
    rebuild();
}

StoichiometryPanel::Step& StoichiometryPanel::currentStep() {
    return pages_[page_].steps[reactionKeys_[reaction_->currentIndex()]];
}

void StoichiometryPanel::rebuild() {
    QSignalBlocker quietLimit(limiting_), quietScale(scale_);
    table_->setRowCount(0);
    components_.clear();
    limiting_->clear();
    const int index = reaction_->currentIndex();
    if (!page_ || index < 0) {
        scale_->setEnabled(false);
        limiting_->setEnabled(false);
        copy_->setEnabled(false);
        message_->setText(tr("Draw a reaction arrow to identify reactants and reagents."));
        return;
    }
    Step& state = currentStep();
    scale_->setValue(state.scale);
    const auto& reaction = reactions_[index];
    QHash<QString, int> occurrences;
    for (const auto& [role, molecules] : {std::pair{tr("Reactant"), &reaction.reactants}, {tr("Reagent"), &reaction.agents}}) {
        for (int i = 0; i < int(molecules->size()); ++i) {
            const Document& mol = (*molecules)[i];
            const QString identity = role + ":" + QString::fromStdString(chem::toSmiles(mol));
            const QString key = identity + ":" + QString::number(occurrences[identity]++);
            const auto properties = chem::properties(mol);
            const bool valid = properties && properties->formula.find('*') == std::string::npos &&
                               std::isfinite(properties->mw) && properties->mw > 0;
            const int row = table_->rowCount();
            table_->insertRow(row);
            components_.push_back({key, valid ? properties->mw : 0});
            const QString name = tr("%1 %2").arg(role).arg(i + 1);
            for (int col = 0; col < table_->columnCount(); ++col) {
                auto* item = new QTableWidgetItem;
                if (col >= 2) item->setTextAlignment(Qt::AlignRight | Qt::AlignVCenter);
                table_->setItem(row, col, item);
            }
            table_->item(row, 0)->setText(name);
            table_->item(row, 1)->setText(valid ? QString::fromStdString(properties->formula) : tr("Unknown structure"));
            table_->item(row, 2)->setText(valid ? QString::number(properties->mw, 'g', 8) : "—");
            if (valid) limiting_->addItem(name + " · " + QString::fromStdString(properties->formula), row);
            Inputs& input = state.inputs[key];
            auto* eq = new NumberBox;
            eq->setDecimals(6);
            eq->setRange(0, 1e6);
            eq->setValue(valid ? input.equivalents : 0);
            eq->setAccessibleName(tr("Equivalents for %1").arg(name));
            auto* density = new NumberBox;
            density->setDecimals(6);
            density->setRange(0, 100);
            density->setSpecialValueText("—");
            density->setValue(input.density);
            density->setEnabled(valid);
            density->setAccessibleName(tr("Density for %1, grams per millilitre").arg(name));
            table_->setCellWidget(row, 3, eq);
            table_->setCellWidget(row, 6, density);
            connect(eq, &QDoubleSpinBox::valueChanged, this, [this, key](double value) {
                currentStep().inputs[key].equivalents = value;
                calculate();
            });
            connect(density, &QDoubleSpinBox::valueChanged, this, [this, key](double value) {
                currentStep().inputs[key].density = value;
                calculate();
            });
        }
    }
    int selected = state.limiting.isEmpty() ? 0 : -1;
    for (int i = 0; i < limiting_->count(); ++i)
        if (components_[limiting_->itemData(i).toInt()].key == state.limiting) selected = i;
    limiting_->setCurrentIndex(selected);
    limiting_->setEnabled(limiting_->count() > 0);
    if (selected >= 0 && selected < limiting_->count()) state.limiting = components_[limiting_->currentData().toInt()].key;
    calculate();
}

void StoichiometryPanel::calculate() {
    if (!page_ || reaction_->currentIndex() < 0) return;
    Step& state = currentStep();
    const bool haveLimit = limiting_->currentIndex() >= 0;
    const int limit = haveLimit ? limiting_->currentData().toInt() : -1;
    scale_->setEnabled(haveLimit);
    copy_->setEnabled(haveLimit);
    bool unknown = false;
    for (int row = 0; row < int(components_.size()); ++row) {
        const auto& component = components_[row];
        Inputs& input = state.inputs[component.key];
        auto* eq = qobject_cast<QDoubleSpinBox*>(table_->cellWidget(row, 3));
        eq->setEnabled(component.mw > 0 && row != limit);
        if (row == limit) {
            QSignalBlocker quiet(eq);
            input.equivalents = 1;
            eq->setValue(1);
        }
        unknown |= component.mw <= 0;
        const double mmol = scale_->value() * input.equivalents;
        const double mass = mmol * component.mw;  // mmol × g/mol = mg; mg ÷ g/mL = µL
        const double volume = input.density > 0 ? mass / input.density : 0;
        for (auto [col, value] : {std::pair{4, mmol}, {5, mass}, {7, volume}})
            table_->item(row, col)->setText(haveLimit && component.mw > 0 && std::isfinite(value) &&
                                           (col != 7 || input.density > 0) ? QString::number(value, 'g', 8) : "—");
    }
    message_->setText(!haveLimit ? tr("Choose a limiting reagent with a complete molecular formula.") :
                     unknown ? tr("Unknown structures have no calculated quantities. Check or complete their labels.") : QString());
}

void StoichiometryPanel::copyText() {
    if (!copy_->isEnabled()) return;
    QStringList lines{tr("Reaction %1; limiting reagent: %2; scale: %3 mmol")
                          .arg(reaction_->currentIndex() + 1).arg(limiting_->currentText()).arg(scale_->value(), 0, 'g', 8)};
    QStringList headers;
    for (int col = 0; col < table_->columnCount(); ++col) headers << table_->horizontalHeaderItem(col)->text();
    lines << headers.join('\t');
    for (int row = 0; row < table_->rowCount(); ++row) {
        QStringList values;
        for (int col = 0; col < table_->columnCount(); ++col) {
            if (col == 3 || col == 6) {
                const double value = qobject_cast<QDoubleSpinBox*>(table_->cellWidget(row, col))->value();
                values << (components_[row].mw > 0 && (col != 6 || value > 0) ? QString::number(value, 'g', 8) : "—");
            } else values << table_->item(row, col)->text();
        }
        lines << values.join('\t');
    }
    QApplication::clipboard()->setText(lines.join('\n'));
}
