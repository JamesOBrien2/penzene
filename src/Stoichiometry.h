#pragma once
#include "Chem.h"

#include <QHash>
#include <QPointer>
#include <QUndoStack>
#include <QWidget>

class QComboBox;
class QDoubleSpinBox;
class QLabel;
class QPushButton;
class QTableWidget;

class StoichiometryPanel : public QWidget {
    Q_OBJECT
public:
    explicit StoichiometryPanel(QWidget* parent = nullptr);
    void setDocument(const Document& doc, QUndoStack* page);
    void reset();  // new file: the first page's undo stack is reused

private:
    struct Inputs { double equivalents = 1, density = 0; };
    struct Step {
        double scale = 1;
        QString limiting;
        QHash<QString, Inputs> inputs;
    };
    struct Page { QString reaction; QHash<QString, Step> steps; };
    struct Component { QString key; double mw = 0; };
    QHash<QUndoStack*, Page> pages_;
    QPointer<QUndoStack> page_;
    std::vector<chem::Reaction> reactions_;
    std::vector<QString> reactionKeys_;
    std::vector<Component> components_;
    QComboBox *reaction_, *limiting_;
    QDoubleSpinBox* scale_;
    QTableWidget* table_;
    QLabel* message_;
    QPushButton* copy_;
    void rebuild();
    Step& currentStep();
    void calculate();
    void copyText();
};
