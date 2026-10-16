#pragma once
#include <QDoubleSpinBox>
#include <QLocale>

// Keep the full editable precision without padding every measurement with zeroes.
class NumberBox : public QDoubleSpinBox {
public:
    using QDoubleSpinBox::QDoubleSpinBox;

protected:
    QString textFromValue(double value) const override {
        QString text = QDoubleSpinBox::textFromValue(value);
        const QString decimal = locale().decimalPoint(), zero = locale().zeroDigit();
        if (text.contains(decimal)) {
            while (text.endsWith(zero)) text.chop(zero.size());
            if (text.endsWith(decimal)) text.chop(decimal.size());
        }
        return text;
    }
};
