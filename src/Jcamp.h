#pragma once
#include <QByteArray>
#include <QPointF>
#include <QString>
#include <vector>

struct JcampSpectrum {
    QString title;
    bool proton = false, peaks = false;
    std::vector<QPointF> points;  // ppm, intensity; NaN intensity marks a missing sample
};
// Processed 1D 1H/13C spectra, including LINK containers and NTUPLES real data.
// An invalid processed block fails the import; the caller keeps its previous spectrum.
std::vector<JcampSpectrum> readJcamp(const QByteArray& data, QString* error = nullptr);
