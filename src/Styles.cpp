#include "Styles.h"
#include <QDir>
#include <QCryptographicHash>
#include <QFile>
#include <QJsonDocument>
#include <QLockFile>
#include <QStandardPaths>
#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>

StyleSettings styleSettings(const DrawingStyle& s) {
    const double k = s.bondLength / kBondLength;
    return {s.name, s.bondLength, s.lineWidth * k, s.boldWidth * k, s.wedgeWidth * k, s.hashSpacing * k,
            s.bondSpacing, s.labelRadius * k, s.labelMargin * k, s.font, int(s.weight), s.fontSize * k};
}

DrawingStyle modelStyle(const StyleSettings& s) {
    const double k = kBondLength / s.bondLength;
    return {s.name, s.bondLength, s.lineWidth * k, s.boldWidth * k, s.wedgeWidth * k, s.hashSpacing * k,
            s.bondSpacing, s.labelRadius * k, s.font, QFont::Weight(s.weight), s.fontSize * k, s.labelMargin * k};
}

QJsonObject styleToJson(const StyleSettings& s) {
    return {{"format", "penzene-drawing-style"}, {"version", 1}, {"name", s.name},
            {"bondLength", s.bondLength}, {"lineWidth", s.lineWidth}, {"boldWidth", s.boldWidth},
            {"wedgeWidth", s.wedgeWidth}, {"hashSpacing", s.hashSpacing}, {"bondSpacing", s.bondSpacing},
            {"labelRadius", s.labelRadius}, {"labelMargin", s.labelMargin},
            {"font", s.font}, {"weight", s.weight}, {"fontSize", s.fontSize}};
}

std::optional<StyleSettings> styleFromJson(const QJsonObject& o) {
    if (o["format"].toString() != "penzene-drawing-style" || o["version"] != 1) return std::nullopt;
    StyleSettings s;
    s.name = o["name"].toString();
    s.font = o["font"].toString();
    if (s.name.isEmpty() || s.name != s.name.trimmed() || s.name.size() > 100 || s.font.trimmed().isEmpty() || s.font.size() > 200) return std::nullopt;
    for (QChar c : s.name) if (c.category() == QChar::Other_Control) return std::nullopt;
    for (auto [key, field, minimum, maximum] : {
        std::tuple{"bondLength", &s.bondLength, 0.1, 1000.0}, {"lineWidth", &s.lineWidth, 0.001, 100.0},
        {"boldWidth", &s.boldWidth, 0.001, 100.0}, {"wedgeWidth", &s.wedgeWidth, 0.001, 100.0},
        {"hashSpacing", &s.hashSpacing, 0.001, 100.0}, {"bondSpacing", &s.bondSpacing, 0.0, 1.0},
        {"labelRadius", &s.labelRadius, 0.0, 1000.0}, {"labelMargin", &s.labelMargin, 0.0, 100.0},
        {"fontSize", &s.fontSize, 1.0, 200.0}}) {
        const QJsonValue v = o[key];
        const double n = v.toDouble(std::numeric_limits<double>::quiet_NaN());
        if (!v.isDouble() || !std::isfinite(n) || n < minimum || n > maximum) return std::nullopt;
        *field = n;
    }
    const double weight = o["weight"].toDouble();
    if (!o["weight"].isDouble() || weight < 1 || weight > 1000 || std::floor(weight) != weight) return std::nullopt;
    s.weight = int(weight);
    return s;
}

QString stylesDirectory() {
    // Independent of the host's Qt application name (CLI and Python may use another).
    return QStandardPaths::writableLocation(QStandardPaths::GenericDataLocation) + "/Penzene/styles";
}

static bool builtInName(const QString& name) {
    const auto& builtins = builtInDrawingStyles();
    return std::any_of(builtins.begin(), builtins.end(), [&](const auto& s) { return s.name.compare(name, Qt::CaseInsensitive) == 0; });
}

static QString stylePath(const QString& name) {
    QDir dir(stylesDirectory());
    for (const auto& path : dir.entryList({"*.json"}, QDir::Files | QDir::NoSymLinks, QDir::Name)) {
        QFile f(dir.filePath(path));
        if (!f.open(QIODevice::ReadOnly) || f.size() > 64 * 1024) continue;
        const auto s = styleFromJson(QJsonDocument::fromJson(f.readAll()).object());
        if (s && !s->name.compare(name, Qt::CaseInsensitive)) return dir.filePath(path);
    }
    // Bounded filenames, even for long Unicode names or names containing path separators.
    const auto hash = QCryptographicHash::hash(name.toCaseFolded().toUtf8(), QCryptographicHash::Sha256).toHex();
    return dir.filePath(QString::fromLatin1(hash) + ".json");
}

std::vector<StyleSettings> userDrawingStyles() {
    std::vector<StyleSettings> out;
    QDir dir(stylesDirectory());
    for (const auto& path : dir.entryList({"*.json"}, QDir::Files | QDir::NoSymLinks, QDir::Name)) {
        QFile f(dir.filePath(path));
        if (!f.open(QIODevice::ReadOnly) || f.size() > 64 * 1024) continue;
        const auto s = styleFromJson(QJsonDocument::fromJson(f.readAll()).object());
        if (!s || builtInName(s->name) || std::any_of(out.begin(), out.end(), [&](const auto& other) {
            return !other.name.compare(s->name, Qt::CaseInsensitive);
        })) continue;
        out.push_back(*s);
    }
    std::sort(out.begin(), out.end(), [](const auto& a, const auto& b) { return a.name < b.name; });
    return out;
}

bool saveDrawingStyle(const StyleSettings& s, const QString& replacing) {
    if (!styleFromJson(styleToJson(s)) || builtInName(s.name) || (!replacing.isEmpty() && s.name != replacing)) return false;
    if (!QDir().mkpath(stylesDirectory())) return false;
    QLockFile lock(stylesDirectory() + "/.lock");
    if (!lock.tryLock()) return false;
    const QString path = stylePath(s.name);
    for (const auto& known : userDrawingStyles())
        if (!known.name.compare(s.name, Qt::CaseInsensitive) && known.name != replacing) return false;
    if (replacing.isEmpty() && QFile::exists(path)) return false;
    if (!replacing.isEmpty() && !QFile::exists(path)) return false;
    return writeWhole(path, QJsonDocument(styleToJson(s)).toJson());
}

bool deleteDrawingStyle(const QString& name) {
    if (builtInName(name)) return false;
    QLockFile lock(stylesDirectory() + "/.lock");
    return lock.tryLock() && QFile::remove(stylePath(name));
}

bool applyDrawingStyle(Document& doc, const QString& name) {
    for (const auto& s : builtInDrawingStyles()) {
        if (s.name != name && !(name.isEmpty() && s.name == "ACS 1996")) continue;
        doc.style = s.name == "ACS 1996" ? QString() : s.name;
        doc.styleSettings.reset();
        doc.labelRatio = 0;
        return true;
    }
    for (const auto& s : userDrawingStyles()) {
        if (s.name != name) continue;
        doc.style = s.name;
        doc.styleSettings = s;
        doc.labelRatio = 0;
        return true;
    }
    return false;
}
