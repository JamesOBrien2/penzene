#include "Jcamp.h"
#include <QHash>
#include <QCoreApplication>
#include <QRegularExpression>
#include <QStringList>
#include <algorithm>
#include <cmath>
#include <limits>
#include <new>

namespace {
constexpr size_t maxSamples = 10'000'000;  // bound corrupt compressed repeat counts before allocating
[[noreturn]] void fail(const char* message) { throw QCoreApplication::translate("Jcamp", message); }
QString compact(QString s) { return s.remove(QRegularExpression("[\\s_-]")).toUpper(); }
double number(QString s) {
    bool ok = false;
    const double value = s.trimmed().toDouble(&ok);
    if (!ok || !std::isfinite(value)) fail(QT_TRANSLATE_NOOP("Jcamp", "Invalid numeric value in the spectrum."));
    return value;
}
size_t count(QString s) {
    const double n = number(s);
    if (n < 1 || n > maxSamples || n != std::floor(n)) fail(QT_TRANSLATE_NOOP("Jcamp", "Invalid spectrum size (maximum 10 million samples)."));
    return size_t(n);
}
enum class Kind { Absolute, Difference, Repeat, Missing };
struct Token { Kind kind; double value = 0; };

std::vector<Token> tokens(const QString& line) {
    std::vector<Token> out;
    // Scientific notation is a whole AFFN word; E/e inside a packed word is SQZ.
    static const QRegularExpression numeric("^[+-]?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)");
    for (const QString& word : line.split(QRegularExpression("[\\s,;]+"), Qt::SkipEmptyParts)) {
        if (out.size() > maxSamples) fail(QT_TRANSLATE_NOOP("Jcamp", "The spectrum contains too many samples."));
        bool ok = false;
        const double plain = word.toDouble(&ok);
        if (ok && std::isfinite(plain)) { out.push_back({Kind::Absolute, plain}); continue; }
        for (qsizetype i = 0; i < word.size();) {
            if (out.size() > maxSamples) fail(QT_TRANSLATE_NOOP("Jcamp", "The spectrum contains too many samples."));
            const QChar c = word[i];
            if (c == '?') { out.push_back({Kind::Missing}); ++i; continue; }
            Kind kind = Kind::Absolute;
            int leading = 0, sign = 1;
            if (c == '@') leading = 0;
            else if (c >= 'A' && c <= 'I') leading = c.unicode() - 'A' + 1;
            else if (c >= 'a' && c <= 'i') leading = c.unicode() - 'a' + 1, sign = -1;
            else if (c == '%') kind = Kind::Difference;
            else if (c >= 'J' && c <= 'R') kind = Kind::Difference, leading = c.unicode() - 'J' + 1;
            else if (c >= 'j' && c <= 'r') kind = Kind::Difference, leading = c.unicode() - 'j' + 1, sign = -1;
            else if ((c >= 'S' && c <= 'Z') || c == 's') kind = Kind::Repeat, leading = c == 's' ? 9 : c.unicode() - 'S' + 1;
            else {
                const auto match = numeric.match(word.mid(i));
                if (!match.hasMatch()) fail(QT_TRANSLATE_NOOP("Jcamp", "Invalid encoded value in the spectrum."));
                out.push_back({Kind::Absolute, number(match.captured())});
                i += match.capturedLength();
                continue;
            }
            ++i;
            double value = leading;
            while (i < word.size() && word[i] >= '0' && word[i] <= '9') value = value * 10 + word[i++].digitValue();
            if (!std::isfinite(value)) fail(QT_TRANSLATE_NOOP("Jcamp", "Encoded spectrum value is too large."));
            out.push_back({kind, sign * value});
        }
    }
    return out;
}

struct Table { QString key, form; QStringList lines; };
struct Block { QHash<QString, QString> headers; std::vector<Table> tables; QString lastKey; int table = -1; };

JcampSpectrum spectrum(const Block& block, const Table& table) {
    const auto& h = block.headers;
    auto scalar = [&](const char* key, double fallback) { return h.contains(key) ? number(h[key]) : fallback; };
    const QString form = compact(table.form);
    const bool incremental = form.contains("++");
    static const QRegularExpression increments("^\\(([A-Z][A-Z0-9]*)\\+\\+\\(([A-Z][A-Z0-9]*)\\.\\.\\2\\)\\)");
    const auto variables = increments.match(form);
    if ((incremental && !variables.hasMatch()) || (!incremental &&
        !form.startsWith("(XY..XY)") && !form.startsWith("(X,Y..X,Y)") && !form.startsWith("(XR..XR)")))
        fail(QT_TRANSLATE_NOOP("Jcamp", "Unsupported spectrum data table."));
    JcampSpectrum out;
    out.title = h.value("TITLE");
    const QString nucleus = compact(h.value(".OBSERVENUCLEUS", h.value("$NUC1"))).remove('^').remove('<').remove('>');
    if (nucleus != "1H" && nucleus != "13C") fail(QT_TRANSLATE_NOOP("Jcamp", "The spectrum must identify its nucleus as 1H or 13C."));
    out.proton = nucleus == "1H";
    out.peaks = table.key == "PEAKTABLE" || form.contains("PEAKS");
    QString units = compact(h.value("XUNITS"));
    double first = scalar("FIRSTX", 0), last = scalar("LASTX", 0);
    double xf = scalar("XFACTOR", 1), yf = scalar("YFACTOR", 1);
    size_t n = h.contains("NPOINTS") ? count(h["NPOINTS"]) : 0;
    if (table.key == "DATATABLE") {
        const auto symbols = compact(h.value("SYMBOL")).split(',');
        const int x = int(symbols.indexOf(incremental ? variables.captured(1) : "X"));
        const int y = int(symbols.indexOf(incremental ? variables.captured(2) : form.startsWith("(XR..XR)") ? "R" : "Y"));
        if (x < 0 || y < 0) fail(QT_TRANSLATE_NOOP("Jcamp", "Missing NTUPLES axis or real intensity variable."));
        auto variable = [&](const char* key, int at) {
            const auto values = h.value(key).split(',');
            if (at >= values.size()) fail(QT_TRANSLATE_NOOP("Jcamp", "Incomplete NTUPLES variable declaration."));
            return values[at].trimmed();
        };
        units = compact(variable("UNITS", x));
        first = number(variable("FIRST", x)), last = number(variable("LAST", x));
        xf = number(variable("FACTOR", x)), yf = number(variable("FACTOR", y));
        n = count(variable("VARDIM", x));
        if (count(variable("VARDIM", y)) != n) fail(QT_TRANSLATE_NOOP("Jcamp", "NTUPLES sample counts disagree."));
    } else if (incremental && (!h.contains("FIRSTX") || !h.contains("LASTX"))) {
        fail(QT_TRANSLATE_NOOP("Jcamp", "Incremental data needs FIRSTX and LASTX."));
    }
    if (incremental && table.key != "DATATABLE" && (variables.captured(1) != "X" ||
        (variables.captured(2) != "Y" && variables.captured(2) != "R")))
        fail(QT_TRANSLATE_NOOP("Jcamp", "Unsupported spectrum data table."));
    if (!xf || !yf) fail(QT_TRANSLATE_NOOP("Jcamp", "Spectrum scale factors must be nonzero."));
    if (units != "PPM" && units != "HZ" && units != "HERTZ") fail(QT_TRANSLATE_NOOP("Jcamp", "The spectrum X axis must use ppm or Hz."));
    const double frequency = h.contains(".OBSERVEFREQUENCY") ? number(h[".OBSERVEFREQUENCY"]) : scalar("$SF", 0);
    const bool hz = units != "PPM";
    if (hz && frequency <= 0) fail(QT_TRANSLATE_NOOP("Jcamp", "A spectrum in Hz needs a positive observe frequency in MHz."));
    if (incremental && (!n || (n > 1 && first == last) || (n == 1 && first != last))) fail(QT_TRANSLATE_NOOP("Jcamp", "Incremental data needs a valid sample count and X range."));
    const double dx = n > 1 ? (last - first) / double(n - 1) : 0;
    double previous = 0, previousDelta = 0;
    Kind previousKind = Kind::Missing;
    bool checkpoint = false;
    auto append = [&](double x, double y) {
        if (out.points.size() >= maxSamples || (n && out.points.size() >= n)) fail(QT_TRANSLATE_NOOP("Jcamp", "The spectrum contains too many samples."));
        if (!std::isfinite(x) || (!std::isnan(y) && !std::isfinite(y))) fail(QT_TRANSLATE_NOOP("Jcamp", "Invalid scaled spectrum value."));
        out.points.emplace_back(x, y);
    };
    for (const QString& line : table.lines) {
        auto values = tokens(line);
        if (values.empty()) continue;
        if (!incremental) {
            if (values.size() % 2) fail(QT_TRANSLATE_NOOP("Jcamp", "An XY data pair is incomplete."));
            for (size_t i = 0; i < values.size(); i += 2) {
                if (values[i].kind != Kind::Absolute || (values[i + 1].kind != Kind::Absolute && values[i + 1].kind != Kind::Missing))
                    fail(QT_TRANSLATE_NOOP("Jcamp", "XY pairs must contain absolute values."));
                append(values[i].value * xf, values[i + 1].kind == Kind::Missing ? std::numeric_limits<double>::quiet_NaN() : values[i + 1].value * yf);
            }
            continue;
        }
        if (values.size() < 2 || values[0].kind != Kind::Absolute) fail(QT_TRANSLATE_NOOP("Jcamp", "An incremental data line needs an X checkpoint and intensities."));
        const size_t at = out.points.size() - (checkpoint ? 1 : 0);
        const double anchor = values[0].value * xf, expected = first + double(at) * dx;
        if (std::abs(anchor - expected) > std::max(std::abs(xf) * 0.51, std::abs(dx) * 0.01) + 1e-7 * std::max(1.0, std::abs(expected)))
            fail(QT_TRANSLATE_NOOP("Jcamp", "Spectrum X checkpoints disagree with the sample sequence."));
        size_t begin = 1;
        if (checkpoint) {
            if (values[1].kind != Kind::Absolute || std::abs(values[1].value - previous) > 1e-8 * std::max(1.0, std::abs(previous)))
                fail(QT_TRANSLATE_NOOP("Jcamp", "Spectrum difference checkpoint is corrupt."));
            previous = values[1].value;
            begin = 2;
        }
        previousKind = Kind::Absolute;
        for (size_t i = begin; i < values.size(); ++i) {
            const Token t = values[i];
            if (t.kind == Kind::Repeat) {
                if (i == 1 || t.value < 1 || t.value != std::floor(t.value) || t.value > maxSamples || values[i - 1].kind == Kind::Repeat)
                    fail(QT_TRANSLATE_NOOP("Jcamp", "Invalid compressed repeat count."));
                const size_t repeats = size_t(t.value) - 1;
                if (repeats > maxSamples - out.points.size() || (n && repeats > n - out.points.size())) fail(QT_TRANSLATE_NOOP("Jcamp", "Compressed spectrum exceeds its declared size."));
                for (size_t k = 0; k < repeats; ++k) {
                    if (previousKind == Kind::Difference) previous += previousDelta;
                    append(first + double(out.points.size()) * dx, previous * yf);
                }
                continue;
            }
            if (t.kind == Kind::Missing) previous = std::numeric_limits<double>::quiet_NaN();
            else if (t.kind == Kind::Difference) {
                if (!std::isfinite(previous) || i == 1) fail(QT_TRANSLATE_NOOP("Jcamp", "A difference needs a preceding absolute intensity."));
                previousDelta = t.value;
                previous += t.value;
            } else previous = t.value;
            previousKind = t.kind;
            append(first + double(out.points.size()) * dx, previous * yf);
        }
        checkpoint = previousKind == Kind::Difference && std::isfinite(previous);
    }
    if (checkpoint) fail(QT_TRANSLATE_NOOP("Jcamp", "Difference-coded data is missing its final checkpoint."));
    if (out.points.empty() || (n && out.points.size() != n)) fail(QT_TRANSLATE_NOOP("Jcamp", "The spectrum sample count does not match its data."));
    double offset = 0;
    if (h.contains(".SHIFTREFERENCE")) {
        const auto ref = h[".SHIFTREFERENCE"].split(',');
        if (ref.size() != 4) fail(QT_TRANSLATE_NOOP("Jcamp", "Invalid chemical shift reference."));
        const size_t point = count(ref[2]);  // JCAMP 5.01 numbers the first data point as 1
        if (point > out.points.size()) fail(QT_TRANSLATE_NOOP("Jcamp", "Chemical shift reference is outside the spectrum."));
        offset = number(ref[3]) - out.points[point - 1].x() / (hz ? frequency : 1);
    } else if (h.contains("$OFFSET")) {
        double high = out.points.front().x();
        if (h.contains("FIRSTX") || table.key == "DATATABLE") high = std::max(high, first);
        if (h.contains("LASTX") || table.key == "DATATABLE") high = std::max(high, last);
        for (const auto& p : out.points) high = std::max(high, p.x());
        offset = number(h["$OFFSET"]) - high / (hz ? frequency : 1);
    }
    else if (hz && h.contains("$REFERENCEPOINT")) offset = -number(h["$REFERENCEPOINT"]) / frequency;
    bool intensity = false;
    for (auto& p : out.points) {
        p.setX(p.x() / (hz ? frequency : 1) + offset);
        if (!std::isfinite(p.x()) || std::abs(p.x()) > 1e6) fail(QT_TRANSLATE_NOOP("Jcamp", "Chemical shifts are outside the supported ppm range."));
        intensity |= std::isfinite(p.y());
    }
    if (!intensity) fail(QT_TRANSLATE_NOOP("Jcamp", "The spectrum has no numeric intensities."));
    return out;
}
}  // namespace

std::vector<JcampSpectrum> readJcamp(const QByteArray& data, QString* error) {
    if (error) error->clear();
    try {
        if (data.size() > 200'000'000) fail(QT_TRANSLATE_NOOP("Jcamp", "The spectrum file is too large."));
        std::vector<Block> stack;
        std::vector<JcampSpectrum> out;
        size_t samples = 0;
        for (QString line : QString::fromUtf8(data).remove(QChar(0xfeff)).split('\n')) {
            line = line.section("$$", 0, 0).trimmed();
            if (line.isEmpty()) continue;
            if (line.startsWith("##")) {
                const qsizetype equal = line.indexOf('=');
                if (equal < 0) fail(QT_TRANSLATE_NOOP("Jcamp", "Malformed spectrum header."));
                const QString key = compact(line.mid(2, equal - 2)), value = line.mid(equal + 1).trimmed();
                if (key == "TITLE") {
                    if (!stack.empty() && compact(stack.back().headers.value("DATATYPE")) != "LINK") fail(QT_TRANSLATE_NOOP("Jcamp", "A spectrum block is missing END."));
                    Block b;
                    if (!stack.empty()) b.headers = stack.back().headers;
                    b.headers["TITLE"] = value;
                    b.lastKey = "TITLE";
                    stack.push_back(std::move(b));
                } else if (key == "END") {
                    if (stack.empty()) fail(QT_TRANSLATE_NOOP("Jcamp", "Spectrum END has no matching TITLE."));
                    const Block& b = stack.back();
                    const QString type = compact(b.headers.value("DATATYPE"));
                    if (type == "NMRSPECTRUM" || type == "NMRPEAKTABLE") {
                        if (number(b.headers.value("NUMDIM", "1")) != 1) fail(QT_TRANSLATE_NOOP("Jcamp", "Only one-dimensional NMR spectra can be overlaid."));
                        if (compact(b.headers.value("VARTYPE")).split(',').count("INDEPENDENT") > 1)
                            fail(QT_TRANSLATE_NOOP("Jcamp", "Only one-dimensional NMR spectra can be overlaid."));
                        bool real = false;
                        for (const Table& t : b.tables) {
                            if (compact(t.form).contains("(I..I)")) continue;
                            if (real) fail(QT_TRANSLATE_NOOP("Jcamp", "A block contains more than one real spectrum."));
                            JcampSpectrum parsed = spectrum(b, t);
                            if (parsed.points.size() > maxSamples - samples) fail(QT_TRANSLATE_NOOP("Jcamp", "The spectrum contains too many samples."));
                            samples += parsed.points.size();
                            out.push_back(std::move(parsed));
                            real = true;
                        }
                        if (b.tables.empty()) fail(QT_TRANSLATE_NOOP("Jcamp", "The NMR block has no real data table."));
                    }
                    stack.pop_back();
                } else {
                    if (stack.empty()) fail(QT_TRANSLATE_NOOP("Jcamp", "Spectrum headers need a TITLE block."));
                    Block& b = stack.back();
                    b.headers[key] = value;
                    b.lastKey = key;
                    b.table = -1;
                    if (key == "XYDATA" || key == "XYPOINTS" || key == "PEAKTABLE" || key == "DATATABLE") {
                        b.tables.push_back({key, value, {}});
                        b.table = int(b.tables.size()) - 1;
                    }
                }
            } else if (!stack.empty()) {
                Block& b = stack.back();
                if (b.table >= 0) b.tables[b.table].lines << line;
                else if (!b.lastKey.isEmpty()) b.headers[b.lastKey] += " " + line;
            } else fail(QT_TRANSLATE_NOOP("Jcamp", "Unexpected text outside a spectrum block."));
        }
        if (!stack.empty()) fail(QT_TRANSLATE_NOOP("Jcamp", "The spectrum file is incomplete (missing END)."));
        if (out.empty()) fail(QT_TRANSLATE_NOOP("Jcamp", "No processed 1H or 13C NMR spectrum was found. FIDs must be processed before import."));
        return out;
    } catch (const QString& message) {
        if (error) *error = message;
        return {};
    } catch (const std::bad_alloc&) {
        if (error) *error = QCoreApplication::translate("Jcamp", "Not enough memory to import the spectrum.");
        return {};
    }
}
