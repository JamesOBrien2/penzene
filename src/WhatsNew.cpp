#include "WhatsNew.h"
#include "Online.h"
#include "Render.h"

#include <QApplication>
#include <QDialog>
#include <QFile>
#include <QGridLayout>
#include <QGuiApplication>
#include <QHBoxLayout>
#include <QLabel>
#include <QPainter>
#include <QPushButton>
#include <QScreen>
#include <QScrollArea>
#include <QSvgRenderer>
#include <QVBoxLayout>

namespace {

struct Tokens {
    QString page, card, border, text, secondary, accent, accentBg;
};

// The app's look: white, one teal accent; its dark twin when the app is dark.
Tokens tokens(const QWidget* w) {
    const Theme& t = theme(w->palette().color(QPalette::Window).lightness() < 128 ? "Dark" : "Light");
    const Chrome c = chrome(t);
    return {t.window.name(), t.surface.name(), c.border.name(), t.text.name(), c.secondary.name(), t.accent.name(),
            c.accentBg.name()};
}

// An SVG from the resources, drawn in `color` (the icons use currentColor).
QPixmap svg(const QString& path, QColor color, int size, qreal ratio) {
    QFile f(path);
    if (!f.open(QIODevice::ReadOnly)) return {};
    QByteArray data = f.readAll();
    data.replace("currentColor", color.name().toUtf8());
    QSvgRenderer renderer(data);
    QPixmap pm(QSize(size, size) * ratio);
    pm.setDevicePixelRatio(ratio);
    pm.fill(Qt::transparent);
    QPainter p(&pm);
    renderer.render(&p, QRectF(0, 0, size, size));
    return pm;
}

}  // namespace

void showWhatsNewDialog(QWidget* parent, const QString& changelog, const QString& version) {
    const auto notes = online::parseReleaseNotes(online::releaseNotes(changelog, version));
    QDialog dialog(parent);
    dialog.setObjectName("whatsNew");
    dialog.setWindowTitle(QObject::tr("What's New"));
    const Tokens t = tokens(&dialog);
    const qreal ratio = dialog.devicePixelRatioF();
    dialog.setStyleSheet(QString(R"(
        QDialog#whatsNew { background: %1; }
        QLabel { color: %4; background: transparent; }
        QLabel#title { font-size: 20px; font-weight: 600; }
        QLabel#highlightTitle { font-size: 14px; font-weight: 600; }
        QLabel#subtitle, QLabel#detail, QLabel#others, QLabel#manyMore { color: %5; }
        QLabel#section { color: %6; font-weight: 600; }
        QFrame#highlight { background: %2; border: 1px solid %3; border-radius: 12px; }
        QLabel#badge { background: %7; border-radius: 21px; }
        QScrollArea { background: transparent; border: none; }
        QScrollArea#moreScroll { border: 1px solid %3; border-radius: 10px; }
        QScrollArea > QWidget > QWidget { background: transparent; }
        QPushButton#continue { background: %6; color: %1; border: none; border-radius: 10px; padding: 7px 22px; font-weight: 600; }
        QPushButton#continue:hover { background: %4; }
    )").arg(t.page, t.card, t.border, t.text, t.secondary, t.accent, t.accentBg));

    auto* layout = new QVBoxLayout(&dialog);
    layout->setContentsMargins(24, 20, 24, 16);
    layout->setSpacing(10);

    auto* logo = new QLabel;
    logo->setPixmap(svg(":/logo.svg", QColor(t.accent), 48, ratio));
    logo->setAlignment(Qt::AlignCenter);
    auto* title = new QLabel(QObject::tr("What's new in Penzene %1").arg(version));
    title->setObjectName("title");
    title->setAlignment(Qt::AlignCenter);
    layout->addWidget(logo);
    layout->addWidget(title);

    if (!notes.highlights.empty()) {
        auto* grid = new QGridLayout;
        grid->setSpacing(10);
        const int n = int(notes.highlights.size());
        const int columns = n > 6 ? 3 : 2;
        for (int k = 0; k < n; ++k) {
            const auto& h = notes.highlights[k];
            auto* card = new QFrame;
            card->setObjectName("highlight");
            auto* row = new QHBoxLayout(card);
            row->setContentsMargins(10, 10, 10, 10);
            row->setSpacing(10);
            auto* badge = new QLabel;
            badge->setObjectName("badge");
            badge->setFixedSize(42, 42);
            badge->setAlignment(Qt::AlignCenter);
            QPixmap icon = svg(":/whatsnew/" + h.icon + ".svg", QColor(t.accent), 24, ratio);
            if (icon.isNull()) icon = svg(":/whatsnew/sparkles.svg", QColor(t.accent), 24, ratio);
            badge->setPixmap(icon);
            row->addWidget(badge, 0, Qt::AlignTop);
            auto* text = new QVBoxLayout;
            text->setSpacing(2);
            auto* name = new QLabel(h.title);
            name->setObjectName("highlightTitle");
            name->setTextFormat(Qt::PlainText);  // changelog text, never markup
            name->setWordWrap(true);
            auto* detail = new QLabel(h.detail);
            detail->setObjectName("detail");
            detail->setTextFormat(Qt::PlainText);
            detail->setWordWrap(true);
            text->addWidget(name);
            text->addWidget(detail);
            text->addStretch();
            row->addLayout(text, 1);
            // The last one out takes the whole last row.
            const bool last = k == n - 1 && n % columns;
            grid->addWidget(card, k / columns, k % columns, 1, last ? columns - k % columns : 1);
        }
        for (int c = 0; c < columns; ++c) grid->setColumnStretch(c, 1);
        layout->addLayout(grid);
    }

    if (!notes.others.isEmpty()) {
        auto* section = new QLabel(QObject::tr("More in this release"));
        section->setObjectName("section");
        layout->addWidget(section);

        auto* list = new QWidget;
        list->setAutoFillBackground(false);
        auto* listLayout = new QVBoxLayout(list);
        listLayout->setContentsMargins(10, 8, 10, 8);
        listLayout->setSpacing(4);
        for (const QString& line : notes.others) {
            auto* other = new QLabel("• " + line);
            other->setObjectName("others");
            other->setTextFormat(Qt::PlainText);
            other->setWordWrap(true);
            listLayout->addWidget(other);
        }
        auto* many = new QLabel(QObject::tr("…and many more small fixes and improvements."));
        many->setObjectName("manyMore");
        many->setTextFormat(Qt::PlainText);
        many->setWordWrap(true);
        listLayout->addSpacing(4);
        listLayout->addWidget(many);
        listLayout->addStretch();

        auto* scroll = new QScrollArea;
        scroll->setObjectName("moreScroll");
        scroll->setFrameShape(QFrame::NoFrame);
        scroll->setWidgetResizable(true);
        scroll->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        scroll->setMinimumHeight(0);
        scroll->setMaximumHeight(220);
        scroll->viewport()->setAutoFillBackground(false);
        scroll->setWidget(list);
        layout->addWidget(scroll, 1);
    }
    if (notes.highlights.empty() && notes.others.isEmpty()) layout->addWidget(new QLabel(QObject::tr("No notes for this version.")));

    auto* footer = new QHBoxLayout;
    auto* all = new QLabel(QString("<a style='color:%1' href='https://github.com/JamesOBrien2/penzene/blob/main/CHANGELOG.md'>%2</a>")
                               .arg(t.accent, QObject::tr("All releases")));
    all->setOpenExternalLinks(true);
    auto* ok = new QPushButton(QObject::tr("Continue"));
    ok->setObjectName("continue");
    ok->setDefault(true);
    QObject::connect(ok, &QPushButton::clicked, &dialog, &QDialog::accept);
    footer->addWidget(all);
    footer->addStretch();
    footer->addWidget(ok);
    layout->addSpacing(4);
    layout->addLayout(footer);

    // A tall changelog scrolls; the window itself stays within the screen so the title bar's
    // close button and Continue are always reachable.
    const QScreen* screen = parent ? parent->screen() : QGuiApplication::primaryScreen();
    const QRect available = screen ? screen->availableGeometry() : QRect(0, 0, 900, 700);
    const int margin = 40;
    const int width = qMin(660, available.width() - 2 * margin);
    dialog.setFixedWidth(width);
    layout->activate();
    int height = dialog.heightForWidth(width);
    if (height <= 0) height = layout->sizeHint().height();
    height = qMin(height, available.height() - 2 * margin);
    dialog.setFixedSize(width, height);
    if (screen) dialog.move(available.center() - QPoint(width / 2, height / 2));
    dialog.exec();
}
