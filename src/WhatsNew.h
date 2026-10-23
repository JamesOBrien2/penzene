#pragma once
class QString;
class QWidget;

// What's New in Penzene <version>: the release's highlights as large icon tiles, with the
// smaller changes in a scroll below (from that version's section of CHANGELOG.md). The window
// stays within the screen so its close button and Continue stay reachable.
void showWhatsNewDialog(QWidget* parent, const QString& changelog, const QString& version);
