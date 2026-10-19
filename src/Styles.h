#pragma once
#include "Render.h"
#include <QJsonObject>

// Saved definitions use real points; the renderer uses model units.
StyleSettings styleSettings(const DrawingStyle& style);
DrawingStyle modelStyle(const StyleSettings& settings);
QJsonObject styleToJson(const StyleSettings& settings);
std::optional<StyleSettings> styleFromJson(const QJsonObject& json);
QString stylesDirectory();
std::vector<StyleSettings> userDrawingStyles();
bool saveDrawingStyle(const StyleSettings& settings, const QString& replacing = {});
bool deleteDrawingStyle(const QString& name);
bool applyDrawingStyle(Document& doc, const QString& name);
const std::vector<DrawingStyle>& builtInDrawingStyles();
