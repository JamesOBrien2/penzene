#pragma once
#include "Chem.h"
#include "Edit.h"
#include "Render.h"

#include <QFont>
#include <QGraphicsView>
#include <QPainterPath>
#include <QPicture>
#include <QPolygonF>
#include <QSet>
#include <optional>
#include <vector>

class QPainter;
class QMenu;
class QUndoStack;

class Canvas : public QGraphicsView {
    Q_OBJECT
public:
    enum class Tool { Select, Rotate3D, Atom, Bond, Wedge, Hash, Chain, Ring, ChargePlus, ChargeMinus, Erase, Arrow, Text, Fill, Colour };

    explicit Canvas(QUndoStack* undo, QWidget* parent = nullptr);

    const Document& document() const { return doc_; }
    // Every edit goes through here so it can be undone.
    void commit(const Document& next, const QString& text);
    void setDocumentSilently(const Document& doc);  // e.g. undo/redo, file open

    const QSet<int>& selection() const { return selectedAtoms_; }
    const QSet<int>& selectedArrows() const { return selectedArrows_; }
    const QSet<int>& selectedTexts() const { return selectedTexts_; }
    void setSelection(QSet<int> atoms, QSet<int> arrows = {}, QSet<int> texts = {});
    Document selectedSubset() const;  // selection (or everything) as a standalone doc
    // Lit without selecting them, e.g. the atoms behind an NMR stick; whoever sets it clears it on documentChanged.
    void setHighlight(QSet<int> atoms);
    const QSet<int>& highlight() const { return highlight_; }
    void deleteSelection();
    void setUndoStack(QUndoStack* undo) { undo_ = undo; }  // each page has its own history
    void setCompoundStart(edit::CompoundCount from) { compoundStart_ = std::move(from); }  // the pages before's numbers (#569)
    void insert(Document fragment, const QString& text);  // centred in view, selected
    void selectAll();
    void rotateSelection(double degrees);
    // The round handle above the selection box that turns it (Shift: 15° steps; Ctrl: square to
    // the page or 45° from it).
    std::optional<QPointF> rotateHandle() const;
    // Out of the page, about the page's x then y axis; stereo is kept (Shift+Alt+drag or arrows).
    void rotate3D(double aboutX, double aboutY);
    // Brackets around the selected atoms (square or round, with a subscript such as "n"); none removes theirs.
    void bracketSelection(bool square, const QString& label);
    void removeBrackets();
    void variableAttachment();  // a bond from the selected atoms to a new carbon, on any one of them
    void arrangeScheme();  // the selection (or everything) as a tidy reaction scheme
    // `t` (rotate, scale, stretch) about the selection's centre; everything if nothing is selected.
    void transformSelection(const QTransform& t, const QString& what);
    void duplicateSelection(QPointF dir);
    void flipSelection(bool horizontal);  // mirror image, as ChemDraw's flip
    void invertStereo();  // wedges become hashes and hashes wedges: the enantiomer, not redrawn
    enum class Align { Left, HCentre, Right, Top, VCentre, Bottom };
    void alignSelection(Align edge);
    void distributeSelection(bool horizontal);
    void centerOnPage();  // the selection (or everything) moved to the middle of the page, spacing kept
    // The right-click menu for whatever is at `scenePos` (public so tests can inspect it).
    QMenu* contextMenuAt(QPointF scenePos);
    void moveHotspot(QPointF dir, bool jump);
    void editLabel(int atom);
    void editAtomProperties(int atom);  // label, charge, map number, lone pairs, radicals, δ
    void expandAbbreviations();  // selection, else hotspot atom, else everything
    void editText(int text, QPointF pos = {});  // text < 0: new text at pos
    int hotspotAtom() const { return hoverAtom_; }
    int hotspotBond() const { return hoverBond_; }
    void setHotspot(int atom, int bond = -1);
    quint64 revision() const { return revision_; }  // counts document changes, for accessibility
    QPointF viewCenter() const;
    void zoomBy(double factor);
    void fitToDocument();
    void fitToSelection();  // everything when nothing is selected

    void setTool(Tool t) { tool_ = t; }
    void setElement(int z) { element_ = z; }
    // The Bond tool's order and style (None, Interaction or Partial).
    void setBondOrder(int order, BondStereo style = BondStereo::None) { bondOrder_ = order, bondStyle_ = style; }
    void setRing(int size, bool aromatic) { ringSize_ = size, ringAromatic_ = aromatic; }
    void setArrow(ArrowKind kind, bool curved, bool dashed = false, OrbitalLook look = OrbitalLook::Outline,
                  bool crossed = false) {
        arrowKind_ = kind, arrowCurved_ = curved, arrowDashed_ = dashed, arrowLook_ = look, arrowCrossed_ = crossed;
    }
    void setTheme(const Theme& t) { theme_ = t, refresh(); }
    void setGuides(bool grid, bool rulers) { grid_ = grid, rulers_ = rulers, viewport()->update(); }  // View menu
    void setFillColor(QColor c) { fillColor_ = c; }
    QColor fillColor() const { return fillColor_; }
    void setColour(QColor c) { colour_ = c; }
    QColor colour() const { return colour_; }
    void colourSelection();  // the current colour on the selected atoms, bonds, arrows and text
    void setArrowHead(double size);  // on the selected arrows, relative to the usual
    // The lone selected curved arrow bowed more (factor > 1) or less; a negative factor flips its side.
    void bendArrow(double factor);
    // A straight arrow, or new text, just right of the selection, or at the hotspot; nothing without either.
    std::optional<QPointF> nextPlace() const;
    void addArrowAfter();
    void addTextAfter();
    void numberCompounds();
    void groupSelection();    // the selected objects, whole molecules, act as one from now on (#410)
    void ungroupSelection();  // the groups the selection touches split up  // a bold number under each selected molecule (or every one) that has none (#504)

signals:
    void documentChanged();
    void selectionChanged();
    void hotspotAtomChanged(int atom);  // under the pointer or the keys; -1 for none
    void toolKey(const QString& key);  // x bond, X chain, j benzene, t text, e arrow, space select

protected:
    void drawBackground(QPainter* p, const QRectF& rect) override;
    void drawForeground(QPainter* p, const QRectF& rect) override;
    void mousePressEvent(QMouseEvent* e) override;
    void mouseMoveEvent(QMouseEvent* e) override;
    void mouseReleaseEvent(QMouseEvent* e) override;
    void mouseDoubleClickEvent(QMouseEvent* e) override;
    void keyPressEvent(QKeyEvent* e) override;
    void wheelEvent(QWheelEvent* e) override;
    void contextMenuEvent(QContextMenuEvent* e) override;
    bool event(QEvent* e) override;          // keys: then announce the hotspot
    bool viewportEvent(QEvent* e) override;  // the mouse: likewise
    void showEvent(QShowEvent* e) override;

private:
    int atomAt(QPointF p) const;
    int bondAt(QPointF p) const;
    int arrowAt(QPointF p) const;
    int textAt(QPointF p) const;
    Arrow draggedArrow() const;
    Arrow curvedArrow(QPointF from, QPointF to, ArrowKind kind) const;
    int draggedRingSize() const;
    void addDraggedRing(Document& doc) const;
    void refresh();
    void drawRulers(QPainter* p);
    void fit(const Document& part);
    std::vector<QPointF> dragPath() const;

    Document doc_;
    bool fitOnShow_ = false;
    QPicture picture_;
    std::vector<QRectF> labels_;  // each atom's drawn label box, as painted into picture_
    Theme theme_;
    bool grid_ = false, rulers_ = false;
    std::vector<QPointF> preview_;
    QUndoStack* undo_;
    edit::CompoundCount compoundStart_;
    void selectGroups();  // a click or marquee on one member takes its whole group
    QPainterPath groupShape(int group) const;  // one smooth outline around a group's members
    Tool tool_ = Tool::Bond;
    int element_ = 6, bondOrder_ = 1, ringSize_ = 6;
    BondStereo bondStyle_ = BondStereo::None;
    bool ringAromatic_ = true;
    ArrowKind arrowKind_ = ArrowKind::Reaction;
    bool arrowCurved_ = false, arrowDashed_ = false, arrowCrossed_ = false;
    OrbitalLook arrowLook_ = OrbitalLook::Outline;
    QColor fillColor_ = QColor(207, 227, 255);
    QColor colour_ = QColor(0xFF, 0x0D, 0x0D);  // CPK oxygen

    QSet<int> selectedAtoms_, selectedArrows_, selectedTexts_;
    int hoverAtom_ = -1, hoverBond_ = -1;
    QSet<int> highlight_;
    quint64 revision_ = 0;
    int announcedAtom_ = -1;
    bool keyHotspot_ = false;  // G or > is picking: the arrow keys move the hotspot, not the selection
    edit::Hotspot arrowMark_;      // where > started a curved arrow
    Document shown_;  // as last drawn: the revision counts real changes
    void announceHotspot();

    // Drag state
    enum class Drag { None, Bond, Chain, Arrow, Ring, Move, Rotate, Rubber, Pan, Scale, Rotate3D, Reshape } drag_ = Drag::None;
    std::vector<int> moleculesOfSelection() const;  // whole molecules; all atoms if none selected
    std::optional<chem::Pose3D> pose_;  // during a 3D rotation drag
    // Scale handles around the selection: corners scale, edges stretch along one axis.
    QRectF selectionBox() const;  // empty unless something with extent is selected
    bool inSelectedBox(QPointF p) const;  // inside the box of a selected molecule, group, arrow or text
    int handleAt(QPointF p) const;  // 0..7 clockwise from the top-left corner, or -1
    int scaleHandle_ = -1;
    // A lone selected arrow is reshaped instead: handles on its ends and at the top of its curve.
    int reshapedArrow() const;  // its index, or -1
    int reshapeHandleAt(QPointF p) const;  // 0 from, 1 to, 2 the curve, or -1
    int reshaping_ = -1, reshape_ = -1;  // during the drag: the arrow, and its handle
    QRectF scaleBox_;
    QPointF pressPos_, curPos_;
    QPointF pressRaw_;  // where the mouse went down, before an orbital snaps pressPos_ to its atom
    QPolygonF lasso_;   // an Alt-drag selection's loop; empty for the rectangle
    bool shift_ = false;  // held during the drag: free bond angle, or move along one axis
    int pressAtom_ = -1;
    Document beforeDrag_;
    QPoint panLast_;
};
