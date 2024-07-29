
#ifndef WIZGRAPHICSVIEW_H
#define WIZGRAPHICSVIEW_H

#include <QGraphicsView>
#include <QGraphicsScene>
#include <QGraphicsRectItem>
#include <QMouseEvent>
#include <QDebug>
#include <enumerator.h>

class WizGraphicsView : public QGraphicsView
{
    Q_OBJECT

public:
    WizGraphicsView(QWidget *parent = nullptr)
        : QGraphicsView(parent), m_rectItem(nullptr)
    {
        m_scene = new QGraphicsScene(this);
        setScene(m_scene);
        m_rect = QRectF(0, 0, 100, 100);
        setMouseTracking(true);
    }

    void setImage(const QImage &image);
    QRectF getRect();
    void setRect(QRectF);
    void setDrawMode(DrawMode);
    void addItem(QRectF);
    void clearRectItemList();

protected:
    void mousePressEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override;
    void mouseReleaseEvent(QMouseEvent *event) override;
    void wheelEvent(QWheelEvent *event) override;
private:
    QGraphicsScene *m_scene;
    QGraphicsRectItem *m_rectItem;  // 仅显示
    QList<QGraphicsRectItem*> m_rectItemsList;
    QPointF m_pressStartPosition;

    bool m_pressed = false; // 鼠标按下状态

    bool m_draggingScene = false;   // 右键拖拽
    bool m_resizingRect = false;    // 调整Rect大小
    bool m_movingRect = false;  // 移动 Rect
    bool m_creatingRect = false; // 创建 Rect
    bool m_movingStampRect = false;     // 移动 Stamp 模式下的 Rect

    QRectF m_rect = QRectF(0, 0, 50, 50);
    QRectF m_sceneRect;
    CursorPosition m_cursorPosition;

    DrawMode m_drawMode = DRAW_MODE_CROP;

    CursorPosition cursorPosition(const QRectF& cropRect, const QPointF& mousePosition);
    void updateCursorIcon(const QPointF& _mousePosition);
    void updateRect(QRectF rect);
    void rectToItem();
    bool isPointNearSide(const int _sideCoordinate, const int _pointCoordinate);
    bool isPointOutsideRect(const QRectF& _rect, const QPointF& _point);
    void correctQRectF(QRectF& rect);
    void correctRectInsideScene(QRectF& rect);
    void correctPointInsideScene(QPointF& point);
    void correctPointInteger(QPointF& point);
    void correctRectInteger(QRectF& rect);
    QRectF resizeRect(QPointF delta);

signals:
    void mouseReleased(QRectF);
};

#endif // WIZGRAPHICSVIEW_H
