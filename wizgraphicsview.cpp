
#include "wizgraphicsview.h"
#include "qapplication.h"
#include "qscrollbar.h"

static const int indentAbs = 10;

void WizGraphicsView::setImage(const QImage &image)
{
    clearRectItemList();
    qreal w = m_rect.width();
    qreal h = m_rect.height();
    qreal l = image.width() / 2 - w / 2;
    qreal t = image.height() / 2 - h / 2;
    m_scene->clear();
    m_scene->addPixmap(QPixmap::fromImage(image));
    m_sceneRect = image.rect();
    m_rectItem = m_scene->addRect(QRectF(l, t, w, h), QPen(Qt::red), QBrush(QColor(69, 170, 242, 0)));
    m_rectItem->setFlag(QGraphicsItem::ItemIsMovable);
    itemToRect();
    centerOn(m_rect.center());
}

QRectF WizGraphicsView::getRect()
{
    return m_rect;
}

void WizGraphicsView::setRect(QRectF rect)
{
    m_rect = rect;
    rectToItem();
    centerOn(m_rect.center());
}

void WizGraphicsView::setDrawMode(DrawMode drawMode)
{
    m_drawMode = drawMode;
    if (m_drawMode == StampMode) {
        setCursor(Qt::CrossCursor);
    }
    else {
        setCursor(Qt::ArrowCursor);
    }
}

void WizGraphicsView::addItem(QRectF rect)
{
    QGraphicsRectItem* rectItem = m_scene->addRect(QRectF(0, 0, rect.width(), rect.height()), QPen(Qt::NoPen), QBrush(QColor(69, 170, 242, 128)));
    rectItem->setPos(rect.x(), rect.y());
    m_rectItemsList.append(rectItem);
    qDebug() << "[AddItem]:\t" << m_rectItemsList.size();
}

void WizGraphicsView::mousePressEvent(QMouseEvent *event)
{
    m_pressed = true;

    // 鼠标右键
    if (event->button() == Qt::RightButton) {

        qDebug() << "[Mouse Action]:\t" << "RightButton Pressed";

        m_draggingScene = true;
        m_pressStartPosition = event->pos();
        QApplication::setOverrideCursor(Qt::ClosedHandCursor);
    }
    else if (event->button() == Qt::LeftButton) {

        qDebug() << "[Mouse Action]:\t" << "LeftButton Pressed";

        if (m_drawMode == CropMode) {

            // 记录鼠标的信息
            m_pressStartPosition = mapToScene(event->pos());

            // 检测创建item时的边界
            correctPointInsideScene(m_pressStartPosition);

            m_cursorPosition = cursorPosition(m_rect, m_pressStartPosition);    // 记录下按下鼠标时鼠标的位置

            // 如果在内部，鼠标为打开的手
            if (cursor().shape() == Qt::OpenHandCursor) {
                m_movingRect = true;
                QApplication::setOverrideCursor(Qt::ClosedHandCursor);
            }

            // 如果在外部，鼠标为十字光标
            if (cursor().shape() == Qt::ArrowCursor) {
                m_creatingRect = true;
                QApplication::setOverrideCursor(Qt::CrossCursor);
            }

            // 如果在矩形框边缘，准备调整大小
            if (cursor().shape() == Qt::SizeFDiagCursor || cursor().shape() == Qt::SizeBDiagCursor ||
                cursor().shape() == Qt::SizeVerCursor || cursor().shape() == Qt::SizeHorCursor) {
                m_resizingRect = true;
            }
        }
        else if (m_drawMode == StampMode) {

            // 将鼠标位置视作中心坐标，生成邮戳模式下的矩形
            QPointF centerPosition = mapToScene(event->pos());
            qreal left = centerPosition.x() - m_rect.width() / 2;
            qreal top = centerPosition.y() - m_rect.height() / 2;
            QRectF stampRect = QRectF(left, top, m_rect.width(), m_rect.height());

            // 边界检测和修正
            correctRectInsideScene(stampRect);

            // 应用矩形
            m_rect = stampRect;
            rectToItem();
        }
    }

    QGraphicsView::mousePressEvent(event);
}

void WizGraphicsView::mouseMoveEvent(QMouseEvent *event)
{
    QPointF mousePos = mapToScene(event->pos());
    correctPointInsideScene(mousePos);

    // 只有鼠标没有被按下的时候，更新鼠标光标
    if (!m_pressed && m_drawMode == CropMode) {
        QPointF mousePos = mapToScene(event->pos());
        updateCursorIcon(mousePos);
    }
    else {
        // 右键拖拽视图
        if (m_draggingScene) {
            // 计算拖拽视图的变化量
            QPointF delta = event->pos() - m_pressStartPosition;
            horizontalScrollBar()->setValue(horizontalScrollBar()->value() - delta.x());
            verticalScrollBar()->setValue(verticalScrollBar()->value() - delta.y());
            m_pressStartPosition = event->pos();
        }

        // 移动矩形区域
        if (m_movingRect) {
            // 使用rect推算是否碰撞到边界
            QPointF moveEndPosition = mapToScene(event->pos());
            QPointF delta = moveEndPosition - m_pressStartPosition;
            QPointF newTopLeft = m_rect.topLeft() + delta;
            QRectF newRect = QRectF(newTopLeft, m_rect.size());

            correctRectInsideScene(newRect);

            m_rect = newRect;
            rectToItem();

            m_pressStartPosition = moveEndPosition;
        }

        // 创建新矩形区域
        if (m_creatingRect) {
            QPointF moveEndPosition = mousePos;
            m_rect = QRectF(m_pressStartPosition, moveEndPosition);
            correctQRectF(m_rect);
            rectToItem();
        }

        // 调整矩形框大小
        if (m_resizingRect) {
            // 计算鼠标移动距离
            QPointF moveEndPosition = mapToScene(event->pos());
            QPointF delta = moveEndPosition - m_pressStartPosition;

            // 根据矩形框的边界位置和鼠标移动距离计算新的矩形框尺寸
            QRectF newRect = m_rect;
            if (m_cursorPosition == CursorPositionTopLeft) {
                newRect.setTopLeft(newRect.topLeft() + delta);
                // 限制最小
                if (newRect.height() <= 1) {
                    qreal bottom = newRect.bottom();
                    newRect.setTop(bottom - 1);
                }
                if (newRect.width() <= 1) {
                    qreal right = newRect.right();
                    newRect.setLeft(right - 1);
                }
            }
            else if (m_cursorPosition == CursorPositionTopRight) {
                newRect.setTopRight(newRect.topRight() + delta);
                // 限制最小
                if (newRect.height() <= 1) {
                    qreal bottom = newRect.bottom();
                    newRect.setTop(bottom - 1);
                }
                if (newRect.width() <= 1) {
                    newRect.setWidth(1);
                }
            }
            else if (m_cursorPosition == CursorPositionBottomLeft) {
                newRect.setBottomLeft(newRect.bottomLeft() + delta);
                // 限制最小
                if (newRect.height() <= 1) {
                    newRect.setHeight(1);
                }
                if (newRect.width() <= 1) {
                    qreal right = newRect.right();
                    newRect.setLeft(right - 1);
                }
            }
            else if (m_cursorPosition == CursorPositionBottomRight) {
                newRect.setBottomRight(newRect.bottomRight() + delta);
                // 限制最小
                if (newRect.height() <= 1) {
                    newRect.setHeight(1);
                }
                if (newRect.width() <= 1) {
                    newRect.setWidth(1);
                }
            }
            else if (m_cursorPosition == CursorPositionTop) {
                newRect.setTop(newRect.top() + delta.y());
                // 限制最小
                if (newRect.height() <= 1) {
                    qreal bottom = newRect.bottom();
                    newRect.setTop(bottom - 1);
                }
            }
            else if (m_cursorPosition == CursorPositionBottom) {
                newRect.setBottom(newRect.bottom() + delta.y());
                // 限制最小
                if (newRect.height() <= 1) {
                    newRect.setHeight(1);
                }
            }
            else if (m_cursorPosition == CursorPositionLeft) {
                newRect.setLeft(newRect.left() + delta.x());
                // 限制最小
                if (newRect.width() <= 1) {
                    qreal right = newRect.right();
                    newRect.setLeft(right - 1);
                }
            }
            else if (m_cursorPosition == CursorPositionRight) {
                newRect.setRight(newRect.right() + delta.x());
                // 限制最小
                if (newRect.width() <= 1) {
                    newRect.setWidth(1);
                }
            }

            // 更新矩形框和场景中的item
            m_rect = newRect;
            rectToItem();
            m_pressStartPosition = moveEndPosition;
        }
    }
}

void WizGraphicsView::mouseReleaseEvent(QMouseEvent *event)
{
    qDebug() << "[Mouse Action]:\t" << "Released";

    if (m_draggingScene) {
        m_draggingScene = false;
    }
    if (m_movingRect) {
        m_movingRect = false;
    }
    if (m_resizingRect) {
        m_resizingRect = false;
    }
    if (m_creatingRect) {
        m_creatingRect = false;
    }

    // 鼠标松开，更新状态
    m_pressed = false;
    QApplication::restoreOverrideCursor();

    // 发送信号
    emit mouseReleased(m_rect);

    QGraphicsView::mouseReleaseEvent(event);
}

void WizGraphicsView::wheelEvent(QWheelEvent *event)
{
    // 获取鼠标滚轮的滚动方向和滚动步数
    const QPoint angleDelta = event->angleDelta();
    const int numDegrees = angleDelta.y() / 8;
    const int numSteps = numDegrees / 15;

    // 获取当前视图的缩放因子
    qreal scaleFactor = transform().m11();

    // 计算视图中心点在场景坐标系中的位置
    QPointF viewCenter = mapToScene(viewport()->rect().center());

    // 根据鼠标位置和视图中心点的差异来调整缩放因子
    constexpr qreal scaleFactorStep = 1.2;
    if (numSteps > 0) {
        scaleFactor *= scaleFactorStep;
    }
    else if (numSteps < 0) {
        scaleFactor /= scaleFactorStep;
    }

    // 限制缩放因子的范围
    constexpr qreal minScaleFactor = 0.2;
    constexpr qreal maxScaleFactor = 10.0;
    scaleFactor = qMin(qMax(scaleFactor, minScaleFactor), maxScaleFactor);

    // 创建变换矩阵并应用缩放因子
    QTransform transform;
    transform.scale(scaleFactor, scaleFactor);
    setTransform(transform);

    // 将视图中心点移动到目标位置
    centerOn(viewCenter);

    // 触发视图重绘
    viewport()->update();
}

CursorPosition WizGraphicsView::cursorPosition(const QRectF& cropRect, const QPointF& mousePosition)
{
    CursorPosition cursorPosition = CursorPositionUndefined;
    if (1) {
        // 左上
        if (isPointNearSide(cropRect.top(), mousePosition.y()) &&
            isPointNearSide(cropRect.left(), mousePosition.x())) {
            cursorPosition = CursorPositionTopLeft;
        }
        // 左下
        else if (isPointNearSide(cropRect.bottom(), mousePosition.y()) &&
                 isPointNearSide(cropRect.left(), mousePosition.x())) {
            cursorPosition = CursorPositionBottomLeft;
        }
        // 右上
        else if (isPointNearSide(cropRect.top(), mousePosition.y()) &&
                 isPointNearSide(cropRect.right(), mousePosition.x())) {
            cursorPosition = CursorPositionTopRight;
        }
        // 右下
        else if (isPointNearSide(cropRect.bottom(), mousePosition.y()) &&
                 isPointNearSide(cropRect.right(), mousePosition.x())) {
            cursorPosition = CursorPositionBottomRight;
        }
        // 左
        else if (isPointNearSide(cropRect.left(), mousePosition.x()) &&
                 !isPointOutsideRect(cropRect, mousePosition)) {
            cursorPosition = CursorPositionLeft;
        }
        // 右
        else if (isPointNearSide(cropRect.right(), mousePosition.x()) &&
                 !isPointOutsideRect(cropRect, mousePosition)) {
            cursorPosition = CursorPositionRight;
        }
        // 上
        else if (isPointNearSide(cropRect.top(), mousePosition.y()) &&
                 !isPointOutsideRect(cropRect, mousePosition)) {
            cursorPosition = CursorPositionTop;
        }
        // 下
        else if (isPointNearSide(cropRect.bottom(), mousePosition.y()) &&
                 !isPointOutsideRect(cropRect, mousePosition)) {
            cursorPosition = CursorPositionBottom;
        }
        // 中间
        else if(cropRect.contains(mousePosition)) {
            cursorPosition = CursorPositionMiddle;
        }
        else {
            cursorPosition = CursorPositionOutside;
        }
    }
    return cursorPosition;
}

void WizGraphicsView::updateCursorIcon(const QPointF & _mousePosition)
{
    QCursor cursorIcon;
    //
    switch (cursorPosition(m_rect, _mousePosition))
    {
    case CursorPositionTopRight:
    case CursorPositionBottomLeft:
        cursorIcon = QCursor(Qt::SizeBDiagCursor);
        break;
    case CursorPositionTopLeft:
    case CursorPositionBottomRight:
        cursorIcon = QCursor(Qt::SizeFDiagCursor);
        break;
    case CursorPositionTop:
    case CursorPositionBottom:
        cursorIcon = QCursor(Qt::SizeVerCursor);
        break;
    case CursorPositionLeft:
    case CursorPositionRight:
        cursorIcon = QCursor(Qt::SizeHorCursor);
        break;
    case CursorPositionMiddle:
        cursorIcon = m_resizingRect ?
                         QCursor(Qt::ClosedHandCursor) :
                         QCursor(Qt::OpenHandCursor);
        break;
    case CursorPositionOutside:
        cursorIcon = QCursor(Qt::ArrowCursor);
        break;
    case CursorPositionUndefined:
    default:
        cursorIcon = QCursor(Qt::ArrowCursor);
        break;
    }
    this->setCursor(cursorIcon);
}

void WizGraphicsView::rectToItem()
{
    QRectF newRect = QRectF(0, 0, m_rect.width(), m_rect.height());
    m_rectItem->setRect(newRect);
    m_rectItem->setPos(m_rect.x(), m_rect.y());
}
void WizGraphicsView::itemToRect()
{
    m_rect = m_rectItem->sceneBoundingRect();
}
bool WizGraphicsView::isPointNearSide(const int sideCoordinate, const int pointCoordinate)
{
    int indent = indentAbs / transform().m11();
    return abs(sideCoordinate - pointCoordinate ) <indent  && abs(pointCoordinate - sideCoordinate) < indent;
}
bool WizGraphicsView::isPointOutsideRect(const QRectF& rect, const QPointF& point)
{
    int indent = indentAbs / transform().m11();
    return point.x() < rect.left() - indent || point.x() > rect.right() + indent || point.y() < rect.top() - indent || point.y() > rect.bottom() + indent;
}

void WizGraphicsView::correctQRectF(QRectF& rect)
{
    int x, y, w, h;

    x = rect.left();
    y = rect.top();
    w = rect.width();
    h = rect.height();

    if (w < 0.0) {
        x += w;
        w = -w;
    }

    if (h < 0.0) {
        y += h;
        h = -h;
    }

    rect = QRectF(x, y, w, h);
}

void WizGraphicsView::correctRectInsideScene(QRectF& rect)
{
    if (rect.left() < m_sceneRect.left()) {
        rect.moveLeft(m_sceneRect.left());
    }
    if (rect.right() > m_sceneRect.right()) {
        rect.moveRight(m_sceneRect.right());
    }
    if (rect.top() < m_sceneRect.top()) {
        rect.moveTop(m_sceneRect.top());
    }
    if (rect.bottom() > m_sceneRect.bottom()) {
        rect.moveBottom(m_sceneRect.bottom());
    }
}

void WizGraphicsView::correctPointInsideScene(QPointF& point)
{
    if (point.x() < m_sceneRect.left())
        point.setX(m_sceneRect.left());
    if (point.x() > m_sceneRect.right())
        point.setX(m_sceneRect.right());
    if (point.y() < m_sceneRect.top())
        point.setY(m_sceneRect.top());
    if (point.y() > m_sceneRect.bottom())
        point.setY(m_sceneRect.bottom());
}

void WizGraphicsView::clearRectItemList()
{
    for (auto it = m_rectItemsList.begin(); it != m_rectItemsList.end(); ++ it) {
        delete *it;
    }
    m_rectItemsList.clear();
}
