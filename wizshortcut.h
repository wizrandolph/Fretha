
#ifndef WIZSHORTCUT_H
#define WIZSHORTCUT_H




#include "qshortcut.h"

class WizShortcut : public QShortcut
{
public:
    WizShortcut(const QKeySequence &key, int pageIndex, QWidget *parent)
        : QShortcut(key, parent), m_pageIndex(pageIndex)
    {
    }

    int getPageIndex() const
    {
        return m_pageIndex;
    }

private:
    int m_pageIndex;
};

#endif // WIZSHORTCUT_H
