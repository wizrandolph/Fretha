
#include "enumerator.h"

Enumerator::Enumerator()
{

}

ChannelName Enumerator::showTypeToChannelName(ShowType type)
{
    switch (type) {
    case AANORM:
        return AA;
    case DANORM:
        return DA;
    case DDNORM:
        return DD;
    default:
        return NaC;
    }
}

QString Enumerator::chartNameToQString(ChartName chartName)
{
    switch (chartName) {
    case EdAfree:
        return "Ed-Afree Chart";
    case EaDfree:
        return "Ea-Dfree Chart";
    case EdAfreeBin:
        return "Binned Ed-Afree Chart";
    case EaDfreeBin:
        return "Binned Ea-Dfree Chart";
    case EdRad:
        return "Ed-Rc Chart";
    case EaRda:
        return "Ea-1/Rc Chart";
    default:
        return "Error";
    }
}
