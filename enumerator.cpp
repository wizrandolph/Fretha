
#include "enumerator.h"

Enumerator::Enumerator()
{

}

ChannelName Enumerator::showTypeToChannelName(ViewType type)
{
    switch (type) {
    case VIEW_TYPE_AANORM:
        return CHANNEL_NAME_AA;
    case VIEW_TYPE_DANORM:
        return CHANNEL_NAME_DA;
    case VIEW_TYPE_DDNORM:
        return CHANNEL_NAME_DD;
    default:
        return CHANNEL_NAME_SIZE;
    }
}

QString Enumerator::chartNameToQString(ChartName chartName)
{
    switch (chartName) {
    case CHART_NAME_ED_AFREE_OURS:
        return "Ed-Afree Chart";
    case CHART_NAME_EA_DFREE_OURS:
        return "Ea-Dfree Chart";
    case CHART_NAME_ED_AFREE_BEN:
        return "Ed-Afree Chart";
    case CHART_NAME_EA_DFREE_BEN:
        return "Ea-Dfree Chart";
    case CHART_NAME_ED_RAD:
        return "Ed-Rc Chart";
    case CHART_NAME_EA_RDA:
        return "Ea-1/Rc Chart";
    default:
        return "Unknown Chart";
    }
}
