
#ifndef ENUMERATOR_H
#define ENUMERATOR_H

#include <QString>

enum ThaResult
{
    KDEFF, ND_NA, EAMAX, EDMAX
};

// Channel Name
enum ChannelName {
    AA, DA, DD, NaC // NaC means not a channel name
};

enum RatioName {
    A, B, C, D, G, K, Y, NaR    // NaP means "Not a ratio"
};

enum CalcProcess {
    SET_PARAM, LOAD_DATA, CORRECT_DATA, CALC_DATA
};

enum CalcResult {
    Ed, Rad, Ea, Rda, Aest, Dest, Afree, Dfree
};
enum ChartName
{
    EdRad,
    EaRda,
    EdAfree,
    EaDfree,
    EdAfreeBin,
    EaDfreeBin
};

enum ShowType
{
    MERGED,
    AANORM,
    DANORM,
    DDNORM,
    RCPSEU,
    EDPSEU
};
enum CursorPosition
{
    CursorPositionUndefined,
    CursorPositionMiddle,
    CursorPositionTop,
    CursorPositionBottom,
    CursorPositionLeft,
    CursorPositionRight,
    CursorPositionTopLeft,
    CursorPositionTopRight,
    CursorPositionBottomLeft,
    CursorPositionBottomRight,
    CursorPositionOutside
};

enum DrawMode {
    CropMode, StampMode
};

class Enumerator
{
public:
    Enumerator();
    static ChannelName showTypeToChannelName(ShowType showType);
    static QString chartNameToQString(ChartName chartName);
};

#endif // ENUMERATOR_H
