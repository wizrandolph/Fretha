
#ifndef ENUMERATOR_H
#define ENUMERATOR_H

#include <QString>

enum RunMode {
    RUNMODE_DEBUG = -1,
    RUNMODE_RELEASE = 0,
    RUNMODE_EXPERIMENT = 1
};

#define RUNMODE 0

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

enum TableHeader {  // 14个
    TABLE_HEADER_IAA,
    TABLE_HEADER_IDA,
    TABLE_HEADER_IDD,
    TABLE_HEADER_ED,
    TABLE_HEADER_RAD,
    TABLE_HEADER_EA,
    TABLE_HEADER_RDA,
    TABLE_HEADER_AEST,
    TABLE_HEADER_DEST,
    TABLE_HEADER_RECTX,
    TABLE_HEADER_RECTY,
    TABLE_HEADER_RECTW,
    TABLE_HEADER_RECTH,
    TABLE_HEADER_VIEW
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
