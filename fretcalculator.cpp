
#include "fretcalculator.h"

FretCalculator::FretCalculator()
{
    // reset calculation process flags to false
    for (bool &flag : calcProcess)
    {
        flag = false;
    }

    // reset data
    for (int i = 0; i < 3; ++ i)
    {
        valFore[i] = 0;
        valBack[i] = 0;
        exposureTime[i] = 0;
    }

    // reset ratios
    for (int i = 0; i < 7; ++ i)
    {
        ratio[i] = 0;
    }
}

bool FretCalculator::checkReady()
{
    bool isReady = true;
    for (bool flag : calcProcess) {
        isReady = flag && isReady;
    }
    return isReady;
}

void FretCalculator::resetCalculator()
{
    // reset calculation process flags to false
    for (bool &flag : calcProcess) {
        flag = false;
    }

    // reset data
    for (int i = 0; i < 3; ++ i) {
        valFore[i] = 0;
        valBack[i] = 0;
        exposureTime[i] = 0;
    }
}
void FretCalculator::resetData()
{
    // reset background data and foreground data
    for (bool &flag : calcProcess) {
        flag = false;
    }
    // reset data
    for (int i = 0; i < 3; ++ i) {
        valFore[i] = 0;
        valBack[i] = 0;
    }
    // reset results
    for (double &t : result) {
        t = 0.0;
    }
}

void FretCalculator::correctData()
{
    if (!calcProcess[LOAD_DATA]) {
        calcProcess[CORRECT_DATA] = false;
        return;
    }
    for (int i = 0; i < 3; ++ i) {
        valCorr[i] = valFore[i] - valBack[i];
        if (valCorr[i] <= 0) {
            qDebug() << "[Correct Data]:\tBackground out of Foreground";
            calcProcess[CORRECT_DATA] = false;
            return;
        }
        valCorr[i] /= exposureTime[i];
    }
    calcProcess[CORRECT_DATA] = true;
}

void FretCalculator::calcEFret()
{
    valFc = valCorr[CHANNEL_NAME_DA]
            - ratio[RATIO_NAME_A] * (valCorr[CHANNEL_NAME_AA] - ratio[RATIO_NAME_C] * valCorr[CHANNEL_NAME_DD])
            - ratio[RATIO_NAME_D] * (valCorr[CHANNEL_NAME_DD] - ratio[RATIO_NAME_B] * valCorr[CHANNEL_NAME_AA]);

    if (valFc <= 0) {
        qDebug() << "[Calc E-FRET]:\tNegative Nonsense Fc";
        calcProcess[CALC_DATA] = false;
        return;
    }

    // ED
    result[Ed] = valFc / (valFc + ratio[RATIO_NAME_G] * valCorr[CHANNEL_NAME_DD]);
    // Rc
    result[Rad] = ratio[RATIO_NAME_K] * valCorr[CHANNEL_NAME_AA] / (valFc / ratio[RATIO_NAME_G] + valCorr[CHANNEL_NAME_DD]);

    calcProcess[CALC_DATA] = true;
}

void FretCalculator::calc3CubeFret()
{
    valFc = valCorr[CHANNEL_NAME_DA]
            - ratio[RATIO_NAME_A] * (valCorr[CHANNEL_NAME_AA] - ratio[RATIO_NAME_C] * valCorr[CHANNEL_NAME_DD])
            - ratio[RATIO_NAME_D] * (valCorr[CHANNEL_NAME_DD] - ratio[RATIO_NAME_B] * valCorr[CHANNEL_NAME_AA]);

    if (valFc <= 0) {
        qDebug() << "[Calc 3 Cube FRET]:\tNegative Nonsense Fc";
        calcProcess[CALC_DATA] = false;
        return;
    }

    // EA
    result[Ea] = valFc * ratio[RATIO_NAME_Y] / (valCorr[CHANNEL_NAME_AA] * ratio[RATIO_NAME_A]);
    // 1/Rc
    result[Rda] = (valFc / ratio[RATIO_NAME_G] + valCorr[CHANNEL_NAME_DD]) / (ratio[RATIO_NAME_K] * valCorr[CHANNEL_NAME_AA]);

    calcProcess[CALC_DATA] = true;
}
void FretCalculator::calc2Hybrid()
{
    //计算估计浓度
    double Ma = ratio[RATIO_NAME_G] * ratio[RATIO_NAME_Y] / ratio[RATIO_NAME_D];
    result[Dest] = ratio[RATIO_NAME_D] * valCorr[CHANNEL_NAME_DD] / (1 - result[Ed]);
    result[Aest] = ratio[RATIO_NAME_A] * (valCorr[CHANNEL_NAME_AA] - ratio[RATIO_NAME_C] * valCorr[CHANNEL_NAME_DD]) / Ma;
}

void FretCalculator::calcData()
{
    if (!calcProcess[CORRECT_DATA]) {
        calcProcess[CALC_DATA] = false;
        return;
    }
    calcProcess[CALC_DATA] = true;
    calcEFret();
    calc3CubeFret();
    calc2Hybrid();
}

bool FretCalculator::isParamSet()
{
    return checkParamSet();
}
bool FretCalculator::isDataLoaded()
{
    return checkDataLoaded();
}
bool FretCalculator::isDataCorrected()
{
    return checkDataCorrected();
}
bool FretCalculator::isDataCalculated()
{
    return calcProcess[CALC_DATA];
}

void FretCalculator::setRatio(RatioName ratioName, double value)
{
    if (ratioName == RATIO_NAME_SIZE) {
        calcProcess[SET_PARAM] = checkParamSet();
        return;
    }
    ratio[ratioName] = value;
    calcProcess[SET_PARAM] = checkParamSet();
}
void FretCalculator::setRatios(double a, double b, double c, double d, double g, double k, double y)
{
    setRatio(RATIO_NAME_A, a);
    setRatio(RATIO_NAME_B, b);
    setRatio(RATIO_NAME_C, c);
    setRatio(RATIO_NAME_D, d);
    setRatio(RATIO_NAME_G, g);
    setRatio(RATIO_NAME_K, k);
    setRatio(RATIO_NAME_Y, y);
}

void FretCalculator::setExposureTime(ChannelName channelName, double value)
{
    if (channelName == CHANNEL_NAME_SIZE) {
        calcProcess[SET_PARAM] = checkParamSet();
        return;
    }
    exposureTime[channelName] = value;
    calcProcess[SET_PARAM] = checkParamSet();
}
void FretCalculator::setExposureTimes(double expsTimeAA, double expsTimeDA, double expsTimeDD)
{
    setExposureTime(CHANNEL_NAME_AA, expsTimeAA);
    setExposureTime(CHANNEL_NAME_DA, expsTimeDA);
    setExposureTime(CHANNEL_NAME_DD, expsTimeDD);
}

void FretCalculator::loadBackGray(double I_AA, double I_DA, double I_DD)
{
    valBack[CHANNEL_NAME_AA] = I_AA;
    valBack[CHANNEL_NAME_DA] = I_DA;
    valBack[CHANNEL_NAME_DD] = I_DD;
    calcProcess[LOAD_DATA] = checkDataLoaded();
}
void FretCalculator::loadForeGray(double I_AA, double I_DA, double I_DD)
{
    valFore[CHANNEL_NAME_AA] = I_AA;
    valFore[CHANNEL_NAME_DA] = I_DA;
    valFore[CHANNEL_NAME_DD] = I_DD;
    calcProcess[LOAD_DATA] = checkDataLoaded();
}
void FretCalculator::loadCorrectedGray(double I_AA, double I_DA, double I_DD)
{
    valCorr[CHANNEL_NAME_AA] = I_AA;
    valCorr[CHANNEL_NAME_DA] = I_DA;
    valCorr[CHANNEL_NAME_DD] = I_DD;
    calcProcess[CORRECT_DATA] = checkDataCorrected();
}

bool FretCalculator::checkDataLoaded()
{
    if (valFore[CHANNEL_NAME_AA] <= 0 || valFore[CHANNEL_NAME_DD] <= 0 || valFore[CHANNEL_NAME_DA] <= 0) {
        return false;
    }
    if (valBack[CHANNEL_NAME_AA] <= 0 || valBack[CHANNEL_NAME_DD] <= 0 || valBack[CHANNEL_NAME_DA] <= 0) {
        return false;
    }
    return true;
}

bool FretCalculator::checkDataCorrected()
{
    if (valCorr[CHANNEL_NAME_AA] <= 0 || valCorr[CHANNEL_NAME_DD] <= 0 || valCorr[CHANNEL_NAME_DA] <= 0) {
        return false;
    }
    return true;
}

bool FretCalculator::checkParamSet()
{
    if (ratio[RATIO_NAME_A] <= 0 || ratio[RATIO_NAME_B] <= 0 || ratio[RATIO_NAME_C] <= 0 || ratio[RATIO_NAME_D] <= 0 || ratio[RATIO_NAME_G] <= 0 || ratio[RATIO_NAME_K] <= 0 || ratio[RATIO_NAME_Y]) {
        return false;
    }
    if (exposureTime[CHANNEL_NAME_AA] <= 0 || exposureTime[CHANNEL_NAME_DA] <= 0 || exposureTime[CHANNEL_NAME_DD] <= 0) {
        return false;
    }
    return true;
}

bool FretCalculator::isSuccess()
{
    return calcProcess[CALC_DATA];
}

void FretCalculator::showParams()
{
    qDebug() << "[RATIO_NAME_A]:\t" << ratio[RATIO_NAME_A];
    qDebug() << "[RATIO_NAME_B]:\t" << ratio[RATIO_NAME_B];
    qDebug() << "[RATIO_NAME_C]:\t" << ratio[RATIO_NAME_C];
    qDebug() << "[RATIO_NAME_D]:\t" << ratio[RATIO_NAME_D];
    qDebug() << "[RATIO_NAME_G]:\t" << ratio[RATIO_NAME_G];
    qDebug() << "[RATIO_NAME_K]:\t" << ratio[RATIO_NAME_K];
    qDebug() << "[RATIO_NAME_Y]:\t" << ratio[RATIO_NAME_Y];
    qDebug() << "[ExpAA]:\t" << exposureTime[CHANNEL_NAME_AA];
    qDebug() << "[ExpDA]:\t" << exposureTime[CHANNEL_NAME_DA];
    qDebug() << "[ExpDD]:\t" << exposureTime[CHANNEL_NAME_DD];
}

double FretCalculator::getResult(CalcResult resultName)
{
    if (result[resultName] < 0) {
        return 0.0;
    }
    return result[resultName];
}

double FretCalculator::getFcValue()
{
    return valFc;
}
