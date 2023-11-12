
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

    //    // reset ratios
    //    for (int i = 0; i < 7; ++ i)
    //    {
    //        ratio[i] = 0;
    //    }
}
void FretCalculator::resetData()
{
    // reset background data and foreground data
    for (bool &flag : calcProcess)
    {
        flag = false;
    }
    // reset data
    for (int i = 0; i < 3; ++ i)
    {
        valFore[i] = 0;
        valBack[i] = 0;
    }
    // reset results
    for (double &t : result)
    {
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
            qDebug() << "扣除背景后荧光强度非正，无法进行计算";
            calcProcess[CORRECT_DATA] = false;
            return;
        }
        valCorr[i] /= exposureTime[i];
    }
    calcProcess[CORRECT_DATA] = true;
}

void FretCalculator::calcEFret()
{
    valFc = valCorr[DA]
            - ratio[A] * (valCorr[AA] - ratio[C] * valCorr[DD])
            - ratio[D] * (valCorr[DD] - ratio[B] * valCorr[AA]);

    if (valFc <= 0)
    {
        qDebug() << "敏化发射强度非正，无法进行计算";
        calcProcess[CALC_DATA] = false;
        return;
    }

    // ED
    result[Ed] = valFc / (valFc + ratio[G] * valCorr[DD]);
    // Rc
    result[Rad] = ratio[K] * valCorr[AA] / (valFc / ratio[G] + valCorr[DD]);

    calcProcess[CALC_DATA] = true;
}

void FretCalculator::calc3CubeFret()
{
    valFc = valCorr[DA]
            - ratio[A] * (valCorr[AA] - ratio[C] * valCorr[DD])
            - ratio[D] * (valCorr[DD] - ratio[B] * valCorr[AA]);

    if (valFc <= 0)
    {
        qDebug() << "敏化发射强度非正，无法进行计算";
        calcProcess[CALC_DATA] = false;
        return;
    }

    // EA
    result[Ea] = valFc * ratio[Y] / (valCorr[AA] * ratio[A]);
    // 1/Rc
    result[Rda] = (valFc / ratio[G] + valCorr[DD]) / (ratio[K] * valCorr[AA]);

    calcProcess[CALC_DATA] = true;
}
void FretCalculator::calc2Hybrid()
{
    //计算估计浓度
    double Ma = ratio[G] * ratio[Y] / ratio[D];
    result[Dest] = ratio[D] * valCorr[DD] / (1 - result[Ed]);
    result[Aest] = ratio[A] * (valCorr[AA] - ratio[C] * valCorr[DD]) / Ma;
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
    if (ratioName == NaR)
    {
        calcProcess[SET_PARAM] = checkParamSet();
        return;
    }
    ratio[ratioName] = value;
    calcProcess[SET_PARAM] = checkParamSet();
}
void FretCalculator::setRatios(double a, double b, double c, double d, double g, double k, double y)
{
    setRatio(A, a);
    setRatio(B, b);
    setRatio(C, c);
    setRatio(D, d);
    setRatio(G, g);
    setRatio(K, k);
    setRatio(Y, y);
}

void FretCalculator::setExposureTime(ChannelName channelName, double value)
{
    if (channelName == NaC)
    {
        calcProcess[SET_PARAM] = checkParamSet();
        return;
    }
    exposureTime[channelName] = value;
    calcProcess[SET_PARAM] = checkParamSet();
}
void FretCalculator::setExposureTimes(double expsTimeAA, double expsTimeDA, double expsTimeDD)
{
    setExposureTime(AA, expsTimeAA);
    setExposureTime(DA, expsTimeDA);
    setExposureTime(DD, expsTimeDD);
}

void FretCalculator::loadBackGray(double I_AA, double I_DA, double I_DD)
{
    valBack[AA] = I_AA;
    valBack[DA] = I_DA;
    valBack[DD] = I_DD;
    calcProcess[LOAD_DATA] = checkDataLoaded();
}
void FretCalculator::loadForeGray(double I_AA, double I_DA, double I_DD)
{
    valFore[AA] = I_AA;
    valFore[DA] = I_DA;
    valFore[DD] = I_DD;
    calcProcess[LOAD_DATA] = checkDataLoaded();
}
void FretCalculator::loadCorrectedGray(double I_AA, double I_DA, double I_DD)
{
    valCorr[AA] = I_AA;
    valCorr[DA] = I_DA;
    valCorr[DD] = I_DD;
    calcProcess[CORRECT_DATA] = checkDataCorrected();
}

bool FretCalculator::checkDataLoaded()
{
    if (valFore[AA] <= 0 || valFore[DD] <= 0 || valFore[DA] <= 0)
    {
        return false;
    }
    if (valBack[AA] <= 0 || valBack[DD] <= 0 || valBack[DA] <= 0)
    {
        return false;
    }
    return true;
}

bool FretCalculator::checkDataCorrected()
{
    if (valCorr[AA] <= 0 || valCorr[DD] <= 0 || valCorr[DA] <= 0)
    {
        return false;
    }
    return true;
}

bool FretCalculator::checkParamSet()
{
    if (ratio[A] <= 0 || ratio[B] <= 0 || ratio[C] <= 0 || ratio[D] <= 0 || ratio[G] <= 0 || ratio[K] <= 0 || ratio[Y])
    {
        return false;
    }
    if (exposureTime[AA] <= 0 || exposureTime[DA] <= 0 || exposureTime[DD] <= 0)
    {
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
    qDebug() << "[A]:\t" << ratio[A];
    qDebug() << "[B]:\t" << ratio[B];
    qDebug() << "[C]:\t" << ratio[C];
    qDebug() << "[D]:\t" << ratio[D];
    qDebug() << "[G]:\t" << ratio[G];
    qDebug() << "[K]:\t" << ratio[K];
    qDebug() << "[Y]:\t" << ratio[Y];
    qDebug() << "[ExpAA]:\t" << exposureTime[AA];
    qDebug() << "[ExpDA]:\t" << exposureTime[DA];
    qDebug() << "[ExpDD]:\t" << exposureTime[DD];
}

double FretCalculator::getResult(CalcResult resultName)
{
    return result[resultName];
}

double FretCalculator::getFcValue()
{
    return valFc;
}
