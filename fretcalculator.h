/* Main work flow is :
// SET_PARAM           ->               LOAD_DATA                  ->           CORRECT_DATA -> CALC_DATA
// Set Ratios -> Set Exposure Times ->  Load ForeGrounds -> Load BackGrounds -> Correct Data -> Calculate Data
*/
#ifndef FRETCALCULATOR_H
#define FRETCALCULATOR_H

#include <QDebug>
#include <enumerator.h>


class FretCalculator
{
public:

    FretCalculator();

    void resetCalculator();
    void resetData();

    bool checkReady();

    void setRatio(RatioName ratioName, double value);
    void setRatios(double a, double b, double c, double d, double G, double K, double Y);

    void setExposureTime(ChannelName channelName, double value);
    void setExposureTimes(double expsTimeAA, double expsTimeDA, double expsTimeDD);

    void loadBackGray(double I_AA, double I_DA, double I_DD);
    void loadForeGray(double I_AA, double I_DA, double I_DD);
    void loadCorrectedGray(double I_AA, double I_DA, double I_DD);

    void correctData();

    void calcEFret();
    void calc3CubeFret();
    void calc2Hybrid();
    void calcData();

    bool isParamSet();
    bool isDataLoaded();
    bool isDataCorrected();
    bool isDataCalculated();

    double getResult(CalcResult result);
    double getFcValue();

    bool isSuccess();

    void showParams();

private:

    bool calcProcess[4];
    double valFore[3];
    double valBack[3];
    double valCorr[3];
    double exposureTime[3];
    double ratio[7];
    double result[6];
    double valFc;

    bool checkParamSet();
    bool checkDataLoaded();
    bool checkDataCorrected();
};

#endif // FRETCALCULATOR_H
