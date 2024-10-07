#ifndef FRETTHASOLVER_H
#define FRETTHASOLVER_H

#include <QDebug>
#include <QString>
#include <QDir>
#include <QThread>
#include <QMap>

#include <vector>

#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "fretcalculator.h"
#include "fretimageprocessor.h"
#include "enumerator.h"



class FretThaSolver : public QThread
{
    Q_OBJECT
public:
    FretThaSolver();

    /* 成员变量 ******************************************************************************************/

    // 路径变量
    QString batchFolderPath;    // 批量数据路径
    QString viewFolderPath; //视野数据路径
    // 计算执行单元
    FretImageProcessor *imageProcessor = new FretImageProcessor();  // FRET 图像处理单元
    FretCalculator *calculator = new FretCalculator();  // FRET 数值计算单元
    // 数据变量
    std::vector<double> vecGrayData[3]; // 三个数组，保存灰度值
    std::vector<double> vecFretData[6]; // 六个数组，保存 FRET 计算结果
    std::vector<double> vecFretDataBin[6]; // 六个数组，保存分箱后的 FRET 计算结果

    double resultsThaOurs[4]; // 优化的双杂交求解的参数结果
    double resultsThaBen[4]; // 分箱以后求解双杂交的参数结果
    double resultsEdRad[4]; // 线性拟合求解的参数结果，注意Kdeff是无效的
    double resultsEaRda[4]; // 线性拟合求解的参数结果，注意Kdeff是无效的
    double m_result_auto_edmax, m_result_auto_eamax, m_result_auto_stoic, m_result_mean_ratio; // 自动求解的结果
    double threshRatio[3] = {3.0, 3.0, 3.0};
    QString m_runFunction;

    Strategy m_strategy = Strategy::MODERATE;


    /***********************清除数据*************************/

    void clearGrayData();
    void clearFretData();
    void clearFretDataBin();

    /************************输出保存数据***************************/
    QString outputData(QString savePath);
    QString outputResults(QString savePath);

    void expandGrayData(double Iaa, double Ida, double Idd);
    void expandFretData(double valEd, double valEa, double valRad, double valRda, double valAest, double valDest);

    void setRatio(RatioName ratioName, double value);
    void setRatios(double a, double b, double c, double d, double G, double K, double Y);
    void setExposureTime(ChannelName channelName, double value);
    void setExposureTimes(double expsTimeAA, double expsTimeDA, double expsTimeDD);


    void setThreshRatio(ChannelName channelName, double value);
    void setThreshRatios(double threshAA, double threshDA, double threshDD);

    /******************************设置文件路径***************************/
    void setBatchPath(QString);
    void loadSourceData(QString folderPath);
    void loadSourceData(QString pathAA, QString pathDA, QString pathDD);

    // 主要活动
    void processOneView(QString viewFolderPath);
    void processOneBatch(QString batchFolderPath);
    void processImageToGrayData();
    void processGrayToFretData();
    // 生成数据
    void generateRoiFromBatch(QString batchFolderPath);
    void generateRoiFromView(QString viewFolderPath, QString folderName);
    void generateRoiFromImage(QString folderName);
    // 导入mask
    void loadRoiFromBatch(QString batchFolderPath);
    void loadRoiFromView(QString viewFolderPath, QString folderName);
    void loadRoiFromImage(QString folderName);


    /*********************************拟合计算函数***********************************/
    void performTwoHybridOurs();  // Dlib求解规划(优化的版本)
    void performTwoHybridBen(double min, double max, double interval, bool bin_flag = false);  // 原版求解规划
    void performTwoHybridDu(double minSlope , double maxSlope, double minAppro, double maxAppro);  // 线性求解
    void performTwoHybridDuAutoRange();


    // 返回最值
    double maxData(CalcResult dataName);
    double maxBinData(CalcResult dataName);
    std::pair<std::vector<double>, std::vector<double>> getDfreeAndAfreeOurs();
    std::pair<std::vector<double>, std::vector<double>> getDfreeAndAfreeBen();
    std::pair<std::vector<double>, std::vector<double>> getDfreeAndAfreeBenBin();


    void setRunFunction(QString str);
    void setStrategy(Strategy strategy);

private:
    /*********************************拟合功能函数********************************/

    double calcSlope(double min, double max, const std::vector<double>& R, const std::vector<double>& E);
    double calcApproach(double min, double max, const std::vector<double>& R, const std::vector<double>& E);
    void fitLangmiurDlib();
    void fitLangmiurMatlab();
    void fitLangmiurMatlabBin(double min, double max, double interval);  // matlab BIN 求解规划
    void lab_fitLangmiurMatlabBinCompare(double min, double max, double interval);

    int binDataByRda(double min, double max, double interval);
    int binDataByRad(double min, double max, double interval);
    int binDataOurs(double min, double max, double interval);

    void autoProcessActivity();
    void autoGenerateActivity();
    void autoImportActivity();

    std::vector<double> getRegionMeans(const cv::Mat& data, const cv::Mat& mask);

    // cv实验
    void generateRoiNew();
    void Segmentation();
    void scorePixel();
    void cellSegmentation();
    void extractRoiByAdaptiveMask(QString viewName);
    void extractRoiByFretImage(QString viewName);
    void extractRoiOrigin(QString viewName);

    // 自动EdRc
    double FindBestFittingRange(CalcResult r_type);

protected:
    void run();
signals:
    //
    void thaFinished();
    // 定义信号发送进度信息
    void progressChanged(int value);
    // 定义发送数据的信号
    void sendData(const QMap<TableHeader, QString> &data);

};

#endif // FRETTHASOLVER_H
