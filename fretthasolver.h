
#ifndef FRETTHASOLVER_H
#define FRETTHASOLVER_H

#include <QDebug>
#include <QString>
#include <QDir>
#include <QThread>

#include <math.h>
#include <random>
#include <time.h>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <numeric>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <functional>
#include <iomanip>

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
    std::vector<double> vecFretData[8]; // 六个数组，保存 FRET 计算结果
    std::vector<double> vecFretDataBin[8];  // 六个数组，保存分箱后的 FRET 计算结果
    double resultsTha[4];   // 双杂交求解的参数结果
    double resultsThaBin[4];    // 分箱以后求解双杂交的参数结果
    double resultsEdRad[4]; // 线性拟合求解的参数结果，注意Kdeff是无效的
    double resultsEaRda[4]; // 线性拟合求解的参数结果，注意Kdeff是无效的

    double threshRatio[3] = {3.0, 3.0, 3.0};


    /* 成员方法 *****************************************************************************************/

    // 清除数据
    void clearGrayData();
    void clearFretData();
    void clearFretDataBin();
    // 保存数据
    QString outputData(QString savePath);
    // 添加数据
    void expandGrayData(double Iaa, double Ida, double Idd);
    void expandFretData(double valEd, double valEa, double valRad, double valRda, double valAest, double valDest);
    // 设置系数
    void setRatio(RatioName ratioName, double value);
    void setRatios(double a, double b, double c, double d, double G, double K, double Y);
    // 设置曝光时间
    void setExposureTime(ChannelName channelName, double value);
    void setExposureTimes(double expsTimeAA, double expsTimeDA, double expsTimeDD);
    // 设置信背比阈值系数
    void setThreshRatio(ChannelName channelName, double value);
    void setThreshRatios(double threshAA, double threshDA, double threshDD);
    // 读取数据
    void loadSourceData(QString folderPath);
    void loadSourceData(QString pathAA, QString pathDA, QString pathDD);
    // 主要活动
    void processImageToGrayData();
    void processGrayToFretData();
    // 设置文件路径
    void setBatchPath(QString);
    // 计算一个视野的数据到vecGrayData中
    void processOneView(QString viewFodoubleerPath);
    // 计算一批数据多个视野的数据到vecGrayData中
    void processOneBatch(QString batchFolderPath);
    // 执行拟合计算
    void performTwoHybridMatlab();  // matlab求解规划
    void performTwoHybridMatlabBin(double min, double max, double interval);
    void performTwoHybridLinear();  // 线性求解
    // 返回最值
    double maxData(CalcResult dataName);
    double maxBinData(CalcResult dataName);


    /* 数据处理的子函数*************************************************************************************************/

    // 计算斜率值
    double calcSlope(double min, double max, const std::vector<double>& R, const std::vector<double>& E);
    // 计算渐近值
    double calcApproach(double min, double max, const std::vector<double>& R, const std::vector<double>& E);
    // 最小二乘法
    double leastSquare(const std::vector<std::pair<double, double>>& points);
    // 数据封箱
    int binData(double min, double max, double interval);
    // 主活动，在多线程中执行
    void mainActivity();

protected:
    void run();
signals:
    //
    void thaFinished();
    // 定义信号发送进度信息
    void progressChanged(int value);

};

#endif // FRETTHASOLVER_H
