
#ifndef FRETIMAGEPROCESSOR_H
#define FRETIMAGEPROCESSOR_H

#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "enumerator.h"

#include <QString>
#include <QDebug>

namespace wiz {
    cv::Mat getMorphologyClose(const cv::Mat& src, int kernelSize);
    cv::Mat getMorphologyOpen(const cv::Mat& src, int kernelSize);

    cv::Mat dilateImage(const cv::Mat& image, int elementSize);
    cv::Mat erodeImage(const cv::Mat& image, int elementSize);

    cv::Mat getPseudoColorImage(const cv::Mat& inputImage8UC1);
    cv::Mat getGreenFireBlue(const cv::Mat& inputImage8UC1);
    cv::Mat getFirePseudocolor(const cv::Mat& input_image);

    cv::Mat normalizeByMinMax(const cv::Mat& image);
    cv::Mat normalizeByZeroMax(const cv::Mat& image);
    cv::Mat normalize(const cv::Mat& image, double min, double max);

    cv::Mat normalizeByMinMax8U(const cv::Mat& image);
    cv::Mat normalizeByZeroMax8U(const cv::Mat& image);
    cv::Mat normalize8U(const cv::Mat& image, double min, double max);

    cv::Mat normalizeByMinMax16U(const cv::Mat& image);
    cv::Mat normalizeByZeroMax16U(const cv::Mat& image);
    cv::Mat normalize16U(const cv::Mat& image, double min, double max);

    cv::Mat applyMaskToImage(const cv::Mat& inputImage, const cv::Mat& inputMask);

    void setNegativeToZero(cv::Mat& inputImage);
    // 三通道图像合并
    cv::Mat mergeChannels(cv::Mat blueChannel, cv::Mat greenChannel, cv::Mat redChannel);
    // 中值滤波
    cv::Mat medianFilter(const cv::Mat& src, int kernelSize);
    cv::Mat medianFilter16U(const cv::Mat& src, int kernelSize);
    // 高斯滤波
    cv::Mat gaussianFilter(const cv::Mat& src, int kernelSize, double sigma);
    cv::Mat gaussianFilter16U(const cv::Mat& src, int kernelSize, double sigma);
    // 均值滤波
    cv::Mat meanFilter(const cv::Mat& src, int kernelSize);
    cv::Mat meanFilter16U(const cv::Mat& src, int kernelSize);
    cv::Mat localStandardDeviation(const cv::Mat& src, int kernelSize);
    cv::Mat localStandardDeviationByChat(const cv::Mat& src, int kernelSize);
    // Otsu算法
    double calcOtsuThreshold(const cv::Mat& src);
    // 极小值点
    cv::Mat findLocalMinima(const cv::Mat& src);
    cv::Mat findLocalMinimaLaplacian(const cv::Mat& src);
    // 边缘检测
    cv::Mat detectEdgesBySobel(const cv::Mat& src, int kernelSize);
    // 获取二值化图像中的最大连通域
    cv::Mat getLargestConnectedComponent(const cv::Mat& binaryImage);
    cv::Mat enhanceImage(const cv::Mat& img, int windowSize);

    // 学习图像算法
    cv::Mat watershedSegmentation(const cv::Mat& src);   // 分水岭算法
    cv::Mat visualizeWatershed(const cv::Mat& markers);  //可视化

    // 直方图计算背景值
    double calcBackgroundGray(cv::Mat mat);
    double calculateAverageGrayValue(const cv::Mat& grayImage, const cv::Mat& mask, cv::Point center, double distanceLimit);

    cv::Mat removeScatter(const cv::Mat& image, int windowSize, double ratio = -1);    // 移除散点
}

class FretImageProcessor
{
public:
    FretImageProcessor();

    void setRatio(RatioName ratioName, double value);
    void setRatios(double a, double b, double c, double d, double G, double K, double Y);

    void setExposureTime(ChannelName channelName, double value);
    void setExposureTimes(double expsTimeAA, double expsTimeDA, double expsTimeDD);

    void loadSourceData(QString folderPath);
    void loadSourceData(QString pathAA, QString pathDA, QString pathDD);

    void preProcessData();
    void correctData(double valBackAA, double valBackDA, double valBackDD);
    void correctData(cv::Mat matBackAA, cv::Mat matBackDA, cv::Mat matBackDD);    // 预留，未实现

    void calcEFret();
    void calc3CubeFret();
    void calc2Hybrid();

    double getBackgroundGray(ChannelName);
    double getOtsuThreshold(ChannelName);

    bool isRectOutOfBounds(const cv::Rect& rect);
    void setRoi(double x, double y, double w, double h);
    double getRoiGrayValue(ChannelName);

    cv::Mat getMaskedResult(CalcResult resultName);
    cv::Mat getMaskedResult8U(CalcResult resultName);
    cv::Mat getMergedImage();
    cv::Mat getPseuImage(CalcResult resultName);
    cv::Mat getNormalizedImage(ChannelName channelName);
    cv::Mat getMatrixCopy(ChannelName channelName);

    bool isParamSet();
    bool isDataLoaded();
    bool isDataCorrected();

    static QStringList findImageFiles(const QString& folderPath, const QString& searchStr);



private:
    QString viewPath;
    bool calcProcess[4];
    cv::Mat matSrc[3];
    cv::Mat matCorr[3];
    double exposureTime[3];
    double ratio[7];
    // Ed, Rad, Ea, Rda, Aest, Dest, Afree, Dfree
    cv::Mat matRst[6];
    cv::Mat matTmp;
    cv::Mat mask;
    cv::Rect roi;
    bool checkParamSet();
    bool checkDataLoaded();
    bool checkDataCorrected();

};

#endif // FRETIMAGEPROCESSOR_H
