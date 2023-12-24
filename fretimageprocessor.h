
#ifndef FRETIMAGEPROCESSOR_H
#define FRETIMAGEPROCESSOR_H

#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "enumerator.h"

#include <QString>
#include <QDebug>

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

    static cv::Mat getPseudoColorImage(const cv::Mat& inputImage8UC1);
    static cv::Mat getGreenFireBlue(const cv::Mat& inputImage8UC1);
    static cv::Mat getFirePseudocolor(const cv::Mat& input_image);

    static cv::Mat normalizeByMinMax(const cv::Mat& image);
    static cv::Mat normalizeByZeroMax(const cv::Mat& image);
    static cv::Mat normalize(const cv::Mat& image, double min, double max);

    static cv::Mat normalizeByMinMax8U(const cv::Mat& image);
    static cv::Mat normalizeByZeroMax8U(const cv::Mat& image);
    static cv::Mat normalize8U(const cv::Mat& image, double min, double max);

    static cv::Mat normalizeByMinMax16U(const cv::Mat& image);
    static cv::Mat normalizeByZeroMax16U(const cv::Mat& image);
    static cv::Mat normalize16U(const cv::Mat& image, double min, double max);

    static cv::Mat applyMaskToImage(const cv::Mat& inputImage, const cv::Mat& inputMask);

    static void setNegativeToZero(cv::Mat& inputImage);
    // 三通道图像合并
    static cv::Mat mergeChannels(cv::Mat blueChannel, cv::Mat greenChannel, cv::Mat redChannel);
    // 中值滤波
    static cv::Mat medianFilter(const cv::Mat& src, int kernelSize);
    static cv::Mat medianFilter16U(const cv::Mat& src, int kernelSize);
    // 高斯滤波
    static cv::Mat gaussianFilter(const cv::Mat& src, int kernelSize, double sigma);
    static cv::Mat gaussianFilter16U(const cv::Mat& src, int kernelSize, double sigma);
    // 均值滤波
    static cv::Mat meanFilter(const cv::Mat& src, int kernelSize);
    static cv::Mat meanFilter16U(const cv::Mat& src, int kernelSize);
    static cv::Mat localStandardDeviation(const cv::Mat& src, int kernelSize);
    static cv::Mat localStandardDeviationByChat(const cv::Mat& src, int kernelSize);
    // Otsu算法
    static double calcOtsuThreshold(const cv::Mat& src);
    // 极小值点
    static cv::Mat findLocalMinima(const cv::Mat& src);
    static cv::Mat findLocalMinimaLaplacian(const cv::Mat& src);
    // 边缘检测
    static cv::Mat detectEdgesBySobel(const cv::Mat& src, int kernelSize);
    // 获取二值化图像中的最大连通域
    static cv::Mat getLargestConnectedComponent(const cv::Mat& binaryImage);
    static cv::Mat enhanceImage(const cv::Mat& img, int windowSize);
    // 直方图计算背景值
    static double calcBackgroundGray(cv::Mat mat);

    static double calculateAverageGrayValue(const cv::Mat& grayImage, const cv::Mat& mask, cv::Point center, double distanceLimit);

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
