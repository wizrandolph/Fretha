
#include "fretimageprocessor.h"
#include "qdir.h"
#include <iostream>

FretImageProcessor::FretImageProcessor()
{

}

void FretImageProcessor::setRatio(RatioName ratioName, double value)
{
    if (ratioName == NaR)
    {
        calcProcess[SET_PARAM] = checkParamSet();
        return;
    }
    ratio[ratioName] = value;
    calcProcess[SET_PARAM] = checkParamSet();
}
void FretImageProcessor::setRatios(double a, double b, double c, double d, double g, double k, double y)
{
    setRatio(A, a);
    setRatio(B, b);
    setRatio(C, c);
    setRatio(D, d);
    setRatio(G, g);
    setRatio(K, k);
    setRatio(Y, y);
}

void FretImageProcessor::setExposureTime(ChannelName channelName, double value)
{
    if (channelName == NaC)
    {
        calcProcess[SET_PARAM] = checkParamSet();
        return;
    }
    exposureTime[channelName] = value;
    calcProcess[SET_PARAM] = checkParamSet();
}
void FretImageProcessor::setExposureTimes(double expsTimeAA,
                                          double expsTimeDA,
                                          double expsTimeDD) {
    setExposureTime(AA, expsTimeAA);
    setExposureTime(DA, expsTimeDA);
    setExposureTime(DD, expsTimeDD);
}

bool FretImageProcessor::checkDataLoaded() {
    if (matSrc[AA].empty() || matSrc[DA].empty() || matSrc[DD].empty()) {
        return false;
    }
    return true;
}

bool FretImageProcessor::checkDataCorrected() {
    if (matCorr[AA].empty() || matCorr[DD].empty() || matCorr[DA].empty())
    {
        return false;
    }
    return true;
}

bool FretImageProcessor::checkParamSet() {
    if (ratio[A] <= 0 || ratio[B] <= 0 || ratio[C] <= 0 || ratio[D] <= 0 || ratio[G] <= 0 || ratio[K] <= 0 || ratio[Y]) {
        return false;
    }
    if (exposureTime[AA] <= 0 || exposureTime[DA] <= 0 || exposureTime[DD] <= 0) {
        return false;
    }
    return true;
}

/**
 * @brief FretImageProcessor::loadSourceData
 * @param folderPath
 */
void FretImageProcessor::loadSourceData(QString folderPath) {
    QStringList pathListAA = findImageFiles(folderPath, "AA.tif");
    QStringList pathListDA = findImageFiles(folderPath, "DA.tif");
    QStringList pathListDD = findImageFiles(folderPath, "DD.tif");
    if (pathListAA.size() != 1 || pathListDA.size() != 1 || pathListDD.size() != 1) {
        calcProcess[LOAD_DATA] = false;
        return;
    }
    QString pathAA = findImageFiles(folderPath, "AA.tif").at(0);
    QString pathDA = findImageFiles(folderPath, "DA.tif").at(0);
    QString pathDD = findImageFiles(folderPath, "DD.tif").at(0);
    loadSourceData(pathAA, pathDA, pathDD);
    viewPath = folderPath;
}

/**
 * @brief FretImageProcessor::loadSourceData
 * @param pathAA
 * @param pathDA
 * @param pathDD
 */
void FretImageProcessor::loadSourceData(QString pathAA,
                                        QString pathDA,
                                        QString pathDD) {
    matSrc[AA] = cv::imread(pathAA.toStdString(), -1);
    matSrc[DA] = cv::imread(pathDA.toStdString(), -1);
    matSrc[DD] = cv::imread(pathDD.toStdString(), -1);
    calcProcess[LOAD_DATA] = checkDataLoaded();
}

void FretImageProcessor::preProcessData() {
    using namespace cv;

    // 【计算背景灰度值】
    double valBack[3];
    for (int i = 0; i < 3; ++ i) {
        valBack[i] = calcBackgroundGray(matSrc[i]);
    }

    // 【进行背景扣除】
    correctData(valBack[AA], valBack[DA], valBack[DD]);

    // 【生成掩膜】
    // Calculate masks for each single channel image
    Mat maskESC[3]; //Masks for Each Single Channel
    double thre[3];
    for (int i = 0; i < 3; ++ i) {
        thre[i] = 2;
    }
    threshold(matSrc[AA], maskESC[AA], valBack[AA] * thre[AA], 255, THRESH_BINARY);
    threshold(matSrc[DA], maskESC[DA], valBack[DA] * thre[DA], 255, THRESH_BINARY);
    threshold(matSrc[DD], maskESC[DD], valBack[DD] * thre[DD], 255, THRESH_BINARY);
    // Convert to 8bit
    for (int i = 0; i < 3; ++ i) {
        maskESC[i].convertTo(maskESC[i], CV_8U);
    }
    // Save to member variable mask
    bitwise_and(maskESC[AA], maskESC[DA], mask);
    bitwise_and(mask, maskESC[DD], mask);
}

void FretImageProcessor::correctData(double valBackAA,
                                     double valBackDA,
                                     double valBackDD) {
    for (int i = 0; i < 3; ++ i) {
        matSrc[i].convertTo(matCorr[i], CV_64F);
    }
    matCorr[AA] = matCorr[AA] - valBackAA;
    matCorr[DA] = matCorr[DA] - valBackDA;
    matCorr[DD] = matCorr[DD] - valBackDD;
    checkDataCorrected();
}

void FretImageProcessor::correctData(cv::Mat matBackAA,
                                     cv::Mat matBackDA,
                                     cv::Mat matBackDD) {
    for (int i = 0; i < 3; ++ i) {
        matSrc[i].convertTo(matCorr[i], CV_64F);
    }
    matCorr[AA] = matCorr[AA] - matBackAA;
    matCorr[DA] = matCorr[DA] - matBackDA;
    matCorr[DD] = matCorr[DD] - matBackDD;
}

void FretImageProcessor::calcEFret() {

    cv::Mat matFc = matCorr[DA]
                    - ratio[A] * (matCorr[AA] - ratio[C] * matCorr[DD])
                    - ratio[D] * (matCorr[DD] - ratio[B] * matCorr[AA]);

    // ED
    matRst[Ed] = matFc / (matFc + ratio[G] * matCorr[DD]);
    // Rc
    matRst[Rad] = ratio[K] * matCorr[AA] / (matFc / ratio[G] + matCorr[DD]);
    if (matRst[Ed].empty() || matRst[Rad].empty()) {
        qDebug() << "E-FRET计算结果为空";
    }

    calcProcess[CALC_DATA] = true;
}

void FretImageProcessor::calc3CubeFret() {
    cv::Mat matFc = matCorr[DA]
                    - ratio[A] * (matCorr[AA] - ratio[C] * matCorr[DD])
                    - ratio[D] * (matCorr[DD] - ratio[B] * matCorr[AA]);

    // EA
    matRst[Ea] = matFc * ratio[Y] / (matCorr[AA] * ratio[A]);
    // 1/Rc
    matRst[Rda] = (matFc / ratio[G] + matCorr[DD]) / (ratio[K] * matCorr[AA]);

    calcProcess[CALC_DATA] = true;
}

/**
 * @brief FretImageProcessor::calc2Hybrid
 * 执行双杂交计算所需的矩阵计算
 * 保存结果在Dest和Aest
 */
void FretImageProcessor::calc2Hybrid() {
    //计算估计浓度
    double Ma = ratio[G] * ratio[Y] / ratio[D];
    matRst[Dest] = ratio[D] * matCorr[DD] / (1 - matRst[Ed]);
    matRst[Aest] = ratio[A] * (matCorr[AA] - ratio[C] * matCorr[DD]) / Ma;
}

/**
 * @brief FretImageProcessor::calcBackgroundGray
 * 计算图像背景灰度值
 * @param mat 图像
 * @return 返回最大频数对应的灰度值
 */
double FretImageProcessor::calcBackgroundGray(cv::Mat mat) {
    using namespace cv;

    //计算直方图
    Mat hist;   //直方图
    float range[2] = {1, 2000};  //统计灰度值范围
    const float* ranges[1] = {range};   //格式需要，指针的指针
    const int bins[1] = {2000};  //宽度，即直方图的柱数，即横轴的分布
    calcHist(&mat, 1, 0, Mat(), hist, 1, bins, ranges);  //计算直方图

    //直方图求峰值
    double minValue, maxValue;
    Point minIdx, maxIdx;
    minMaxLoc(hist, &minValue, &maxValue, &minIdx, &maxIdx);

    double valBack = maxIdx.y;   //直方图峰值

    return valBack;
}

/**
 * @brief FretImageProcessor::getBackgroundGray
 * 通过统计图像的直方图，将其中频数最高的灰度值视作背景灰度值
 * @param channelName 通道名
 * @return 灰度值
 */
double FretImageProcessor::getBackgroundGray(ChannelName channelName) {
    return calcBackgroundGray(matSrc[channelName]);
}

/**
 * @brief FretImageProcessor::getOtsuThreshold
 * 大津阈值法计算前后景的区分阈值
 * @param channelName 通道名
 * @return 灰度值
 */
double FretImageProcessor::getOtsuThreshold(ChannelName channelName) {
    return calcOtsuThreshold(matSrc[channelName]);
}

bool FretImageProcessor::isRectOutOfBounds(const cv::Rect& rect)
{
    cv::Size imageSize = matSrc[0].size();
    return (rect.x < 0 || rect.y < 0 || rect.x + rect.width > imageSize.width || rect.y + rect.height > imageSize.height);
}

/**
 * @brief FretImageProcessor::setRoi
 * 设置ROI的区域信息
 * @param x
 * @param y
 * @param w
 * @param h
 */
void FretImageProcessor::setRoi(double x, double y, double w, double h)
{
    cv::Rect cvRect(x, y, w, h);
    if (isRectOutOfBounds(cvRect)) {
        roi = cv::Rect(0, 0, 0, 0);
    }
    else {
        roi = cvRect;
    }
}

/**
 * @brief FretImageProcessor::getRoiGrayValue
 * @param channelName
 * @return 统计ROI区域内的平均灰度值
 */
double FretImageProcessor::getRoiGrayValue(ChannelName channelName)
{
    cv::Scalar mean = cv::mean(matSrc[channelName](roi));
    double averageGrayCV = mean.val[0];
    return averageGrayCV;
}

/**
 * @brief FretImageProcessor::getMaskedResult
 * 使用掩膜的内容
 * @param resultName
 * @return 掩膜合并的结果
 */
cv::Mat FretImageProcessor::getMaskedResult(CalcResult resultName)
{
    cv::Mat mat = applyMaskToImage(matRst[resultName], mask);
    qDebug() << "getMaskedResult";
    return mat;
}

cv::Mat FretImageProcessor::getMaskedResult8U(CalcResult resultName)
{
    // Convert the input image to CV_8UC1 format
    cv::Mat outputImage;
    double minVal, maxVal;
    cv::Mat matMasked = applyMaskToImage(matRst[resultName], mask);
    setNegativeToZero(matMasked);
    cv::minMaxLoc(matMasked, &minVal, &maxVal);
    qDebug() << minVal << maxVal;
    cv::convertScaleAbs(matMasked, outputImage, 255.0 / (maxVal - minVal), -minVal * 255.0 / (maxVal - minVal));
    return outputImage;
}

void FretImageProcessor::setNegativeToZero(cv::Mat& inputImage)
{
    // Iterate over all pixels in the input image
    for (int i = 0; i < inputImage.rows; i++) {
        for (int j = 0; j < inputImage.cols; j++) {
            double pixelValue = inputImage.at<double>(i, j);
            if (pixelValue < 0) {
                inputImage.at<double>(i, j) = 0;
            }
        }
    }
}

cv::Mat FretImageProcessor::getPseudoColorImage(const cv::Mat& inputImage8UC1)
{
    using namespace cv;
    Mat outputImage;

    if (inputImage8UC1.type() != CV_8UC1 || inputImage8UC1.empty()) {
        return outputImage;
    }
    applyColorMap(inputImage8UC1, outputImage, cv::COLORMAP_JET);

    return outputImage;
}

cv::Mat FretImageProcessor::getGreenFireBlue(const cv::Mat& input_image)
{
    // 创建调色板
    cv::Mat colors(256, 1, CV_8UC3);
    for (int i = 0; i < 256; i++)
    {
        colors.at<cv::Vec3b>(i, 0)[0] = 0;                                 // 蓝色通道
        colors.at<cv::Vec3b>(i, 0)[1] = static_cast<int>(pow(i, 0.7));    // 绿色通道
        colors.at<cv::Vec3b>(i, 0)[2] = 255;                               // 红色通道
    }

    // 应用颜色映射
    cv::Mat pseudocolor_image;
    cv::applyColorMap(input_image, pseudocolor_image, colors);

    return pseudocolor_image;
}

cv::Mat FretImageProcessor::getFirePseudocolor(const cv::Mat& input_image) {
    // 创建灰度图像
    cv::Mat grayscale_image;
    cv::cvtColor(input_image, grayscale_image, cv::COLOR_BGR2GRAY);

    // 创建自定义颜色映射表
    cv::Mat color_table(256, 1, CV_8UC3);
    for (int i = 0; i < 256; i++) {
        uchar* color_ptr = color_table.ptr<uchar>(i);
        if (i < 64) {
            color_ptr[0] = 0;
            color_ptr[1] = 0;
            color_ptr[2] = i * 4;
        } else if (i < 128) {
            color_ptr[0] = 0;
            color_ptr[1] = (i - 64) * 4;
            color_ptr[2] = 255;
        } else if (i < 192) {
            color_ptr[0] = (i - 128) * 4;
            color_ptr[1] = 255;
            color_ptr[2] = 255 - (i - 128) * 4;
        } else {
            color_ptr[0] = 255;
            color_ptr[1] = 255 - (i - 192) * 4;
            color_ptr[2] = 0;
        }
    }

    // 应用颜色映射表
    cv::Mat pseudocolor_image;
    cv::LUT(grayscale_image, color_table, pseudocolor_image);

    return pseudocolor_image;
}

/**
 * @brief FretImageProcessor::applyMaskToImage
 * @param image
 * @param mask
 * @return
 */
cv::Mat FretImageProcessor::applyMaskToImage(const cv::Mat& image,
                                             const cv::Mat& mask)
{
    CV_Assert(image.type() == CV_8UC1 || image.type() == CV_16UC1 || image.type() == CV_64FC1);
    CV_Assert(mask.type() == CV_8UC1);

    cv::Mat outputImage;
    image.copyTo(outputImage);

    for (int i = 0; i < mask.rows; i++) {
        for (int j = 0; j < mask.cols; j++) {
            if (mask.at<uchar>(i, j) == 0) {
                if (image.type() == CV_8UC1) {
                    outputImage.at<uchar>(i, j) = 0;
                } else if (image.type() == CV_16UC1) {
                    outputImage.at<ushort>(i, j) = 0;
                } else if (image.type() == CV_64FC1) {
                    outputImage.at<double>(i, j) = 0.0;
                }
            }
        }
    }

    return outputImage;
}

cv::Mat FretImageProcessor::medianFilter(const cv::Mat& src,
                                         int kernelSize) {
    cv::Mat dst;
    cv::medianBlur(src, dst, kernelSize);
    return dst;
}

cv::Mat FretImageProcessor::medianFilter16U(const cv::Mat& src,
                                            int kernelSize) {
    cv::Mat dst = cv::Mat::zeros(src.size(), CV_16UC1);

    // 检查输入图像的类型和通道数
    if (src.type() != CV_16U || src.channels() != 1) {
        throw std::invalid_argument("Input image must be 16-bit single-channel.");
    }

    // 应用中值滤波
    cv::medianBlur(src, dst, kernelSize);

    return dst;
}

cv::Mat FretImageProcessor::gaussianFilter(const cv::Mat& src,
                                           int kernelSize,
                                           double sigma) {
    cv::Mat dst;
    cv::GaussianBlur(src, dst, cv::Size(kernelSize, kernelSize), sigma);
    return dst;
}

cv::Mat FretImageProcessor::gaussianFilter16U(const cv::Mat& src, int kernelSize, double sigma) {
    cv::Mat dst = cv::Mat::zeros(src.size(), CV_16UC1);

    // 检查输入图像的类型和通道数
    if (src.type() != CV_16U || src.channels() != 1) {
        throw std::invalid_argument("Input image must be 16-bit single-channel.");
    }
    // 将输入图像转换为32位浮点数格式
    cv::Mat imageFloat;
    src.convertTo(imageFloat, CV_32F);
    // 应用高斯模糊
    cv::GaussianBlur(imageFloat, imageFloat, cv::Size(kernelSize, kernelSize), sigma);
    // 将结果图像转换回16位图像
    imageFloat.convertTo(dst, CV_16U);

    return dst;
}

cv::Mat FretImageProcessor::meanFilter(const cv::Mat& src,
                                       int kernelSize) {
    cv::Mat dst;
    cv::blur(src, dst, cv::Size(kernelSize, kernelSize));
    return dst;
}

cv::Mat FretImageProcessor::meanFilter16U(const cv::Mat& src,
                                          int kernelSize) {
    cv::Mat dst = cv::Mat::zeros(src.size(), CV_16UC1);

    // 检查输入图像的类型和通道数
    if (src.type() != CV_16U || src.channels() != 1) {
        throw std::invalid_argument("Input image must be 16-bit single-channel.");
    }

    // 将输入图像转换为32位浮点数格式
    cv::Mat imageFloat;
    src.convertTo(imageFloat, CV_32F);

    // 应用平均滤波
    cv::blur(imageFloat, imageFloat, cv::Size(kernelSize, kernelSize));

    // 将结果图像转换回16位图像
    imageFloat.convertTo(dst, CV_16U);

    return dst;
}

/**
 * @brief FretImageProcessor::standardDeviation
 * @param src
 * @param kernelSize
 * @return
 */
cv::Mat FretImageProcessor::localStandardDeviation(const cv::Mat& src,
                                                   int kernelSize) {

    cv::Mat converted;
    if (src.type() != CV_64FC1) {
        src.convertTo(converted, CV_64F);
    } else {
        src.copyTo(converted, CV_64F);
    }

    int kernelRadius = kernelSize / 2;
    cv::Mat standardDeviationImage(src.size(), CV_64FC1);  // 存储每个像素点的标准差

    for (int y = kernelRadius; y < src.rows - kernelRadius; ++y) {
        for (int x = kernelRadius; x < src.cols - kernelRadius; ++x) {
            cv::Rect roi(x - kernelRadius, y - kernelRadius, kernelSize, kernelSize);

            cv::Mat roiImage = converted(roi);
            cv::Scalar mean, stddev;
            cv::meanStdDev(roiImage, mean, stddev);

            double standardDeviation = stddev.val[0];
            // double average = mean.val[0];

            standardDeviationImage.at<double>(y, x) = standardDeviation ;
        }
    }

    return standardDeviationImage;
}

cv::Mat FretImageProcessor::localStandardDeviationByChat(const cv::Mat& src, int kernelSize)
{

    // 使用均值滤波计算图像的平均值
    cv::Mat meanImage;
    cv::blur(src, meanImage, cv::Size(kernelSize, kernelSize));

    // 计算图像的方差和标准差
    cv::Mat squaredDiff, meanSquaredDiff;
    cv::subtract(src, meanImage, squaredDiff);
    cv::multiply(squaredDiff, squaredDiff, meanSquaredDiff);
    cv::Mat variance, stdDeviation;
    //    cv::blur(meanSquaredDiff, variance, cv::Size(kernelSize, kernelSize));
    cv::sqrt(meanSquaredDiff, stdDeviation);
    return stdDeviation;
}

double FretImageProcessor::calcOtsuThreshold(const cv::Mat& image) {
    // 计算图像直方图
    int histSize = 65536;  // 16-bit图像的灰度级数
    float range[] = {0, 65536};  // 灰度范围
    const float* histRange = {range};
    cv::Mat hist;
    cv::calcHist(&image, 1, 0, cv::Mat(), hist, 1, &histSize, &histRange, true, false);

    // 统计像素总数
    int totalPixels = image.rows * image.cols;

    // 计算每个灰度级的像素数和灰度级的累积概率
    std::vector<double> pixelCount(histSize, 0.0);
    std::vector<double> cumulativeProbability(histSize, 0.0);
    for (int i = 0; i < histSize; i++) {
        pixelCount[i] = hist.at<float>(i);
        cumulativeProbability[i] = pixelCount[i] / totalPixels;
    }

    // 计算全局平均灰度
    double globalMean = 0.0;
    for (int i = 0; i < histSize; i++) {
        globalMean += i * cumulativeProbability[i];
    }

    // 计算类间方差和阈值
    double maxVariance = 0.0;
    double threshold = 0.0;
    double sum1 = 0.0;
    double sum2 = 0.0;
    double weight1 = 0.0;
    double mean1 = 0.0;
    double mean2 = 0.0;

    for (int i = 0; i < histSize; i++) {
        weight1 += cumulativeProbability[i];
        if (weight1 != 0) {
            sum1 += i * cumulativeProbability[i];
            mean1 = sum1 / weight1;
        }
        if (weight1 != 1) {
            sum2 = globalMean - sum1;
            mean2 = sum2 / (1 - weight1);
        }

        double currentVariance = weight1 * (1 - weight1) * (mean1 - mean2) * (mean1 - mean2);
        if (currentVariance > maxVariance) {
            maxVariance = currentVariance;
            threshold = i;
        }
    }

    return threshold;
}

cv::Mat FretImageProcessor::findLocalMinima(const cv::Mat& src)
{

    cv::Mat kernels[4];
    cv::Mat middleMatrixs[4];
    kernels[0] = (cv::Mat_<float>(5, 5)
                      << 1, 0, 0, 0, 0,
                  1, 1, 0, 0, 0,
                  1, 1, -8, 0, 0,
                  1, 1, 0, 0, 0,
                  1, 0, 0, 0, 0);
    kernels[1] = (cv::Mat_<float>(5, 5)
                      << 0, 0, 0, 0, 1,
                  0, 0, 0, 1, 1,
                  0, 0, -8, 1, 1,
                  0, 0, 0, 1, 1,
                  0, 0, 0, 0, 1);
    kernels[2] = (cv::Mat_<float>(5, 5)
                      << 1, 1, 1, 1, 1,
                  0, 1, 1, 1, 0,
                  0, 0, -8, 0, 0,
                  0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0);
    kernels[3] = (cv::Mat_<float>(5, 5)
                      << 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0,
                  0, 0, -8, 0, 0,
                  0, 1, 1, 1, 0,
                  1, 1, 1, 1, 1);

    for (int i = 0; i < 4; ++ i) {
        cv::filter2D(src, middleMatrixs[i], -1, kernels[i]);
    }

    cv::Mat dst(src.size(), CV_8U);

    for (int i = 0; i < middleMatrixs[0].rows; ++ i)
        for (int j = 0; j < middleMatrixs[0].cols; ++ j) {
            bool flag = false;    // 检测是否有0
            int sum = 0;
            for (int k = 0; k < 4; ++ k) {
                sum += middleMatrixs[k].at<uchar>(i, j);
                if (middleMatrixs[k].at<uchar>(i, j) == 0)
                    flag = true; // 找到了
            }
            if (flag) {
                dst.at<uchar>(i, j) = 0;
            } else {
                dst.at<uchar>(i, j) = sum / 4;
            }
        }

    return dst;
}

cv::Mat FretImageProcessor::findLocalMinimaLaplacian(const cv::Mat& src)
{

    cv::Mat kernel;
    kernel = (cv::Mat_<float>(5, 5)
                      <<  1, 1, 1, 1, 1,
                          1, 1, 1, 1, 1,
                          1, 1, -24, 1, 1,
                          1, 1, 1, 1, 1,
                          1, 1, 1, 1, 1);
    cv::Mat dst;
    cv::filter2D(src, dst, -1, kernel);

    return dst;
}

cv::Mat FretImageProcessor::detectEdgesBySobel(const cv::Mat& src, int kernelSize) {
    // 检查是否是多通道图像
    cv::Mat gray;
    if (src.channels() == 1) {
        // 如果是单通道图像，则直接使用作为灰度图像处理
        gray = src.clone();
    } else {
        // 如果是多通道图像，则将其转换为灰度图像
        cv::cvtColor(src, gray, cv::COLOR_BGR2GRAY);
    }

    // 高斯滤波，降噪
    cv::Mat blurred;
    cv::GaussianBlur(gray, blurred, cv::Size(kernelSize, kernelSize), 0);

    // 计算水平和垂直梯度
    cv::Mat gradX, gradY;
    cv::Sobel(blurred, gradX, CV_16S, 1, 0);
    cv::Sobel(blurred, gradY, CV_16S, 0, 1);

    // 转换为绝对值并归一化梯度图像
    cv::Mat absGradX, absGradY;
    cv::convertScaleAbs(gradX, absGradX);
    cv::convertScaleAbs(gradY, absGradY);

    // 合并梯度分量
    cv::Mat edges;
    cv::addWeighted(absGradX, 0.5, absGradY, 0.5, 0, edges);

    return edges;
}

cv::Mat FretImageProcessor::mergeChannels(cv::Mat blueChannel,
                                          cv::Mat greenChannel,
                                          cv::Mat redChannel) {
    // 创建一个空的三通道图像
    cv::Mat mergedImage(blueChannel.rows, blueChannel.cols, CV_8UC3);

    // 将三个单通道图像合并到一个三通道图像中
    std::vector<cv::Mat> channels;
    channels.push_back(blueChannel);
    channels.push_back(greenChannel);
    channels.push_back(redChannel);
    cv::merge(channels, mergedImage);

    return mergedImage;
}

/**
 * @brief FretImageProcessor::normalizeByMinMax
 Perform normalization by automatically setting the range from minimum gray to the maximum gray
 * @param image 输入的整数
 * @return normalized 8-bit image
 */
cv::Mat FretImageProcessor::normalizeByMinMax(const cv::Mat& image) {
    // 寻找输入图像的最小值和最大值
    double minValue, maxValue;
    cv::minMaxLoc(image, &minValue, &maxValue);

    return normalize(image, minValue, maxValue);
}

/**
 * @brief FretImageProcessor::normalizeByZeroMax
 * Perform normalization by automatically setting the range from 0 to the maximum gray
 * @param image
 * @return normalized 8-bit image
 */
cv::Mat FretImageProcessor::normalizeByZeroMax(const cv::Mat& image) {
    // 寻找输入图像的最小值和最大值
    double minValue, maxValue;
    cv::minMaxLoc(image, &minValue, &maxValue);

    return normalize(image, 0, maxValue);
}


cv::Mat FretImageProcessor::normalize(const cv::Mat& image,
                                      double min,
                                      double max) {
    cv::Mat outputImage;
    double scale = 1.0 / (max - min);
    image.convertTo(outputImage, CV_64FC1, scale, -min * scale);

    return outputImage;
}

cv::Mat FretImageProcessor::normalizeByMinMax16U(const cv::Mat& image) {
    // 寻找输入图像的最小值和最大值
    double minValue, maxValue;
    cv::minMaxLoc(image, &minValue, &maxValue);

    return normalize16U(image, minValue, maxValue);
}

cv::Mat FretImageProcessor::normalizeByZeroMax16U(const cv::Mat& image) {
    // 寻找输入图像的最小值和最大值
    double minValue, maxValue;
    cv::minMaxLoc(image, &minValue, &maxValue);

    return normalize16U(image, 0, maxValue);
}

cv::Mat FretImageProcessor::normalize16U(const cv::Mat& image,
                                         double min,
                                         double max) {
    cv::Mat outputImage;
    double scale = 65535.0 / (max - min);
    image.convertTo(outputImage, CV_16UC1, scale, -min * scale);

    return outputImage;
}

cv::Mat FretImageProcessor::normalizeByMinMax8U(const cv::Mat& image) {
    // 寻找输入图像的最小值和最大值
    double minValue, maxValue;
    cv::minMaxLoc(image, &minValue, &maxValue);

    return normalize8U(image, minValue, maxValue);
}

cv::Mat FretImageProcessor::normalizeByZeroMax8U(const cv::Mat& image) {
    // 寻找输入图像的最小值和最大值
    double minValue, maxValue;
    cv::minMaxLoc(image, &minValue, &maxValue);

    return normalize8U(image, 0, maxValue);
}

cv::Mat FretImageProcessor::normalize8U(const cv::Mat& image,
                                        double min,
                                        double max) {
    cv::Mat outputImage;
    double scale = 255.0 / (max - min);
    image.convertTo(outputImage, CV_8UC1, scale, -min * scale);

    return outputImage;
}

QStringList FretImageProcessor::findImageFiles(const QString& folderPath,
                                               const QString& searchStr) {
    QStringList imageFiles;
    QDir folder(folderPath);
    QStringList filters;
    filters << "*.tif";
    folder.setNameFilters(filters);
    folder.setFilter(QDir::Files | QDir::NoSymLinks);

    QFileInfoList fileList = folder.entryInfoList();
    for (int i = 0; i < fileList.size(); ++i) {
        QFileInfo fileInfo = fileList.at(i);
        QString fileName = fileInfo.fileName();
        if (fileName.contains(searchStr)) {
            imageFiles.append(fileInfo.absoluteFilePath());
        }
    }

    return imageFiles;
}

cv::Mat FretImageProcessor::getMergedImage() {
    cv::Mat mergedImage;
    if (calcProcess[LOAD_DATA]) {
        cv::Mat matNorm[3];
        for (int i = 0; i < 3; ++ i) {
            matNorm[i] = normalizeByZeroMax8U(matSrc[i]);
        }
        mergedImage = mergeChannels(matNorm[AA], matNorm[DA], matNorm[DD]);
    }
    return mergedImage;
}

cv::Mat FretImageProcessor::getPseuImage(CalcResult resultName) {
    cv::Mat pseuImage;
    if (calcProcess[CALC_DATA]) {
        pseuImage = getPseudoColorImage(getMaskedResult8U(resultName));
        // pseuImage = getFirePseudocolor(getMaskedResult8U(resultName));
    }
    return pseuImage;
}

cv::Mat FretImageProcessor::getNormalizedImage(ChannelName channelName) {
    cv::Mat normalizedImage;
    if (calcProcess[CALC_DATA]) {
        normalizedImage = normalizeByMinMax8U(matSrc[channelName]);
    }
    return normalizedImage;
}

cv::Mat FretImageProcessor::getMatrixCopy(ChannelName channelName)
{
    cv::Mat rst = matSrc[channelName].clone();
    return rst;
}

bool FretImageProcessor::isParamSet()
{
    return checkParamSet();
}

bool FretImageProcessor::isDataLoaded()
{
    return checkDataLoaded();
}

bool FretImageProcessor::isDataCorrected()
{
    return checkDataCorrected();
}
