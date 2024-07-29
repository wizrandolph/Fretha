
#include "fretimageprocessor.h"
#include "qdir.h"
#include <iostream>

namespace wiz {

    cv::Mat skel(cv::Mat img, cv::Mat element)
    {
        cv::Mat skel(img.size(), CV_8UC1, cv::Scalar(0));
        cv::Mat temp;
        cv::Mat eroded;

        bool done;
        do
        {
            cv::erode(img, eroded, element);
            cv::dilate(eroded, temp, element);
            cv::subtract(img, temp, temp);
            cv::bitwise_or(skel, temp, skel);
            eroded.copyTo(img);

            done = (cv::countNonZero(img) == 0);
        } while (!done);

        return skel;
    }
    //
    cv::Mat computeGradientSobel(
        const cv::Mat& src,
        int ddepth,
        int scale,
        int delta,
        int kernel_size
    ) {
        cv::Mat grad_x, grad_y;
        cv::Mat abs_grad_x, abs_grad_y;

        // Gradient X
        cv::Sobel(src, grad_x, ddepth, 1, 0, kernel_size, scale, delta, cv::BORDER_DEFAULT);
        cv::convertScaleAbs(grad_x, abs_grad_x);

        // Gradient Y
        cv::Sobel(src, grad_y, ddepth, 0, 1, kernel_size, scale, delta, cv::BORDER_DEFAULT);
        cv::convertScaleAbs(grad_y, abs_grad_y);

        cv::Mat grad;
        cv::addWeighted(abs_grad_x, 0.5, abs_grad_y, 0.5, 0, grad);

        return grad;
    }

    cv::Mat computeGradientLaplacian(
        const cv::Mat& src,
        int ddepth,
        int scale,
        int delta,
        int kernel_size
    ) {
        cv::Mat dst, abs_dst;

        cv::Laplacian(src, dst, ddepth, kernel_size, scale, delta, cv::BORDER_DEFAULT);
        cv::convertScaleAbs(dst, abs_dst);

        return abs_dst;
    }

    cv::Mat computeGradientCanny(
        const cv::Mat& src,
        double lowThreshold,
        int ratio,
        int kernel_size
    ) {
        cv::Mat edges;

        cv::Canny(src, edges, lowThreshold, lowThreshold*ratio, kernel_size);

        return edges;
    }

    cv::Mat getMorphologyClose(const cv::Mat& src, int kernelSize) {
        // 创建一个kernelSize x kernelSize的结构元素
        cv::Mat element = cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(kernelSize, kernelSize));

        // 进行闭运算
        cv::Mat close;
        cv::morphologyEx(src, close, cv::MORPH_CLOSE, element);

        return close;
    }
    cv::Mat getMorphologyOpen(const cv::Mat& src, int kernelSize) {
        // 创建一个kernelSize x kernelSize的结构元素
        cv::Mat element = cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(kernelSize, kernelSize));

        // 进行开运算
        cv::Mat open;
        cv::morphologyEx(src, open, cv::MORPH_OPEN, element);

        return open;
    }

    cv::Mat dilateImage(const cv::Mat& image, int elementSize) {
        // 创建一个结构元素
        cv::Mat element = cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(elementSize, elementSize));

        // 膨胀图像
        cv::Mat dilated;
        cv::dilate(image, dilated, element);

        return dilated;
    }

    cv::Mat erodeImage(const cv::Mat& image, int elementSize) {
        // 创建一个结构元素
        cv::Mat element = cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(elementSize, elementSize));

        // 膨胀图像
        cv::Mat eroded;
        cv::erode(image, eroded, element);

        return eroded;
    }

    /**
     * @brief FretImageProcessor::calcBackgroundGray
     * 计算图像背景灰度值
     * @param mat 图像
     * @return 返回最大频数对应的灰度值
     */
    double calcBackgroundGray(cv::Mat mat) {
        using namespace cv;

        //计算直方图
        Mat hist;   //直方图
        float range[2] = {1, 65535};  //统计灰度值范围
        const float* ranges[1] = {range};   //格式需要，指针的指针
        const int bins[1] = {65535};  //宽度，即直方图的柱数，即横轴的分布
        calcHist(&mat, 1, 0, Mat(), hist, 1, bins, ranges);  //计算直方图

        //直方图求峰值
        double minValue, maxValue;
        Point minIdx, maxIdx;
        minMaxLoc(hist, &minValue, &maxValue, &minIdx, &maxIdx);

        double valBack = maxIdx.y;   //直方图峰值

        return valBack;
    }

    cv::Mat getLargestConnectedComponent(const cv::Mat& binaryImage) {
        // 进行连通组件标记
        cv::Mat labels, stats, centroids;
        int numLabels = cv::connectedComponentsWithStats(binaryImage, labels, stats, centroids);

        // 找到最大连通域的索引和面积
        int maxArea = -1;
        int maxAreaIndex = -1;
        for (int i = 1; i < numLabels; i++) {
            int area = stats.at<int>(i, cv::CC_STAT_AREA);
            if (area > maxArea) {
                maxArea = area;
                maxAreaIndex = i;
            }
        }

        // 创建一个新的二值化图像，只包含最大连通域的区域
        cv::Mat largestComponent;
        cv::compare(labels, maxAreaIndex, largestComponent, cv::CMP_EQ);

        return largestComponent;
    }

    cv::Mat enhanceImage(const cv::Mat& img, int windowSize)
    {
        cv::Mat result = cv::Mat::zeros(img.size(), img.type());
        int halfSize = windowSize / 2;

        for(int i = 0; i < img.rows; i++)
        {
            for(int j = 0; j < img.cols; j++)
            {
                int rStart = std::max(i - halfSize, 0);
                int rEnd = std::min(i + halfSize, img.rows - 1);
                int cStart = std::max(j - halfSize, 0);
                int cEnd = std::min(j + halfSize, img.cols - 1);

                cv::Mat window = img(cv::Range(rStart, rEnd + 1), cv::Range(cStart, cEnd + 1));
                double minVal, maxVal;
                cv::minMaxLoc(window, &minVal, &maxVal);
                result.at<ushort>(i, j) = static_cast<ushort>((static_cast<double>(img.at<ushort>(i, j) - minVal) / (maxVal - minVal)) * 65535);
            }
        }

        return result;
    }

    cv::Mat watershedSegmentation(const cv::Mat& src) {
        // Convert to 8-bit
        cv::Mat src8bit = normalizeByMinMax8U(src);

        // Convert to 3-channel image
        cv::Mat srcColor;
        cv::cvtColor(src8bit, srcColor, cv::COLOR_GRAY2BGR);

        // Step 1: Convert to binary image
        cv::Mat binary;
        cv::threshold(src8bit, binary, 0, 255, cv::THRESH_BINARY_INV | cv::THRESH_OTSU);

        // Step 2: Remove noise
        cv::Mat kernel = cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(3, 3));
        cv::Mat opening;
        cv::morphologyEx(binary, opening, cv::MORPH_OPEN, kernel, cv::Point(-1, -1), 2);

        // Step 3: Compute distance transform
        cv::Mat distTransform;
        cv::distanceTransform(opening, distTransform, cv::DIST_L2, 5);
        cv::normalize(distTransform, distTransform, 0, 1., cv::NORM_MINMAX);

        // Step 4: Find sure foreground
        cv::Mat sureForeground;
        cv::threshold(distTransform, sureForeground, 0.5, 1., cv::THRESH_BINARY);

        // Step 5: Find unknown region
        cv::Mat sureForeground8U;
        sureForeground.convertTo(sureForeground8U, CV_8U, 255);
        cv::Mat unknown;
        cv::subtract(opening, sureForeground8U, unknown);

        // Step 6: Label markers
        cv::Mat markers;
        cv::connectedComponents(sureForeground8U, markers);
        markers = markers + 1;
        markers.setTo(0, unknown == 255);

        // Step 7: Apply watershed
        cv::watershed(srcColor, markers);

        return markers;
    }

    cv::Mat visualizeWatershed(const cv::Mat& markers) {
        // Create an output image
        cv::Mat dst = cv::Mat::zeros(markers.size(), CV_8UC3);

        // Create a random color for each label
        std::vector<cv::Vec3b> colors;
        int numLabels = markers.rows * markers.cols;
        for (int i = 0; i < numLabels; i++) {
            colors.push_back(cv::Vec3b(rand()&255, rand()&255, rand()&255));
        }

        // Fill the output image with the color of the corresponding label
        for (int i = 0; i < markers.rows; i++) {
            for (int j = 0; j < markers.cols; j++) {
                int index = markers.at<int>(i, j);
                if (index > 0 && index <= numLabels) {
                    dst.at<cv::Vec3b>(i, j) = colors[index - 1];
                } else {
                    dst.at<cv::Vec3b>(i, j) = cv::Vec3b(0, 0, 0);
                }
            }
        }

        return dst;
    }

    double calculateAverageGrayValue(const cv::Mat& grayImage, const cv::Mat& mask, cv::Point center, double distanceLimit) {
        // 确保灰度图像和掩膜尺寸相同
        if (grayImage.size() != mask.size()) {
            std::cerr << "Error: Gray image and mask should have the same size!" << std::endl;
            return 0.0;
        }

        int imageHeight = grayImage.rows;
        int imageWidth = grayImage.cols;
        double totalGrayValue = 0.0;
        int count = 0;

        // 遍历指定坐标周围的像素
        for (int i = 0; i < imageHeight; i++) {
            for (int j = 0; j < imageWidth; j++) {
                cv::Point currentPixel(j, i);
                double distance = cv::norm(currentPixel - center);

                // 检查距离是否小于给定上限且掩膜中像素值为255
                if (distance <= distanceLimit && mask.at<uchar>(i, j) == 255) {
                    totalGrayValue += grayImage.at<ushort>(i, j);
                    count++;
                }
            }
        }

        // 计算平均灰度值
        double averageGrayValue = (count > 0) ? (totalGrayValue / count) : 0.0;

        return averageGrayValue;
    }

    void setNegativeToZero(cv::Mat& inputImage)
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

    cv::Mat getPseudoColorImage(const cv::Mat& inputImage8UC1)
    {
        using namespace cv;
        Mat outputImage;

        if (inputImage8UC1.type() != CV_8UC1 || inputImage8UC1.empty()) {
            return outputImage;
        }
        applyColorMap(inputImage8UC1, outputImage, cv::COLORMAP_JET);

        return outputImage;
    }

    cv::Mat getGreenFireBlue(const cv::Mat& input_image)
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

    cv::Mat getFirePseudocolor(const cv::Mat& input_image) {
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
    cv::Mat applyMaskToImage(const cv::Mat& image, const cv::Mat& mask)
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

    cv::Mat medianFilter(const cv::Mat& src,
                                             int kernelSize) {
        cv::Mat dst;
        cv::medianBlur(src, dst, kernelSize);
        return dst;
    }

    cv::Mat medianFilter16U(const cv::Mat& src, int kernelSize) {
        cv::Mat dst = cv::Mat::zeros(src.size(), CV_16UC1);

        // 检查输入图像的类型和通道数
        if (src.type() != CV_16U || src.channels() != 1) {
            throw std::invalid_argument("Input image must be 16-bit single-channel.");
        }

        // 应用中值滤波
        cv::medianBlur(src, dst, kernelSize);

        return dst;
    }

    cv::Mat gaussianFilter(
    const cv::Mat& src,
    int kernelSize,
    double sigma) {

        cv::Mat dst;
        cv::GaussianBlur(src, dst, cv::Size(kernelSize, kernelSize), sigma);
        return dst;
    }

    cv::Mat gaussianFilter16U(const cv::Mat& src, int kernelSize, double sigma) {
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

    cv::Mat meanFilter(const cv::Mat& src,
                                           int kernelSize) {
        cv::Mat dst;
        cv::blur(src, dst, cv::Size(kernelSize, kernelSize));
        return dst;
    }

    cv::Mat meanFilter16U(const cv::Mat& src,
                                              int kernelSize) {
        cv::Mat dst = cv::Mat::zeros(src.size(), CV_16UC1);

        // 检查输入图像的类型和通道数
        if (src.type() != CV_16U || src.channels() != 1) {
            throw std::invalid_argument("Input image must be 16-bit single-channel.");
        }

        // 将输入图像转换为32位浮点数格式
        cv::Mat imageFloat;
        src.convertTo(imageFloat, CV_64F);

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
    cv::Mat localStandardDeviation(const cv::Mat& src,
                                   int kernelSize,
                                   bool normalize,
                                   cv::BorderTypes border) {

        cv::Mat converted;
        if (src.type() != CV_64FC1) {
            src.convertTo(converted, CV_64FC1);
        } else {
            src.copyTo(converted);
        }

        cv::Mat meanImg, meanImg_pow;
        cv::boxFilter(converted, meanImg, CV_64FC1, cv::Size(kernelSize, kernelSize));
        cv::pow(meanImg, 2, meanImg_pow);

        // 计算原图的平方图像及该图像的均值图像
        cv::Mat originalImg_pow, originalImg_pow_mean;
        cv::pow(converted, 2, originalImg_pow);
        cv::boxFilter(originalImg_pow, originalImg_pow_mean, CV_64FC1, cv::Size(kernelSize, kernelSize));

        // 作差然后开根号，得出标准差图像
        cv::Mat varianceImg, stdDevImg;    //方差图像

        cv::subtract(originalImg_pow_mean, meanImg_pow, varianceImg);
        cv::sqrt(varianceImg, stdDevImg);
        stdDevImg.convertTo(stdDevImg, CV_64FC1);

        // 计算标准差除以平均值的结果
        cv::Mat result;
        if (normalize) {
            cv::divide(stdDevImg, meanImg, result);
        } else {
            result = stdDevImg;
        }

        // 将边缘区域置为非特征
        int borderSize = kernelSize;
        cv::Mat nonEdgeRegion = result(cv::Rect(borderSize, borderSize, result.cols - 2 * borderSize, result.rows - 2 * borderSize));
        double minVal, maxVal;
        cv::minMaxLoc(nonEdgeRegion, &minVal, &maxVal);
        result.rowRange(0, borderSize).setTo(minVal);  // top edge
        result.rowRange(result.rows - borderSize, result.rows).setTo(minVal);  // bottom edge
        result.colRange(0, borderSize).setTo(minVal);  // left edge
        result.colRange(result.cols - borderSize, result.cols).setTo(minVal);  // right edge

        return result;
    }

    cv::Mat localStandardDeviationByChat(const cv::Mat& src, int kernelSize) {
        cv::Mat converted;
        if (src.type() != CV_64FC1) {
            src.convertTo(converted, CV_64FC1);
        } else {
            src.copyTo(converted, CV_64FC1);
        }

        cv::Mat meanImg, meanImg_pow;
        cv::boxFilter(converted, meanImg, CV_64FC1, cv::Size(kernelSize, kernelSize));
        cv::pow(meanImg, 2, meanImg_pow);

        // 计算原图的平方图像及该图像的均值图像
        cv::Mat originalImg_pow, originalImg_pow_mean;
        cv::pow(converted, 2, originalImg_pow);
        cv::boxFilter(originalImg_pow, originalImg_pow_mean, CV_64FC1, cv::Size(kernelSize, kernelSize));

        // 作差然后开根号，得出标准差图像
        cv::Mat varianceImg, stdDevImg;    //方差图像

        cv::subtract(originalImg_pow_mean, meanImg_pow, varianceImg);
        // varianceImg.setTo(0, varianceImg < 0);  // 将负值设置为0
        cv::sqrt(varianceImg, stdDevImg);
        stdDevImg.convertTo(stdDevImg, CV_64FC1);

        return stdDevImg;
    }

    double calcOtsuThreshold(const cv::Mat& image) {
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

    cv::Mat findLocalMinima(const cv::Mat& src)
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
                    sum += middleMatrixs[RATIO_NAME_K].at<uchar>(i, j);
                    if (middleMatrixs[RATIO_NAME_K].at<uchar>(i, j) == 0)
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

    cv::Mat findLocalMinimaLaplacian(const cv::Mat& src)
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

    cv::Mat detectEdgesBySobel(const cv::Mat& src, int kernelSize) {
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

    cv::Mat mergeChannels(cv::Mat blueChannel,
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
    cv::Mat normalizeByMinMax(const cv::Mat& image) {
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
    cv::Mat normalizeByZeroMax(const cv::Mat& image) {
        // 寻找输入图像的最小值和最大值
        double minValue, maxValue;
        cv::minMaxLoc(image, &minValue, &maxValue);

        return normalize(image, 0, maxValue);
    }


    cv::Mat normalize(const cv::Mat& image,
                                          double min,
                                          double max) {
        cv::Mat outputImage;
        double scale = 1.0 / (max - min);
        image.convertTo(outputImage, CV_64FC1, scale, -min * scale);

        return outputImage;
    }

    cv::Mat normalizeByMinMax16U(const cv::Mat& image) {
        // 寻找输入图像的最小值和最大值
        double minValue, maxValue;
        cv::minMaxLoc(image, &minValue, &maxValue);

        return normalize16U(image, minValue, maxValue);
    }

    cv::Mat normalizeByZeroMax16U(const cv::Mat& image) {
        // 寻找输入图像的最小值和最大值
        double minValue, maxValue;
        cv::minMaxLoc(image, &minValue, &maxValue);

        return normalize16U(image, 0, maxValue);
    }

    cv::Mat normalize16U(const cv::Mat& image,
                                             double min,
                                             double max) {
        cv::Mat outputImage;
        double scale = 65535.0 / (max - min);
        image.convertTo(outputImage, CV_16UC1, scale, -min * scale);

        return outputImage;
    }

    cv::Mat normalizeByMinMax8U(const cv::Mat& image) {
        // 寻找输入图像的最小值和最大值
        double minValue, maxValue;
        cv::minMaxLoc(image, &minValue, &maxValue);

        return normalize8U(image, minValue, maxValue);
    }

    cv::Mat normalizeByZeroMax8U(const cv::Mat& image) {
        // 寻找输入图像的最小值和最大值
        double minValue, maxValue;
        cv::minMaxLoc(image, &minValue, &maxValue);

        return normalize8U(image, 0, maxValue);
    }

    cv::Mat normalize8U(const cv::Mat& image,
                                            double min,
                                            double max) {
        cv::Mat outputImage;
        double scale = 255.0 / (max - min);
        image.convertTo(outputImage, CV_8UC1, scale, -min * scale);

        return outputImage;
    }

    cv::Mat removeScatter(const cv::Mat& image, int windowSize, double ratio) {
        int height = image.rows;
        int width = image.cols;
        int halfWindow = windowSize / 2;

        // 如果ratio是负数，那么自动计算ratio
        if (ratio < 0) {
            // 计算白色和黑色像素的数量
            int whitePixels = cv::countNonZero(image);
            int blackPixels = image.total() - whitePixels;
            ratio = static_cast<double>(whitePixels) / (whitePixels + blackPixels);
        }

        double threshold = 255.0 * (windowSize * windowSize) * ratio;

        cv::Mat dilatedImage(image.size(), CV_8UC1, cv::Scalar(0));

        // 创建积分图像
        cv::Mat integralImage;
        cv::integral(image, integralImage, CV_64F);

        for (int i = halfWindow; i < height - halfWindow; i++) {
            for (int j = halfWindow; j < width - halfWindow; j++) {
                // 使用积分图像计算窗口内的白色像素数量
                double total = integralImage.at<double>(i + halfWindow + 1, j + halfWindow + 1)
                               - integralImage.at<double>(i - halfWindow, j + halfWindow + 1)
                               - integralImage.at<double>(i + halfWindow + 1, j - halfWindow)
                               + integralImage.at<double>(i - halfWindow, j - halfWindow);
                if (total > threshold) {
                    dilatedImage.at<uchar>(i, j) = 255;
                }
            }
        }

        return dilatedImage;
    }

    cv::Mat fillHoles(const cv::Mat& binImg) {
        // 克隆输入图像
        cv::Mat temp = binImg.clone();

        // 反转图像
        cv::bitwise_not(temp, temp);

        // Floodfill从点(0, 0)开始
        cv::floodFill(temp, cv::Point(0,0), cv::Scalar(255));

        // 反转图像
        cv::bitwise_not(temp, temp);

        // 合并填充的洞和原始图像
        cv::Mat filledImg = (binImg | temp);

        return filledImg;
    }

    void adaptiveThreshold16(const cv::Mat& src, cv::Mat& dst, double maxValue, int adaptiveMethod, int thresholdType, int blockSize, double C) {
        cv::Mat mean;
        if (adaptiveMethod == cv::ADAPTIVE_THRESH_MEAN_C) {
            cv::blur(src, mean, cv::Size(blockSize, blockSize));
        } else if (adaptiveMethod == cv::ADAPTIVE_THRESH_GAUSSIAN_C) {
            mean = gaussianFilter16U(src, blockSize, (blockSize - 1) / 6.0);
        } else {
            throw std::invalid_argument("Invalid adaptiveMethod.");
        }

        // 计算阈值图像
        cv::Mat thresh;
        cv::subtract(mean, cv::Scalar(C), thresh);

        // 进行阈值处理
        if (thresholdType == cv::THRESH_BINARY) {
            cv::compare(src, thresh, dst, cv::CMP_GT);
        } else if (thresholdType == cv::THRESH_BINARY_INV) {
            cv::compare(src, thresh, dst, cv::CMP_LE);
        } else {
            throw std::invalid_argument("Only THRESH_BINARY and THRESH_BINARY_INV are supported.");
        }

        // 将结果乘以最大值
        dst.convertTo(dst, dst.type(), maxValue);
    }

} // namespace wiz end

FretImageProcessor::FretImageProcessor()
{

}

void FretImageProcessor::setRatio(RatioName ratioName, double value)
{
    if (ratioName == RATIO_NAME_SIZE)
    {
        calcProcess[SET_PARAM] = checkParamSet();
        return;
    }
    ratio[ratioName] = value;
    calcProcess[SET_PARAM] = checkParamSet();
}
void FretImageProcessor::setRatios(double a, double b, double c, double d, double g, double k, double y)
{
    setRatio(RATIO_NAME_A, a);
    setRatio(RATIO_NAME_B, b);
    setRatio(RATIO_NAME_C, c);
    setRatio(RATIO_NAME_D, d);
    setRatio(RATIO_NAME_G, g);
    setRatio(RATIO_NAME_K, k);
    setRatio(RATIO_NAME_Y, y);
}

void FretImageProcessor::setExposureTime(ChannelName channelName, double value)
{
    if (channelName == CHANNEL_NAME_SIZE)
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
    setExposureTime(CHANNEL_NAME_AA, expsTimeAA);
    setExposureTime(CHANNEL_NAME_DA, expsTimeDA);
    setExposureTime(CHANNEL_NAME_DD, expsTimeDD);
}

bool FretImageProcessor::checkDataLoaded() {
    if (matSrc[CHANNEL_NAME_AA].empty() || matSrc[CHANNEL_NAME_DA].empty() || matSrc[CHANNEL_NAME_DD].empty()) {
        return false;
    }
    return true;
}

bool FretImageProcessor::checkDataCorrected() {
    if (matCorr[CHANNEL_NAME_AA].empty() || matCorr[CHANNEL_NAME_DD].empty() || matCorr[CHANNEL_NAME_DA].empty())
    {
        return false;
    }
    return true;
}

bool FretImageProcessor::checkParamSet() {
    if (ratio[RATIO_NAME_A] <= 0 || ratio[RATIO_NAME_B] <= 0 || ratio[RATIO_NAME_C] <= 0 || ratio[RATIO_NAME_D] <= 0 || ratio[RATIO_NAME_G] <= 0 || ratio[RATIO_NAME_K] <= 0 || ratio[RATIO_NAME_Y]) {
        return false;
    }
    if (exposureTime[CHANNEL_NAME_AA] <= 0 || exposureTime[CHANNEL_NAME_DA] <= 0 || exposureTime[CHANNEL_NAME_DD] <= 0) {
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
    matSrc[CHANNEL_NAME_AA] = cv::imread(pathAA.toStdString(), -1);
    matSrc[CHANNEL_NAME_DA] = cv::imread(pathDA.toStdString(), -1);
    matSrc[CHANNEL_NAME_DD] = cv::imread(pathDD.toStdString(), -1);
    calcProcess[LOAD_DATA] = checkDataLoaded();
}

void FretImageProcessor::preProcessData() {
    using namespace cv;

    // 计算背景灰度值
    double valBack[3];
    for (int i = 0; i < 3; ++ i) {
        valBack[i] = wiz::calcBackgroundGray(matSrc[i]);
    }
    // 进行背景扣除
    correctData(valBack[CHANNEL_NAME_AA], valBack[CHANNEL_NAME_DA], valBack[CHANNEL_NAME_DD]);
    // 生成掩膜
    Mat maskESC[3]; //Masks for Each Single Channel
    double thre[3];
    for (int i = 0; i < 3; ++ i) {
        thre[i] = 2;
    }
    threshold(matSrc[CHANNEL_NAME_AA], maskESC[CHANNEL_NAME_AA], valBack[CHANNEL_NAME_AA] * thre[CHANNEL_NAME_AA], 255, THRESH_BINARY);
    threshold(matSrc[CHANNEL_NAME_DA], maskESC[CHANNEL_NAME_DA], valBack[CHANNEL_NAME_DA] * thre[CHANNEL_NAME_DA], 255, THRESH_BINARY);
    threshold(matSrc[CHANNEL_NAME_DD], maskESC[CHANNEL_NAME_DD], valBack[CHANNEL_NAME_DD] * thre[CHANNEL_NAME_DD], 255, THRESH_BINARY);
    for (int i = 0; i < 3; ++ i) {
        maskESC[i].convertTo(maskESC[i], CV_8U);
    }

    for (int i = 0; i < 3; ++ i) {
        matCorr[i] = matCorr[i] / exposureTime[i];
    }

    bitwise_and(maskESC[CHANNEL_NAME_AA], maskESC[CHANNEL_NAME_DA], mask);
    bitwise_and(mask, maskESC[CHANNEL_NAME_DD], mask);
}

void FretImageProcessor::correctData(double valBackAA,
                                     double valBackDA,
                                     double valBackDD) {
    for (int i = 0; i < 3; ++ i) {
        matSrc[i].convertTo(matCorr[i], CV_64F);
    }
    matCorr[CHANNEL_NAME_AA] = matCorr[CHANNEL_NAME_AA] - valBackAA;
    matCorr[CHANNEL_NAME_DA] = matCorr[CHANNEL_NAME_DA] - valBackDA;
    matCorr[CHANNEL_NAME_DD] = matCorr[CHANNEL_NAME_DD] - valBackDD;
    checkDataCorrected();
}

void FretImageProcessor::correctData(cv::Mat matBackAA,
                                     cv::Mat matBackDA,
                                     cv::Mat matBackDD) {
    for (int i = 0; i < 3; ++ i) {
        matSrc[i].convertTo(matCorr[i], CV_64F);
    }
    matCorr[CHANNEL_NAME_AA] = matCorr[CHANNEL_NAME_AA] - matBackAA;
    matCorr[CHANNEL_NAME_DA] = matCorr[CHANNEL_NAME_DA] - matBackDA;
    matCorr[CHANNEL_NAME_DD] = matCorr[CHANNEL_NAME_DD] - matBackDD;
}

void FretImageProcessor::calcEFret() {

    cv::Mat matFc = matCorr[CHANNEL_NAME_DA]
                    - ratio[RATIO_NAME_A] * (matCorr[CHANNEL_NAME_AA] - ratio[RATIO_NAME_C] * matCorr[CHANNEL_NAME_DD])
                    - ratio[RATIO_NAME_D] * (matCorr[CHANNEL_NAME_DD] - ratio[RATIO_NAME_B] * matCorr[CHANNEL_NAME_AA]);

    // ED
    matRst[Ed] = matFc / (matFc + ratio[RATIO_NAME_G] * matCorr[CHANNEL_NAME_DD]);
    // Rc
    matRst[Rad] = ratio[RATIO_NAME_K] * matCorr[CHANNEL_NAME_AA] / (matFc / ratio[RATIO_NAME_G] + matCorr[CHANNEL_NAME_DD]);
    if (matRst[Ed].empty() || matRst[Rad].empty()) {
        qDebug() << "E-FRET计算结果为空";
    }

    calcProcess[CALC_DATA] = true;
}

void FretImageProcessor::calc3CubeFret() {
    cv::Mat matFc = matCorr[CHANNEL_NAME_DA]
                    - ratio[RATIO_NAME_A] * (matCorr[CHANNEL_NAME_AA] - ratio[RATIO_NAME_C] * matCorr[CHANNEL_NAME_DD])
                    - ratio[RATIO_NAME_D] * (matCorr[CHANNEL_NAME_DD] - ratio[RATIO_NAME_B] * matCorr[CHANNEL_NAME_AA]);

    // EA
    matRst[Ea] = matFc * ratio[RATIO_NAME_Y] / (matCorr[CHANNEL_NAME_AA] * ratio[RATIO_NAME_A]);
    // 1/Rc
    matRst[Rda] = (matFc / ratio[RATIO_NAME_G] + matCorr[CHANNEL_NAME_DD]) / (ratio[RATIO_NAME_K] * matCorr[CHANNEL_NAME_AA]);

    calcProcess[CALC_DATA] = true;
}

/**
 * @brief FretImageProcessor::calc2Hybrid
 * 执行双杂交计算所需的矩阵计算
 * 保存结果在Dest和Aest
 */
void FretImageProcessor::calc2Hybrid() {
    // 计算估计浓度
    double Ma = ratio[RATIO_NAME_G] * ratio[RATIO_NAME_Y] / ratio[RATIO_NAME_D];
    matRst[Dest] = ratio[RATIO_NAME_D] * matCorr[CHANNEL_NAME_DD] / (1 - matRst[Ed]);
    matRst[Aest] = ratio[RATIO_NAME_A] * (matCorr[CHANNEL_NAME_AA] - ratio[RATIO_NAME_C] * matCorr[CHANNEL_NAME_DD]) / Ma;
}












/**
 * @brief FretImageProcessor::getBackgroundGray
 * 通过统计图像的直方图，将其中频数最高的灰度值视作背景灰度值
 * @param channelName 通道名
 * @return 灰度值
 */
double FretImageProcessor::getBackgroundGray(ChannelName channelName) {
    return wiz::calcBackgroundGray(matSrc[channelName]);
}

/**
 * @brief FretImageProcessor::getOtsuThreshold
 * 大津阈值法计算前后景的区分阈值
 * @param channelName 通道名
 * @return 灰度值
 */
double FretImageProcessor::getOtsuThreshold(ChannelName channelName) {
    return wiz::calcOtsuThreshold(matSrc[channelName]);
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

cv::Mat FretImageProcessor::getRawResult(CalcResult resultName) {

    cv::Mat result = matRst[resultName].clone();
    wiz::setNegativeToZero(result);

    return result;
}

/**
 * @brief FretImageProcessor::getMaskedResult
 * 使用掩膜的内容
 * @param resultName
 * @return 掩膜合并的结果
 */
cv::Mat FretImageProcessor::getMaskedResult(CalcResult resultName)
{
    cv::Mat mat = wiz::applyMaskToImage(matRst[resultName], mask);
    qDebug() << "getMaskedResult";
    return mat;
}

cv::Mat FretImageProcessor::getMaskedResult8U(CalcResult resultName)
{
    // Convert the input image to CV_8UC1 format
    cv::Mat outputImage;
    double minVal, maxVal;
    cv::Mat matMasked = wiz::applyMaskToImage(matRst[resultName], mask);
    wiz::setNegativeToZero(matMasked);
    cv::minMaxLoc(matMasked, &minVal, &maxVal);
    qDebug() << minVal << maxVal;
    cv::convertScaleAbs(matMasked, outputImage, 255.0 / (maxVal - minVal), -minVal * 255.0 / (maxVal - minVal));
    return outputImage;
}




QStringList FretImageProcessor::findImageFiles(const QString& folderPath,
                                               const QString& searchStr) {
    QStringList imageFiles;
    QDir folder(folderPath);
    QStringList filters;
    filters << "*.tif*";
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
        preProcessData();
        cv::Mat matNorm[3];
        for (int i = 0; i < 3; ++ i) {
            matNorm[i] = wiz::normalizeByZeroMax8U(matSrc[i]);
        }
        mergedImage = wiz::mergeChannels(matNorm[CHANNEL_NAME_AA], matNorm[CHANNEL_NAME_DA], matNorm[CHANNEL_NAME_DD]);
    }
    return mergedImage;
}

cv::Mat FretImageProcessor::getPseuImage(CalcResult resultName) {
    cv::Mat pseuImage;
    if (calcProcess[CALC_DATA]) {
        pseuImage = wiz::getPseudoColorImage(getMaskedResult8U(resultName));
        // pseuImage = getFirePseudocolor(getMaskedResult8U(resultName));
    }
    return pseuImage;
}

cv::Mat FretImageProcessor::getNormalizedImage(ChannelName channelName) {
    cv::Mat normalizedImage;
    if (calcProcess[CALC_DATA]) {
        normalizedImage = wiz::normalizeByMinMax8U(matSrc[channelName]);
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
