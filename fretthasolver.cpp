
#include "fretthasolver.h"
#include "TwoHybridSolver/TwoHybridSolver.h"
#include "TwoHybridSolver/TwoHybridSolver_terminate.h"
#include "TwoHybridSolver/coder_array.h"

#include <dlib/optimization.h>
#include <vector>
#include <cmath>
#include <QFile>
#include <QTextStream>
#include <QStringList>
#include <QDateTime>

struct WizPoint {
    int x;
    int y;
};

// 线性回归 最小二乘
double leastSquare(const std::vector<std::pair<double, double>>& points)
{
    double sumXY = 0, sumXX = 0;
    for (uint i = 0; i < points.size(); ++ i) {
        sumXY += points[i].first * points[i].second;
        sumXX += points[i].first * points[i].first;
    }
    return sumXY / sumXX;
}

// 基于dlib的数据拟合算法
typedef dlib::matrix<double,0,1> column_vector;

double gaussian_pdf(double x, double mean, double stddev) {
    return std::exp(-0.5 * std::pow((x - mean) / stddev, 2)) / (stddev * std::sqrt(2 * M_PI));
}
struct MyData {
    std::vector<double> AEST;
    std::vector<double> DEST;
    std::vector<double> EACORR;
    std::vector<double> EDCORR;
    std::vector<double> weights;  // 新增权重数组
};

double HuberLoss(double x) {
    double delta = 0.01;
    if (std::abs(x) <= delta) {
        return 0.5 * x * x;
    } else {
        return delta * (std::abs(x) - 0.5 * delta);
    }
}

double Error(double x) {
    // Huber loss
    return HuberLoss(x);
}

double Loss(MyData& data, const column_vector& X) {
    std::vector<double> Dfree(data.AEST.size());
    std::vector<double> Afree(data.AEST.size());
    std::vector<double> EAPRED(data.AEST.size());
    std::vector<double> EDPRED(data.AEST.size());

    for (uint i = 0; i < data.AEST.size(); ++i) {
        Dfree[i] = ((data.DEST[i]-X(0)-data.AEST[i]*X(1))+std::sqrt(std::pow(data.DEST[i]-X(0)-data.AEST[i]*X(1), 2)+4*X(0)*data.DEST[i]))/2;
        Afree[i] = data.AEST[i]-(data.DEST[i]-Dfree[i])/X(1);
        EAPRED[i] = X(2)*Dfree[i]/(Dfree[i]+X(0));
        EDPRED[i] = X(3)*Afree[i]/(Afree[i]+X(0)/X(1));
    }

    // 计算误差
    std::vector<double> EaError(data.AEST.size());
    std::vector<double> EdError(data.AEST.size());
    for (uint i = 0; i < data.AEST.size(); ++i) {
        EaError[i] = Error(data.EACORR[i] - EAPRED[i]);
        EdError[i] = Error(data.EDCORR[i] - EDPRED[i]);
    }

    double SUMEAERROR = 0.0;
    double SUMEDERROR = 0.0;
    for (uint i = 0; i < data.AEST.size(); ++i) {
        double weight = data.weights[i];
        // weight = 1.0;
        SUMEAERROR += EaError[i] * weight;
        SUMEDERROR += EdError[i] * weight;
    }

    // std::cout << "[Current Loss]:\t"
    //           << "Total Error: " << SUMEAERROR + SUMEDERROR << " "
    //           << "Data Size: " << data.AEST.size() << " "
    //           << "Mean Loss: " << (std::accumulate(EaError.begin(), EaError.end(), 0.0) + std::accumulate(EdError.begin(), EdError.end(), 0.0)) / 2 / data.AEST.size() << "\t";
    // for (uint i = 0; i < 4; ++ i) {
    //     std::cout << X(i) << " ";
    // }
    // std::cout << std::endl;

    return SUMEAERROR + SUMEDERROR;
}

double Objective(const column_vector& x, void *my_func_data) {
    auto data = static_cast<MyData*>(my_func_data);
    return Loss(*data, x);
}

column_vector TwoHybridSolver(const std::vector<double>& AAEST,
                              const std::vector<double>& DDEST,
                              const std::vector<double>& EACORRR,
                              const std::vector<double>& EDCORRR,
                              bool weight_flag) {
    column_vector starting_point(4);
    starting_point = 1, 1, 0.5, 0.5;

    // data weight calculation
    std::vector<double> da_ratios(AAEST.size());
    for (uint i = 0; i < AAEST.size(); ++i) {
        da_ratios[i] = std::log(DDEST[i] / AAEST[i]);
    }
    double mean_da_ratio = std::accumulate(da_ratios.begin(), da_ratios.end(), 0.0) / da_ratios.size();
    double stddev_da_ratio = std::sqrt(std::accumulate(
                                           da_ratios.begin(),
                                           da_ratios.end(),
                                           0.0,
                                           [mean_da_ratio](double sum, double ratio) {
                                           return sum + std::pow(ratio - mean_da_ratio, 2);
                                       }) / da_ratios.size());
    std::cout << "[Mean DA Ratio]:\t" << mean_da_ratio << std::endl;
    std::cout << "[Stddev DA Ratio]:\t" << stddev_da_ratio << std::endl;
    std::vector<double> weights(AAEST.size());
    for (uint i = 0; i < AAEST.size(); ++i) {
        if (weight_flag) {
            if (std::abs(da_ratios[i] - mean_da_ratio) > 3 * stddev_da_ratio)
                weights[i] = 1;
            else
                weights[i] = 1 / gaussian_pdf(da_ratios[i], mean_da_ratio, stddev_da_ratio);
        } else {
            weights[i] = 1;
        }
    } // end data weight calculation

    MyData data = {AAEST, DDEST, EACORRR, EDCORRR, weights};  // 添加数据和权重

    try {
        auto objective_wrapper = [&data](const column_vector& x) {
            return Objective(x, &data);
        };

        dlib::find_min_using_approximate_derivatives(dlib::bfgs_search_strategy(),
                                                     dlib::objective_delta_stop_strategy(1e-7),
                                                     objective_wrapper,
                                                     starting_point,
                                                     -1,
                                                     0.01);
    } catch (std::exception& e) {
        std::cout << e.what() << std::endl;
    }

    return starting_point;
} // end 拟合算法

// 质心算法
double calculateDistance(const WizPoint& p1, const WizPoint& p2) {
    int dx = p2.x - p1.x;
    int dy = p2.y - p1.y;
    return std::sqrt(dx * dx + dy * dy);
}

WizPoint getCentroid(const std::vector<WizPoint>& pixels) {
    double minDistanceSum = std::numeric_limits<double>::max();
    WizPoint centroid;

    for (const auto& currentPixel : pixels) {
        double distanceSum = 0.0;

        for (const auto& otherPixel : pixels) {
            if (&currentPixel != &otherPixel) {
                distanceSum += calculateDistance(currentPixel, otherPixel);
            }
        }

        if (distanceSum < minDistanceSum) {
            minDistanceSum = distanceSum;
            centroid = currentPixel;
        }
    }

    return centroid;
}

void markVisited(cv::Mat& mask, int x, int y, std::vector<WizPoint>& component) {
    WizPoint point;
    point.x = x;
    point.y = y;
    component.push_back(point);
    mask.at<uchar>(y, x) = 0;

    if (x > 0 && mask.at<uchar>(y, x - 1) == 255)
        markVisited(mask, x - 1, y, component); // 左
    if (x < mask.cols - 1 && mask.at<uchar>(y, x + 1) == 255)
        markVisited(mask, x + 1, y, component); // 右
    if (y > 0 && mask.at<uchar>(y - 1, x) == 255)
        markVisited(mask, x, y - 1, component); // 上
    if (y < mask.rows - 1 && mask.at<uchar>(y + 1, x) == 255)
        markVisited(mask, x, y + 1, component); // 下
}

std::vector<WizPoint> getConnectedComponentCentroids(const cv::Mat& mask) {
    std::vector<WizPoint> centroids;
    cv::Mat tempMask = mask.clone();

    for (int y = 0; y < tempMask.rows; ++y) {
        for (int x = 0; x < tempMask.cols; ++x) {
            if (tempMask.at<uchar>(y, x) == 255) {
                std::vector<WizPoint> component;
                markVisited(tempMask, x, y, component);

                if (component.size() > 50) {
                    WizPoint centroid = getCentroid(component);
                    centroids.push_back(centroid);
                }
            }
        }
    }

    return centroids;
} // End 质心算法

// 取最亮点
WizPoint getGreatestPoint(std::vector<WizPoint> component,
                          const cv::Mat& mat) {

    double maxVal = std::numeric_limits<double>::lowest();
    if (component.empty()) {
        return {-1, -1};
    }
    WizPoint ret = component[0];
    for (const WizPoint& pixel : component) {
        // 如果在边缘8个像素距离内，跳过
        if (pixel.x < 8 || pixel.x > mat.cols - 8 || pixel.y < 8 || pixel.y > mat.rows - 8) {
            continue;
        }

        double val = mat.at<double>(pixel.y, pixel.x);
        if (val > maxVal) {
            maxVal = val;
            ret = pixel;
        }
    }
    return ret;
}

std::vector<WizPoint> getConnectedComponentMaxValue(const cv::Mat& mask, const cv::Mat& mat) {
    std::vector<WizPoint> centroids;
    cv::Mat tempMask = mask.clone();

    for (int y = 0; y < tempMask.rows; ++y) {
        for (int x = 0; x < tempMask.cols; ++x) {
            if (tempMask.at<uchar>(y, x) == 255) {
                std::vector<WizPoint> component;
                markVisited(tempMask, x, y, component);

                if (component.size() > 20) {
                    WizPoint centroid = getGreatestPoint(component, mat);
                    centroids.push_back(centroid);
                    if (centroids.size() > 61) {  // 限制最大数量
                        return centroids;
                    }
                }
            }
        }
    }

    return centroids;
} // End 最值点取点




FretThaSolver::FretThaSolver()
{
    qRegisterMetaType<QMap<TableHeader, QString>>("QMap<TableHeader, QString>");
}

/**
 * @brief FretThaSolver::clearGrayData 清除灰度数据
 */
void FretThaSolver::clearGrayData()
{
    for (auto& it : vecGrayData) {
        it.clear();
    }
}

/**
 * @brief FretThaSolver::clearFretData 清除 FRET 数据
 */
void FretThaSolver::clearFretData()
{
    for (auto& it : vecFretData) {
        it.clear();
    }
}

/**
 * @brief FretThaSolver::clearFretDataBin 清除 FRET 分箱数据
 */
void FretThaSolver::clearFretDataBin()
{
    for (auto& it : vecFretDataBin) {
        it.clear();
    }
}

/**
 * @brief FretThaSolver::outputData 保存后台数据到CSV文件
 * @param savePath 目标文件夹的路径
 * @return 保存文件的完整路径
 */
QString FretThaSolver::outputData(QString savePath)
{
    // 开始
    auto data = getDfreeAndAfreeOurs();
    std::vector<double> Dfree = data.first;
    std::vector<double> Afree = data.second;

    qDebug() << "[Write Data]:\tStart";
    // 检查文件，若有同名文件会覆盖
    QFile file(savePath + "/FretThaData.csv");
    if (file.exists()) file.remove();
    file.open(QIODevice::Append);
    QTextStream out(&file);
    out << "Iaa,Ida,Idd,Rc,Ed,1/Rc,Ea,Aest,Dest,Afree,Dfree\n";
    uint i = 0;
    for (; i < vecGrayData[CHANNEL_NAME_AA].size(); ++ i){
        out << vecGrayData[CHANNEL_NAME_AA][i] << ","
            << vecGrayData[CHANNEL_NAME_DA][i] << ","
            << vecGrayData[CHANNEL_NAME_DD][i] << ","
            << vecFretData[Rad][i] << ","
            << vecFretData[Ed][i] << ","
            << vecFretData[Rda][i] << ","
            << vecFretData[Ea][i] << ","
            << vecFretData[Aest][i] << ","
            << vecFretData[Dest][i] << ","
            << Afree[i] << ","
            << Dfree[i]
            <<"\n";
    }
    file.close();

    qDebug() << "[Writie Data]:\t" << i << "pieces of data are saved";
    return savePath + "/FretThaData.csv";
}

QString FretThaSolver::outputResults(QString savePath)
{
    // 开始
    qDebug() << "[Write Result]:\tStart";
    QString res = savePath + "/FretThaResults.csv";
    // 检查文件，若有同名文件会覆盖
    QFile file(res);
    if (file.exists()) file.remove();
    file.open(QIODevice::Append);
    QTextStream out(&file);
    out << "Methods,Eamax,Edmax,Kdeff,nD/nA\n";

    out << "Ours,"
        << resultsThaOurs[THA_RESULT_EAMAX] << ","
        << resultsThaOurs[THA_RESULT_EDMAX] << ","
        << resultsThaOurs[THA_RESULT_KDEFF] << ","
        << resultsThaOurs[THA_RESULT_STOIC] << "\n";

    out << "Ben,"
        << resultsThaBen[THA_RESULT_EAMAX] << ","
        << resultsThaBen[THA_RESULT_EDMAX] << ","
        << resultsThaBen[THA_RESULT_KDEFF] << ","
        << resultsThaBen[THA_RESULT_STOIC] << "\n";

    out << "Du-EdRc,"
        << resultsEdRad[THA_RESULT_EAMAX] << ","
        << resultsEdRad[THA_RESULT_EDMAX] << ","
        << resultsEdRad[THA_RESULT_KDEFF] << ","
        << resultsEdRad[THA_RESULT_STOIC] << "\n";

    out << "Du-EA-1/RC,"
        << resultsEaRda[THA_RESULT_EAMAX] << ","
        << resultsEaRda[THA_RESULT_EDMAX] << ","
        << resultsEaRda[THA_RESULT_KDEFF] << ","
        << resultsEaRda[THA_RESULT_STOIC] << "\n";

    file.close();
    return res;
}

/**
 * @brief FretThaSolver::expandGrayData 向数组中记录一条灰度数据
 * @param Iaa AA 通道灰度值
 * @param Ida DA 通道灰度值
 * @param Idd DD 通道灰度值
 */
void FretThaSolver::expandGrayData(double Iaa, double Ida, double Idd)
{
    vecGrayData[CHANNEL_NAME_AA].push_back(Iaa);
    vecGrayData[CHANNEL_NAME_DA].push_back(Ida);
    vecGrayData[CHANNEL_NAME_DD].push_back(Idd);
}

/**
 * @brief FretThaSolver::expandFretData 向 FRET 数组中记录一条数据
 * @param valEd
 * @param valEa
 * @param valRad
 * @param valRda
 * @param valAest
 * @param valDest
 */
void FretThaSolver::expandFretData(double valEd, double valEa, double valRad, double valRda, double valAest, double valDest)
{
    // 进行数据范围筛选，没有意义的数据会被跳过
    if (valEd <= 0 || valEa <= 0 || valRad <= 0 || valRda <= 0 || valAest <= 0 || valDest <= 0) {
        return ;
    }
    vecFretData[Ed].push_back(valEd);
    vecFretData[Ea].push_back(valEa);
    vecFretData[Rad].push_back(valRad);
    vecFretData[Rda].push_back(valRda);
    vecFretData[Aest].push_back(valAest);
    vecFretData[Dest].push_back(valDest);
}
/**
 * @brief FretThaSolver::processGrayToFretData 完成灰度值记录后，进行 FRET 计算，存入 FRET 数据中
 */
void FretThaSolver::processGrayToFretData()
{
    for (uint i = 0; i < vecGrayData[0].size(); ++ i) {
        double Iaa = vecGrayData[CHANNEL_NAME_AA][i];
        double Ida = vecGrayData[CHANNEL_NAME_DA][i];
        double Idd = vecGrayData[CHANNEL_NAME_DD][i];
        // 使用calculator进行计算
        calculator->loadCorrectedGray(Iaa, Ida, Idd);
        calculator->calcData();
        if (!calculator->isDataCalculated()) continue;
        // 获取数据
        double valEa = calculator->getResult(Ea);
        double valEd = calculator->getResult(Ed);
        double valRad = calculator->getResult(Rad);
        double valRda = calculator->getResult(Rda);
        double valAest = calculator->getResult(Aest);
        double valDest = calculator->getResult(Dest);
        // 添加数据
        expandFretData(valEd, valEa, valRad, valRda, valAest, valDest);
    }
    qDebug() << "[Final Data Size]:\t" << vecFretData[0].size();
}
/**
 * @brief FretThaSolver::setBatchPath 设置样本文件夹的路径
 * @param batchFolderPath
 */
void FretThaSolver::setBatchPath(QString batchFolderPath)
{
    this->batchFolderPath = batchFolderPath;
}
/**
 * @brief FretThaSolver::setRatio 设置某项系数的值
 * @param ratioName 系数名
 * @param value 系数的值
 */
void FretThaSolver::setRatio(RatioName ratioName, double value)
{
    imageProcessor->setRatio(ratioName, value);
    calculator->setRatio(ratioName, value);
}
/**
 * @brief FretThaSolver::setRatios 批量设置系数
 * @param a
 * @param b
 * @param c
 * @param d
 * @param G
 * @param K
 * @param Y
 */
void FretThaSolver::setRatios(double a, double b, double c, double d, double G, double K, double Y)
{
    imageProcessor->setRatios(a, b, c, d, G, K, Y);
    calculator->setRatios(a, b, c, d, G, K, Y);
}
/**
 * @brief FretThaSolver::setExposureTime 设置某一通道的曝光时间
 * @param channelName 通道名
 * @param value 曝光时间数值
 */
void FretThaSolver::setExposureTime(ChannelName channelName, double value)
{
    imageProcessor->setExposureTime(channelName, value);
    calculator->setExposureTime(channelName, value);
}
/**
 * @brief FretThaSolver::setExposureTimes 批量设置曝光时间
 * @param expsTimeAA
 * @param expsTimeDA
 * @param expsTimeDD
 */
void FretThaSolver::setExposureTimes(double expsTimeAA, double expsTimeDA, double expsTimeDD)
{
    imageProcessor->setExposureTimes(expsTimeAA, expsTimeDA, expsTimeDD);
    calculator->setExposureTimes(expsTimeAA, expsTimeDA, expsTimeDD);
}
void FretThaSolver::setThreshRatio(ChannelName channelName, double value)
{
    threshRatio[channelName] = value;
}
void FretThaSolver::setThreshRatios(double threshAA, double threshDA, double threshDD)
{
    setThreshRatio(CHANNEL_NAME_AA, threshAA);
    setThreshRatio(CHANNEL_NAME_DA, threshDA);
    setThreshRatio(CHANNEL_NAME_DD, threshDD);
}
/**
 * @brief FretThaSolver::loadSourceData 加载数据
 * @param fodoubleerPath
 */
void FretThaSolver::loadSourceData(QString folderPath)
{
    imageProcessor->loadSourceData(folderPath);
}

void FretThaSolver::loadSourceData(QString pathAA, QString pathDA, QString pathDD)
{
    imageProcessor->loadSourceData(pathAA, pathDA, pathDD);
}

void FretThaSolver::processImageToGrayData()
{
    using namespace cv;

    // 赋值图片数据
    Mat matSrc[3];
    for (int i = 0; i < 3; ++ i) {
        matSrc[i] = imageProcessor->getMatrixCopy(ChannelName(i));
    }

    // 对数据进行均值滤波，并且计算得到背景灰度值
    Mat matFiltered[3];
    double grayBackground[3];
    for (int i = 0; i < 3; ++ i) {
        matFiltered[i] = wiz::meanFilter16U(matSrc[i], 5);
        grayBackground[i] = wiz::calcBackgroundGray(matFiltered[i]);
    }

    // 根据SBR阈值进行滤波
    double threshRatio[3] = {3.0, 3.0, 3.0};
    Mat maskSingleChannel[3];
    Mat maskSbr, maskSbr8U;    // 0-1
    for (int i = 0; i < 3; ++ i) {
        threshold(matFiltered[i], maskSingleChannel[i], grayBackground[i] * threshRatio[i], 255, THRESH_BINARY);
        maskSingleChannel[i].convertTo(maskSingleChannel[i], CV_8U);
    }
    bitwise_and(maskSingleChannel[CHANNEL_NAME_AA], maskSingleChannel[CHANNEL_NAME_DA], maskSbr8U);
    bitwise_and(maskSbr8U, maskSingleChannel[CHANNEL_NAME_DD], maskSbr8U);
    maskSbr = wiz::normalizeByMinMax(maskSbr8U);
    //if (RUNMODE == RUNMODE_EXPERIMENT) cv::imwrite((viewFolderPath + "/maskSBRatio.tif").toStdString(), maskSbr);

    // 按照我们设计的算法生成评价结果
    Mat maskStd;
    Mat scoreStd(matSrc[CHANNEL_NAME_AA].size(), CV_64FC1, Scalar(0.0));
    for (int i = 0; i < 3; ++ i) {

        Mat matStd = wiz::localStandardDeviation(matFiltered[i], 5);
        Mat matStd8U = wiz::normalizeByMinMax8U(matStd);
        //if (RUNMODE == RUNMODE_EXPERIMENT) cv::imwrite((viewFolderPath + "/LocalStandardDeviation" + QString::number(i) + ".tif").toStdString(), matStd8U);

        Mat matStdMinima = wiz::findLocalMinima(matStd8U);
        Mat matStdMinima8U = wiz::normalizeByMinMax8U(matStdMinima);
        //if (RUNMODE == RUNMODE_EXPERIMENT) cv::imwrite((viewFolderPath + "/StdLocalMinima" + QString::number(i) + ".tif").toStdString(), matStdMinima8U);

        Mat matStdMinimaLap = wiz::findLocalMinimaLaplacian(matStd8U);
        Mat matStdMinimaLap8U = wiz::normalizeByMinMax8U(matStdMinimaLap);
        //if (RUNMODE == RUNMODE_EXPERIMENT) cv::imwrite((viewFolderPath + "/StdLocalMinimaLaplacian" + QString::number(i) + ".tif").toStdString(), matStdMinimaLap8U);

        scoreStd = scoreStd + wiz::normalizeByMinMax(matStdMinima);
    }

    Mat scoreStd16U = wiz::normalizeByMinMax16U(scoreStd);
    //if (RUNMODE == RUNMODE_EXPERIMENT) cv::imwrite((viewFolderPath + "/StdMinimaFinal.tif").toStdString(), scoreStd16U);
    cv::threshold(scoreStd16U, maskStd, 0, 255, cv::THRESH_BINARY | cv::THRESH_OTSU);
    maskStd.convertTo(maskStd, CV_8U);
    //if (RUNMODE == RUNMODE_EXPERIMENT) cv::imwrite((viewFolderPath + "/StdMinimaMaskOtsu.tif").toStdString(), maskStd);

    // 合成最后的结果
    Mat maskUlt, scoreMat;
    scoreMat = maskSbr.mul(scoreStd);
    bitwise_and(maskSbr8U, maskStd, maskUlt);
    //if (RUNMODE == RUNMODE_EXPERIMENT) cv::imwrite((viewFolderPath + "/FinalJudge.tif").toStdString(), maskUlt);

    // 计算质心
    std::vector<WizPoint> centroids = getConnectedComponentCentroids(maskUlt);

    for (const auto& centroid : centroids) {
        int x = centroid.x;
        int y = centroid.y;
        double valGrayAA = matFiltered[CHANNEL_NAME_AA].at<ushort>(y, x) - grayBackground[CHANNEL_NAME_AA];
        double valGrayDA = matFiltered[CHANNEL_NAME_DA].at<ushort>(y, x) - grayBackground[CHANNEL_NAME_DA];
        double valGrayDD = matFiltered[CHANNEL_NAME_DD].at<ushort>(y, x) - grayBackground[CHANNEL_NAME_DD];

        expandGrayData(valGrayAA, valGrayDA, valGrayDD);
    }
    qDebug() << "[Expand Gray Data]:\t" << centroids.size();
}

void FretThaSolver::processOneView(QString viewFolderPath)
{
    qDebug() << "[Processing View]:\t" << viewFolderPath;
    this->viewFolderPath = viewFolderPath;
    // 使用图像处理类读入数据,检查数据完整性
    imageProcessor->loadSourceData(viewFolderPath);
    if (imageProcessor->isDataLoaded()) {
        // 处理数据
        processImageToGrayData();   // 从图片中取得灰度值
    }
    else {
        qDebug() << "[Expand Gray Data]:\tFailed";
    }
}

void FretThaSolver::processOneBatch(QString batchFolderPath)
{
    qDebug() << "[Processing Batch]:\t" << batchFolderPath;
    this->batchFolderPath = batchFolderPath;

    // 清空数据存储
    clearGrayData();
    clearFretData();

    // 检查目录是否存在
    QDir batchDir(batchFolderPath);
    if (!batchDir.exists()) return;

    // 对目录下的每个子文件夹进行处理
    QFileInfoList folderlist = batchDir.entryInfoList(QDir::Dirs | QDir::NoDotAndDotDot);
    for (int i = 0; i < folderlist.size(); ++ i) {
        // 处理一个视野的数据，扩展进灰度数据中
        processOneView(folderlist.at(i).absoluteFilePath());
        // 发送进度条
        int progressValue = (i + 1) * 100 / folderlist.size();
        emit progressChanged(progressValue);
    }
}

void FretThaSolver::performTwoHybridOurs() {
    fitLangmiurDlib();
}

void FretThaSolver::performTwoHybridBen(double min, double max, double interval, bool bin_flag)
{
    lab_fitLangmiurMatlabBinCompare(min, max, interval);
    if (bin_flag) fitLangmiurMatlabBin(min, max, interval);
    else fitLangmiurMatlab();
}

void FretThaSolver::fitLangmiurDlib() {
    column_vector result = TwoHybridSolver(vecFretData[Aest], vecFretData[Dest], vecFretData[Ea], vecFretData[Ed], false);
    resultsThaOurs[THA_RESULT_KDEFF] = result(0);
    resultsThaOurs[THA_RESULT_STOIC] = result(1);
    resultsThaOurs[THA_RESULT_EAMAX] = result(2);
    resultsThaOurs[THA_RESULT_EDMAX] = result(3);

    qDebug() << "[Fit Langmiur With Dlib]:";
    qDebug() << "Ed,max:" << resultsThaOurs[THA_RESULT_EDMAX]
             << "\tEa,max:" << resultsThaOurs[THA_RESULT_EAMAX]
             << "\tKd,EFF:" << resultsThaOurs[THA_RESULT_KDEFF]
             << "\tNd/Na:" << resultsThaOurs[THA_RESULT_STOIC];
}

void FretThaSolver::fitLangmiurMatlab() {
    // 进行规划求解
    coder::array<double, 2U> arrayAest, arrayDest, arrayEacorr, arrayEdcorr;
    uint n = vecFretData[Ed].size();

    // 检查大小匹配
    if (n == vecFretData[Dest].size() &&
        n == vecFretData[Ea].size() &&
        n == vecFretData[Aest].size()) {
        arrayAest.set_size(1, n);
        arrayDest.set_size(1, n);
        arrayEacorr.set_size(1, n);
        arrayEdcorr.set_size(1, n);
    }
    else {
        qDebug() << "[Fit With Matlab]:\t Data Size Not Matched";
        return;
    }

    // 读取数据
    for (uint i = 0; i < n; ++ i) {
        arrayAest[i] = vecFretData[Aest][i];
        arrayDest[i] = vecFretData[Dest][i];
        arrayEacorr[i] = vecFretData[Ea][i];
        arrayEdcorr[i] = vecFretData[Ed][i];
    }

    // 进行计算
    double loss;
    TwoHybridSolver(arrayAest, arrayDest, arrayEacorr, arrayEdcorr, resultsThaBen, &loss);
    TwoHybridSolver_terminate();

    qDebug() << "[Fit With Matlab]:";
    qDebug() << "Ed,max:" << resultsThaBen[THA_RESULT_EDMAX]
             << "\tEa,max:" << resultsThaBen[THA_RESULT_EAMAX]
             << "\tKd,EFF:" << resultsThaBen[THA_RESULT_KDEFF]
             << "\tNd/Na:" << resultsThaBen[THA_RESULT_STOIC]
             << "\tloss:" << loss;
}

void FretThaSolver::fitLangmiurMatlabBin(double min, double max, double interval)  // matlab求解规划()
{
    // 执行数据封箱
    binDataByRda(min, max, interval);

    // 进行规划求解
    coder::array<double, 2U> arrayAest, arrayDest, arrayEacorr, arrayEdcorr;
    uint n = vecFretDataBin[Ed].size();

    // 检查大小匹配
    if (n == vecFretDataBin[Dest].size() &&
        n == vecFretDataBin[Ea].size() &&
        n == vecFretDataBin[Aest].size()) {
        arrayAest.set_size(1, n);
        arrayDest.set_size(1, n);
        arrayEacorr.set_size(1, n);
        arrayEdcorr.set_size(1, n);
    } else {
        qDebug() << "[Matlab Fmincon Fit With Bin]:\t Data Size Not Matched";
        return;
    }

    // 读取数据
    for (uint i = 0; i < n; ++ i) {
        arrayAest[i] = vecFretDataBin[Aest][i];
        arrayDest[i] = vecFretDataBin[Dest][i];
        arrayEacorr[i] = vecFretDataBin[Ea][i];
        arrayEdcorr[i] = vecFretDataBin[Ed][i];
    }

    // 进行计算
    double loss;
    TwoHybridSolver(arrayAest, arrayDest, arrayEacorr, arrayEdcorr, resultsThaBen, &loss);
    TwoHybridSolver_terminate();

    qDebug() << "[Matlab Fmincon Fit With Bin]:";
    qDebug() << "Ed,max:" << resultsThaBen[THA_RESULT_EDMAX]
             << "\tEa,max:" << resultsThaBen[THA_RESULT_EAMAX]
             << "\tKd,EFF:" << resultsThaBen[THA_RESULT_KDEFF]
             << "\tNd/Na:" << resultsThaBen[THA_RESULT_STOIC]
             << "\tloss:" << loss;
}

void FretThaSolver::lab_fitLangmiurMatlabBinCompare(double min, double max, double interval)  // matlab求解规划()
{
    qDebug();
    qDebug() << "********************************LAB********************************";

    // Rda
    binDataByRda(min, max, interval);
    coder::array<double, 2U> arrayAest, arrayDest, arrayEacorr, arrayEdcorr;
    uint n = vecFretDataBin[Ed].size();
    if (n == vecFretDataBin[Dest].size() &&
        n == vecFretDataBin[Ea].size() &&
        n == vecFretDataBin[Aest].size()) {
        arrayAest.set_size(1, n);
        arrayDest.set_size(1, n);
        arrayEacorr.set_size(1, n);
        arrayEdcorr.set_size(1, n);
    } else {
        qDebug() << "[Matlab Fmincon Fit With Bin]:\t Data Size Not Matched";
        return;
    }
    for (uint i = 0; i < n; ++ i) {
        arrayAest[i] = vecFretDataBin[Aest][i];
        arrayDest[i] = vecFretDataBin[Dest][i];
        arrayEacorr[i] = vecFretDataBin[Ea][i];
        arrayEdcorr[i] = vecFretDataBin[Ed][i];
    }
    double loss;
    TwoHybridSolver(arrayAest, arrayDest, arrayEacorr, arrayEdcorr, resultsThaBen, &loss);
    TwoHybridSolver_terminate();
    qDebug() << "[Lab Bin Rda]:";
    qDebug() << "Ed,max:" << resultsThaBen[THA_RESULT_EDMAX]
             << "\tEa,max:" << resultsThaBen[THA_RESULT_EAMAX]
             << "\tKd,EFF:" << resultsThaBen[THA_RESULT_KDEFF]
             << "\tNd/Na:" << resultsThaBen[THA_RESULT_STOIC]
             << "\tloss:" << loss;

    auto origin = getDfreeAndAfreeBenBin();
    std::vector<double> Dfree = origin.first;
    std::vector<double> Afree = origin.second;
    // 将vecFretDataBin中Ed和Afree、Ea和Dfree的数据成对保存
    std::pair<std::vector<double>, std::vector<double>> EdAfreeOrigin;
    std::pair<std::vector<double>, std::vector<double>> EaDfreeOrigin;
    for (uint i = 0; i < n; ++ i) {
        EdAfreeOrigin.first.push_back(vecFretDataBin[Ed][i]);
        EdAfreeOrigin.second.push_back(Afree[i]);
        EaDfreeOrigin.first.push_back(vecFretDataBin[Ea][i]);
        EaDfreeOrigin.second.push_back(Dfree[i]);
    }


    // Ours
    binDataOurs(-2.5, 2.5, 0.05);
    arrayAest.clear();
    arrayDest.clear();
    arrayEacorr.clear();
    arrayEdcorr.clear();
    n = vecFretDataBin[Ed].size();
    if (n == vecFretDataBin[Dest].size() &&
        n == vecFretDataBin[Ea].size() &&
        n == vecFretDataBin[Aest].size()) {
        arrayAest.set_size(1, n);
        arrayDest.set_size(1, n);
        arrayEacorr.set_size(1, n);
        arrayEdcorr.set_size(1, n);
    } else {
        qDebug() << "[Matlab Fmincon Fit With Bin]:\t Data Size Not Matched";
        return;
    }
    for (uint i = 0; i < n; ++ i) {
        arrayAest[i] = vecFretDataBin[Aest][i];
        arrayDest[i] = vecFretDataBin[Dest][i];
        arrayEacorr[i] = vecFretDataBin[Ea][i];
        arrayEdcorr[i] = vecFretDataBin[Ed][i];
    }
    TwoHybridSolver(arrayAest, arrayDest, arrayEacorr, arrayEdcorr, resultsThaBen, &loss);
    TwoHybridSolver_terminate();
    qDebug() << "[Lab Bin Ours]:";
    qDebug() << "Ed,max:" << resultsThaBen[THA_RESULT_EDMAX]
             << "\tEa,max:" << resultsThaBen[THA_RESULT_EAMAX]
             << "\tKd,EFF:" << resultsThaBen[THA_RESULT_KDEFF]
             << "\tNd/Na:" << resultsThaBen[THA_RESULT_STOIC]
             << "\tloss:" << loss;

    auto ours = getDfreeAndAfreeBenBin();
    Dfree = ours.first;
    Afree = ours.second;
    // 将vecFretDataBin中Ed和Afree、Ea和Dfree的数据成对保存
    std::pair<std::vector<double>, std::vector<double>> EdAfreeOurs;
    std::pair<std::vector<double>, std::vector<double>> EaDfreeOurs;
    for (uint i = 0; i < n; ++ i) {
        EdAfreeOurs.first.push_back(vecFretDataBin[Ed][i]);
        EdAfreeOurs.second.push_back(Afree[i]);
        EaDfreeOurs.first.push_back(vecFretDataBin[Ea][i]);
        EaDfreeOurs.second.push_back(Dfree[i]);
    }
    qDebug() << "********************************END********************************";
    qDebug();
}

void FretThaSolver::performTwoHybridDu(double minSlope , double maxSlope, double minAppro, double maxAppro)
{
    // 经验值范围应是0-0.5，3-5
    resultsEdRad[THA_RESULT_EAMAX] = calcSlope(minSlope, maxSlope, vecFretData[Rad], vecFretData[Ed]);
    resultsEdRad[THA_RESULT_EDMAX] = calcApproach(minAppro, maxAppro, vecFretData[Rad], vecFretData[Ed]);
    resultsEdRad[THA_RESULT_STOIC] = resultsEdRad[THA_RESULT_EAMAX] / resultsEdRad[THA_RESULT_EDMAX];

    resultsEaRda[THA_RESULT_EDMAX] = calcSlope(minSlope, maxSlope, vecFretData[Rda], vecFretData[Ea]);
    resultsEaRda[THA_RESULT_EAMAX] = calcApproach(minAppro, maxAppro, vecFretData[Rda], vecFretData[Ea]);
    resultsEaRda[THA_RESULT_STOIC] = resultsEaRda[THA_RESULT_EAMAX] / resultsEaRda[THA_RESULT_EDMAX];

    // Debug
    qDebug() << "[Linear Fit]:";
    qDebug() << "[Donor View Results]:";
    qDebug() << "Ed,max:" << resultsEdRad[THA_RESULT_EDMAX]
             << "\tEa,max:" << resultsEdRad[THA_RESULT_EAMAX]
             << "\tNd/Na:" << resultsEdRad[THA_RESULT_STOIC];
    qDebug() << "[Accept View Results]:";
    qDebug() << "Ed,max:" << resultsEaRda[THA_RESULT_EDMAX]
             << "\tEa,max:" << resultsEaRda[THA_RESULT_EAMAX]
             << "\tNd/Na:" << resultsEaRda[THA_RESULT_STOIC];
}

std::pair<std::vector<double>, std::vector<double>> sortArrays(const std::vector<double>& a, const std::vector<double>& b) {
    // 创建一个索引向量，初识状态下索引i的值就是i
    std::vector<size_t> indices(a.size());
    std::iota(indices.begin(), indices.end(), 0);

    // 根据a的值对索引向量进行排序
    std::sort(indices.begin(), indices.end(), [&a](size_t i1, size_t i2) { return a[i1] < a[i2]; });

    // 根据索引向量对a和b进行排序
    std::vector<double> aSorted(a.size());
    std::vector<double> bSorted(b.size());
    for (size_t i = 0; i < indices.size(); ++i) {
        aSorted[i] = a[indices[i]];
        bSorted[i] = b[indices[i]];
    }

    return {aSorted, bSorted};
}

int findTrueValueDynamicEfficient(const std::vector<double>& arr) {
    if (arr.size() < 2) return arr.empty() ? -1 : arr[0]; // 如果数组为空或只有一个元素，直接返回

    double true_value = arr[0], true_index = 0;
    std::vector<double> diffs; // 存储相邻元素之间的差异

    for (size_t i = 1; i < arr.size(); ++i) {
        // 更新diffs数组
        diffs.push_back(std::abs(arr[i] - arr[i - 1]));

        // 计算当前阈值：基于diffs的平均值
        double threshold = std::accumulate(diffs.begin(), diffs.end(), 0.0) / diffs.size() * 0.5;

        // 如果是最后一个元素，就不需要比较差异了
        if (i == arr.size() - 1) break;

        // 使用动态阈值判断稳定性
        if (std::abs(arr[i] - arr[i - 1]) < threshold) {
            true_value = std::max(true_value, arr[i]);
            true_index = i;
        }
    }

    qDebug() << true_index << "/" << arr.size();

    return true_index; // 如果没有找到稳定点，返回遍历过程中的最大值
}

void saveToCSV(const std::vector<double>& vec, const std::string& filename) {
    std::ofstream file(filename); // 创建文件流，准备写入文件

    if (file.is_open()) {
        // 遍历vector中的每个元素
        for (size_t i = 0; i < vec.size(); ++i) {
            // 写入元素到文件，用逗号分隔
            file << vec[i];
            if (i < vec.size() - 1) {
                file << "\n"; // 除了最后一个元素外，其他元素后面添加逗号
            }
        }
        file.close(); // 关闭文件流
        std::cout << "数据已保存到CSV文件。" << std::endl;
    } else {
        std::cerr << "无法打开文件。" << std::endl;
    }
}

double FretThaSolver::FindBestFittingRange(CalcResult r_type) {
    using namespace std;
    CalcResult e_type;
    if (r_type == Rad)
    {
        e_type = Ed;
    } else {
        e_type = Ea;
    }
    int size = vecFretData[r_type].size();
    vector<double> slopes(size, 0.0);
    auto arrays = sortArrays(vecFretData[r_type], vecFretData[e_type]);
    vector<double> sorted_r = arrays.first, sorted_e = arrays.second;
    for (int i = 0; i < size; ++ i) {
        slopes[i] = calcSlope(0, sorted_r[i], sorted_r, sorted_e);
    }

    // 保存到CSV文件，命名自动加上时间
    QString filename = QString("D:\\Files\\FRET\\Data\\slopes_%1.csv").arg(QDateTime::currentDateTime().toString("yyyyMMddhhmmss"));
    qDebug() << filename;
    saveToCSV(slopes, filename.toStdString());

    return sorted_r[findTrueValueDynamicEfficient(slopes)];
}

double calcSlopeIndex(std::vector<double> x, std::vector<double> y, int startIndex, int endIndex)
{
    std::vector<std::pair<double, double>> points;

    for (int i = startIndex; i <= endIndex; ++ i) {
        points.push_back({x[i], y[i]});
    }
    return leastSquare(points);
}

void FretThaSolver::performTwoHybridDuAutoRange()
{
    if (RUNMODE != RUNMODE_DEBUG) {
        return ;
    }

    std::vector<double> log_rad(vecFretData[Rad]);
    for (uint i = 0; i < log_rad.size(); ++i) {
        log_rad[i] = std::log(log_rad[i]);
    }
    double mean_log_ad_ratio = std::accumulate(log_rad.begin(), log_rad.end(), 0.0) / log_rad.size();
    double mean_ad_ratio = std::exp(mean_log_ad_ratio);
    double rad_bound = mean_ad_ratio * 0.5;
    double rda_bound = 0.5 * 1.0 / mean_ad_ratio;
    double eamax = calcSlope(0, rad_bound, vecFretData[Rad], vecFretData[Ed]);
    double edmax = calcSlope(0, rda_bound, vecFretData[Rda], vecFretData[Ea]);
    qDebug() << "[Rc Bound & 1/Rc Bound]:\t\t" << rad_bound << rda_bound;
    qDebug() << "[Auto EdRc Result]:\t\t" << eamax << edmax << eamax / edmax;

    m_result_auto_edmax = edmax;
    m_result_auto_eamax = eamax;
    m_result_auto_stoic = eamax / edmax;
    m_result_mean_ratio = mean_ad_ratio;
}

double FretThaSolver::maxData(CalcResult dataName)
{
    auto maxElem = std::max_element(vecFretData[dataName].begin(), vecFretData[dataName].end());
    if (maxElem != vecFretData[dataName].end()) {
        return *maxElem;
    }
    else return 1;
}

double FretThaSolver::maxBinData(CalcResult dataName)
{
    auto maxElem = std::max_element(vecFretDataBin[dataName].begin(), vecFretDataBin[dataName].end());
    if (maxElem != vecFretDataBin[dataName].end()) {
        return *maxElem;
    }
    else return 1;
}

std::pair<std::vector<double>, std::vector<double>> FretThaSolver::getDfreeAndAfreeOurs() {
    std::vector<double> Afree, Dfree;
    for (uint i = 0; i < vecFretData[Aest].size(); ++ i) {
        double valDfree = ( ((vecFretData[Dest][i] - resultsThaOurs[THA_RESULT_KDEFF] - vecFretData[Aest][i] * resultsThaOurs[THA_RESULT_STOIC]))
                           + sqrt(pow((vecFretData[Dest][i] - resultsThaOurs[THA_RESULT_KDEFF] - vecFretData[Aest][i] * resultsThaOurs[THA_RESULT_STOIC]), 2) +
                                  4 * resultsThaOurs[THA_RESULT_KDEFF] * vecFretData[Dest][i])
                           ) / 2;
        double valAfree = vecFretData[Aest][i] - (vecFretData[Dest][i] - valDfree) / resultsThaOurs[THA_RESULT_STOIC];
        Afree.push_back(valAfree);
        Dfree.push_back(valDfree);
    }

    return {Dfree, Afree};
}


std::pair<std::vector<double>, std::vector<double>> FretThaSolver::getDfreeAndAfreeBen() {
    std::vector<double> Afree, Dfree;
    for (uint i = 0; i < vecFretData[Aest].size(); ++ i) {
        double valDfree = ( ((vecFretData[Dest][i] - resultsThaBen[THA_RESULT_KDEFF] - vecFretData[Aest][i] * resultsThaBen[THA_RESULT_STOIC]))
                           + sqrt(pow((vecFretData[Dest][i] - resultsThaBen[THA_RESULT_KDEFF] - vecFretData[Aest][i] * resultsThaBen[THA_RESULT_STOIC]), 2) +
                                  4 * resultsThaBen[THA_RESULT_KDEFF] * vecFretData[Dest][i])
                           ) / 2;
        double valAfree = vecFretData[Aest][i] - (vecFretData[Dest][i] - valDfree) / resultsThaBen[THA_RESULT_STOIC];
        Afree.push_back(valAfree);
        Dfree.push_back(valDfree);
    }

    return {Dfree, Afree};
}

std::pair<std::vector<double>, std::vector<double>> FretThaSolver::getDfreeAndAfreeBenBin()
{
    std::vector<double> Afree, Dfree;
    for (uint i = 0; i < vecFretDataBin[Aest].size(); ++ i) {
        double valDfree = ( ((vecFretDataBin[Dest][i] - resultsThaBen[THA_RESULT_KDEFF] - vecFretDataBin[Aest][i] * resultsThaBen[THA_RESULT_STOIC]))
                           + sqrt(pow((vecFretDataBin[Dest][i] - resultsThaBen[THA_RESULT_KDEFF] - vecFretDataBin[Aest][i] * resultsThaBen[THA_RESULT_STOIC]), 2) +
                                  4 * resultsThaBen[THA_RESULT_KDEFF] * vecFretDataBin[Dest][i]) ) / 2;
        double valAfree = vecFretDataBin[Aest][i] - (vecFretDataBin[Dest][i] - valDfree) / resultsThaBen[THA_RESULT_STOIC];
        Afree.push_back(valAfree);
        Dfree.push_back(valDfree);
    }

    return {Dfree, Afree};
}

void FretThaSolver::generateRoiFromBatch(QString batchFolderPath)
{
    qDebug() << "[Processing Batch]:\t" << batchFolderPath;
    this->batchFolderPath = batchFolderPath;

    // 检查目录是否存在
    QDir batchDir(batchFolderPath);
    if (!batchDir.exists()) return;

    // 对目录下的每个子文件夹进行处理
    QFileInfoList folderlist = batchDir.entryInfoList(QDir::Dirs | QDir::NoDotAndDotDot);
    for (int i = 0; i < folderlist.size(); ++ i) {
        // 删除目标文件夹下除AA.tif、DA.tif、DD.tif以外的所有文件
        QString viewFolderPath = folderlist.at(i).absoluteFilePath();
        QDir viewDir(viewFolderPath);
        QStringList filters;
        filters << "*.tif";
        viewDir.setNameFilters(filters);
        QFileInfoList filelist = viewDir.entryInfoList(QDir::Files);
        for (int j = 0; j < filelist.size(); ++ j) {
            QString fileName = filelist.at(j).fileName();
            if (fileName != "AA.tif" && fileName != "DA.tif" && fileName != "DD.tif") {
                QFile::remove(filelist.at(j).absoluteFilePath());
            }
        }
        // 处理一个视野的数据，扩展进灰度数据中
        generateRoiFromView(folderlist.at(i).absoluteFilePath(), folderlist.at(i).fileName());
        // 发送进度条
        int progressValue = (i + 1) * 100 / folderlist.size();
        qDebug() << "[Auto Generate Progress]" << progressValue << "%";
        emit progressChanged(progressValue);
    }
}

void FretThaSolver::loadRoiFromBatch(QString batchFolderPath)
{
    qDebug() << "[Processing Batch]:\t" << batchFolderPath;
    this->batchFolderPath = batchFolderPath;

    // 检查目录是否存在
    QDir batchDir(batchFolderPath);
    if (!batchDir.exists()) return;

    // 对目录下的每个子文件夹进行处理
    QFileInfoList folderlist = batchDir.entryInfoList(QDir::Dirs | QDir::NoDotAndDotDot);
    for (int i = 0; i < folderlist.size(); ++ i) {

        // 处理一个视野的数据，扩展进灰度数据中
        loadRoiFromView(folderlist.at(i).absoluteFilePath(), folderlist.at(i).fileName());
        // 发送进度条
        int progressValue = (i + 1) * 100 / folderlist.size();
        qDebug() << "[Auto Generate Progress]" << progressValue << "%";
        emit progressChanged(progressValue);
    }
}

void FretThaSolver::generateRoiFromView(QString viewFolderPath, QString viewName)
{
    qDebug() << "[Processing View]:\t" << viewFolderPath;
    this->viewFolderPath = viewFolderPath;
    // 使用图像处理类读入数据,检查数据完整性
    imageProcessor->loadSourceData(viewFolderPath);
    if (imageProcessor->isDataLoaded()) {
        imageProcessor->preProcessData();
        // imageProcessor->calcEFret();
        // 处理数据
        generateRoiFromImage(viewName);   // 从图片中取得灰度值
    }
    else {
        qDebug() << "[Expand Gray Data]:\tFailed";
    }
}
void FretThaSolver::loadRoiFromView(QString viewFolderPath, QString folderName)
{
    using namespace cv;
    qDebug() << "[Processing View]:\t" << viewFolderPath;
    this->viewFolderPath = viewFolderPath;
    imageProcessor->loadSourceData(viewFolderPath);
    if (imageProcessor->isDataLoaded()) {
        // 读数据
        std::vector<Mat> datas;
        Mat matSrc[3];
        double valBg[3];
        for (int i = 0; i < 3; ++ i) {
            matSrc[i] = imageProcessor->getMatrixCopy(ChannelName(i));
            datas.push_back(imageProcessor->getMatrixCopy(ChannelName(i)));
            valBg[i] = wiz::calcBackgroundGray(matSrc[i]);
        }
        QString maskPath = "D:/mask/" + folderName + ".png";

        Mat mask = imread(maskPath.toStdString(), IMREAD_GRAYSCALE);
        if (mask.empty()) {
            qDebug() << "mask empty";
            return;
        }
        for (int i = 0; i < 3; ++ i) {
            mask = wiz::getMorphologyClose(mask, 3); // 进行形态学操作去除小的噪点
            mask = wiz::getMorphologyOpen(mask, 3);
        }

        std::vector<WizPoint> centroids = getConnectedComponentMaxValue(mask,
                                                                        matSrc[CHANNEL_NAME_AA] / valBg[CHANNEL_NAME_AA] + matSrc[CHANNEL_NAME_DA] / valBg[CHANNEL_NAME_DA] + matSrc[CHANNEL_NAME_DD] / valBg[CHANNEL_NAME_DD]);

        int roiSize = 5;
        for (const WizPoint& centroid : centroids) {
            int x = centroid.x;
            int y = centroid.y;

            int rectx = x - roiSize / 2;
            int recty = y - roiSize / 2;
            int rectw = roiSize;
            int recth = roiSize;

            int n = matSrc[0].rows;
            if (!(0 <= rectx && 0 <= rectw && rectx + rectw <= n && 0 <= recty && 0 <= recth && recty + recth <= n)) continue;

            // 计算MatSrc在rect区域的平均灰度值
            double sbr = 0.0;
            cv::Rect roi(rectx, recty, rectw, recth);
            double valFg[3];
            for (int i = 0; i < 3; ++ i) {
                Mat matRoi = matSrc[i](roi);
                Scalar mean = cv::mean(matRoi);
                valFg[i] = mean.val[0];

                sbr += valFg[i] / valBg[i];
            }
            if (sbr < 9) continue;

            double valGrayAA = valFg[CHANNEL_NAME_AA] - valBg[CHANNEL_NAME_AA];
            double valGrayDA = valFg[CHANNEL_NAME_DA] - valBg[CHANNEL_NAME_DA];
            double valGrayDD = valFg[CHANNEL_NAME_DD] - valBg[CHANNEL_NAME_DD];

            // 使用calculator进行计算
            calculator->loadCorrectedGray(valGrayAA, valGrayDA, valGrayDD);
            calculator->calcData();
            if (!calculator->isDataCalculated()) continue;
            // 获取数据
            double valEa = calculator->getResult(Ea);
            double valEd = calculator->getResult(Ed);
            double valRad = calculator->getResult(Rad);
            double valRda = calculator->getResult(Rda);
            double valAest = calculator->getResult(Aest);
            double valDest = calculator->getResult(Dest);

            if (valEa > 0 && valEd > 0 && valRad > 0 && valRda > 0 && valAest > 0 && valDest > 0) {
                QMap<TableHeader, QString> map;
                map.insert(TABLE_HEADER_IAA, QString::number(valGrayAA));
                map.insert(TABLE_HEADER_IDA, QString::number(valGrayDA));
                map.insert(TABLE_HEADER_IDD, QString::number(valGrayDD));
                map.insert(TABLE_HEADER_EA, QString::number(valEa));
                map.insert(TABLE_HEADER_ED, QString::number(valEd));
                map.insert(TABLE_HEADER_RAD, QString::number(valRad));
                map.insert(TABLE_HEADER_RDA, QString::number(valRda));
                map.insert(TABLE_HEADER_AEST, QString::number(valAest));
                map.insert(TABLE_HEADER_DEST, QString::number(valDest));
                map.insert(TABLE_HEADER_RECTX, QString::number(rectx));
                map.insert(TABLE_HEADER_RECTY, QString::number(recty));
                map.insert(TABLE_HEADER_RECTW, QString::number(rectw));
                map.insert(TABLE_HEADER_RECTH, QString::number(recth));
                map.insert(TABLE_HEADER_VIEW, folderName);
                if (rectx >= 0 && recty >= 0) emit sendData(map);
            }
        }



        /* 提取ROI灰度值

        // 使用connectedComponents找到所有连通域
        cv::Mat labels;
        int num_labels = cv::connectedComponents(mask, labels, 8, CV_32S);

        // 创建一个向量来存储每个连通域的灰度平均值，对于每个data矩阵
        std::vector<std::vector<double>> avgs(datas.size(), std::vector<double>(num_labels, 0.0));
        std::vector<int> counts(num_labels, 0);

        // 遍历每个像素，计算每个连通域的灰度总和和像素数量，对于每个data矩阵
        for (int y = 0; y < labels.rows; ++y) {
            for (int x = 0; x < labels.cols; ++x) {
                int label = labels.at<int>(y, x);
                for (size_t i = 0; i < datas.size(); ++i) {
                    avgs[i][label] += datas[i].at<ushort>(y, x);
                }
                counts[label]++;
            }
        }

        // 计算每个连通域的灰度平均值，对于每个data矩阵
        for (size_t i = 0; i < datas.size(); ++i) {
            for (int j = 0; j < num_labels; ++j) {
                avgs[i][j] /= counts[j];
            }
        }

        for (int i = 0; i < num_labels; ++ i) {
            double valGrayAA = avgs[0][i];
            double valGrayDA = avgs[1][i];
            double valGrayDD = avgs[2][i];


        }*/


    }
    else {
        qDebug() << "[Expand Gray Data]:\tFailed";
    }
}

void FretThaSolver::Segmentation()
{
    using namespace cv;
    Mat  matSrc[3], matMean[3], matBlur[3];
    for (int i = 0; i < 3; ++ i) {
        matSrc[i] = imageProcessor->getMatrixCopy(ChannelName(i));
        Mat markers = wiz::watershedSegmentation(matSrc[i]);
        Mat colored = wiz::visualizeWatershed(markers);
        cv::imwrite(
            (viewFolderPath + "/Watershed" + QString::number(i)+ ".tif").toStdString(),
            colored
        );
    }
}

namespace biocv {
    using namespace std;
    using namespace cv;

    // 计算两个连通域之间的距离
    double calculateDistance(const Moments& m1, const Moments& m2) {
        double dx = m1.m10 / m1.m00 - m2.m10 / m2.m00;
        double dy = m1.m01 / m1.m00 - m2.m01 / m2.m00;
        return sqrt(dx*dx + dy*dy);
    }

    // 合并连通域
    void mergeConnectedComponents(Mat& labels, double threshold) {
        int num_labels = labels.cols;

        // 合并连通域
        for (int i = 1; i < num_labels; ++i) {
            for (int j = i + 1; j < num_labels; ++j) {
                // 计算连通域之间的距离
                double distance = calculateDistance(moments(labels == i), moments(labels == j));
                if (distance < threshold) {
                    // 合并连通域
                    labels.setTo(i, labels == j);
                }
            }
        }

        // 更新标记，确保每个连通域都具有唯一的标记
        Mat mergedLabels = labels.clone();
        for (int i = 1; i < num_labels; ++i) {
            mergedLabels.setTo(i, labels == i);
        }

        labels = mergedLabels.clone(); // 更新原始标记
    }
}

void FretThaSolver::extractRoiByAdaptiveMask(QString viewName) {
    QElapsedTimer timer;
    timer.start();
    using namespace cv;

    // 设置参数
    int kernelSize = 21;
    int roiSize = 5;

    // 读数据
    Mat matSrc[3];
    for (int i = 0; i < 3; ++ i) {
        matSrc[i] = imageProcessor->getMatrixCopy(ChannelName(i));
    }

    // 基于灰度值的自适应分割掩膜
    double valBg[3];
    Mat matMean[3], matBlur[3], maskGrayVal[3], abs_mask[3];
    for (int i = 0; i < 3; ++ i) {
        matSrc[i] = imageProcessor->getMatrixCopy(ChannelName(i));
        matBlur[i] = wiz::medianFilter16U(matSrc[i], roiSize);
        valBg[i] = wiz::calcBackgroundGray(matBlur[i]);
        matMean[i] = wiz::meanFilter16U(matSrc[i], roiSize);
        cv::threshold(matBlur[i], abs_mask[i], valBg[i]*3, 255, cv::THRESH_BINARY);
        abs_mask[i].convertTo(abs_mask[i], CV_8UC1);
        wiz::adaptiveThreshold16(matBlur[i], maskGrayVal[i], 255, cv::ADAPTIVE_THRESH_MEAN_C, cv::THRESH_BINARY, kernelSize, -valBg[i]);
        // 删除四周的边缘
        for (int j = 0; j < kernelSize; ++ j) {
            maskGrayVal[i].row(j).setTo(0);
            maskGrayVal[i].row(maskGrayVal[i].rows - j - 1).setTo(0);
            maskGrayVal[i].col(j).setTo(0);
            maskGrayVal[i].col(maskGrayVal[i].cols - j - 1).setTo(0);
        }
    }

    Mat merged_abs_mask = abs_mask[0] & abs_mask[1] & abs_mask[2];
    Mat maskAdaptiveMerge(matSrc[0].size(), CV_8UC1, Scalar(0));

    for (int i = 0; i < CHANNEL_NAME_SIZE; ++ i) {
        maskAdaptiveMerge = maskAdaptiveMerge | maskGrayVal[i];
        cv::imwrite(
            (viewFolderPath + "/maskAdaptive" + QString::number(i) + ".tif").toStdString(),
            maskGrayVal[i]
        );
    }
    cv::imwrite(
        (viewFolderPath + "/maskAdaptiveMerge.tif").toStdString(),
        maskAdaptiveMerge
    );
    // cv::imwrite(
    //     (viewFolderPath + "/merged_abs_mask.tif").toStdString(),
    //     merged_abs_mask
    //     );

    maskAdaptiveMerge = maskAdaptiveMerge | merged_abs_mask;

    // 基于标准差的自适应分割掩膜
    Mat matStd(matSrc[0].size(), CV_64FC1, Scalar(0.0));
    Mat maskStdDev[3];
    for (int i = 0; i < 3; ++ i) {
        Mat localStd = wiz::localStandardDeviation(matBlur[i], roiSize, true);
        matStd += localStd;
        Mat localStd16U = wiz::normalizeByMinMax16U(localStd);
        // cv::imwrite(
        //     (viewFolderPath + "/localStd" + QString::number(i) + ".tif").toStdString(),
        //     localStd16U
        // );
        double backgroundStd = wiz::calcBackgroundGray(localStd16U);
        // 自适应阈值法
        wiz::adaptiveThreshold16(localStd16U, maskStdDev[i], 255, cv::ADAPTIVE_THRESH_MEAN_C, cv::THRESH_BINARY_INV, kernelSize, backgroundStd);
        // // 绝对阈值法
        // cv::threshold(localStd16U, maskStdDev[i], backgroundStd, 255, cv::THRESH_BINARY_INV);
        // maskStdDev[i].convertTo(maskStdDev[i], CV_8UC1);
    }
    Mat maskStdDevMerge(matSrc[0].size(), CV_8UC1, Scalar(255));
    for (int i = 0; i < CHANNEL_NAME_SIZE; ++ i) {
        maskStdDevMerge = maskStdDevMerge & maskStdDev[i];
        cv::imwrite(
            (viewFolderPath + "/maskStdDev" + QString::number(i) + ".tif").toStdString(),
            maskStdDev[i]);
    }
    cv::imwrite(
        (viewFolderPath + "/maskStdDevMerge.tif").toStdString(),
        maskStdDevMerge);

    // 生成最终掩膜
    Mat maskFinal = maskAdaptiveMerge & maskStdDevMerge;

    for (int i = 0; i < 3; ++ i) {
        maskFinal = wiz::getMorphologyClose(maskFinal, 3); // 进行形态学操作去除小的噪点
        maskFinal = wiz::getMorphologyOpen(maskFinal, 3);
    }

    cv::imwrite(
        (viewFolderPath + "/maskFinal.tif").toStdString(),
        maskFinal);

    std::vector<WizPoint> centroids = getConnectedComponentMaxValue(maskFinal,
                                                                    matMean[CHANNEL_NAME_AA] / valBg[CHANNEL_NAME_AA] + matMean[CHANNEL_NAME_DA] / valBg[CHANNEL_NAME_DA] + matMean[CHANNEL_NAME_DD] / valBg[CHANNEL_NAME_DD]);

    // int kernel_size = 5;
    for (const WizPoint& centroid : centroids) {
        int x = centroid.x;
        int y = centroid.y;

        int rectx = x - roiSize / 2;
        int recty = y - roiSize / 2;
        int rectw = roiSize;
        int recth = roiSize;
        if (rectx < 0 || recty < 0 || rectx + rectw >= matSrc[0].cols || recty + recth >= matSrc[0].rows) continue;

        // 计算MatSrc在rect区域的平均灰度值
        double sbr = 0.0;
        cv::Rect roi(rectx, recty, rectw, recth);
        double valFg[3];
        for (int i = 0; i < 3; ++ i) {
            Mat matRoi = matSrc[i](roi);
            Scalar mean = cv::mean(matRoi);
            valFg[i] = mean.val[0];

            sbr += valFg[i] / valBg[i];
        }
        if (sbr < 6) continue;

        double valGrayAA = valFg[CHANNEL_NAME_AA] - valBg[CHANNEL_NAME_AA];
        double valGrayDA = valFg[CHANNEL_NAME_DA] - valBg[CHANNEL_NAME_DA];
        double valGrayDD = valFg[CHANNEL_NAME_DD] - valBg[CHANNEL_NAME_DD];

        // 使用calculator进行计算
        calculator->loadCorrectedGray(valGrayAA, valGrayDA, valGrayDD);
        calculator->calcData();
        if (!calculator->isDataCalculated()) continue;
        // 获取数据
        double valEa = calculator->getResult(Ea);
        double valEd = calculator->getResult(Ed);
        double valRad = calculator->getResult(Rad);
        double valRda = calculator->getResult(Rda);
        double valAest = calculator->getResult(Aest);
        double valDest = calculator->getResult(Dest);

        // if (true) {
        if (valEa > 0 && valEd > 0 && valRad > 0 && valRda > 0 && valAest > 0 && valDest > 0) {
            QMap<TableHeader, QString> map;
            map.insert(TABLE_HEADER_IAA, QString::number(valGrayAA));
            map.insert(TABLE_HEADER_IDA, QString::number(valGrayDA));
            map.insert(TABLE_HEADER_IDD, QString::number(valGrayDD));
            map.insert(TABLE_HEADER_EA, QString::number(valEa));
            map.insert(TABLE_HEADER_ED, QString::number(valEd));
            map.insert(TABLE_HEADER_RAD, QString::number(valRad));
            map.insert(TABLE_HEADER_RDA, QString::number(valRda));
            map.insert(TABLE_HEADER_AEST, QString::number(valAest));
            map.insert(TABLE_HEADER_DEST, QString::number(valDest));
            map.insert(TABLE_HEADER_RECTX, QString::number(rectx));
            map.insert(TABLE_HEADER_RECTY, QString::number(recty));
            map.insert(TABLE_HEADER_RECTW, QString::number(rectw));
            map.insert(TABLE_HEADER_RECTH, QString::number(recth));
            map.insert(TABLE_HEADER_VIEW, viewName);
            if (rectx >= 0 && recty >= 0) emit sendData(map);
        }
    }

    qDebug() << "[Func Time Cost]:\t" << timer.elapsed();
}

void FretThaSolver::extractRoiByFretImage(QString viewName) {
    QElapsedTimer timer;
    timer.start();
    using namespace cv;

    // 设置参数
    int kernelSize = 21;

    // 读数据
    Mat matSrc[3];
    for (int i = 0; i < 3; ++ i) {
        matSrc[i] = imageProcessor->getMatrixCopy(ChannelName(i));
    }
    double valBg[3];
    for (int i = 0; i < 3; ++ i) {
        valBg[i] = wiz::calcBackgroundGray(matSrc[i]);
    }
    Mat matEd = imageProcessor->getRawResult(Ed), maskEd;
    Mat localStd = wiz::localStandardDeviation(matEd, kernelSize, false);
    Mat localStd16U = wiz::normalizeByMinMax16U(localStd);
    wiz::adaptiveThreshold16(localStd16U, maskEd, 255, cv::ADAPTIVE_THRESH_GAUSSIAN_C, cv::THRESH_BINARY_INV, kernelSize, 0);
    cv::imwrite(
        (viewFolderPath + "/maskEd.tif").toStdString(),
        maskEd
    );


    // 基于灰度值的自适应分割掩膜
    Mat matMean[3], matBlur[3], maskGrayVal[3];
    for (int i = 0; i < 3; ++ i) {
        matSrc[i] = imageProcessor->getMatrixCopy(ChannelName(i));
        matBlur[i] = wiz::medianFilter16U(matSrc[i], 5);
        matMean[i] = wiz::meanFilter16U(matSrc[i], 5);
        wiz::adaptiveThreshold16(matBlur[i], maskGrayVal[i], 255, cv::ADAPTIVE_THRESH_MEAN_C, cv::THRESH_BINARY, kernelSize, 0);
        // 删除四周的边缘
        for (int j = 0; j < kernelSize; ++ j) {
            maskGrayVal[i].row(j).setTo(0);
            maskGrayVal[i].row(maskGrayVal[i].rows - j - 1).setTo(0);
            maskGrayVal[i].col(j).setTo(0);
            maskGrayVal[i].col(maskGrayVal[i].cols - j - 1).setTo(0);
        }
    }
    Mat maskAdaptiveMerge(matSrc[0].size(), CV_8UC1, Scalar(0));
    for (int i = 0; i < CHANNEL_NAME_SIZE; ++ i) {
        maskAdaptiveMerge = maskAdaptiveMerge | maskGrayVal[i];
        cv::imwrite(
            (viewFolderPath + "/maskAdaptive" + QString::number(i) + ".tif").toStdString(),
            maskGrayVal[i]);
    }
    cv::imwrite(
        (viewFolderPath + "/maskAdaptiveMerge.tif").toStdString(),
        maskAdaptiveMerge
    );

    // 生成最终掩膜
    Mat maskFinal = maskAdaptiveMerge & maskEd;

    for (int i = 0; i < 3; ++ i) {
        maskFinal = wiz::getMorphologyClose(maskFinal, 3); // 进行形态学操作去除小的噪点
        maskFinal = wiz::getMorphologyOpen(maskFinal, 3);
    }

    cv::imwrite(
        (viewFolderPath + "/maskFinal.tif").toStdString(),
        maskFinal);

    // 计算质心
    // std::vector<WizPoint> centroids = getConnectedComponentMaxMean(maskFinal, matMean[CHANNEL_NAME_DA]);
    std::vector<WizPoint> centroids = getConnectedComponentCentroids(maskFinal);

    int kernel_size = 5;
    for (const auto& centroid : centroids) {
        int x = centroid.x;
        int y = centroid.y;

        int rectx = x - kernel_size / 2;
        int recty = y - kernel_size / 2;
        int rectw = kernel_size;
        int recth = kernel_size;

        // 计算MatSrc在rect区域的平均灰度值
        double sbr = 0.0;
        cv::Rect roi(rectx, recty, rectw, recth);
        double valFg[3];
        for (int i = 0; i < 3; ++ i) {
            Mat matRoi = matSrc[i](roi);
            Scalar mean = cv::mean(matRoi);
            valFg[i] = mean.val[0];

            if (valFg[i] / valBg[i] < 2) {
                sbr += -1e7;
                continue;
            }
        }
        if (sbr < 0) continue;

        double valGrayAA = valFg[CHANNEL_NAME_AA] - valBg[CHANNEL_NAME_AA];
        double valGrayDA = valFg[CHANNEL_NAME_DA] - valBg[CHANNEL_NAME_DA];
        double valGrayDD = valFg[CHANNEL_NAME_DD] - valBg[CHANNEL_NAME_DD];

        // 使用calculator进行计算
        calculator->loadCorrectedGray(valGrayAA, valGrayDA, valGrayDD);
        calculator->calcData();
        if (!calculator->isDataCalculated()) continue;
        // 获取数据
        double valEa = calculator->getResult(Ea);
        double valEd = calculator->getResult(Ed);
        double valRad = calculator->getResult(Rad);
        double valRda = calculator->getResult(Rda);
        double valAest = calculator->getResult(Aest);
        double valDest = calculator->getResult(Dest);

        if (valEa > 0 && valEd > 0 && valRad > 0 && valRda > 0 && valAest > 0 && valDest > 0) {
            QMap<TableHeader, QString> map;
            map.insert(TABLE_HEADER_IAA, QString::number(valGrayAA));
            map.insert(TABLE_HEADER_IDA, QString::number(valGrayDA));
            map.insert(TABLE_HEADER_IDD, QString::number(valGrayDD));
            map.insert(TABLE_HEADER_EA, QString::number(valEa));
            map.insert(TABLE_HEADER_ED, QString::number(valEd));
            map.insert(TABLE_HEADER_RAD, QString::number(valRad));
            map.insert(TABLE_HEADER_RDA, QString::number(valRda));
            map.insert(TABLE_HEADER_AEST, QString::number(valAest));
            map.insert(TABLE_HEADER_DEST, QString::number(valDest));
            map.insert(TABLE_HEADER_RECTX, QString::number(rectx));
            map.insert(TABLE_HEADER_RECTY, QString::number(recty));
            map.insert(TABLE_HEADER_RECTW, QString::number(rectw));
            map.insert(TABLE_HEADER_RECTH, QString::number(recth));
            map.insert(TABLE_HEADER_VIEW, viewName);
            if (rectx >= 0 && recty >= 0) emit sendData(map);
        }
    }

    qDebug() << "[Func Time Cost]:\t" << timer.elapsed();
}

void FretThaSolver::scorePixel() {
    QElapsedTimer timer;
    timer.start();
    using namespace cv;

    // 读数据
    Mat matSrc[3];
    for (int i = 0; i < 3; ++ i) {
        matSrc[i] = imageProcessor->getMatrixCopy(ChannelName(i));
    }

    // 产生中间数据
    Mat matMean[3], matBlur[3];
    for (int i = 0; i < 3; ++ i) {
        matSrc[i] = imageProcessor->getMatrixCopy(ChannelName(i));
        matBlur[i] = wiz::medianFilter16U(matSrc[i], 5);
        matMean[i] = wiz::meanFilter16U(matSrc[i], 5);
    }

    // 生成灰度值分数
    Mat scoreBright[3];
    for (int i = 0; i < 3; ++ i) {
        scoreBright[i] = wiz::normalizeByMinMax(matMean[i]);
        cv::exp(scoreBright[i], scoreBright[i]);
        scoreBright[i] = scoreBright[i] - 1.0;
        scoreBright[i] = wiz::normalizeByMinMax(scoreBright[i]);
    }

    // 生成标准差分数
    Mat floatStd(matSrc[0].size(), CV_64FC1, Scalar(0.0));
    Mat maskStd(matSrc[0].size(), CV_8UC1, Scalar(255));
    Mat scoreStdDev[3];
    for (int i = 0; i < 3; ++ i) {
        Mat localStd = wiz::localStandardDeviation(matBlur[i], 5, true);

        Mat matStd = wiz::normalizeByMinMax(localStd);
        Mat matStd16U = wiz::normalizeByMinMax16U(localStd);

        scoreStdDev[i] = 1.0 - matStd;
        cv::exp(scoreStdDev[i], scoreStdDev[i]);
        scoreStdDev[i] = scoreStdDev[i] - 1.0;
        scoreStdDev[i] = wiz::normalizeByMinMax(scoreStdDev[i]);
        floatStd += scoreStdDev[i];

        Mat mask;
        double thresh = wiz::calcBackgroundGray(matStd16U);
        Mat matStdBlured16U = wiz::medianFilter16U(matStd16U, 5);
        cv::threshold(matStdBlured16U, mask, thresh * 2, 255, THRESH_BINARY_INV);

        mask.convertTo(mask, CV_8UC1);
        bitwise_and(maskStd, mask, maskStd);
    }

    for (int i = 0; i < 3; ++ i) {
        maskStd = wiz::getMorphologyClose(maskStd, 11);
        maskStd = wiz::erodeImage(maskStd, 11);
        maskStd = wiz::getMorphologyOpen(maskStd, 11);
        maskStd = wiz::dilateImage(maskStd, 11);
    }

    Mat maskBg = 255 - wiz::getLargestConnectedComponent(maskStd);
    cv::imwrite(
        (viewFolderPath + "/maskCellSegment" + ".tif").toStdString(),
        maskBg
    );

    // Mat product;
    Mat sum(matSrc[CHANNEL_NAME_AA].size(), CV_64FC1, Scalar(0.0));
    Mat sumBright(matSrc[CHANNEL_NAME_AA].size(), CV_64FC1, Scalar(0.0));
    Mat sumStdDev(matSrc[CHANNEL_NAME_AA].size(), CV_64FC1, Scalar(0.0));
    for (int i = 0; i < 3; ++ i) {
        sumBright += scoreBright[i];
        sumStdDev += scoreStdDev[i];
    }
    add(sumBright, sumStdDev, sum);

    // 生成最后结果
    Mat result, resultDisp;
    result = wiz::normalizeByMinMax(sum);
    resultDisp = wiz::normalizeByMinMax8U(result);
    Mat maskAdaptive;
    cv::adaptiveThreshold(resultDisp, maskAdaptive, 255, cv::ADAPTIVE_THRESH_MEAN_C, cv::THRESH_BINARY, 11, -2);

    for (int i = 0; i < 1; ++ i) {
        int m = 1;
        while (m --) {
            maskAdaptive = wiz::erodeImage(maskAdaptive, 3);
            maskAdaptive = wiz::dilateImage(maskAdaptive, 3);
        }
        maskAdaptive = wiz::erodeImage(maskAdaptive, 3);
    }

    // 基于opencv获取maskAdaptive中的所有连通域，去除面积过小的连通域
    Mat maskAdaptive8U;
    maskAdaptive.convertTo(maskAdaptive8U, CV_8UC1);
    std::vector<std::vector<cv::Point>> contours;
    std::vector<cv::Vec4i> hierarchy;
    cv::findContours(maskAdaptive8U, contours, hierarchy, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_SIMPLE);
    for (uint i = 0; i < contours.size(); ++ i) {
        double area = cv::contourArea(contours[i]);
        if (area < 50) {
            cv::drawContours(maskAdaptive8U, contours, i, cv::Scalar(0), -1);
        }
    }

    // 计算质心
    std::vector<WizPoint> centroids = getConnectedComponentCentroids(maskAdaptive8U);

    qDebug() << "[Func Time Cost]:\t\t" << timer.elapsed();
    cv::imwrite((viewFolderPath + "/FinalScore" + ".tif").toStdString(), resultDisp);
    cv::imwrite((viewFolderPath + "/FinalMask" + ".tif").toStdString(), maskAdaptive);
    cv::imwrite((viewFolderPath + "/FinalMask8U" + ".tif").toStdString(), maskAdaptive8U);

}

void FretThaSolver::extractRoiOrigin(QString viewName) {
    using namespace cv;

    // 赋值图片数据
    Mat matSrc[3];
    for (int i = 0; i < 3; ++ i) {
        matSrc[i] = imageProcessor->getMatrixCopy(ChannelName(i));
    }

    // 对数据进行均值滤波，并且计算得到背景灰度值
    Mat matFiltered[3];
    Mat matData[3];
    double grayBackground[3];
    for (int i = 0; i < 3; ++ i) {
        matFiltered[i] = wiz::gaussianFilter16U(matSrc[i], 5, 1);
        matData[i] = wiz::meanFilter16U(matFiltered[i], 5);
        grayBackground[i] = wiz::calcBackgroundGray(matFiltered[i]);
    }

    // 根据SBR阈值进行滤波
    double threshRatio[3] = {3.0, 3.0, 3.0};
    Mat maskSingleChannel[3];
    Mat maskSbr, maskSbr8U;    // 0-1
    for (int i = 0; i < 3; ++ i) {
        threshold(matData[i], maskSingleChannel[i], grayBackground[i] * threshRatio[i], 255, THRESH_BINARY);
        maskSingleChannel[i].convertTo(maskSingleChannel[i], CV_8U);
    }
    bitwise_and(maskSingleChannel[CHANNEL_NAME_AA], maskSingleChannel[CHANNEL_NAME_DA], maskSbr8U);
    bitwise_and(maskSbr8U, maskSingleChannel[CHANNEL_NAME_DD], maskSbr8U);
    maskSbr = wiz::normalizeByMinMax(maskSbr8U);
    //if (RUNMODE == RUNMODE_EXPERIMENT) cv::imwrite((viewFolderPath + "/maskSBRatio.tif").toStdString(), maskSbr);

    // 按照我们设计的算法生成评价结果
    Mat maskStd;
    Mat scoreStd(matSrc[CHANNEL_NAME_AA].size(), CV_64FC1, Scalar(0.0));
    Mat maskBg(matSrc[CHANNEL_NAME_AA].size(), CV_8UC1, cv::Scalar(0));

    for (int i = 0; i < 3; ++ i) {

        Mat matStd = wiz::localStandardDeviation(matFiltered[i], 5);
        Mat matStd8U = wiz::normalizeByMinMax8U(matStd);

        Mat mask;
        threshold(matStd8U, mask, 4, 255, cv::THRESH_BINARY);
        Mat kernel = getStructuringElement(MORPH_RECT, Size(3, 3));  // 创建一个 3x3 的矩形内核
        cv::Mat dilatedImage;
        mask = wiz::removeScatter(mask, 10, 0.5);
        if (RUNMODE == RUNMODE_EXPERIMENT) cv::imwrite((viewFolderPath + "/mask" + QString::number(i) + ".tif").toStdString(), mask);
        dilate(mask, dilatedImage, kernel);

        bitwise_or(dilatedImage, maskBg, maskBg);   // 存到maskBg中

        Mat matStdMinima = wiz::findLocalMinima(matStd8U);
        // Mat matStdMinima8U = FretImageProcessor::normalizeByMinMax8U(matStdMinima);

        scoreStd = scoreStd + wiz::normalizeByMinMax(matStdMinima);
    }

    // 定义膨胀操作的内核
    int m = 1;
    while (m --) {
        maskBg = wiz::removeScatter(maskBg, 5);
        Mat kernel = getStructuringElement(MORPH_RECT, Size(3, 3));  // 创建一个 3x3 的矩形内核
        dilate(maskBg, maskBg, kernel);
    }
    Mat maskCell = maskBg.clone();

    if (RUNMODE == RUNMODE_EXPERIMENT) cv::imwrite((viewFolderPath + "/BackgroundMask" + ".tif").toStdString(), maskBg);
    maskBg = wiz::getLargestConnectedComponent(255 - maskBg);
    if (RUNMODE == RUNMODE_EXPERIMENT) cv::imwrite((viewFolderPath + "/BackgroundAreaMask" + ".tif").toStdString(), maskBg);

    Mat scoreStd16U = wiz::normalizeByMinMax16U(scoreStd);
    cv::threshold(scoreStd16U, maskStd, 0, 255, cv::THRESH_BINARY | cv::THRESH_OTSU);
    maskStd.convertTo(maskStd, CV_8U);

    // 合成最后的结果
    Mat maskUlt, scoreMat;
    scoreMat = maskSbr.mul(scoreStd);
    bitwise_and(maskSbr8U, maskStd, maskUlt);
    //if (RUNMODE == RUNMODE_EXPERIMENT) cv::imwrite((viewFolderPath + "/FinalJudge.tif").toStdString(), maskUlt);

    // 图像增强
    if (RUNMODE == RUNMODE_EXPERIMENT) {

        //    给王婧臻
        //    QDir dir("D:/test");
        //    if (dir.exists()) {
        //        QString path;
        //        path = "D:/test/Ed" + viewName + ".tif";
        //        imageProcessor->preProcessData();
        //        imageProcessor->calcEFret();
        //        if (RUNMODE == RUNMODE_EXPERIMENT) cv::imwrite(path.toStdString(), imageProcessor->getPseuImage(Ed));
        //        path = "D:/test/Rc" + viewName + ".tif";
        //        if (RUNMODE == RUNMODE_EXPERIMENT) cv::imwrite(path.toStdString(), imageProcessor->getPseuImage(Rad));
        //        path = "D:/test/Merge" + viewName + ".tif";
        //        if (RUNMODE == RUNMODE_EXPERIMENT) cv::imwrite(path.toStdString(), imageProcessor->getMergedImage());
        //        path = "D:/test/AA" + viewName + ".tif";
        //        if (RUNMODE == RUNMODE_EXPERIMENT) cv::imwrite(path.toStdString(), imageProcessor->getNormalizedImage(AA));
        //        path = "D:/test/DA" + viewName + ".tif";
        //        if (RUNMODE == RUNMODE_EXPERIMENT) cv::imwrite(path.toStdString(), imageProcessor->getNormalizedImage(DA));
        //        path = "D:/test/DD" + viewName + ".tif";
        //        if (RUNMODE == RUNMODE_EXPERIMENT) cv::imwrite(path.toStdString(), imageProcessor->getNormalizedImage(DD));
        //    }

        Mat matEnh[3];
        for (int i = 0; i < 3; ++ i) {
            matEnh[i] = wiz::enhanceImage(matSrc[i], 64);
            matEnh[i] = wiz::gaussianFilter16U(matEnh[i], 5, 1);
            imwrite((viewFolderPath + "/Enhanced" + QString::number(i) + ".tif").toStdString(),
                    matEnh[i]);
            matEnh[i] = wiz::applyMaskToImage(matEnh[i], maskCell);
            imwrite((viewFolderPath + "/EnhancedMasked" + QString::number(i) + ".tif").toStdString(),
                    matEnh[i]);
        }
    }

    // 计算质心
    std::vector<WizPoint> centroids = getConnectedComponentCentroids(maskUlt);

    int kernel_size = 5;

    for (const auto& centroid : centroids) {
        int x = centroid.x;
        int y = centroid.y;
        double valBg[3];
        for (int i = 0; i < 3; ++ i) {
            valBg[i] = wiz::calculateAverageGrayValue(matData[CHANNEL_NAME_AA], maskBg, cv::Point(x, y), 1024);
        }
        double valGrayAA = matData[CHANNEL_NAME_AA].at<ushort>(y, x) - valBg[CHANNEL_NAME_AA];
        double valGrayDA = matData[CHANNEL_NAME_DA].at<ushort>(y, x) - valBg[CHANNEL_NAME_DA];
        double valGrayDD = matData[CHANNEL_NAME_DD].at<ushort>(y, x) - valBg[CHANNEL_NAME_DD];
        // 使用calculator进行计算
        calculator->loadCorrectedGray(valGrayAA, valGrayDA, valGrayDD);
        calculator->calcData();
        if (!calculator->isDataCalculated()) continue;
        // 获取数据
        double valEa = calculator->getResult(Ea);
        double valEd = calculator->getResult(Ed);
        double valRad = calculator->getResult(Rad);
        double valRda = calculator->getResult(Rda);
        double valAest = calculator->getResult(Aest);
        double valDest = calculator->getResult(Dest);
        int rectx = x - kernel_size / 2;
        int recty = y - kernel_size / 2;
        int rectw = kernel_size;
        int recth = kernel_size;

        QMap<TableHeader, QString> map;
        map.insert(TABLE_HEADER_IAA, QString::number(valGrayAA));
        map.insert(TABLE_HEADER_IDA, QString::number(valGrayDA));
        map.insert(TABLE_HEADER_IDD, QString::number(valGrayDD));
        map.insert(TABLE_HEADER_EA, QString::number(valEa));
        map.insert(TABLE_HEADER_ED, QString::number(valEd));
        map.insert(TABLE_HEADER_RAD, QString::number(valRad));
        map.insert(TABLE_HEADER_RDA, QString::number(valRda));
        map.insert(TABLE_HEADER_AEST, QString::number(valAest));
        map.insert(TABLE_HEADER_DEST, QString::number(valDest));
        map.insert(TABLE_HEADER_RECTX, QString::number(rectx));
        map.insert(TABLE_HEADER_RECTY, QString::number(recty));
        map.insert(TABLE_HEADER_RECTW, QString::number(rectw));
        map.insert(TABLE_HEADER_RECTH, QString::number(recth));
        map.insert(TABLE_HEADER_VIEW, viewName);
        if (rectx >= 0 && recty >= 0) emit sendData(map);
    }
}

void FretThaSolver::generateRoiNew() {
    using namespace cv;

    QElapsedTimer timer;
    timer.start();

    // 取数据，并得到中间数据
    Mat matSrc[3];
    for (int i = 0; i < 3; ++ i) {
        matSrc[i] = imageProcessor->getMatrixCopy(ChannelName(i));
    }
    Mat matMean[3], matBlur[3];
    Mat maskAdaptive(matSrc[0].size(), CV_8UC1, Scalar(255));

    for (int i = 0; i < 3; ++ i) {
        matSrc[i] = imageProcessor->getMatrixCopy(ChannelName(i));
        matBlur[i] = wiz::gaussianFilter16U(matSrc[i], 3, 1);
        matMean[i] = wiz::meanFilter16U(matBlur[i], 3);

        Mat mat8U, dst;
        mat8U = wiz::normalizeByZeroMax8U(matBlur[i]);
        cv::adaptiveThreshold(mat8U, dst, 255, cv::ADAPTIVE_THRESH_MEAN_C, cv::THRESH_BINARY, 51, 2);
        cv::imwrite(
            (viewFolderPath + "/maskAdaptive" + QString::number(i) + ".tif").toStdString(),
            dst
        );

        dst.convertTo(dst, CV_8UC1);

        bitwise_and(maskAdaptive, dst, maskAdaptive);

    }
    cv::imwrite(
        (viewFolderPath + "/maskAdaptiveMerge" + ".tif").toStdString(),
        maskAdaptive
    );

    for (int i = 0; i < 10; ++ i) {
        maskAdaptive = wiz::getMorphologyClose(maskAdaptive, 5);
        maskAdaptive = wiz::erodeImage(maskAdaptive, 5);
        maskAdaptive = wiz::getMorphologyOpen(maskAdaptive, 5);

        Mat maskBg = wiz::getLargestConnectedComponent(maskAdaptive);

        cv::imwrite(
            (viewFolderPath + "/maskAdaptiveMerge" + QString::number(i) + ".tif").toStdString(),
            maskBg
        );
    }

    qDebug() << "[Function Time Cost]:\t\t" << timer.elapsed();
}

void FretThaSolver::generateRoiFromImage(QString viewName)
{
    if (RUNMODE == RUNMODE_EXPERIMENT) {
        cellSegmentation();

        // Segmentation();
    } else if (RUNMODE == RUNMODE_RELEASE || RUNMODE == RUNMODE_DEBUG) {
        extractRoiByAdaptiveMask(viewName);
    }

}


int FretThaSolver::binDataByRda(double min, double max, double interval)
{
    using namespace std;
    qDebug() << "[Bin Data]:\tFrom"
             << min
             << "To" << max
             << "By" << interval;

    // 清空bin数据
    clearFretDataBin();

    // 定义数据
    int num = (max - min) / interval;
    std::vector<std::vector<double>> Easum(num), Edsum(num), Destsum(num), Aestsum(num);
    // 对所有数据进行存入
    for (int i = 0; i < (int)vecFretData[Ea].size(); ++ i)
    {
        // 按照Rad划分区间并存入数据
        // 原本
        double Rda = (vecFretData[Dest][i] / vecFretData[Aest][i]);
        int idx = (Rda - min) / interval;
        if (idx < 0 || idx >= num) continue;

        Easum[idx].push_back(vecFretData[Ea][i]);
        Edsum[idx].push_back(vecFretData[Ed][i]);
        Destsum[idx].push_back(vecFretData[Dest][i]);
        Aestsum[idx].push_back(vecFretData[Aest][i]);
    }
    //按区间内均值生成bin的数据
    int cnt = 0;
    for (int i = 0; i < num; ++ i) {
        if(!Easum[i].size()) continue;
        //计算均值
        double sum = accumulate(begin(Easum[i]), end(Easum[i]), 0.0);
        double newEa = sum / Easum[i].size();
        sum = accumulate(begin(Edsum[i]), end(Edsum[i]), 0.0);
        double newEd = sum / Edsum[i].size();
        sum = accumulate(begin(Aestsum[i]), end(Aestsum[i]), 0.0);
        double newAest = sum / Aestsum[i].size();
        sum = accumulate(begin(Destsum[i]), end(Destsum[i]), 0.0);
        double newDest = sum / Destsum[i].size();

        vecFretDataBin[Ea].push_back(newEa);
        vecFretDataBin[Ed].push_back(newEd);
        vecFretDataBin[Aest].push_back(newAest);
        vecFretDataBin[Dest].push_back(newDest);
        cnt ++ ;
    }

    return cnt;
}

int FretThaSolver::binDataByRad(double min, double max, double interval)
{
    using namespace std;
    qDebug() << "[Bin Data]:\tFrom"
             << min
             << "To" << max
             << "By" << interval;

    // 清空bin数据
    clearFretDataBin();

    // 定义数据
    int num = (max - min) / interval;
    std::vector<std::vector<double>> Easum(num), Edsum(num), Destsum(num), Aestsum(num);
    // 对所有数据进行存入
    for (int i = 0; i < (int)vecFretData[Ea].size(); ++ i)
    {
        // 按照Rad划分区间并存入数据
        double Rad = (vecFretData[Aest][i] / vecFretData[Dest][i]);
        int idx = (Rad - min) / interval;
        if (idx < 0 || idx >= num) continue;

        Easum[idx].push_back(vecFretData[Ea][i]);
        Edsum[idx].push_back(vecFretData[Ed][i]);
        Destsum[idx].push_back(vecFretData[Dest][i]);
        Aestsum[idx].push_back(vecFretData[Aest][i]);
    }
    //按区间内均值生成bin的数据
    int cnt = 0;
    for (int i = 0; i < num; ++ i) {
        if(!Easum[i].size()) continue;
        //计算均值
        double sum = accumulate(begin(Easum[i]), end(Easum[i]), 0.0);
        double newEa = sum / Easum[i].size();
        sum = accumulate(begin(Edsum[i]), end(Edsum[i]), 0.0);
        double newEd = sum / Edsum[i].size();
        sum = accumulate(begin(Aestsum[i]), end(Aestsum[i]), 0.0);
        double newAest = sum / Aestsum[i].size();
        sum = accumulate(begin(Destsum[i]), end(Destsum[i]), 0.0);
        double newDest = sum / Destsum[i].size();

        vecFretDataBin[Ea].push_back(newEa);
        vecFretDataBin[Ed].push_back(newEd);
        vecFretDataBin[Aest].push_back(newAest);
        vecFretDataBin[Dest].push_back(newDest);
        cnt ++ ;
    }

    return cnt;
}

int FretThaSolver::binDataOurs(double min, double max, double interval)
{
    using namespace std;

    // data weight calculation
    int size = vecFretData[Aest].size();
    std::vector<double> da_ratios(size);
    for (int i = 0; i < size; ++i) {
        da_ratios[i] = std::log(vecFretData[Dest][i] / vecFretData[Aest][i]);
    }
    double mean_da_ratio = std::accumulate(da_ratios.begin(), da_ratios.end(), 0.0) / da_ratios.size();
    double stddev_da_ratio = std::sqrt(std::accumulate(
                                           da_ratios.begin(),
                                           da_ratios.end(),
                                           0.0,
                                           [mean_da_ratio](double sum, double ratio) {
                                               return sum + std::pow(ratio - mean_da_ratio, 2);
                                           }) / da_ratios.size());
    // std::cout << "[Mean DA Ratio]:\t" << mean_da_ratio << std::endl;
    // std::cout << "[Stddev DA Ratio]:\t" << stddev_da_ratio << std::endl;

    min = mean_da_ratio - 2 * stddev_da_ratio;
    max = mean_da_ratio + 2 * stddev_da_ratio;
    interval = (max - min) / 25;

    qDebug() << "[Bin Data]:\tFrom"
             << min
             << "To" << max
             << "By" << interval;

    // 清空bin数据
    clearFretDataBin();

    // 定义数据
    int num = (max - min) / interval;
    std::vector<std::vector<double>> Easum(num), Edsum(num), Destsum(num), Aestsum(num);
    // 对所有数据进行存入
    for (int i = 0; i < (int)vecFretData[Ea].size(); ++ i)
    {
        // 按照Rad划分区间并存入数据
        double Rda = log(vecFretData[Dest][i] / vecFretData[Aest][i]);
        int idx = (Rda - min) / interval;
        if (idx < 0 || idx >= num) continue;

        Easum[idx].push_back(vecFretData[Ea][i]);
        Edsum[idx].push_back(vecFretData[Ed][i]);
        Destsum[idx].push_back(vecFretData[Dest][i]);
        Aestsum[idx].push_back(vecFretData[Aest][i]);
    }
    //按区间内均值生成bin的数据
    int cnt = 0;
    for (int i = 0; i < num; ++ i) {
        if(!Easum[i].size()) continue;
        //计算均值
        double sum = accumulate(begin(Easum[i]), end(Easum[i]), 0.0);
        double newEa = sum / Easum[i].size();
        sum = accumulate(begin(Edsum[i]), end(Edsum[i]), 0.0);
        double newEd = sum / Edsum[i].size();
        sum = accumulate(begin(Aestsum[i]), end(Aestsum[i]), 0.0);
        double newAest = sum / Aestsum[i].size();
        sum = accumulate(begin(Destsum[i]), end(Destsum[i]), 0.0);
        double newDest = sum / Destsum[i].size();

        vecFretDataBin[Ea].push_back(newEa);
        vecFretDataBin[Ed].push_back(newEd);
        vecFretDataBin[Aest].push_back(newAest);
        vecFretDataBin[Dest].push_back(newDest);
        cnt ++ ;
    }

    return cnt;
}

void FretThaSolver::autoProcessActivity()
{
    qDebug() << "[Auto Process Data]:\t Start";
    void clearGrayData();
    void clearFretData();
    void clearFretDataBin();
    processOneBatch(batchFolderPath);
    // 计算得到FRET数据
    processGrayToFretData();
    // 进行计算
    performTwoHybridBen(0, 5, 0.1);
    performTwoHybridDu(0, 0.5, 2, 5);
    performTwoHybridOurs();
}

void FretThaSolver::autoGenerateActivity()
{
    qDebug() << "[Auto Generate Data]:\t Start";
    generateRoiFromBatch(batchFolderPath);
}

void FretThaSolver::autoImportActivity()
{
    qDebug() << "[Auto Import Mask]:\t Start";
    loadRoiFromBatch(batchFolderPath);
}

void FretThaSolver::setRunFunction(QString str)
{
    m_runFunction = str;
}

std::vector<double> FretThaSolver::getRegionMeans(const cv::Mat& data, const cv::Mat& mask)
{
    // 确保输入图像是8位单通道图像
    CV_Assert(mask.type() == CV_8UC1);
    cv::Mat converted;
    data.convertTo(converted, CV_64F);

    // 存储结果的向量
    std::vector<double> means;

    // 找到掩膜中的所有连通域
    cv::Mat labels;
    int num_labels = cv::connectedComponents(mask, labels, 8, CV_32S);

    // 对每个连通域，计算数据图中对应区域的灰度均值
    for (int i = 1; i < num_labels; ++i) {
        cv::Mat region = (labels == i);
        cv::Scalar mean = cv::mean(converted, region);
        means.push_back(mean[0]);
    }

    return means;
}

double FretThaSolver::calcApproach(double min, double max, const std::vector<double>& R, const std::vector<double>& E)
{
    double sum = 0;
    int cnt = 0;
    for (uint i = 0; i < R.size(); ++ i) {
        if (R[i] < min || R[i] > max) continue;
        sum += E[i];
        cnt ++;
    }
    double approach = sum / (double)cnt;

    return approach;
}

double FretThaSolver::calcSlope(double min, double max, const std::vector<double>& R, const std::vector<double>& E)
{
    double slope;
    std::vector<std::pair<double, double>> points;
    for (uint i = 0; i < R.size(); ++ i) {
        if (R[i] < min || R[i] > max) continue;
        points.push_back({R[i], E[i]});
    }
    slope = leastSquare(points);

    return slope;
}

void FretThaSolver::run()
{
    if (m_runFunction == "AutoGen") {
        autoGenerateActivity();
    } else if (m_runFunction == "ImportMask") {
        autoImportActivity();
    }

    emit thaFinished();
}
