
#include "fretthasolver.h"
#include "TwoHybridSolver/TwoHybridSolver.h"
#include "TwoHybridSolver/TwoHybridSolver_terminate.h"
#include "TwoHybridSolver/coder_array.h"

/*******************质心算法********************/

// 定义保存坐标的结构体
struct Point2D {
    int x;
    int y;
};

// 连通域对象
struct ConnectedComponent {
    std::vector<Point2D> pixels;
    Point2D centroid;
};

// 获取连通域的质心坐标
Point2D getCentroid(const std::vector<Point2D>& pixels) {
    Point2D centroid;
    int totalPixels = pixels.size();
    int sumX = 0;
    int sumY = 0;
    for (const auto& point : pixels) {
        sumX += point.x;
        sumY += point.y;
    }
    centroid.x = sumX / totalPixels;
    centroid.y = sumY / totalPixels;

    return centroid;
}

// 标记访问过的像素
void markVisited(cv::Mat& mask, int x, int y, std::vector<Point2D>& component) {
    Point2D point;
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

// 获取连通域的质心坐标向量
std::vector<Point2D> getConnectedComponentCentroids(const cv::Mat& mask) {
    std::vector<Point2D> centroids;
    cv::Mat tempMask = mask.clone();

    for (int y = 0; y < tempMask.rows; ++y) {
        for (int x = 0; x < tempMask.cols; ++x) {
            if (tempMask.at<uchar>(y, x) == 255) {
                std::vector<Point2D> component;
                markVisited(tempMask, x, y, component);

                if (component.size() > 5) {
                    Point2D centroid = getCentroid(component);
                    centroids.push_back(centroid);
                }
            }
        }
    }

    return centroids;
}

/*******************质心算法结尾********************/

FretThaSolver::FretThaSolver()
{

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
    qDebug() << "[Write Data]:\tStart";
    // 检查文件，若有同名文件会覆盖
    QFile file(savePath + "/FretThaData.csv");
    if (file.exists()) file.remove();
    file.open(QIODevice::Append);
    QTextStream out(&file);
    out << "Iaa,Ida,Idd,Rc,Ed,1/Rc,Ea,Aest,Dest,Afree,Dfree\n";
    uint i = 0;
    for (; i < vecGrayData[Ed].size(); ++ i){
        out << vecGrayData[AA][i] << ","
            << vecGrayData[DA][i] << ","
            << vecGrayData[DD][i] << ","
            << vecFretData[Rad][i] << ","
            << vecFretData[Ed][i] << ","
            << vecFretData[Rda][i] << ","
            << vecFretData[Ea][i] << ","
            << vecFretData[Aest][i] << ","
            << vecFretData[Dest][i] << ","
            << vecFretData[Afree][i] << ","
            << vecFretData[Dfree][i] << ","
            <<"\n";
    }
    file.close();

    qDebug() << "[Writie Data]:\t" << i << "pieces of data are saved";
    return savePath + "/FretThaData.csv";
}

/**
 * @brief FretThaSolver::expandGrayData 向数组中记录一条灰度数据
 * @param Iaa AA 通道灰度值
 * @param Ida DA 通道灰度值
 * @param Idd DD 通道灰度值
 */
void FretThaSolver::expandGrayData(double Iaa, double Ida, double Idd)
{
    vecGrayData[AA].push_back(Iaa);
    vecGrayData[DA].push_back(Ida);
    vecGrayData[DD].push_back(Idd);
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
        double Iaa = vecGrayData[AA][i];
        double Ida = vecGrayData[DA][i];
        double Idd = vecGrayData[DD][i];
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
    setThreshRatio(AA, threshAA);
    setThreshRatio(DA, threshDA);
    setThreshRatio(DD, threshDD);
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
        matFiltered[i] = FretImageProcessor::meanFilter16U(matSrc[i], 5);
        grayBackground[i] = FretImageProcessor::calcBackgroundGray(matFiltered[i]);
    }

    // 根据SBR阈值进行滤波
    double threshRatio[3] = {3.0, 3.0, 3.0};
    Mat maskSingleChannel[3];
    Mat maskSbr, maskSbr8U;    // 0-1
    for (int i = 0; i < 3; ++ i) {
        threshold(matFiltered[i], maskSingleChannel[i], grayBackground[i] * threshRatio[i], 255, THRESH_BINARY);
        maskSingleChannel[i].convertTo(maskSingleChannel[i], CV_8U);
    }
    bitwise_and(maskSingleChannel[AA], maskSingleChannel[DA], maskSbr8U);
    bitwise_and(maskSbr8U, maskSingleChannel[DD], maskSbr8U);
    maskSbr = FretImageProcessor::normalizeByMinMax(maskSbr8U);
    // 按照我们设计的算法生成评价结果
    Mat maskStd;
    Mat scoreStd(matSrc[AA].size(), CV_64FC1, Scalar(0.0));
    for (int i = 0; i < 3; ++ i) {
        Mat matStd = FretImageProcessor::localStandardDeviation(matFiltered[i], 5);
        Mat matStd8U = FretImageProcessor::normalizeByMinMax8U(matStd);
        Mat matStdMinima = FretImageProcessor::findLocalMinima(matStd8U);
        scoreStd = scoreStd + FretImageProcessor::normalizeByMinMax(matStdMinima);
    }

    Mat scoreStd16U = FretImageProcessor::normalizeByMinMax16U(scoreStd);
    cv::imwrite((viewFolderPath + "/StdMinima.tif").toStdString(), maskStd);
    cv::threshold(scoreStd16U, maskStd, 0, 255, cv::THRESH_BINARY | cv::THRESH_OTSU);
    maskStd.convertTo(maskStd, CV_8U);
    cv::imwrite((viewFolderPath + "/StdMinimaMaskOtsu.tif").toStdString(), maskStd);

    // 合成最后的结果
    Mat maskUlt, scoreMat;
    scoreMat = maskSbr.mul(scoreStd);
    bitwise_and(maskSbr8U, maskStd, maskUlt);
    cv::imwrite((viewFolderPath + "/FinalJudge.tif").toStdString(), maskUlt);

    // 计算质心
    std::vector<Point2D> centroids = getConnectedComponentCentroids(maskUlt);

    for (const auto& centroid : centroids) {
        int x = centroid.x;
        int y = centroid.y;
        double valGrayAA = matFiltered[AA].at<ushort>(y, x) - grayBackground[AA];
        double valGrayDA = matFiltered[DA].at<ushort>(y, x) - grayBackground[DA];
        double valGrayDD = matFiltered[DD].at<ushort>(y, x) - grayBackground[DD];

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

// 执行拟合计算
void FretThaSolver::performTwoHybridMatlab()  // matlab求解规划()
{
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
        qDebug() << "fitNonlinear failed:\t size not match";
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
    TwoHybridSolver(arrayAest, arrayDest, arrayEacorr, arrayEdcorr, resultsTha, &loss);
    TwoHybridSolver_terminate();

    // 按照求得的结果生成Afree和Dfree
    for (uint i = 0; i < n; ++ i) {
        double vadoublefree = ( ((vecFretData[Dest][i] - resultsTha[KDEFF] - vecFretData[Aest][i] * resultsTha[ND_NA]))
                               + sqrt(pow((vecFretData[Dest][i] - resultsTha[KDEFF] - vecFretData[Aest][i] * resultsTha[ND_NA]), 2) +
                                      4 * resultsTha[KDEFF] * vecFretData[Dest][i])
                               ) / 2;
        double valAfree = vecFretData[Aest][i] - (vecFretData[Dest][i] - vadoublefree) / resultsTha[ND_NA];
        vecFretData[Afree].push_back(valAfree);
        vecFretData[Dfree].push_back(vadoublefree);
    }
    qDebug() << "****************Fret Two Hybrid Matlab Solution****************";
    qDebug() << "Ed,max:\t" << resultsTha[EDMAX]
             << "Ea,max:\t" << resultsTha[EAMAX]
             << "Kd,EFF:\t" << resultsTha[KDEFF]
             << "Nd/Na:\t" << resultsTha[ND_NA]
             << "loss:\t" << loss;
    qDebug() << "Solve succeeded";
}

void FretThaSolver::performTwoHybridMatlabBin(double min, double max, double interval)  // matlab求解规划()
{
    // 执行数据封箱
    binData(min, max, interval);

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
    }
    else {
        qDebug() << "fitNonlinear failed:\t size not match";
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
    TwoHybridSolver(arrayAest, arrayDest, arrayEacorr, arrayEdcorr, resultsThaBin, &loss);
    TwoHybridSolver_terminate();

    // 按照求得的结果生成Afree和Dfree
    for (uint i = 0; i < n; ++ i) {
        double valDfree = ( ((vecFretDataBin[Dest][i] - resultsThaBin[KDEFF] - vecFretDataBin[Aest][i] * resultsThaBin[ND_NA]))
                           + sqrt(pow((vecFretDataBin[Dest][i] - resultsThaBin[KDEFF] - vecFretDataBin[Aest][i] * resultsThaBin[ND_NA]), 2) +
                                  4 * resultsThaBin[KDEFF] * vecFretDataBin[Dest][i])
                           ) / 2;
        double valAfree = vecFretDataBin[Aest][i] - (vecFretDataBin[Dest][i] - valDfree) / resultsThaBin[ND_NA];
        vecFretDataBin[Afree].push_back(valAfree);
        vecFretDataBin[Dfree].push_back(valDfree);
    }
    qDebug() << "****************Fret Two Hybrid Nonlinear Solution****************";
    qDebug() << "Ed,max:\t" << resultsThaBin[EDMAX]
             << "Ea,max:\t" << resultsThaBin[EAMAX]
             << "Kd,EFF:\t" << resultsThaBin[KDEFF]
             << "Nd/Na:\t" << resultsThaBin[ND_NA]
             << "loss:\t" << loss;
    qDebug() << "fitNonlinear succeeded:\t loss =" << loss;
}

void FretThaSolver::performTwoHybridLinear()
{
    // 经验值范围应是0-0.5，3-5
    resultsEdRad[EAMAX] = calcSlope(0, 1, vecFretData[Rad], vecFretData[Ed]);
    resultsEdRad[EDMAX] = calcApproach(1, 5, vecFretData[Rad], vecFretData[Ed]);
    resultsEdRad[ND_NA] = resultsEdRad[EAMAX] / resultsEdRad[EDMAX];

    resultsEaRda[EDMAX] = calcSlope(0, 1, vecFretData[Rda], vecFretData[Ea]);
    resultsEaRda[EAMAX] = calcApproach(1, 5, vecFretData[Rda], vecFretData[Ea]);
    resultsEaRda[ND_NA] = resultsEaRda[EAMAX] / resultsEaRda[EDMAX];


    // Debug
    qDebug() << "****************Fret Two Hybrid Linear Solution****************";
    qDebug() << "[Donor View Results]";
    qDebug() << "Ed,max:\t" << resultsEdRad[EDMAX] << "\t"
             << "Ea,max:\t" << resultsEdRad[EAMAX] << "\t"
             << "Nd/Na:\t" << resultsEdRad[ND_NA];
    qDebug() << "[Accept View Results]";
    qDebug() << "Ed,max:\t" << resultsEaRda[EDMAX] << "\t"
             << "Ea,max:\t" << resultsEaRda[EAMAX] << "\t"
             << "Nd/Na:\t" << resultsEaRda[ND_NA];
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

int FretThaSolver::binData(double min, double max, double interval)
{
    using namespace std;
    qDebug() << "*********************数据合并预处理Data*******************";
    qDebug() << (double)min
             << "-" << (double)max
             << "，" << (double)interval;
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

void FretThaSolver::mainActivity()
{
    void clearGrayData();
    void clearFretData();
    void clearFretDataBin();
    processOneBatch(batchFolderPath);
    // 计算得到FRET数据
    processGrayToFretData();
    // 进行计算
    performTwoHybridMatlab();
    performTwoHybridLinear();
    performTwoHybridMatlabBin(0, 5, 0.1);
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
double FretThaSolver::leastSquare(const std::vector<std::pair<double, double>>& points)
{
    double sumXY = 0, sumXX = 0;
    for (uint i = 0; i < points.size(); ++ i) {
        sumXY += points[i].first * points[i].second;
        sumXX += points[i].first * points[i].first;
    }
    return sumXY / sumXX;
}
void FretThaSolver::run()
{
    mainActivity();
    emit thaFinished();
}
