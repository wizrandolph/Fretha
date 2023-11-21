#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent) , ui(new Ui::MainWindow) {

    for (int i = 0; i < 6; ++ i) {
        m_chartView[i] = new QChartView(this);
    }

    ui->setupUi(this);
    ui->stackedWidget->setCurrentIndex(0);

    m_setting = new QSettings(QApplication::applicationDirPath() + "/config.ini", QSettings::IniFormat);

    loadLastRatio();

    // 连接信号和槽
    connectSignalsAndSlots();

    // 初始化界面
    initUICharts();  // 将类成员中的 QChartView 初始化后加入到界面中
    initUI();

    // 【用于测试】
    // 用来测试, 省去输入路径的步骤
    // ui->lineEditPath->setText("D:\\Files\\FRET\\Data\\MyData\\THA\\Template\\lite");
}

static bool checkChannelImages(const QString& folderPath)
{
    QDir dir(folderPath);

    if (!dir.exists())
    {
        qDebug() << "Invalid folder path.";
        return false;
    }

    bool hasDA = false;
    bool hasDD = false;
    bool hasAA = false;

    QStringList files = dir.entryList(QDir::Files);
    for (const QString& file : files)
    {
        if (file.endsWith("DA.tif"))
        {
            hasDA = true;
        }
        else if (file.endsWith("DD.tif"))
        {
            hasDD = true;
        }
        else if (file.endsWith("AA.tif"))
        {
            hasAA = true;
        }
    }

    return hasDA && hasDD && hasAA;
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButtonBrowse_clicked()
{
    QString dirStr = QFileDialog::getExistingDirectory(this, tr("选择目录"), "/home", QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
    // 选中目录后，将路径显示到编辑框中
    ui->lineEditPath->setText(dirStr);
}

void MainWindow::on_pushButtonAuto_clicked()
{
    QString strPath = ui->lineEditPath->text();

    QDir dir(strPath);
    if(!dir.exists() || strPath.isEmpty()) {
        showAlertDialog("当前路径不存在！");
        return;
    }

    m_savePath = strPath + "\\THA_Result";
    fretThaSolver->setBatchPath(strPath);
    showProgressDialog("请等待数据处理...");
    connect(fretThaSolver, SIGNAL(thaFinished()), this, SLOT(solveFinished()));
    fretThaSolver->start();
    m_dialog->exec();


    //solver->setBatchPath(strPath);
    //connect(solver, SIGNAL(finish()), this, SLOT(ready()));
    //solver->start();

}

void MainWindow::solveFinished()
{
    ui->pushButtonBin->setChecked(false);
    ui->pushButtonEdRc->setChecked(false);
    ui->pushButtonTHA->setChecked(true);
    ui->stackedWidget->setCurrentIndex(1);
    ui->stackedWidgetResult->setCurrentIndex(0);

    updateThaCharts();
    updateThaBinCharts();
    updateEdRcCharts();
    updateResultInterface();

    m_dialog->close();
}
void MainWindow::updateStatusBar(QRectF rect)
{
    int x = rect.x();
    int y = rect.y();
    int w = rect.width();
    int h = rect.height();
    fretImageProcessor->setRoi(x, y, w, h);

    double foreGround[3];
    double backGround[3];

    /*获取灰度值*/
    for (int i = 0; i < 3; ++ i)
    {
        backGround[i] = fretImageProcessor->getBackgroundGray((ChannelName)i);
        foreGround[i] = fretImageProcessor->getRoiGrayValue((ChannelName)i);
    }

    fretCalculator->resetData();
    fretCalculator->loadBackGray(backGround[0], backGround[1], backGround[2]);
    fretCalculator->loadForeGray(foreGround[0], foreGround[1], foreGround[2]);

    fretCalculator->correctData();

    // 在执行计算后更新当前的进度是否正常
    fretCalculator->calcEFret();
    fretCalculator->calc3CubeFret();
    fretCalculator->calc2Hybrid();

    double valEd = fretCalculator->getResult(Ed);
    double valRad = fretCalculator->getResult(Rad);
    double valFc = fretCalculator->getFcValue();

    ui->lineEditStatusGrayAA->setText(QString::number(foreGround[AA]));
    ui->lineEditStatusGrayDA->setText(QString::number(foreGround[DA]));
    ui->lineEditStatusGrayDD->setText(QString::number(foreGround[DD]));
    ui->lineEditBackGroundAA->setText(QString::number(backGround[AA]));
    ui->lineEditBackGroundDA->setText(QString::number(backGround[DA]));
    ui->lineEditBackGroundDD->setText(QString::number(backGround[DD]));
    ui->lineEditStatusFretEd->setText(QString::number(valEd));
    ui->lineEditStatusFretFc->setText(QString::number(valFc));
    ui->lineEditStatusFretRc->setText(QString::number(valRad));
    ui->lineEditStatusSbrAA->setText(QString::number(foreGround[AA] / backGround[AA]));
    ui->lineEditStatusSbrDA->setText(QString::number(foreGround[DA] / backGround[DA]));
    ui->lineEditStatusSbrDD->setText(QString::number(foreGround[DD] / backGround[DD]));
}

void MainWindow::on_pushButtonSave_clicked()
{
    // 使用弹窗获取保存文件夹路径
    QString folderPath = getSaveFolderPath(this);
    QString filePath = fretThaSolver->outputData(folderPath);   // 保存后台数据
    m_savePath = folderPath;

    // 保存图表文件
    QString filename[6] = {"Ed-Rad图.png",
                           "Ea-Rda图.png",
                           "Ed-Afree图.png",
                           "Ea-Dfree图.png",
                           "Ed-Afree图Bin.png",
                           "Ea-Dfree图Bin.png"};
    for (int i = 0; i < 6; i ++)
    {
        QPixmap pixmap = m_chartView[i]->grab();
        QImage image = pixmap.toImage();
        image.save(folderPath + "\\" + filename[i]);
    }

    // 弹出保存完成提示框
    showAlertDialog("数据已保存到：\n" + filePath);
}

void MainWindow::on_pushButtonTHA_clicked()
{
    ui->pushButtonTHA->setChecked(true);
    ui->pushButtonBin->setChecked(false);
    ui->pushButtonEdRc->setChecked(false);
    ui->stackedWidgetResult->setCurrentIndex(0);
}

void MainWindow::on_pushButtonBin_clicked()
{
    ui->pushButtonEdRc->setChecked(false);
    ui->pushButtonTHA->setChecked(false);
    ui->pushButtonBin->setChecked(true);
    ui->stackedWidgetResult->setCurrentIndex(1);
}

void MainWindow::on_pushButtonEdRc_clicked()
{
    ui->pushButtonTHA->setChecked(false);
    ui->pushButtonBin->setChecked(false);
    ui->pushButtonEdRc->setChecked(true);
    ui->stackedWidgetResult->setCurrentIndex(2);
}

void MainWindow::on_pushButtonBinExec_clicked()
{
    double min = ui->doubleSpinBoxMin->value(), max = ui->doubleSpinBoxMax->value(), itv = ui->doubleSpinBoxItv->value();
    fretThaSolver->performTwoHybridMatlabBin(min, max, itv);
    updateThaBinCharts();
}

void MainWindow::on_pushButtonHome_clicked()
{
    ui->stackedWidget->setCurrentIndex(0);
}



void findAllButtons(QLayout *layout, QList<QPushButton *> &buttonList)
{
    for (int i = 0; i < layout->count(); ++i)
    {
        QLayoutItem *item = layout->itemAt(i);
        if (QLayout *childLayout = item->layout())
        {
            findAllButtons(childLayout, buttonList);  // 子布局中继续查找
        }
        else if (QPushButton *button = qobject_cast<QPushButton *>(item->widget()))
        {
            buttonList.append(button);  // 找到一个 QPushButton
        }
    }
}

void MainWindow::on_pushButtonRatio_clicked()
{
    ui->lineEditA->setText(QString::number(ui->lineEditA->text().toDouble()));
    ui->lineEditB->setText(QString::number(ui->lineEditB->text().toDouble()));
    ui->lineEditC->setText(QString::number(ui->lineEditC->text().toDouble()));
    ui->lineEditD->setText(QString::number(ui->lineEditD->text().toDouble()));
    ui->lineEditG->setText(QString::number(ui->lineEditG->text().toDouble()));
    ui->lineEditK->setText(QString::number(ui->lineEditK->text().toDouble()));
    ui->lineEditY->setText(QString::number(ui->lineEditY->text().toDouble()));
    ui->lineEditExpsAA->setText(QString::number(ui->lineEditExpsAA->text().toDouble()));
    ui->lineEditExpsDA->setText(QString::number(ui->lineEditExpsDA->text().toDouble()));
    ui->lineEditExpsDD->setText(QString::number(ui->lineEditExpsDD->text().toDouble()));

    //将十个控件放进数组，方便循环
    QLineEdit *lineEdits[10] = {ui->lineEditA, ui->lineEditB, ui->lineEditC, ui->lineEditD, ui->lineEditG, ui->lineEditK, ui->lineEditY, ui->lineEditExpsAA, ui->lineEditExpsDA, ui->lineEditExpsDD};
    //创建对应的十个控件的默认值数组
    double defaultValues[10] = {0.186046857, 0.002056633, 0.001193742, 0.781544552, 4.6658, 0.645755859, 0.061747, 100, 100, 100};
    //创建对应的十个系数的名称数组
    QString names[10] = {"Ratio/A", "Ratio/B", "Ratio/C", "Ratio/D", "Ratio/G", "Ratio/K", "Ratio/Y", "ExposureTime/AA", "ExposureTime/DA", "ExposureTime/DD"};

    //循环判断是否有控件的值为空或者不是数字，如果是则将其值设置为默认值，如果不是则将其值设置为用户输入的值，并写入到配置文件中
    for (int i = 0; i < 10; ++ i)
    {
        if (lineEdits[i]->text().isEmpty() || !lineEdits[i]->text().toDouble())
        {
            lineEdits[i]->setText(QString::number(defaultValues[i]));
        }
        else
        {
            defaultValues[i] = lineEdits[i]->text().toDouble();
            m_setting->setValue(names[i], defaultValues[i]);
        }

    }
    loadLastRatio();

    ui->stackedWidget->setCurrentIndex(0);
}

void MainWindow::on_pushButtonSetRatio_clicked()
{
    ui->stackedWidget->setCurrentIndex(2);
    loadLastRatio();
}

void MainWindow::on_pushButtonManu_clicked()
{

    // 检查目录是否正确
    m_globalPath = ui->lineEditPath->text();
    qDebug() << m_globalPath;
    if (!checkPathLegal(m_globalPath)) return;

    initViewTableModel(m_globalPath);
    initRecordTableModel();
    updateStatusBar(ui->graphicsView->getRect());

//    // 设置快捷键
//    setupShortcuts();

    ui->stackedWidget->setCurrentIndex(3);    // 页面跳转
    showMaximized();    // 全屏
}

void MainWindow::loadLastRatio()
{
    // 内嵌在程序的默认参数
    double ratio[7] = {0.186046857, 0.002056633, 0.001193742, 0.781544552, 4.6658, 0.645755859, 0.061747};
    QString ratioName[7] = {"A", "B", "C", "D", "G", "K", "Y"};
    double exps[3] = {100, 100, 100};
    QString expsName[3] = {"AA", "DA", "DD"};
    double thre[3] = {3.0, 3.0, 3.0};
    QString threName[3] = {"threAA", "threDA", "threDD"};

    for (int i = 0; i < 7; ++ i)
    {
        // 读取配置，如果存在那么更新ratio数组，如果不存在那么按照ratio写入配置文件
        if (m_setting->contains("/Ratio/" + ratioName[i]))
            ratio[i] = m_setting->value("/Ratio/" + ratioName[i]).toDouble();
        else
            m_setting->setValue("/Ratio/" + ratioName[i], ratio[i]);
    }

    // 将ratio数组写入ui
    ui->lineEditA->setText(QString::number(ratio[0]));
    ui->lineEditB->setText(QString::number(ratio[1]));
    ui->lineEditC->setText(QString::number(ratio[2]));
    ui->lineEditD->setText(QString::number(ratio[3]));
    ui->lineEditG->setText(QString::number(ratio[4]));
    ui->lineEditK->setText(QString::number(ratio[5]));
    ui->lineEditY->setText(QString::number(ratio[6]));

    for (int i = 0; i < 3; ++ i)
    {
        // 读取配置，如果存在那么更新ratio数组，如果不存在那么按照ratio写入配置文件
        if (m_setting->contains("/ExposureTime/" + expsName[i]))
            exps[i] = m_setting->value("/ExposureTime/" + expsName[i]).toDouble();
        else
            m_setting->setValue("/ExposureTime/" + expsName[i], exps[i]);
    }

    ui->lineEditExpsAA->setText(QString::number(exps[0]));
    ui->lineEditExpsDA->setText(QString::number(exps[1]));
    ui->lineEditExpsDD->setText(QString::number(exps[2]));

    for (int i = 0; i < 3; ++ i)
    {
        // 读取配置，如果存在那么更新ratio数组，如果不存在那么按照ratio写入配置文件
        if (m_setting->contains("/AutoParams/" + threName[i]))
            thre[i] = m_setting->value("/AutoParams/" + threName[i]).toDouble();
        else
            m_setting->setValue("/AutoParams/" + threName[i], thre[i]);
    }

    /*
     * 添加数据到控件 */

    // 同步FRET计算器中的参数设置
    fretCalculator->setExposureTimes(exps[0], exps[1], exps[2]);
    fretCalculator->setRatios(ratio[0],
                              ratio[1],
                              ratio[2],
                              ratio[3],
                              ratio[4],
                              ratio[5],
                              ratio[6]);

    // 同步FRET图像处理器中的参数设置
    fretImageProcessor->setExposureTimes(exps[0], exps[1], exps[2]);
    fretImageProcessor->setRatios(ratio[0],
                                  ratio[1],
                                  ratio[2],
                                  ratio[3],
                                  ratio[4],
                                  ratio[5],
                                  ratio[6]);
    // 同步双杂交求解器中的参数设置
    fretThaSolver->setExposureTimes(exps[0], exps[1], exps[2]);
    fretThaSolver->setRatios( ratio[0],
                              ratio[1],
                              ratio[2],
                              ratio[3],
                              ratio[4],
                              ratio[5],
                              ratio[6]);
    fretThaSolver->setThreshRatios(thre[0], thre[1], thre[2]);

    fretCalculator->showParams();

}

// 一个简单的由自定义的channel_name映射到QString的函数
QString getChannelNameString(ShowType type)
{
    QString channelName;
    switch (type) {
        case AANORM: {
            channelName = "AA.tif";
            break;
        }
        case DANORM: {
            channelName = "DA.tif";
            break;
        }
        case DDNORM: {
            channelName = "DD.tif";
            break;
        }
        default: {
            break;
        }
    }
    return channelName;
}

void MainWindow::updateGraphicsView()
{
    m_image2Show = getImageToShow();
    ui->graphicsView->setImage(m_image2Show);
}

bool MainWindow::checkPathLegal(QString path)
{
    QDir dir(path);
    if(!dir.exists() || path.isEmpty())
    {
        showAlertDialog("当前路径不存在！");
        return false;
    }
    return true;
}

void MainWindow::on_pushButtonAdd_clicked()
{
    QStandardItemModel *model = dynamic_cast<QStandardItemModel*>(ui->tableRecord->model());
    if (model == nullptr) {
        qDebug() << "[Add data]:\tModel Not Found";
        return;
    }

    // 获取ROI
    QRectF rectf = ui->graphicsView->getRect();
    ui->graphicsView->addItem(rectf);   // 记录到图形场景中
    double x, y, w, h;
    x = rectf.x();
    y = rectf.y();
    w = rectf.width();
    h = rectf.height();
    int cnt = model->rowCount();
    qDebug() << "[Add Data]:\t" << x << y << w << h;

    double foreGround[3];

    // 获取灰度值
    for (int i = 0; i < 3; ++ i)
    {
        fretImageProcessor->setRoi(x, y, w, h);
        foreGround[i] = fretImageProcessor->getRoiGrayValue((ChannelName)i);
    }
    // 用来判断数据是否正常

    double valEd = fretCalculator->getResult(Ed);
    double valRad = fretCalculator->getResult(Rad);
    double valEa = fretCalculator->getResult(Ea);
    double valRda = fretCalculator->getResult(Rda);
    double valAest = fretCalculator->getResult(Aest);
    double valDest = fretCalculator->getResult(Dest);

    int columnIndex = 0;

    // 如果能够完成数据的计算，那么会向结果列表中添加当前数据
    QModelIndex index = model->index(cnt, 0);

    model->setItem(cnt, columnIndex ++, new QStandardItem(QString::number(foreGround[AA])));
    model->setItem(cnt, columnIndex ++, new QStandardItem(QString::number(foreGround[DA])));
    model->setItem(cnt, columnIndex ++, new QStandardItem(QString::number(foreGround[DD])));
    model->setItem(cnt, columnIndex ++, new QStandardItem(QString::number(valEd)));
    model->setItem(cnt, columnIndex ++, new QStandardItem(QString::number(valRad)));
    model->setItem(cnt, columnIndex ++, new QStandardItem(QString::number(valEa)));
    model->setItem(cnt, columnIndex ++, new QStandardItem(QString::number(valRda)));
    model->setItem(cnt, columnIndex ++, new QStandardItem(QString::number(valAest)));
    model->setItem(cnt, columnIndex ++, new QStandardItem(QString::number(valDest)));
    model->setItem(cnt, columnIndex ++, new QStandardItem(QString::number(x)));
    model->setItem(cnt, columnIndex ++, new QStandardItem(QString::number(y)));
    model->setItem(cnt, columnIndex ++, new QStandardItem(QString::number(w)));
    model->setItem(cnt, columnIndex ++, new QStandardItem(QString::number(h)));
    model->setItem(cnt, columnIndex ++, new QStandardItem(currentViewName));
    model->setRowCount(cnt + 1);

    ui->tableRecord->setCurrentIndex(index);


    // 设置表格格式
    QHeaderView *headerR = ui->tableRecord->horizontalHeader(); // 获取水平表头对象
    // 将表头的调整模式设置为ResizeToContents
    headerR->setSectionResizeMode(QHeaderView::ResizeToContents);
    // 自动调整所有列的宽度以适应内容
    ui->tableRecord->resizeColumnsToContents();
    ui->tableRecord->scrollToBottom();
}

void MainWindow::exportToCSV(const QTableView* tableView, const QString& filePath)
{
    QFile file(filePath);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        return;
    }

    QTextStream stream(&file);

    // 写入表头
    QStringList headerList;
    for (int section = 0; section < tableView->horizontalHeader()->count(); ++section) {
        QString headerText = tableView->model()->headerData(section, Qt::Horizontal, Qt::DisplayRole).toString();
        headerList.append(headerText);
    }
    stream << " , " << headerList.join(",") << "\n";

    // 写入数据
    for (int row = 0; row < tableView->model()->rowCount(); ++row) {
        QStringList dataList;

        // 写入行序号
        dataList.append(QString::number(row + 1));

        // 写入每一行的数据
        for (int column = 0; column < tableView->model()->columnCount(); ++column) {
            QString cellText = tableView->model()->data(tableView->model()->index(row, column)).toString();
            dataList.append(cellText);
        }

        stream << dataList.join(",") << "\n";
    }

    file.close();
}

void MainWindow::on_pushButtonExport_clicked()
{
    QString filePath = getSaveCsvFilePath(this);
    if (!filePath.isEmpty()) {
        exportToCSV(ui->tableRecord, filePath);
    }
}

QString MainWindow::getSaveCsvFilePath(QWidget* parent)
{
    QString filePath;

    // 弹出保存路径对话框
    filePath = QFileDialog::getSaveFileName(parent, "保存文件", QDir::homePath(), "CSV 文件 (*.csv)");

    return filePath;
}

QString MainWindow::getSaveFolderPath(QWidget* parent)
{
    QString folderPath;

    // 弹出选择文件夹对话框
    folderPath = QFileDialog::getExistingDirectory(parent, "选择文件夹", QDir::homePath());

    return folderPath;
}


QString MainWindow::getOpenCsvFilePath()
{
    QString filePath = QFileDialog::getOpenFileName(nullptr, "Select CSV file", "", "CSV files (*.csv)");
    return filePath;
}

void MainWindow::initRecordTableModel()
{
    // 设置表头
    int columnCount = 14;
    QStringList tableHeader = {"IAA",
                               "IDA",
                               "IDD",
                               "Ed",
                               "Rc",
                               "Ea",
                               "1/Rc",
                               "Aest",
                               "Dest",
                               "Left",
                               "Top",
                               "Width",
                               "Height",
                               "View"};
    if (columnCount != tableHeader.size()) {
        qDebug() << "表头数量不匹配";
        return;
    }
    QStandardItemModel *modelR = new QStandardItemModel(this);
    modelR->clear();
    modelR->setColumnCount(columnCount);
    ui->tableRecord->setModel(modelR);

    for (int i = 0; i < columnCount; ++ i) {
        modelR->setHeaderData(i, Qt::Horizontal, tableHeader[i]);
    }

    // 设置表格的点击行为和其他特性
    ui->tableRecord->setEditTriggers(QAbstractItemView::NoEditTriggers);    // 设置表格只能点击不能进入编辑状态
    ui->tableRecord->setSelectionBehavior(QAbstractItemView::SelectRows);    // 设置表格只能选择一整行
    ui->tableRecord->verticalHeader()->setHidden(true); // 隐藏默认序号
    ui->tableRecord->resizeColumnsToContents(); // 自适应列宽
}

void MainWindow::initViewTableModel(QString batchPath)
{
    QDir dir(batchPath);

    // 【显示视野列表】
    // 创建一个标准项模型
    QStandardItemModel *model = new QStandardItemModel;
    model->setColumnCount(3);
    ui->tableView->setModel(model);
    // 设置表头
    model->setHeaderData(0, Qt::Horizontal, QString("目录名"));
    model->setHeaderData(1, Qt::Horizontal, QString("目录类型"));
    model->setHeaderData(2, Qt::Horizontal, QString("状态"));
    // 获取子文件夹，并将其显示到界面中
    QStringList subfolders = dir.entryList(QDir::Dirs | QDir::NoDotAndDotDot);
    int cnt = 0;
    bool flag = false;
    for (const QString& subfolder : subfolders)
    {
        // 添加数据
        QString completeness;
        if (checkChannelImages(m_globalPath + "/" + subfolder)) completeness = "FRET";
        else completeness = "Others";
        model->setItem(cnt, 0, new QStandardItem(subfolder));
        model->setItem(cnt, 1, new QStandardItem(completeness));
        model->setItem(cnt, 2, new QStandardItem(""));

        // 同时寻找第一个为FRET视野的子文件夹选中
        if (!flag && completeness == "FRET")
        {
            m_showType = MERGED;

            QModelIndex index = model->index(cnt, 0);
            changeView(index);
            flag = true;
        }
        cnt ++ ;
    }

    // 设置表格的点击行为

    ui->tableView->setEditTriggers(QAbstractItemView::NoEditTriggers);    // 设置表格只能点击不能进入编辑状态
    // 设置选择模式为单选模式，且只能选择整行
    ui->tableView->setSelectionMode(QAbstractItemView::SingleSelection);
    ui->tableView->setSelectionBehavior(QAbstractItemView::SelectRows);
    ui->tableView->verticalHeader()->setHidden(true);   // 隐藏默认序号列

    // 自动调整列宽
    QHeaderView *header = ui->tableView->horizontalHeader();    // 获得表头
    header->setSectionResizeMode(QHeaderView::ResizeToContents);    // 将表头的调整模式设置为ResizeToContents
    ui->tableView->resizeColumnsToContents();   // 自动调整所有列的宽度以适应内容
}

void MainWindow::showAlertDialog(QString info)
{
    m_dialog = new QDialog(this);
    // 创建一个 QLabel 对象
    QLabel* label = new QLabel(info, m_dialog);
    label->setStyleSheet("QLabel {"
                         "border: 2px solid #555555;"
                         "border-radius: 10px;"
                         "font-size: 20px;"
                         "font-family: MiSans;"
                         "color: #ffffff;"
                         "background-color: #333333;"
                         "}");
    label->setAlignment(Qt::AlignCenter);
    label->setWordWrap(true);
    // 创建一个 QVBoxLayout 布局管理器
    QVBoxLayout* layout = new QVBoxLayout();
    // 将 QLabel 和 QPushButton 添加到布局管理器中
    layout->addWidget(label);
    // 设置 QDialog 的布局管理器
    m_dialog->setLayout(layout);
    // 将 QDialog 设置为模态对话框并显示它
    m_dialog->setModal(true);
    m_dialog->setWindowFlags(Qt::Dialog | Qt::WindowTitleHint | Qt::WindowCloseButtonHint);
    m_dialog->setSizeGripEnabled(true);
    m_dialog->setFixedSize(800, 600);
    m_dialog->setStyleSheet("QDialog {background-color: #121212;}");
    m_dialog->exec();
}

void MainWindow::showProgressDialog(QString info)
{
    m_dialog = new QDialog(this);
    // 创建一个 QLabel 对象
    QLabel* label = new QLabel(info, m_dialog);
    label->setAlignment(Qt::AlignCenter);
    label->setWordWrap(true);
    label->setStyleSheet("QLabel {"
                            "border: 2px solid #555555;"
                            "border-radius: 10px;"
                            "font-size: 20px;"
                            "font-family: MiSans;"
                            "color: #ffffff;"
                            "background-color: #333333;"
                        "}");

    // 创建一个 QProgressBar 对象
    QProgressBar* progressBar = new QProgressBar(m_dialog);
    progressBar->setRange(0, 100); // 设置进度条的范围

    // 创建一个 QVBoxLayout 布局管理器
    QVBoxLayout* layout = new QVBoxLayout();
    // 将 QLabel 和 QProgressBar 添加到布局管理器中
    layout->addWidget(label);
    layout->addWidget(progressBar);

    // 设置 QDialog 的布局管理器
    m_dialog->setLayout(layout);
    // 将 QDialog 设置为模态对话框并显示它
    m_dialog->setModal(true);
    m_dialog->setWindowFlags(Qt::Dialog | Qt::WindowTitleHint | Qt::WindowCloseButtonHint);
    m_dialog->setSizeGripEnabled(true);
    m_dialog->setFixedSize(800, 600);
    m_dialog->setStyleSheet("QDialog {background-color: #121212;}");

    // 连接信号与槽
    connect(fretThaSolver, &FretThaSolver::progressChanged, progressBar, &QProgressBar::setValue);
}

void MainWindow::on_pushButtonCalc_clicked()
{
    // 清除数据
    fretThaSolver->clearFretData();
    fretThaSolver->clearFretDataBin();

    QStandardItemModel *model = dynamic_cast<QStandardItemModel*>(ui->tableRecord->model());
    if (model == nullptr) return;

    for (int i = 0; i < model->rowCount(); ++ i)
    {
        QModelIndex index = model->index(i, 3);
        QVariant data = model->data(index);
        double valEd = data.toDouble();
        index = model->index(i, 4);
        data = model->data(index);
        double valRad = data.toDouble();
        index = model->index(i, 5);
        data = model->data(index);
        double valEa= data.toDouble();
        index = model->index(i, 6);
        data = model->data(index);
        double valRda = data.toDouble();
        index = model->index(i, 7);
        data = model->data(index);
        double valAest = data.toDouble();
        index = model->index(i, 8);
        data = model->data(index);
        double valDest= data.toDouble();

        fretThaSolver->expandFretData(valEd, valEa, valRad, valRda, valAest, valDest);
    }

    m_savePath = m_globalPath + "\\THA_Result";

    fretThaSolver->performTwoHybridMatlab();
    fretThaSolver->performTwoHybridMatlabBin(0, 5, 0.1);
    fretThaSolver->performTwoHybridLinear();

    ui->pushButtonBin->setChecked(false);
    ui->pushButtonEdRc->setChecked(false);
    ui->pushButtonTHA->setChecked(true);

    updateThaCharts();
    updateThaBinCharts();
    updateEdRcCharts();
    updateResultInterface();

    ui->stackedWidget->setCurrentIndex(1);
    ui->stackedWidgetResult->setCurrentIndex(0);

}

void MainWindow::on_pushButtonBack_clicked()
{
    QStandardItemModel *model = dynamic_cast<QStandardItemModel*>(ui->tableRecord->model());
    if (model == nullptr) {
        qDebug() << "[Exit Manually Draw]:\tModel Not Found";
        return;
    } else {
        // 清空数据
        model->clear();
    }

    fretCalculator->resetData();

    // 页面跳转
    qDebug() << ui->stackedWidget->currentIndex();
    ui->stackedWidget->setCurrentIndex(0);
    qDebug() << ui->stackedWidget->currentIndex();

}

void MainWindow::on_pushButtonImportScreen_clicked()
{
    QStandardItemModel *model = dynamic_cast<QStandardItemModel*>(ui->tableRecord->model());
    if (model == nullptr) {
        qDebug() << "[Exit Manually Draw]:\tModel Not Found";
        return;
    }

    // 非测试
    QString csvPath = getOpenCsvFilePath();


    initRecordTableModel();
    int cnt = model->rowCount();

    QFile file(csvPath);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        return;

    QTextStream in(&file);
    bool isHeader = true;
    while (!in.atEnd()) {
        QString line = in.readLine();
        QStringList fields = line.split(",");
        if (fields.size() != 15) {
            return;
        }
        if (fields.size() >= 13 && !isHeader) {
            int columnIndex = 0;
            int i = 1;
            for (; i < fields.size(); ++ i)
            {
                model->setItem(cnt, columnIndex ++, new QStandardItem(fields.at(i)));
            }
            model->setRowCount(++ cnt);
        }
        isHeader = false;
    }
    // 设置表格只能点击不能进入编辑状态
    ui->tableRecord->setEditTriggers(QAbstractItemView::NoEditTriggers);
    // 设置只能选择整行
    ui->tableRecord->setSelectionBehavior(QAbstractItemView::SelectRows);
    // 自动调整所有列的宽度以适应内容
    ui->tableRecord->resizeColumnsToContents();
    file.close();
}

void MainWindow::on_pushButtonDelete_clicked()
{
    QModelIndexList selectedRows = ui->tableRecord->selectionModel()->selectedRows();

    if (selectedRows.size() > 0) {

        int row = selectedRows.at(0).row();

        QModelIndex viewIndex = ui->tableRecord->model()->index(row, 13);
        QString view = viewIndex.data().toString();

        ui->graphicsView->clearRectItemList();
        ui->tableRecord->model()->removeRow(row);
        updateRectRecorded(view);
        if (row - 1 >= 0) {
            QModelIndex index = ui->tableRecord->model()->index(row - 1, 0);
            ui->tableRecord->setCurrentIndex(index);
        }
    }
}

void MainWindow::comboBoxViewTypeChanged(int index)
{

    switch (index) {
        case 0: {
            m_showType = MERGED;
            break;
        }
        case 1: {
            m_showType = AANORM;
            break;
        }
        case 2: {
            m_showType = DANORM;
            break;
        }
        case 3: {
            m_showType = DDNORM;
            break;
        }
        case 4: {
            m_showType = EDPSEU;
            break;
        }
        case 5: {
            m_showType = RCPSEU;
            break;
        }
        default: {
            break;
        }
    }

    updateGraphicsView();

}

void MainWindow::comboBoxDrawModeChanged(int index)
{
    switch (index) {
        case 0: {
            ui->graphicsView->setDrawMode(CropMode);
            break;
        }
        case 1: {
            ui->graphicsView->setDrawMode(StampMode);
            break;
        }
        default: {
            break;
        }
    }
}

QImage MainWindow::getImageToShow()
{
    cv::Mat cvImage;
    QImage qtImage;

    switch (m_showType){
        case MERGED: {
            cvImage= fretImageProcessor->getMergedImage();
            break;
        }
        case EDPSEU: {
            fretImageProcessor->preProcessData();
            fretImageProcessor->calcEFret();
            cvImage= fretImageProcessor->getPseuImage(Ed);
            break;
        }
        case RCPSEU: {
            fretImageProcessor->preProcessData();
            fretImageProcessor->calcEFret();
            cvImage = fretImageProcessor->getPseuImage(Rad);
            break;
        }
        default: {
            ChannelName channelName = Enumerator::showTypeToChannelName(m_showType);
            if (channelName != NaC) {
                cvImage = fretImageProcessor->getNormalizedImage(channelName);
            }
            break;
        }
    }

    if (cvImage.empty()) {
        qDebug() << "当前用于显示的图片出现问腿";
        return m_image2Show;
    }

    QString tmpFilePath = m_tmpDirPath + "\\tmp.jpg";
    cv::imwrite(tmpFilePath.toStdString(), cvImage);
    qtImage = QImage(tmpFilePath);

    // 删除临时文件
    QFile file(tmpFilePath);
    QFile::remove(tmpFilePath);

    return qtImage;
}

/**
 * @brief MainWindow::changeView 切换视野
 * @param index 视野表格中的索引
 */
void MainWindow::changeView(QModelIndex &index)
{
    // 获取点击的单元格的索引、视野名
    int row = index.row();
    QModelIndex secondColumnIndex = ui->tableView->model()->index(row, 1);
    QString viewName = ui->tableView->model()->index(row, 0).data().toString();

    // 如果当前显示的是相同视野，那么不会产生操作
    if (viewName == m_currentViewPath) {
        return;
    }

    if (ui->tableView->model()->data(secondColumnIndex).toString() != "Unknown") {

        // 修改该行的第三列数据为"Yes"
        QModelIndex thirdColumnIndex = ui->tableView->model()->index(row, 2);
        ui->tableView->model()->setData(thirdColumnIndex, "√");
        // 将其他行修改为空
        for (int i = 0; i < ui->tableView->model()->rowCount(); i++) {
            if (i != row) {
                QModelIndex otherThirdColumnIndex = ui->tableView->model()->index(i, 2);
                ui->tableView->model()->setData(otherThirdColumnIndex, "");
            }
        }

        // 切换为视野的图片内容
        m_currentViewPath = m_globalPath + "\\" + viewName;
        currentViewName = viewName;
        qDebug() << "[Change View]:\t" << m_currentViewPath;
        fretImageProcessor->loadSourceData(m_currentViewPath);

        // 更新绘图控件、表格、状态栏
        updateGraphicsView();
        updateRectRecorded(viewName);
        updateStatusBar(ui->graphicsView->getRect());
    }
}

/**
 * @brief MainWindow::trackRect 跟踪定位ROI
 * @param index
 */
void MainWindow::trackRect(QModelIndex &index)
{
    // 获取点击的单元格中的数据
    int row = index.row();
    QModelIndex viewNameIndex = ui->tableRecord->model()->index(row, 13);
    QModelIndex leftIndex = ui->tableRecord->model()->index(row, 9);
    QModelIndex topIndex = ui->tableRecord->model()->index(row, 10);
    QModelIndex widthIndex = ui->tableRecord->model()->index(row, 11);
    QModelIndex heightIndex = ui->tableRecord->model()->index(row, 12);
    QString viewName = viewNameIndex.data().toString();    // 记录表格中的视野名

    int targetRow = -1;
    QAbstractItemModel *model = ui->tableView->model();

    // 找到视野表格中同视野名的索引
    for (int i = 0; i < model->rowCount(); ++ i) {
        QModelIndex index = model->index(i, 0);  // 第一列的索引为(行, 列)
        QVariant data = model->data(index);

        if (data.toString() == viewName) {
            targetRow = i;
            break;
        }
    }

    if (targetRow == -1) {
        showAlertDialog("[Track Data]:\tView Not Found");
        return;
    }
    ui->tableView->selectRow(targetRow);

    QModelIndex targetIndex = model->index(targetRow, 0); // 目标行索引
    changeView(targetIndex);
    updateRectRecorded(viewName);

    // 更新ROI到绘图界面
    double left = (leftIndex.data().toString()).toDouble();
    double top = (topIndex.data().toString()).toDouble();
    double width = (widthIndex.data().toString()).toDouble();
    double height = (heightIndex.data().toString()).toDouble();
    QRectF roi(left, top, width, height);
    qDebug() << "[Track Data]:\t" << roi;
    ui->graphicsView->setRect(roi);
    updateStatusBar(roi);

}

/**
 * @brief MainWindow::initUICharts 初始化界面的中的表格
 */
void MainWindow::initUICharts()
{
    // 初始化所有的图表
    initAllCharts();
    for (int i = 0; i < 6; i++) {
        assert(m_chartView[i] != nullptr);
    }
    // 向布局中添加图表
    ui->horizontalLayoutChartP1->insertWidget(0, m_chartView[EdRad]);
    ui->horizontalLayoutChartP1->insertWidget(1, m_chartView[EaRda]);
    ui->horizontalLayoutChartP2->insertWidget(0, m_chartView[EdAfree]);
    ui->horizontalLayoutChartP2->insertWidget(1, m_chartView[EaDfree]);
    ui->horizontalLayoutChartP3->insertWidget(0, m_chartView[EdAfreeBin]);
    ui->horizontalLayoutChartP3->insertWidget(1, m_chartView[EaDfreeBin]);
}

/**
 * @brief MainWindow::initUI 初始化界面
 */
void MainWindow::initUI()
{
    ui->pushButtonEdRc->setCheckable(true);
    ui->pushButtonTHA->setCheckable(true);
    ui->pushButtonBin->setCheckable(true);
}

void MainWindow::initAllCharts()
{
    ChartName chartNames[6] = {
        EdAfree, EaDfree, EdAfreeBin, EaDfreeBin, EdRad, EaRda
    };
    for (auto it : chartNames) {
        initializeChart(m_chartView[it], Enumerator::chartNameToQString(it));
    }
}

/**
 * @brief MainWindow::updateResultInterface 更新结果到界面
 */
void MainWindow::updateResultInterface()
{
    // 更新THA的结果
    ui->lineEditEamaxP2->setText(QString::number(fretThaSolver->resultsTha[EAMAX]));
    ui->lineEditEdmaxP2->setText(QString::number(fretThaSolver->resultsTha[EDMAX]));
    ui->lineEditNdNaP2->setText(QString::number(fretThaSolver->resultsTha[ND_NA]));
    ui->lineEditKdeffP2->setText(QString::number(fretThaSolver->resultsTha[KDEFF]));
    // 更新BIN的结果
    ui->lineEditEamaxP3->setText(QString::number(fretThaSolver->resultsThaBin[EAMAX]));
    ui->lineEditEdmaxP3->setText(QString::number(fretThaSolver->resultsThaBin[EDMAX]));
    ui->lineEditNdNaP3->setText(QString::number(fretThaSolver->resultsThaBin[ND_NA]));
    ui->lineEditKdeffP3->setText(QString::number(fretThaSolver->resultsThaBin[KDEFF]));
    // 更新EdRc结果
    ui->lineEditEamaxEdP1->setText(QString::number(fretThaSolver->resultsEdRad[EAMAX]));
    ui->lineEditEdmaxEdP1->setText(QString::number(fretThaSolver->resultsEdRad[EDMAX]));
    ui->lineEditNNEdP1->setText(QString::number(fretThaSolver->resultsEdRad[ND_NA]));
    ui->lineEditEamaxEaP1->setText(QString::number(fretThaSolver->resultsEaRda[EAMAX]));
    ui->lineEditEdmaxEaP1->setText(QString::number(fretThaSolver->resultsEaRda[EDMAX]));
    ui->lineEditNNEaP1->setText(QString::number(fretThaSolver->resultsEaRda[ND_NA]));
    ui->lineEditEamaxAvgP1->setText(QString::number((fretThaSolver->resultsEdRad[EAMAX] + fretThaSolver->resultsEaRda[EAMAX]) / 2));
    ui->lineEditEdmaxAvgP1->setText(QString::number((fretThaSolver->resultsEdRad[EDMAX] + fretThaSolver->resultsEaRda[EDMAX]) / 2));
    ui->lineEditNNAvgP1->setText(QString::number((fretThaSolver->resultsEdRad[ND_NA] + fretThaSolver->resultsEaRda[ND_NA]) / 2));
}

/**
 * @brief MainWindow::updateThaCharts 更新双杂交图
 */
void MainWindow::updateThaCharts()
{
    // 清空图标中的数据
    clearChart(m_chartView[EdAfree]);
    clearChart(m_chartView[EaDfree]);

    // 绘制散点图
    drawScatter(m_chartView[EdAfree], "Ed-Afree Scatter", fretThaSolver->vecFretData[Afree], fretThaSolver->vecFretData[Ed]);
    drawScatter(m_chartView[EaDfree], "Ea-Dfree Scatter", fretThaSolver->vecFretData[Dfree], fretThaSolver->vecFretData[Ea]);

    // 生成线条的点
    std::vector<double> vecAfree, vecDfree, vecEa, vecEd;
    double maxAfree = fretThaSolver->maxData(Afree),
        maxDfree = fretThaSolver->maxData(Dfree),
        maxEa = fretThaSolver->maxData(Ea),
        maxEd = fretThaSolver->maxData(Ed);

    // 生成模型上的点
    int intervals = 100;    // 设置曲线精细度，数值越大，线条越平滑
    for (int i = 0; i < intervals; i++) {
        double preAfree = maxAfree / intervals * i;
        double preEd = fretThaSolver->resultsTha[EDMAX] * preAfree / (preAfree + fretThaSolver->resultsTha[KDEFF] * fretThaSolver->resultsTha[ND_NA]);
        double preDfree = maxDfree / intervals * i;
        double preEa = fretThaSolver->resultsTha[EAMAX] * preDfree / (preDfree + fretThaSolver->resultsTha[KDEFF]);
        vecAfree.push_back(preAfree);
        vecEd.push_back(preEd);
        vecDfree.push_back(preDfree);
        vecEa.push_back(preEa);
    }

    // 绘制曲线
    drawLine(m_chartView[EdAfree], "Ed-Afree Curve", vecAfree, vecEd);
    drawLine(m_chartView[EaDfree], "Ea-Dfree Curve", vecDfree, vecEa);

    // 设置坐标轴范围
    setupChartAxis(m_chartView[EdAfree], "Afree", 0, maxAfree * 1.1, "Ed", 0, maxEd * 1.1);
    setupChartAxis(m_chartView[EaDfree], "Dfree", 0, maxDfree * 1.1, "Ea", 0, maxEa * 1.1);
}
/**
 * @brief MainWindow::updateThaBinCharts 更新数据分箱界面的图表
 */
void MainWindow::updateThaBinCharts()
{
    // 清空图标中的数据
    clearChart(m_chartView[EdAfreeBin]);
    clearChart(m_chartView[EaDfreeBin]);

    // 绘制散点图
    drawScatter(m_chartView[EdAfreeBin], "Ed-Afree Scatter", fretThaSolver->vecFretDataBin[Afree], fretThaSolver->vecFretDataBin[Ed]);
    drawScatter(m_chartView[EaDfreeBin], "Ea-Dfree Scatter", fretThaSolver->vecFretDataBin[Dfree], fretThaSolver->vecFretDataBin[Ea]);

    // 生成线条的点
    std::vector<double> vecAfree, vecDfree, vecEa, vecEd;
    double maxAfree = fretThaSolver->maxBinData(Afree),
        maxDfree = fretThaSolver->maxBinData(Dfree),
        maxEa = fretThaSolver->maxBinData(Ea),
        maxEd = fretThaSolver->maxBinData(Ed);

    // 生成模型上的点
    int intervals = 100;    // 设置曲线精细度，数值越大，线条越平滑
    for (int i = 0; i < intervals; i++) {
        double preAfree = maxAfree / intervals * i;
        double preEd = fretThaSolver->resultsThaBin[EDMAX] * preAfree / (preAfree + fretThaSolver->resultsThaBin[KDEFF] * fretThaSolver->resultsThaBin[ND_NA]);
        double preDfree = maxDfree / intervals * i;
        double preEa = fretThaSolver->resultsThaBin[EAMAX] * preDfree / (preDfree + fretThaSolver->resultsThaBin[KDEFF]);
        vecAfree.push_back(preAfree);
        vecEd.push_back(preEd);
        vecDfree.push_back(preDfree);
        vecEa.push_back(preEa);
    }

    // 绘制曲线
    drawLine(m_chartView[EdAfreeBin], "Ed-Afree Curve", vecAfree, vecEd);
    drawLine(m_chartView[EaDfreeBin], "Ea-Dfree Curve", vecDfree, vecEa);

    // 设置坐标轴范围
    setupChartAxis(m_chartView[EdAfreeBin], "Afree", 0, maxAfree * 1.1, "Ed", 0, maxEd * 1.1);
    setupChartAxis(m_chartView[EaDfreeBin], "Dfree", 0, maxDfree * 1.1, "Ea", 0, maxEa * 1.1);
}
/**
 * @brief MainWindow::updateEdRcCharts 更新Ed-RC图
 */
void MainWindow::updateEdRcCharts()
{
    // 清空图标中的数据
    clearChart(m_chartView[EdRad]);
    clearChart(m_chartView[EaRda]);

    // 绘制散点图
    drawScatter(m_chartView[EdRad], "Ed-Rc Scatter", fretThaSolver->vecFretData[Rad], fretThaSolver->vecFretData[Ed]);
    drawScatter(m_chartView[EaRda], "Ea-1/Rc Scatter", fretThaSolver->vecFretData[Rda], fretThaSolver->vecFretData[Ea]);

    // 生成线条的点
    std::vector<double> vecEdSlope, vecRadSlope, vecEdAppro, vecRadAprro,
        vecEaSlope, vecRdaSlope, vecEaAppro, vecRdaAprro;

    double maxRad = fretThaSolver->maxData(Rad),
        maxRda = fretThaSolver->maxData(Rda),
        maxEa = fretThaSolver->maxData(Ea),
        maxEd = fretThaSolver->maxData(Ed);

    // 生成模型上的点
    int intervals = 10;    // 设置曲线精细度，数值越大，线条越平滑
    for (int i = 0; i < intervals; ++ i) {
        // Ed Rc(Rad)
        double preRadSlope = 1.0 / intervals * i;
        double preEdSlope = preRadSlope * fretThaSolver->resultsEdRad[EAMAX];
        double preRadAppro = i * 3.0 / intervals + 1.5;
        double preEdAppro = fretThaSolver->resultsEdRad[EDMAX];
        vecEdSlope.push_back(preEdSlope);
        vecRadSlope.push_back(preRadSlope);
        vecEdAppro.push_back(preEdAppro);
        vecRadAprro.push_back(preRadAppro);

        // Ea 1/Rc(Rda)
        double preRdaSlope = 1.0 / intervals * i;
        double preEaSlope = preRdaSlope * fretThaSolver->resultsEaRda[EDMAX];
        double preRdaAppro = i * 3.0 / intervals + 1.5;
        double preEaAppro = fretThaSolver->resultsEaRda[EAMAX];
        vecEaSlope.push_back(preEaSlope);
        vecRdaSlope.push_back(preRdaSlope);
        vecEaAppro.push_back(preEaAppro);
        vecRdaAprro.push_back(preRdaAppro);
    }

    // 绘制曲线
    drawLine(m_chartView[EdRad], "Ed-Rc Approach", vecRadAprro, vecEdAppro);
    drawLine(m_chartView[EdRad], "Ed-Rc Slope", vecRadSlope, vecEdSlope);
    drawLine(m_chartView[EaRda], "Ea-1/Rc Approach", vecRdaAprro, vecEaAppro);
    drawLine(m_chartView[EaRda], "Ea-1/Rc Slope", vecRdaSlope, vecEaSlope);

    // 设置坐标轴
    setupChartAxis(m_chartView[EdRad], "Rc", 0, maxRad * 1.1, "Ed", 0, maxEd * 1.1);
    setupChartAxis(m_chartView[EaRda], "1/Rc", 0, maxRda * 1.1, "Ea", 0, maxEa * 1.1);
}
/**
 * @brief MainWindow::initializeChart 初始化图表
 * @param chartView ChartView指针
 * @param chartName 图表的名字
 */
void MainWindow::initializeChart(QChartView* chartView, QString chartName)
{
    // 创建图表并设置标题、背景等属性
    QChart* chart = new QChart();
    chart->setTitle(chartName);

    // 创建X轴和Y轴，并进行初始化设置
    QValueAxis* xAxis = new QValueAxis();
    QValueAxis* yAxis = new QValueAxis();
    xAxis->setTitleText("X Axis");
    yAxis->setTitleText("Y Axis");
    chart->addAxis(xAxis, Qt::AlignBottom);
    chart->addAxis(yAxis, Qt::AlignLeft);
    chart->setTheme(QtCharts::QChart::ChartThemeDark);
    chart->setMargins(QMargins(20, 20, 20, 20));    // 设置外边界全部为0
    chart->layout()->setContentsMargins(5, 5, 9, 9);    // 设置内边界全部为0

    // 设置图表视图的属性
    chartView->setChart(chart);
    chartView->setStyleSheet("background: transparent;");   // 使背景为透明
    chartView->setRenderHint(QPainter::Antialiasing); // 启用抗锯齿
}

/**
 * @brief MainWindow::clearChart 清除图表中的数据和坐标轴
 * @param chartView ChartView指针
 */
void MainWindow::clearChart(QChartView* chartView)
{
    // 获取图表
    QChart* chart = chartView->chart();

    // 移除所有的系列
    QList<QAbstractSeries*> seriesList = chart->series();
    foreach (QAbstractSeries* series, seriesList) {
        chart->removeSeries(series);
        delete series;
    }

    // 更新图表视图
    chartView->update();
}

/**
 * @brief MainWindow::drawScatter 向图表中添加一个散点系列
 * @param chartView 目标图表
 * @param seriesName 散点系列的名字
 * @param xData 横坐标数据
 * @param yData 纵坐标数据
 */
void MainWindow::drawScatter(QChartView* chartView, QString seriesName, const std::vector<double> xData, const std::vector<double> yData)
{
    // 获取图表
    QChart* chart = chartView->chart();

    // 创建散点系列
    QScatterSeries* series = new QScatterSeries();
    series->setName(seriesName);

    // 向散点系列中添加数据
    for (size_t i = 0; i < xData.size(); i++) {
        series->append(xData[i], yData[i]);
    }

    // 设置散点样式
    series->setMarkerShape(QScatterSeries::MarkerShapeRectangle);
    series->setMarkerSize(5);

    // 添加散点系列到图表
    chart->addSeries(series);

    // 更新图表视图的显示
    chartView->repaint(); // 标记视图需要重新绘制
    chartView->update(); // 更新视图

    QAbstractAxis *axisX = chart->axes(Qt::Horizontal).at(0);
    QAbstractAxis *axisY = chart->axes(Qt::Vertical).at(0);
    series->attachAxis(axisX);
    series->attachAxis(axisY);
}

/**
 * @brief MainWindow::drawScatter 向图表中添加一个线条系列
 * @param chartView 目标图表
 * @param seriesName 线条系列的名字
 * @param xData 横坐标数据
 * @param yData 纵坐标数据
 */
void MainWindow::drawLine(QChartView* chartView, QString seriesName, const std::vector<double> xData, const std::vector<double> yData)
{
    // 获取图表
    QChart* chart = chartView->chart();

    // 创建散点系列
    QLineSeries* series = new QLineSeries();
    series->setName(seriesName);

    // 向散点系列中添加数据
    for (size_t i = 0; i < xData.size(); i++) {
        series->append(xData[i], yData[i]);
        // qDebug() << xData[i] << yData[i];
    }

    // 添加散点系列到图表
    chart->addSeries(series);

    // 将线条系列与坐标轴关联
    QAbstractAxis *axisX = chart->axes(Qt::Horizontal).at(0);
    QAbstractAxis *axisY = chart->axes(Qt::Vertical).at(0);
    series->attachAxis(axisX);
    series->attachAxis(axisY);

    // 设置线条样式
    QPen pen = series->pen();
    pen.setWidth(5);
    series->setPen(pen);
}

/**
 * @brief MainWindow::setupChartAxis 设置坐标轴的相关参数
 * @param chartView ChartView指针
 * @param xAxisLabel X轴的名称
 * @param xAxisMin X轴的最小值
 * @param xAxisMax X轴的最大值
 * @param yAxisLabel Y轴的名称
 * @param yAxisMin Y轴的最小值
 * @param yAxisMax Y轴的最大值
 */
void MainWindow::setupChartAxis(QChartView *chartView, const QString &xAxisLabel, qreal xAxisMin, qreal xAxisMax, const QString &yAxisLabel, qreal yAxisMin, qreal yAxisMax)
{
    // 获取指向图表的指针
    QChart *chart = chartView->chart();
    QAbstractAxis *axisX = chart->axes(Qt::Horizontal).at(0);
    QAbstractAxis *axisY = chart->axes(Qt::Vertical).at(0);

    if (axisX && axisY) {
        // 设置坐标轴范围
        axisX->setTitleText(xAxisLabel);
        axisX->setRange(xAxisMin, xAxisMax);
        axisY->setTitleText(yAxisLabel);
        axisY->setRange(yAxisMin, yAxisMax);

        // 更新图表显示
        chart->update();
    }

    // 设置绘图区域的圆角
    QRectF plotArea = chart->plotArea();
    QPainterPath path;
    path.addRoundedRect(plotArea, 10, 10);
}

/**
 * @brief MainWindow::updateRecord 更新视图中已记录的矩形控件
 * @param viewPath 视野文件夹的名字
 */
void MainWindow::updateRectRecorded(QString viewPath)
{
    QAbstractItemModel *model = ui->tableRecord->model();
    if (model != nullptr) {
        for (int row = 0; row < model->rowCount(); ++ row) {
            QString viewName = model->index(row, 13).data().toString();
            if (viewName == viewPath) {

                QModelIndex leftIndex = model->index(row, 9);
                QModelIndex topIndex = model->index(row, 10);
                QModelIndex widthIndex = model->index(row, 11);
                QModelIndex heightIndex = model->index(row, 12);

                int left = (leftIndex.data().toString()).toDouble();
                int top = (topIndex.data().toString()).toDouble();
                int width = (widthIndex.data().toString()).toDouble();
                int height = (heightIndex.data().toString()).toDouble();

                ui->graphicsView->addItem(QRectF(left, top, width, height));
            }
        }
    }
}

/**
 * @brief MainWindow::setupShortcuts 设置快捷键，仅在圈点界面有效
 */
void MainWindow::setupShortcuts()
{
    // 添加数据按键
    WizShortcut *shortcutAdd = new WizShortcut(QKeySequence(Qt::Key_A), 3, ui->stackedWidget);
    QObject::connect(shortcutAdd, &WizShortcut::activated, this, [this, shortcutAdd]{
        if (m_currentPageIndex == shortcutAdd->getPageIndex()) {
            // 当前活动页面与按钮所在页面相同时执行操作
            ui->pushButtonAdd->click();
        }
    });

    // 删除数据按键
    WizShortcut *shortcutRemove = new WizShortcut(QKeySequence(Qt::Key_D), 3, ui->stackedWidget);
    QObject::connect(shortcutRemove, &WizShortcut::activated, this, [this, shortcutRemove]{
        if (m_currentPageIndex == shortcutRemove->getPageIndex()) {
            // 当前活动页面与按钮所在页面相同时执行操作
            ui->pushButtonDelete->click();
        }
    });
}

/**
 * @brief MainWindow::connectSignalsAndSlots 连接信号和槽
 */
void MainWindow::connectSignalsAndSlots()
{
    // 设置tableView的信号与槽函数
    connect(ui->tableView, &QTableView::clicked, [=](const QModelIndex& index) {
        QModelIndex newIndex = index;
        changeView(newIndex);
    });
    // 设置tableViewResult的信号与槽函数
    connect(ui->tableRecord, &QTableView::clicked, [=](const QModelIndex& index) {
        QModelIndex newIndex = index;
        trackRect(newIndex);
    });

    // 设置comboBoxView的信号与槽函数
    connect(ui->comboBoxViewType, SIGNAL(currentIndexChanged(int)), this, SLOT(comboBoxViewTypeChanged(int)));
    connect(ui->comboBoxDrawMode, SIGNAL(currentIndexChanged(int)), this, SLOT(comboBoxDrawModeChanged(int)));

    // 绘图完成
    connect(ui->graphicsView, SIGNAL(mouseReleased(QRectF)), this, SLOT(updateStatusBar(QRectF)));

    // 记录页面索引
    QObject::connect(ui->stackedWidget, &QStackedWidget::currentChanged, this, [this](int index){
        m_currentPageIndex = index;
    });

    // 设置快捷键
    setupShortcuts();
}
