#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QtCharts>
#include <QFileDialog>
#include <QMap>

#include <enumerator.h>
#include <fretcalculator.h>
#include <fretimageprocessor.h>
#include <fretthasolver.h>
#include <wizgraphicsview.h>
#include <wizshortcut.h>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private:

    QSettings *m_setting;   // 配置文件
    QString m_tmpDirPath = QCoreApplication::applicationDirPath();    // 程序运行的目录
    QDialog *m_dialog;    // 警告对话框
    QString m_savePath;   // 保存目录，一般为 strPath + "/THAResult"
    QString m_globalPath;    // 总文件夹的完整绝对路径，其下包含多个FRET共定位视野
    QString m_currentViewName;    // 视野名字
    int m_currentPageIndex;

    FretCalculator *fretCalculator = new FretCalculator();
    FretImageProcessor *fretImageProcessor = new FretImageProcessor();
    FretThaSolver *fretThaSolver = new FretThaSolver();

    QImage m_image2Show;
    QString m_currentViewPath; // 当前显示图片的视野的完整绝对路径
    ViewType m_viewType;  // 当前要显示的通道

    //图表
    QChartView *m_chartView[6];

    Ui::MainWindow *ui;


    void loadLastRatio();
    void initSignalsAndSlots();
    void initUICharts();
    void initUI();
    void initAllCharts();

    bool checkPathLegal(QString path);

    void showAlertDialog(QString);
    void showProgressDialog(QString);

    void updateGraphicsView();
    void exportToCSV(const QTableView* tableView, const QString& filePath);

    QString getSaveFolderPath(QWidget* parent);
    QString getSaveCsvFilePath(QWidget* parent);
    QString getOpenCsvFilePath();

    void initRecordTableModel();
    void initViewTableModel(QString batchPath);

    QImage getImageToShow();

    void changeView(QModelIndex &index);
    void trackRect(QModelIndex &index);

    void updateResultValue();

    void updateChartsOurs();
    void updateChartsBen();
    void updateChartsDu();

    void initializeChart(QChartView* chartView, QString chartName);
    void clearChart(QChartView* chartView);
    void drawScatter(QChartView* chartView, QString seriesName, const std::vector<double> xData, const std::vector<double> yData);
    void drawLine(QChartView* chartView, QString seriesName, const std::vector<double> xData, const std::vector<double> yData);
    void setupChartAxis(QChartView *chartView, const QString &xAxisLabel, qreal xAxisMin, qreal xAxisMax, const QString &yAxisLabel, qreal yAxisMin, qreal yAxisMax);

    void updateRectRecorded(QString viewPath);

    void updateResultsDu();
    void updateResultsOurs();
    void updateResultsBen();

    void setupShortcuts();


private slots:
    void on_pushButtonBrowse_clicked();

    void on_pushButtonAuto_clicked();
    void on_pushButtonManu_clicked();
    void on_pushButtonRatio_clicked();
    void on_pushButtonSetRatio_clicked();

    void on_pushButtonDu_clicked();
    void on_pushButtonOurs_clicked();
    void on_pushButtonBen_clicked();

    void on_pushButtonBinExec_clicked();
    void on_pushButtonUpdateEdRc_clicked();

    void on_pushButtonHome_clicked();

    void on_pushButtonSave_clicked();

    void on_pushButtonAdd_clicked();
    void on_pushButtonDelete_clicked();

    void on_pushButtonExport_clicked();
    void on_pushButtonCalc_clicked();
    void on_pushButtonBack_clicked();
    void on_pushButtonImportScreen_clicked();

    void on_pushButtonAutoGen_clicked();
    void on_pushButtonImportMask_clicked();

    void on_pushButtonClear_clicked();

    void on_ratioButtonLoose_clicked();
    void on_ratioButtonModerate_clicked();
    void on_ratioButtonStrict_clicked();

    void solveFinished();
    void updateStatusBar(QRectF);
    void receiveData(const QMap<TableHeader, QString> &data);

    void comboBoxViewTypeChanged(int index);
    void comboBoxDrawModeChanged(int index);
    void checkBoxBinChanged(int state);
};
#endif // MAINWINDOW_H
