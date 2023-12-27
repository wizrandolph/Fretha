
#include "mainwindow.h"
#include "enumerator.h"
#include <QApplication>

#include <QLocale>
#include <QTranslator>
#include <iostream>
#include <cstdlib>
#include <QFile>
#include <QString>
#include <QTextStream>
#include <QMutex>
#include <QDateTime>
using namespace std;

void setDefaultFont(QApplication *app)
{
    int loadedFontId = QFontDatabase::addApplicationFont(":/font/MiSans.ttf");
    QStringList loadedFontFamilies = QFontDatabase::applicationFontFamilies(loadedFontId);
    for (int i = 0 ;i < loadedFontFamilies.size() ; i++) {
        QString sansCNFamily = loadedFontFamilies.at(i);
        if (i == 0 ) {
            QFont defaultFont = (app)->font();

            //设置字体大小格式为pixel，字体占用固定像素。相同像素的大小屏，不会发生字体遮挡现象
            defaultFont.setFamily(sansCNFamily);
            defaultFont.setPixelSize(16);
            (app)->setFont(defaultFont);
        }

        qDebug()<<"[Default Font]:\t" << sansCNFamily ;
    }
}

void myMessageOutput(QtMsgType type, const QMessageLogContext &context, const QString &msg)
{
    //加锁，防止多用户同时读写
    static QMutex mutex;
    mutex.lock();

    QString strMsg("");
    switch(type)
    {
    case QtDebugMsg:
        strMsg = QString("Debug:");
        break;
    case QtWarningMsg:
        strMsg = QString("Warning:");
        break;
    case QtCriticalMsg:
        strMsg = QString("Critical:");
        break;
    case QtFatalMsg:
        strMsg = QString("Fatal:");
        break;
    default:
        break;
    }

    // 设置输出信息格式
    QByteArray localMsg = msg.toLocal8Bit();
    QString strDateTime = QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss");
    QString strMessage = QString("%1     File:%2    Function:%3    Line:%4    %5 %6")
                             .arg(strDateTime).arg(context.file).arg(context.function).arg(context.line).arg(strMsg).arg(localMsg.constData());

    // 输出信息至文件中（读写、追加形式）
    QFile file(QApplication::applicationDirPath() + "/log.txt");
    file.open(QIODevice::ReadWrite | QIODevice::Append);
    QTextStream stream(&file);
    stream << strMessage << "\n";
    file.flush();
    file.close();

    // 解锁
    mutex.unlock();
}

int main(int argc, char *argv[])
{

    if (RUNMODE == RUNMODE_RELEASE)
        qInstallMessageHandler(myMessageOutput);

    QApplication a(argc, argv);

    QTranslator translator;
    const QStringList uiLanguages = QLocale::system().uiLanguages();
    for (const QString &locale : uiLanguages) {
        const QString baseName = "Fretha_" + QLocale(locale).name();
        if (translator.load(":/i18n/" + baseName)) {
            a.installTranslator(&translator);
            break;
        }
    }

    setDefaultFont(&a);
    QString styleSheet;
    QFile styleFile(":/style/darkblue.qss");
    if (styleFile.open(QFile::ReadOnly)) {
        // 读取样式表内容
        styleSheet = styleFile.readAll();
        styleFile.close();
        qDebug() << "[Apply Stylesheet] : Succeeded";
    }
    else {
        qDebug() << "[Apply Stylesheet] : Failed";
    }

    a.setStyleSheet(styleSheet);
    MainWindow w;
    w.show();
    return a.exec();
}
