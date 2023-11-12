
#include "mainwindow.h"

#include <QApplication>

#include <QLocale>
#include <QTranslator>

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

int main(int argc, char *argv[])
{
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
