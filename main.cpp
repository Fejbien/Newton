#include "mainwindow.h"

#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    QPixmap pixmap( 32, 32 );
    pixmap.fill( Qt::transparent );
    w.setWindowIcon(QIcon(pixmap));
    w.show();
    return a.exec();
}
