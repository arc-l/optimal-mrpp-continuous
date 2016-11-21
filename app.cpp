/*
 *	Main app class
 *
 *  Created on: Jan 30, 2015
 *  Author: Jingjin Yu
 */

#include "helper_functions.h"
#include "main_window.h"

// Qt headers
#include <QtGui>
#include <QString>
#include <QFileDialog>
#include <QGraphicsLineItem>
#include <QApplication>

// GraphicsView items and event filters (input classes)
#include <CGAL/Qt/GraphicsViewPolylineInput.h>
#include <CGAL/Qt/PolygonGraphicsItem.h>
#include <CGAL/Qt/PolygonWithHolesGraphicsItem.h>
#include <CGAL/Qt/LineGraphicsItem.h>

// Base window class
#include <CGAL/Qt/DemosMainWindow.h>
#include <CGAL/Qt/Converter.h>


int main(int argc, char* argv[]){
	QApplication a(argc, argv, true, 0);
	MainWindow mainWindow(&a, (argc > 1 ? QString(argv[1]) : QString("F:\\temp\\mr")), "Multi-Robot Path Planning", "./resources/app.ico", false);
	mainWindow.show();


	return a.exec();
}
