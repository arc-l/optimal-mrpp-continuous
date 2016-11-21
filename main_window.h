/*
 *  The main window class. 
 *
 *  Created on: Jan 30, 2015
 *  Author: Jingjin Yu
 */

#ifndef _MAIN_WINDOW_H_
#define _MAIN_WINDOW_H_

#include <QMainWindow>
#include <QMenuBar>
#include <QToolBar>
#include <QAction>
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QLineEdit>
#include <QGraphicsEllipseItem>
#include <QTimer>
#include <QApplication>
#include <CGAL/Qt/DemosMainWindow.h>

#include "helper_functions.h"
#include "types.h"
#include "graphics_item.h"
#include "roadmap.h"

class MainWindow : public CGAL::Qt::DemosMainWindow
{
    Q_OBJECT

public:
    MainWindow( QApplication *pApp, QString& fileFolder, 
		QString title = "Main Window", QString iconPath = "", bool addDefaultContent = false);

    void setupMenuAndToolbar();

public slots:
    void openEnvironment();

	void overlayLattice();
	void locateBoundingLatticeCycle();

    void createRandomProblem();
    void solveProblem();
	void animate();

	void fitView();

	void animatePath();

protected:
	QApplication * m_pApp;
	QString	m_fileFolder;

    QToolBar * m_pMainToolBar;

    QAction * m_pOpenEnvAction;

    QAction * m_pOverlayLatticAction;
    QAction * m_pLocateBoundaryAction;
    QAction * m_pConnectGraphAction;

    QAction * m_pCreateAction;
	QAction * m_pSolveAction;
	QAction * m_pPlayAction;

	QLabel	* m_pNumberofRobotLabel;
	QLineEdit *	m_pLineEdit;

	QLabel	*	m_pMinDistLabel;
	QLineEdit *	m_pMinDistLineEdit;

    QAction * m_pFitScreenAction;

	QString	m_envFileName;

    QGraphicsView* m_pView;
    QGraphicsScene m_scene;

	// Internal building state
	int	m_state;
	static const int BUILD_STATE_UNKNOWN = -1;
	static const int BUILD_STATE_MAP_OPEN = 0;
	static const int BUILD_STATE_LATTICE_OVERLAYED = 1;
	static const int BUILD_STATE_TRIMMED = 2;
	static const int BUILD_STATE_CONNECTION_FIXED = 3;

	// The environment
	double	m_radius;
    Polygon_2 m_boundingPoly;
    Polygon_2 m_boundingPoly2;
    Polygon2_list m_envPolyList;
    Polygon2_list m_envObsPolyList;

    AdvancedGraphicsItem<Polygon_2> * m_pBoundingPolyAGI;
    AdvancedGraphicsItem<Polygon_2> * m_pBoundingPolyAGI2;
    std::vector<AdvancedGraphicsItem<Polygon_2> *>  m_PolyAGIList;

	// Bounding rect
	QRectF m_boundingRect;

	// Roadmap
	Roadmap m_roadmap;

	// Number of robots
	int			m_numRobots;

	// Timer for animation
	QTimer *	m_pTimer;

	// Start and goal locations 
	vector<pair<double, double> >	m_svVec;
	vector<pair<double, double> >	m_gvVec;
	vector<pair<int, int> >			m_sgVec;
	
	// Paths
	map<int, vector<int> >			m_pathMap;

	// Current 
	map<int, pair<double, double> >	m_currentStartMap;
	map<int, pair<double, double> >	m_currentGoalMap;

	// Vector of robot graphics item pointers
	vector<QGraphicsEllipseItem*>	m_pRobotItemVec;
	vector<QGraphicsSimpleTextItem*>	m_pTextItemVec;
	
	static const int ANI_STATE_FROM_START = 0;
	static const int ANI_STATE_GRID_PATH = 10;
	static const int ANI_STATE_TO_GOAL = 20;
	static const int ANI_STATE_COMPLETE = 30;

	// Shortest path related
	double m_shotestPathLength;
	vector<QGraphicsLineItem*> m_pPathLineItems;

	int			m_aniState;			
	int			m_pathStep;
	double		m_aniStatePercent;

private:
	// Rendering
	void drawBasicEnvironment();

} ;

#endif //_MAIN_WINDOW_H_
