/*
 *	Main window class implemenation. 
 *
 *  Created on: Jan 30, 2015
 *  Author: Jingjin Yu
 */

#include "main_window.h"
#include "helper_functions.h"

#include <QLabel>
#include <QMenuBar>
#include <QMenu>
#include <QToolBar>
#include <QToolButton>
#include <QFileDialog>
#include <QGraphicsView>
#include <QMainWindow>
#if _MSC_VER == 1600
#include <qgl.h>
#endif

#include <CGAL/Qt/DemosMainWindow.h>

#include <iostream>
#include <ios>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/thread/thread.hpp> 

MainWindow::MainWindow( QApplication *pApp, QString& fileFolder, 
	QString title, QString iconPath, bool addDefaultContent)
:m_pApp(pApp), m_fileFolder(fileFolder), m_state(BUILD_STATE_UNKNOWN),CGAL::Qt::DemosMainWindow()
{
    m_pBoundingPolyAGI = 0;
	m_pBoundingPolyAGI2 = 0;

    // Set title
    setWindowTitle(title);

    // Set icon if the icon path is given
    if(iconPath.length() > 0){
        this->setWindowIcon(QIcon(iconPath));
    }

    // Check whether to display some default content
    if(addDefaultContent){
        resize(600, 400);
        QLabel * label = new QLabel(tr("Central Widget"));
        setCentralWidget(label);
        label->setAlignment(Qt::AlignCenter);
		m_pView->setMouseTracking(true);
    }
    else{
        m_pView = new QGraphicsView(this);
#if _MSC_VER == 1600
		m_pView->setViewport(new QGLWidget(QGLFormat(QGL::SampleBuffers)));
#endif

        m_scene.setItemIndexMethod(QGraphicsScene::NoIndex);
        m_scene.setSceneRect(-100, -100, 100, 100);
        m_pView->setScene(&m_scene);
        m_pView->setMouseTracking(true);

        this->view = m_pView;

        resize(620, 710);
        this->view->resize(620, 710);
        this->view->scale(1, 1);
        this->addNavigation(this->view);

        this->setupStatusBar();
        setCentralWidget(m_pView);
	}

	m_pView->setRenderHint(QPainter::Antialiasing);

	// Setup timer
	m_pTimer = new QTimer(this);
	connect(m_pTimer, SIGNAL(timeout()), this, SLOT(animatePath()));

	// Setup toolbar
    setupMenuAndToolbar();
}

void MainWindow::setupMenuAndToolbar(){
	int iconSize = 40;
    m_pMainToolBar = new QToolBar(tr("Main Toolbar"), this);
    m_pMainToolBar->setIconSize(QSize(iconSize, iconSize));
    addToolBar(m_pMainToolBar);
    // m_pMainToolBar->setVisible(false);


    // Create the environment menu
    QMenu* fileMenu = menuBar()->addMenu(tr("&Environment"));

    // Open action
    m_pOpenEnvAction = new QAction(tr("&Open Environment"), this);
    m_pOpenEnvAction->setStatusTip(tr("Open a new envrionment"));
    m_pOpenEnvAction->setIcon(QIcon("./resources/open.png"));
    connect(m_pOpenEnvAction, SIGNAL(triggered()), this, SLOT(openEnvironment()));
    fileMenu->addAction(m_pOpenEnvAction);
    m_pMainToolBar->addAction(m_pOpenEnvAction);

	m_pMainToolBar->addSeparator();

    // Connect graph action
    m_pOverlayLatticAction = new QAction(tr("&Overlay Lattice Graph"), this);
    m_pOverlayLatticAction->setStatusTip(tr("Overlay covering lattice graph"));
    m_pOverlayLatticAction->setIcon(QIcon("./resources/hexgrid.png"));
    connect(m_pOverlayLatticAction, SIGNAL(triggered()), this, SLOT(overlayLattice()));
	fileMenu->addSeparator();
    fileMenu->addAction(m_pOverlayLatticAction);
    m_pMainToolBar->addAction(m_pOverlayLatticAction);
	m_pOverlayLatticAction->setDisabled(true); // Disabled at beginning 

    m_pLocateBoundaryAction = new QAction(tr("&Trim Lattice Graph"), this);
    m_pLocateBoundaryAction->setStatusTip(tr("Remove lattice parts outside the free space"));
    m_pLocateBoundaryAction->setIcon(QIcon("./resources/hexcycle.png"));
    connect(m_pLocateBoundaryAction, SIGNAL(triggered()), this, SLOT(locateBoundingLatticeCycle()));
    fileMenu->addAction(m_pLocateBoundaryAction);
    m_pMainToolBar->addAction(m_pLocateBoundaryAction);
	m_pLocateBoundaryAction->setDisabled(true); // Disabled at beginning 

	QWidget* empty3 = new QWidget();
	empty3->setFixedWidth(4);
	m_pMainToolBar->addWidget(empty3);

	m_pNumberofRobotLabel = new QLabel(this);
	m_pNumberofRobotLabel->setText(" Robots: ");
	m_pNumberofRobotLabel->setFixedHeight(iconSize-2);
	m_pMainToolBar->addWidget(m_pNumberofRobotLabel);

	m_pLineEdit = new QLineEdit(this);
	m_pLineEdit->setFixedWidth(30);
	m_pLineEdit->setFixedHeight(iconSize-2);
	m_pLineEdit->setText(QString("25"));
	m_pMainToolBar->addWidget(m_pLineEdit);

	QWidget* empty1 = new QWidget();
	empty1->setFixedWidth(4);
	m_pMainToolBar->addWidget(empty1);

	m_pMinDistLabel = new QLabel(this);
	m_pMinDistLabel->setText(" Spacing: ");
	m_pMinDistLabel->setFixedHeight(iconSize-2);
	m_pMainToolBar->addWidget(m_pMinDistLabel);

	m_pMinDistLineEdit = new QLineEdit(this);
	m_pMinDistLineEdit->setFixedWidth(30);
	m_pMinDistLineEdit->setFixedHeight(iconSize-2);
	m_pMinDistLineEdit->setText(QString("3"));
	m_pMainToolBar->addWidget(m_pMinDistLineEdit);

	QWidget* empty2 = new QWidget();
	empty2->setFixedWidth(4);
	m_pMainToolBar->addWidget(empty2);

	m_pCreateAction = new QAction(tr("&Create Random Problem"), this);
    m_pCreateAction->setStatusTip(tr("Create a random instance"));
    m_pCreateAction->setIcon(QIcon("./resources/sg.png"));
    connect(m_pCreateAction, SIGNAL(triggered()), this, SLOT(createRandomProblem()));
    fileMenu->addAction(m_pCreateAction);
    m_pMainToolBar->addAction(m_pCreateAction);
	m_pCreateAction->setDisabled(true); // Disabled at beginning 

    m_pSolveAction = new QAction(tr("&Solve"), this);
    m_pSolveAction->setStatusTip(tr("Solve a randomly created instance"));
    m_pSolveAction->setIcon(QIcon("./resources/solution.png"));
    connect(m_pSolveAction, SIGNAL(triggered()), this, SLOT(solveProblem()));
    fileMenu->addAction(m_pSolveAction);
    m_pMainToolBar->addAction(m_pSolveAction);
	m_pSolveAction->setDisabled(true); // Disabled at beginning 

    m_pPlayAction = new QAction(tr("&Animate"), this);
    m_pPlayAction->setStatusTip(tr("Animate the solution"));
    m_pPlayAction->setIcon(QIcon("./resources/play.png"));
    connect(m_pPlayAction, SIGNAL(triggered()), this, SLOT(animate()));
    fileMenu->addAction(m_pPlayAction);
    m_pMainToolBar->addAction(m_pPlayAction);
	m_pPlayAction->setDisabled(true); // Disabled at beginning 

    // Fit to screen action
	m_pFitScreenAction = new QAction(tr("&Fit content to view"), this);
    m_pFitScreenAction->setStatusTip(tr("Fit content to view"));
    m_pFitScreenAction->setIcon(QIcon("./resources/fitview.png"));
    connect(m_pFitScreenAction, SIGNAL(triggered()), this, SLOT(fitView()));
	fileMenu->addSeparator();
    fileMenu->addAction(m_pFitScreenAction);
    m_pMainToolBar->addAction(m_pFitScreenAction);

}

void MainWindow::openEnvironment(){
	// Open file dialog for file selection
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Environment Description File"), "", tr("Files (*.*)"));
    if(fileName == 0 || fileName.length() < 0) return;

	QStringList qsl = fileName.split("/");
	m_envFileName = qsl[qsl.size()-1].split(".")[0];
	
	// Clean up
    m_scene.clear();
    m_boundingPoly.clear();
	m_boundingPoly2.clear();
	m_envPolyList.clear();
	m_envObsPolyList.clear();
	m_PolyAGIList.clear(); /* The element pointers are gone when the scene object is cleared */

	// Reading in the environment, first the raidus of the (disc) robot
    std::ifstream ifs(qPrintable(fileName));
    double radius;
    ifs >> radius;

	// Then the number of obstacle polygons
    int numberOfPolygons;
    ifs >> numberOfPolygons;

	// Then read in all obstacle polygons and compute the configuration space for a single disc
    for(int i = 0; i < numberOfPolygons; i ++){
        // Load polygon
        Polygon_2 tp;
        ifs >> tp;
        m_envPolyList.push_back(tp);

		// Add raw obstacle to scene and set fill to be true with fill transparency
        AdvancedGraphicsItem<Polygon_2>* pAGI = new AdvancedGraphicsItem<Polygon_2>(&(m_envPolyList.back()));
        m_PolyAGIList.push_back(pAGI);
        // m_scene.addItem(pAGI);
        pAGI->m_bShowEdge = false;
        pAGI->m_bShowVertices = false;
        pAGI->m_bFill = true;
        pAGI->m_fillBrush = QColor(16,16,16,192);

        // Computing the Minkowski sum
        Polygon_2 ep = growPolygonByRadius(tp, radius);
        m_envObsPolyList.push_back(ep);

		// Add the configuration space obstacle to the scene
        pAGI = new AdvancedGraphicsItem<Polygon_2>(&(m_envObsPolyList.back()));
        m_PolyAGIList.push_back(pAGI);
        // m_scene.addItem(pAGI);
        pAGI->m_bFill = false;
        pAGI->m_edgePen = QPen(Qt::gray, 0.25, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
        pAGI->m_vertexPen = QPen(Qt::black, 0.5, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
    }

	// Then read the bounding rectangle (configuration space)
    ifs >> m_boundingPoly;

	// Add bounding walls
	CGAL::Cartesian_converter<K,ICK> K_ICK_converter;

	// Add to scence
    m_pBoundingPolyAGI = new AdvancedGraphicsItem<Polygon_2>(&m_boundingPoly);
    m_pBoundingPolyAGI->m_bFill = false;
    m_pBoundingPolyAGI->m_bShowVertices = false;
    m_pBoundingPolyAGI->m_bShowEdge = true;
	m_pBoundingPolyAGI->m_edgePen = QPen(Qt::gray, 0.25, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
    m_PolyAGIList.push_back(m_pBoundingPolyAGI);
    // m_scene.addItem(m_pBoundingPolyAGI);

	// Then compute the outside bounding rectangle
	double x1, y1, x2, y2;
	x1 = K_ICK_converter(m_boundingPoly[0].x()); y1 = K_ICK_converter(m_boundingPoly[0].y());
	x2 = K_ICK_converter(m_boundingPoly[2].x()); y2 = K_ICK_converter(m_boundingPoly[2].y());
	m_boundingPoly2.push_back(Point_2(x1 - radius, y1 - radius));
	m_boundingPoly2.push_back(Point_2(x1 - radius, y2 + radius));
	m_boundingPoly2.push_back(Point_2(x2 + radius, y2 + radius));
	m_boundingPoly2.push_back(Point_2(x2 + radius, y1 - radius));

	// Add to scene
    m_pBoundingPolyAGI2 = new AdvancedGraphicsItem<Polygon_2>(&m_boundingPoly2);
    m_pBoundingPolyAGI2->m_bFill = false;
    m_pBoundingPolyAGI2->m_bShowVertices = false;
    m_pBoundingPolyAGI2->m_bShowEdge = true;
	m_pBoundingPolyAGI2->m_edgePen = QPen(Qt::black, 0.25, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
    m_PolyAGIList.push_back(m_pBoundingPolyAGI2);
    // m_scene.addItem(m_pBoundingPolyAGI2);

	m_radius = radius;

	drawBasicEnvironment();

	// Do roadmap building setup
	m_roadmap.buildRoadmap(&m_envObsPolyList, &m_boundingPoly, radius);
	// m_roadmap.addToScene(m_scene);

	m_boundingRect = QRectF(x1 - radius*4, y1 - radius*4, x2 + radius*4, y2 + radius*4);

    // Fit to the view
    fitView();

	// Enable the actions for building the graph
	m_state = BUILD_STATE_MAP_OPEN;
	m_pOverlayLatticAction->setEnabled(true);
    m_pLocateBoundaryAction->setEnabled(true);
	m_pCreateAction->setEnabled(false);
	m_pSolveAction->setEnabled(false);
	m_pPlayAction->setEnabled(false);

}

void MainWindow::drawBasicEnvironment(){
	m_PolyAGIList.clear();
	for(Polygon2_list::iterator pit = m_envPolyList.begin(); pit != m_envPolyList.end(); pit++){
		// Add raw obstacle to scene and set fill to be true with fill transparency
        AdvancedGraphicsItem<Polygon_2>* pAGI = new AdvancedGraphicsItem<Polygon_2>(&(*pit));
        m_PolyAGIList.push_back(pAGI);
        pAGI->m_bShowEdge = false;
        pAGI->m_bShowVertices = false;
        pAGI->m_bFill = true;
        pAGI->m_fillBrush = QColor(16,16,16,192);
    }

	for(Polygon2_list::iterator pit = m_envObsPolyList.begin(); pit != m_envObsPolyList.end(); pit++){
		// Add the configuration space obstacle to the scene
        AdvancedGraphicsItem<Polygon_2>* pAGI = new AdvancedGraphicsItem<Polygon_2>(&(*pit));
        m_PolyAGIList.push_back(pAGI);
        pAGI->m_bFill = false;
        pAGI->m_edgePen = QPen(Qt::gray, 0.25, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
        pAGI->m_vertexPen = QPen(Qt::black, 0.5, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
    }

	// Add to scence
    m_pBoundingPolyAGI = new AdvancedGraphicsItem<Polygon_2>(&m_boundingPoly);
    m_pBoundingPolyAGI->m_bFill = false;
    m_pBoundingPolyAGI->m_bShowVertices = false;
    m_pBoundingPolyAGI->m_bShowEdge = true;
	m_pBoundingPolyAGI->m_edgePen = QPen(Qt::gray, 0.25, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
    m_PolyAGIList.push_back(m_pBoundingPolyAGI);

	// Add to scene
    m_pBoundingPolyAGI2 = new AdvancedGraphicsItem<Polygon_2>(&m_boundingPoly2);
    m_pBoundingPolyAGI2->m_bFill = false;
    m_pBoundingPolyAGI2->m_bShowVertices = false;
    m_pBoundingPolyAGI2->m_bShowEdge = true;
	m_pBoundingPolyAGI2->m_edgePen = QPen(Qt::black, 0.25, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
    m_PolyAGIList.push_back(m_pBoundingPolyAGI2);

	for(std::vector<AdvancedGraphicsItem<Polygon_2> *>::iterator pagiItor = m_PolyAGIList.begin();
		pagiItor != m_PolyAGIList.end(); pagiItor ++){
        AdvancedGraphicsItem<Polygon_2>* pAGI = *pagiItor;
        m_scene.addItem(pAGI);
    }
}


void MainWindow::overlayLattice(){
	m_scene.clear();
	m_roadmap.buildHexgaonLattice();
	drawBasicEnvironment();
	m_roadmap.drawHexagonLattice(m_scene, true);

	// Set build state and enable the actions for building the graph
	m_state = BUILD_STATE_LATTICE_OVERLAYED;
	m_pOverlayLatticAction->setEnabled(false);
}

void MainWindow::locateBoundingLatticeCycle(){
	// Clean up
	m_scene.clear();

	// Build depending on the build state
	switch(m_state){
	case BUILD_STATE_MAP_OPEN:
		m_roadmap.buildHexgaonLattice();
		m_pOverlayLatticAction->setEnabled(false);
	case BUILD_STATE_LATTICE_OVERLAYED: 
		m_roadmap.removeExcessEdges();
		m_pLocateBoundaryAction->setEnabled(false);
	default:;
	}

	// Draw environment
	drawBasicEnvironment();

	// Draw the lattice (now truncated)
	m_roadmap.drawHexagonLattice(m_scene, true);
#ifdef _DEBUG
	m_roadmap.drawBoundingCycle(m_scene);
	m_roadmap.drawVertexIds(m_scene);
#endif

	// Set build state
	m_state = BUILD_STATE_TRIMMED;

	// Enable random instance creation
	m_pCreateAction->setEnabled(true); 
}

// Create and solve a random problem over current graph
void MainWindow::createRandomProblem() {
	// Stop any ongoing timer
	m_pTimer->stop();

	// Temporarily disable the action
	m_pCreateAction->setEnabled(false);
	m_pApp->processEvents();

	// Do clean up
	m_pRobotItemVec.clear();
	m_pTextItemVec.clear();
	m_svVec.clear();
	m_gvVec.clear();
	m_pathMap.clear();
	m_pPathLineItems.clear();
	
	// Get number of robots
	m_numRobots = m_pLineEdit->text().toInt();
	double spacing = m_pMinDistLineEdit->text().toDouble();

	srand (time(NULL));

	// Setup problem 
	vector<pair<double, double> > vvg;
	vector<int> svv;
	vector<int> gvv;
	while(true){
		try{
			// Create random start and goal vertices
			m_svVec.clear(); m_gvVec.clear(); svv.clear(); gvv.clear();
			m_roadmap.createRandomStartGoalPairs(m_numRobots, spacing, m_svVec, m_gvVec);

			// Snap to graph
			m_sgVec.clear();
			m_roadmap.snapToGraph(m_svVec, svv);
			m_roadmap.snapToGraph(m_gvVec, gvv);
			break;
		}
		catch(...){
		}
	}

	for(int r = 0;  r < m_numRobots; r ++){
		m_sgVec.push_back(pair<int, int>(svv[r], gvv[r]));
	}
#ifdef _DEBUG
	for(int r = 0; r < m_numRobots; r ++){
		cout << setw(3) << svv[r];
	}
	cout << endl;
	for(int r = 0; r < m_numRobots; r ++){
		cout << setw(3) << gvv[r];
	}
	cout << endl;
#endif
	
	// Draw start and end robot locations
	m_scene.clear();
	drawBasicEnvironment();
#ifdef _DEBUG
	m_roadmap.drawBoundingCycle(m_scene);
	m_roadmap.drawVertexIds(m_scene);
#endif

	// Draw the robots
	QFont font;
	QPainterPath path;
	font.setPointSizeF(m_radius/1.5);
	font.setBold(true);

	for(int r = 0; r < m_numRobots; r ++){
		// Robots at current (start) locations
		QGraphicsEllipseItem *prs = m_scene.addEllipse(m_svVec[r].first - m_radius, m_svVec[r].second - m_radius, m_radius*2-0.2, m_radius*2-0.2, 
			QPen(Qt::black, 0.2, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin),
			QBrush(Qt::blue));
		prs->setZValue(10);
		m_pRobotItemVec.push_back(prs);

		// Labels for robots at current (start) locations
		QGraphicsSimpleTextItem *ti = m_scene.addSimpleText(QString::number(r+1), font);
		ti->setPos(m_svVec[r].first - m_radius/2*(3.5/(m_radius+1)) + (r < 9 ? m_radius / 6 *(2.2/(m_radius-0.3)) :(r > 99 ? -m_radius / 6 *(2.2/(m_radius-0.3)) :0 )) , m_svVec[r].second - m_radius/2*(2./(m_radius-0.5)));
		ti->setPen(QPen(QColor(Qt::white), 0.15, Qt::SolidLine, Qt::RoundCap,Qt::RoundJoin));
		ti->setZValue(15);
		m_pTextItemVec.push_back(ti);

		// Robots at goal locations
		QGraphicsEllipseItem *prg = m_scene.addEllipse(m_gvVec[r].first - m_radius, m_gvVec[r].second - m_radius, m_radius*2-0.2, m_radius*2-0.2, 
			QPen(QColor(127, 127, 127, 64), 0.2, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin),
			QBrush(QColor(255, 127, 127, 64)));
		prg->setZValue(5);

		// Labels for robots at goal locations
		ti = m_scene.addSimpleText(QString::number(r+1), font);
		ti->setPos(m_gvVec[r].first - m_radius/2*(3.5/(m_radius+1)) + (r < 9 ? m_radius / 6 *(2.2/(m_radius-0.3)) :(r > 99 ? -m_radius / 6 *(2.2/(m_radius-0.3)) :0 )), m_gvVec[r].second - m_radius/2*(2./(m_radius-0.5)));
		ti->setPen(QPen(QColor(Qt::red), 0.15, Qt::SolidLine, Qt::RoundCap,Qt::RoundJoin));
		ti->setZValue(6);
	}

	// Always build visibility graph and compute shortest possible makespan
	m_roadmap.buildVisibilityGraph();
	double maxLength = 0.;
	for(int r = 0; r < m_numRobots; r ++){
		vector<pair<double, double> > path;
		double length = m_roadmap.computeShortestPath(m_svVec[r].first, m_svVec[r].second, m_gvVec[r].first, m_gvVec[r].second, path);
		if(length > maxLength){ maxLength = length; }
	}
	m_shotestPathLength = maxLength;
#ifdef _DEBUG
	cout << "Min possible makespan: " << m_shotestPathLength << endl;
#endif 

	// Update action enable/disable status
	m_pCreateAction->setEnabled(true);
	m_pSolveAction->setEnabled(true);
	m_pPlayAction->setEnabled(false);
}


void MainWindow::solveProblem() {
	// Setup animation. Note that the solution will be computed the first time animatePath()
	// is called
	m_pCreateAction->setEnabled(false);
	m_pSolveAction->setEnabled(false);
	m_roadmap.solveProblem(m_sgVec, m_pathMap, m_fileFolder.toStdString(), m_envFileName.toStdString());
	m_pPlayAction->setEnabled(true);
	m_pCreateAction->setEnabled(true);
}

void MainWindow::animate(){
	// Remove all paths line items from scence
	for(vector<QGraphicsLineItem*>::iterator lii = m_pPathLineItems.begin(); 
		lii != m_pPathLineItems.end(); lii ++)
	{
		m_scene.removeItem(*lii);
	}
	m_pPathLineItems.clear();

	// Disable play action
	m_pPlayAction->setEnabled(false);
	m_pApp->processEvents();

	// Sleep for half a second
	boost::this_thread::sleep(boost::posix_time::milliseconds(500));

	// Put all start to the initial location
	m_aniState = ANI_STATE_FROM_START;
	m_pathStep = 0;
	m_aniStatePercent = 0.0;
	for(int r = 0; r < m_numRobots; r++){
		m_currentStartMap[r] = m_svVec[r];
		m_currentGoalMap[r] = m_roadmap.getVertexLocationFromID(m_pathMap[r][0]);
		double x = m_currentStartMap[r].first*(1-m_aniStatePercent) + m_currentGoalMap[r].first*m_aniStatePercent;
		double y = m_currentStartMap[r].second*(1-m_aniStatePercent) + m_currentGoalMap[r].second*m_aniStatePercent;
		m_pRobotItemVec[r]->setPos(x - m_svVec[r].first, y - m_svVec[r].second);
		m_pTextItemVec[r]->setPos(x - m_radius/2*(3.5/(m_radius+1)) + (r < 9 ? m_radius / 6*(2.2/(m_radius-0.3)):0), y - m_radius/2*(2./(m_radius-0.5)));
	}
	m_pApp->processEvents();

	// Sleep for half a second
	boost::this_thread::sleep(boost::posix_time::milliseconds(2000));

	// Start animation
	m_pTimer->start(35);
}

void MainWindow::animatePath(){
	// If we are at the beginning of a stage, set the temp starts and goals
	if(m_aniStatePercent == 0.0){
		for(int r = 0; r < m_numRobots; r ++){
			switch(m_aniState){
			// Snapping from start to the graph
			case ANI_STATE_FROM_START:
				m_currentStartMap[r] = m_svVec[r];
				m_currentGoalMap[r] = m_roadmap.getVertexLocationFromID(m_pathMap[r][0]);
				break;
			case ANI_STATE_GRID_PATH:
				m_currentStartMap[r] = m_roadmap.getVertexLocationFromID(m_pathMap[r][m_pathStep]);
				m_currentGoalMap[r] = m_roadmap.getVertexLocationFromID(m_pathMap[r][m_pathStep + 1]);
				break;
			case ANI_STATE_TO_GOAL:
				m_currentStartMap[r] = m_roadmap.getVertexLocationFromID(m_pathMap[r][m_pathStep]);
				m_currentGoalMap[r] = m_gvVec[r];
				break;
			case ANI_STATE_COMPLETE:
			default:
				return;
			}
		}
	}

	// Animate
	m_aniStatePercent += 0.100001;
	if(m_aniState == ANI_STATE_FROM_START || m_aniState == ANI_STATE_TO_GOAL) m_aniStatePercent += 0.100001;
	for(int r = 0; r < m_numRobots; r++){
		double x = m_currentStartMap[r].first*(1-m_aniStatePercent) + m_currentGoalMap[r].first*m_aniStatePercent;
		double y = m_currentStartMap[r].second*(1-m_aniStatePercent) + m_currentGoalMap[r].second*m_aniStatePercent;
		m_pRobotItemVec[r]->setPos(x - m_svVec[r].first, y - m_svVec[r].second);
		m_pTextItemVec[r]->setPos(x - m_radius/2*(3.5/(m_radius+1)) + (r < 9 ? m_radius / 6*(2.2/(m_radius-0.3)):0), y - m_radius/2*(2./(m_radius-0.5)));
	}

	// Check current stage perecentage
	if(m_aniStatePercent > 1.0){
		m_aniStatePercent = 0.0;
		switch(m_aniState){
		// Snapping from start to the graph
		case ANI_STATE_FROM_START:
			m_aniState = ANI_STATE_GRID_PATH;
			break;
		case ANI_STATE_GRID_PATH:
			m_pathStep ++;
			if(m_pathStep == m_pathMap[0].size() - 1){
				m_aniState = ANI_STATE_TO_GOAL;
			}
			break;
		case ANI_STATE_TO_GOAL:
			m_aniState = ANI_STATE_COMPLETE;
			m_pTimer->stop();
			m_pPlayAction->setEnabled(true);
			break;
		}
	}
}

// Fit the content to the current view port
void MainWindow::fitView(){
    if(m_pBoundingPolyAGI2 != 0){
        m_pView->setSceneRect(m_pBoundingPolyAGI2->boundingRect());
        m_pView->fitInView(m_boundingRect, Qt::KeepAspectRatio);
    }
}
