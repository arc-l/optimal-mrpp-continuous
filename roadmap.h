/*
 *  The core class (header) for roadmap building. Roadmap building is done in four steps:
 *	1. Overlay a regular lattice (currently only hexagonal lattice) 
 *	2. Remove lattice eges that intersect with configuration space obstacles, and compute the 
 *	   smallest cycle on the remaining lattice that encloses the C-space obstacles. 
 *  3. When a cycle computed about does not lie completely in the free C-space, compute 
 *	   shortest paths in C-space that complete the cycle
 *	4. Post-processing the computed shortest paths for inclusion into the lattice graph
 *
 *  Created on: Jan 30, 2015
 *  Author: Jingjin Yu
 */

#ifndef _O_ROADMAP_H_
#define _O_ROADMAP_H_

#include "types.h"
#include "graph.h"
#include "shortest_path/visilibity.hpp"

#include <QColor>
#include <QPen>
#include <QBrush>
#include <QGraphicsScene>

#include <vector>
using namespace std;
using namespace VisiLibity;

class Roadmap{
private:
	double							m_edgeLength;			// Edge length of the hexgaon
	Polygon_2*						m_pBoundingRect;		// Bounding polygon
	Polygon2_list					m_obstaclePolyList;		// Obstacle polygon list pointer

	std::list<Point_2>				m_pointList;			// List of all vertices
	std::map<int, Point_2>			m_vidPointMap;			// Vertex id to point map
	std::map<Point_2, int>			m_pointVidMap;			// Point to vertex id map

	Graph							m_graph;				// The lattice graph
	Graph							m_finalGraph;			// The final graph for doing planning
	std::map<int, int>				m_vidFGMap;				// Vertex id MAP from final graph to graph
	std::map<int, int>				m_vidGFMap;				// Vertex id MAP from graph to final graph

	std::set<std::pair<int, int> >  m_edgeToBeRemovedSet;	// Edges to be removed
	std::set<int>					m_vertexToBeRemovedSet;	// Vertices to be removed (single degree)

	Graph							m_boundaryBoundingCycle;// Cycle just inside the bounding rect
	std::vector<Graph*>				m_obsBoundingCycleVec;	// List of cycles around obstacle

	typedef std::map<std::pair<int, int>, std::list<Point_2> > Path_Map;
	Path_Map m_connectingPathMap;							// Look up map for "bridging" paths

	// Temporary local variables
	double bottomLeftX, bottomLeftY, width, height;			// Bounding rectangle paramaters
	double sqrt3;											// sqrt(3)
	double xs, ys;											// Starting x, y for the raw lattice 
	int	   n_w, n_h;										// Raw lattice columns and rows

	
	// Converter
	CGAL::Cartesian_converter<K,ICK> K_ICK_converter; 

	// Visibility graph
	double						m_epsilon;
	Environment *				m_pEnvironment;
	Visibility_Graph *			m_pVisibilityGraph;

	double							m_radius;

public:
	Roadmap(){m_epsilon = 0.000001; m_pEnvironment = 0; m_pVisibilityGraph = 0;}

	// Build roadmap
	void buildRoadmap(Polygon2_list* pObsList, Polygon_2* pbondingRect, double radius);

	// Retrive all vertices and edges
	std::list<Point_2> & getAllVertices(){return m_pointList;};

	// Paint method
	void addToScene(QGraphicsScene& scene, bool drawEdge = true, 
		QPen edgePen = QPen(Qt::red, 0.025, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin),
		bool drawVertex = true, 
		QPen vertexPen = QPen(Qt::red, 0.05, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));

	// Build a hexagonal lattice that covers the entire bounding rectangle. The rectangle should 
	// contain the configuraiton space and not the free space
	void buildHexgaonLattice();
	void drawHexagonLattice(QGraphicsScene& scene, bool drawVetex = false);

	// Remove edges that do not fully belong to the configuration space
	void removeExcessEdges();
	void drawBoundingCycle(QGraphicsScene& scene);

	void drawVertexIds(QGraphicsScene& scene);

	// Build the final graph for graph-based computation
	void buildFinalGraph();

	// Generate random start and goal pairs
	void createRandomStartGoalPairs(int numRobots, double spacing, vector<pair<double, double> >& starts, 
		vector<pair<double, double> >& goals);

	// Solve a problem
	bool solveProblem(vector<pair<int, int> >& sgVec, map<int, vector<int> >& paths, 
		string& fileFolder, string& fileNameExtra = string(""));

	// Post process path to reduce ossilation 
	void improvePaths(map<int, vector<int> >& paths);

	// Retrive the location of a vertex given the int id in m_finalGraph
	pair<double, double> getVertexLocationFromID(int vid);

	// Build visibility graph
	void buildVisibilityGraph();

	// Compute shortest path between points 
	double computeShortestPath(Point_2& p1, Point_2& p2, std::list<Point_2> & path);
	double computeShortestPath(double x1, double y1, double x2, double y2, std::list<Point_2> & path);
	double computeShortestPath(double x1, double y1, double x2, double y2, vector<pair<double, double> >& path);

private:
	// Refresh edge list
	void refreshEdgeList();

	// Check whether a point is inside the configuration space
	bool isPointInCSpace(Point_2 &p);

	// Get the (col, row) of the hexagon a point belongs to 
	std::pair<int, int> locateHexagonForPoint(ICPoint_2 &p);

	// Check whether and edge is in an edge set
	bool edgeInSet(int v1, int v2, std::set<std::pair<int, int> >& edgeSet);

	// Add edge to graph if one is not already in
	void addEdgeIfNotThere(int v1, int v2, Graph & graph);

	// Get the next hexagon bordering the current at an edge
	std::pair<int, int> getBorderHexagon(int col, int row, int v1, int v2, int pIndex[]);

	// Get the next vertex in a sequence of vertices
	int getNextNodeInSequence(int i1, int i2, int pIndex[]);

	// Compute a hexgon on the lattice given a (col, row)
	void populateHexagon(int col, int row, Point_2* p, Segment_2* e, Polygon_2 &hex, int* pIndex);

	// Compute all intersecting edges of the full lattice with the polygon obstacle
	void getIntersectingEdges(Polygon_2 & poly, Graph& boundingCycle, bool outerBoundary = false);

	// Check whether a point belongs to some edge not in the obstacle
	int pointBelongToEdgeOutideObstacle(Polygon_2 & poly, int pIndex, std::set<std::pair<int,int> > &tempEdgeToRemoveSet, bool outerBoundary = false);

	// Check whether a point is outside an obstacle
	bool pointOutsideObstacle(Polygon_2 & poly, int pIndex, bool outerBoundary = false);

	// Generate random start and goal pairs
	void createRandomCoords(int numRobots, double spacing, vector<pair<double, double> >& coords);

public:

	// Snapping a set of robots to grid points
	void snapToGraph(vector<pair<double, double> >&coords, vector<int>&snapped);

};

#endif //_OPTMLY_ROADMAP_H_
