/*
 *  The core class (implementation) for roadmap building
 *
 *  Created on: Jan 30, 2015
 *  Author: Jingjin Yu
 */

#include "roadmap.h"
#include "helper_functions.h"

#include <utility>   
#include <algorithm> 
#include <vector>

#include <CGAL/Qt/Converter.h>
#include <CGAL/Boolean_set_operations_2.h>

#include <QGraphicsSimpleTextItem>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/lookup_edge.hpp>

/*------------------------------------------------------------------------------------------
			see resources/lattice-indexing.pptx for the processing logic
--------------------------------------------------------------------------------------------*/

static const QColor BASIC_QCOLORS8[] = {Qt::red, Qt::blue, Qt::green, Qt::magenta, Qt::cyan, Qt::gray, Qt::black, Qt::yellow};
static int colorCounter = 0;

void Roadmap::buildRoadmap(Polygon2_list* pObsList, Polygon_2 *pBoundingRect, double radius){
	// Populate some internal variables for use across calls
	m_radius = radius;
	m_edgeLength = radius/0.43;
	m_obstaclePolyList = *pObsList;
	m_pBoundingRect = pBoundingRect;
	if(m_pVisibilityGraph != 0){
		delete m_pVisibilityGraph;
		m_pVisibilityGraph = 0;
	}
	if(m_pEnvironment != 0){
		delete m_pEnvironment;
		m_pEnvironment = 0;
	}

	// Some basic setup
	ICPoint_2 bottomLeft = K_ICK_converter((*m_pBoundingRect)[0]);
	ICPoint_2 topRight = K_ICK_converter((*m_pBoundingRect)[2]);

	bottomLeftX = bottomLeft.x();
	bottomLeftY = bottomLeft.y();
	width = topRight.x() - bottomLeft.x();
	height = topRight.y() - bottomLeft.y();;
	sqrt3 = sqrt(3.0);

	// Compute number of columns and rows
	n_w = (int)(ceil(width/(m_edgeLength*3/2))) + 3;
	n_h = (int)(ceil(height/(m_edgeLength*sqrt3))) + 3;

	// The lattice start x, y
	xs = bottomLeftX - (3/2)*m_edgeLength*1.35;
	ys = bottomLeftY - sqrt3*m_edgeLength*1.4;

	// Clean up from previous build
	m_pointList.clear();
	m_vidPointMap.clear();
	m_pointVidMap.clear();
	m_graph.clear();
	m_finalGraph.clear();
	m_vidFGMap.clear();
	m_vidGFMap.clear();
	m_edgeToBeRemovedSet.clear();
	m_vertexToBeRemovedSet.clear();
	m_boundaryBoundingCycle.clear();
	m_connectingPathMap.clear();
	for(std::vector<Graph*>::iterator git = m_obsBoundingCycleVec.begin(); git != m_obsBoundingCycleVec.end(); git++){
		delete (*git);
	}
	m_obsBoundingCycleVec.clear();

	/*
	// Build the roadmap, first obtain a lattice that cover the outer boundary
	buildHexgaonLattice();

	// Remove extra edges, at the same time, find smallest cycles enclosing the obstacles
	removeExcessEdges();

	// Preserve connectivity
	checkAndFixConnectivity();
*/
}

void Roadmap::removeExcessEdges(){
	// We do this in several steps. First, we go through each polygon obstacle boundary 
	// and delete all edges of the lattice that intersect with these boundaries. Then, 
	// we find the connected components of the remaining lattice graph. For each component
	// we only need to test one vertex to know whether it belongs to the configuration space
	// or not. We keep all components that belong to the configuration space

	// =====================================================================================
	// Compute the set of edges that falls on obstacle boundaries, at the same time also
	// compute the smallest cycle in the full lattice that encloses the obstacle. First do
	// it for the bounding polygon
	getIntersectingEdges(*m_pBoundingRect, m_boundaryBoundingCycle, true);

	// Then for all obstacles 
	for(Polygon2_list::iterator obsit = m_obstaclePolyList.begin(); obsit != m_obstaclePolyList.end(); obsit++){
		Graph* pg = new Graph();
		m_obsBoundingCycleVec.push_back(pg);
		getIntersectingEdges(*obsit, *pg);
	}

	// =====================================================================================
	// Remove edges that do not belong to the graph
	Graph g;
	std::set<std::pair<int, int> > edgeSet = m_graph.getEdgeSet();
	for(std::set<std::pair<int, int> >::iterator eit = edgeSet.begin(); eit != edgeSet.end(); eit++){
		int fv = (*eit).first;
		int sv = (*eit).second;
		if(!edgeInSet(fv, sv, m_edgeToBeRemovedSet)){
			g.addEdge(fv, sv);
		}
	}
	m_graph = g;

	// =====================================================================================
	// Go through all vertices and find all vertices and edges inside the configuration
	// space. To do so, iterate over all lattices points (the full lattice minus the edges that
	// crosses obstacle boundaries) and for each point, if it has not been checked, see whether
	// the point is inside the configuration space. Then we do a BFS from the point and visit
	// all points/edges connected to the point. We add the edges to our final graph if and only
	// if the starting vertex is in the c-space. This way, we need to do c-space membership
	// check geometrically only very limited number of times, usually around the number of
	// obstacles in the c-space.

	g.clear();
	std::set<int> visitedVertices;
	for(std::map<int, Point_2>::iterator vit = m_vidPointMap.begin(); vit != m_vidPointMap.end(); vit++){
		int vid = vit->first;
		if(visitedVertices.find(vid) == visitedVertices.end()){
			visitedVertices.insert(vid);
			// Test whether the vertex is inside the configuration space
			bool inCSpace = isPointInCSpace(vit->second);

			// Do BFS
			std::list<int> tempQueue;
			tempQueue.push_back(vid);
			while(tempQueue.size() > 0){
				int current = tempQueue.front();
				tempQueue.pop_front();

				std::set<int> neighborSet = m_graph.getNeighborSet(current);

				for(std::set<int>::iterator vit = neighborSet.begin(); vit != neighborSet.end(); vit++){
					// Retrieve first and second vertices
					int et = *vit;
					if(visitedVertices.find(et) == visitedVertices.end()){
						visitedVertices.insert(et);
						tempQueue.push_back(et);
					}
					// Add edge as needed
					if(inCSpace) g.addEdge(current, et);
				}
			}	
		}
	}
	m_graph.clear();
	m_graph = g;

	// Remove isolated vertices
	for(std::set<int>::iterator vit = m_vertexToBeRemovedSet.begin(); vit != m_vertexToBeRemovedSet.end(); vit++){
		// std::cout << "Removing vertex with id: " << *vit << std::endl;
		m_graph.removeVertex(*vit);
	}

	// Remove single degree edges
	std::set<int> vSetCopy = m_graph.getVertexSet();
	for(std::set<int>::iterator vit = vSetCopy.begin(); vit != vSetCopy.end(); vit++){
		if(m_graph.hasVertex(*vit)){
			if(m_graph.getNeighborSet(*vit).size() == 0){
				m_graph.removeVertex(*vit);
				// std::cout << "Removing vertex with id: " << *vit << std::endl;
			}
		}
	}

	buildFinalGraph();
}

void Roadmap::buildFinalGraph(){
	// Construct m_finalGraph
	std::set<std::pair<int, int>>& eSet = m_graph.getEdgeSet();
	int vIDCount = 0;
	for(std::set<std::pair<int, int> >::iterator eit = eSet.begin(); eit != eSet.end(); eit++){
		int fv = eit->first;
		int sv = eit->second;

		// Check whether fv was added
		if(m_vidGFMap.find(fv) == m_vidGFMap.end()){
			m_vidGFMap[fv] = vIDCount;
			m_vidFGMap[vIDCount] = fv;
			fv = vIDCount++;
		}
		else{
			fv = m_vidGFMap[fv];
		}

		// Check whether sv was added
		if(m_vidGFMap.find(sv) == m_vidGFMap.end()){
			m_vidGFMap[sv] = vIDCount;
			m_vidFGMap[vIDCount] =sv;
			sv = vIDCount++;
		}
		else{
			sv = m_vidGFMap[sv];
		}
		m_finalGraph.addEdge(fv, sv);
	}
}

void Roadmap::drawBoundingCycle(QGraphicsScene& scene){
	// Draw the cycles
	QPen regularPen = QPen(QColor(0, 255, 0, 127), 0.5, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
	QPen illegalPen = QPen(QColor(255, 0, 0, 127), 0.5, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
	for(std::vector<Graph*>::iterator git = m_obsBoundingCycleVec.begin(); git != m_obsBoundingCycleVec.end(); git++){
		Graph& g = *(*git);
		std::set<std::pair<int, int> > edgeSet = g.getEdgeSet();
		for(std::set<std::pair<int, int> >::iterator eit = edgeSet.begin(); eit != edgeSet.end(); eit++){
			int v1 = (*eit).first;
			int v2 = (*eit).second;
			ICPoint_2 p1 = K_ICK_converter(m_vidPointMap[v1]);
			ICPoint_2 p2 = K_ICK_converter(m_vidPointMap[v2]);
			if(m_graph.hasEdge(v1, v2)){
				scene.addLine(p1.x(), p1.y(), p2.x(), p2.y(), regularPen);
			}
			else{
				scene.addLine(p1.x(), p1.y(), p2.x(), p2.y(), illegalPen);
			}
		}
	}
}

bool Roadmap::isPointInCSpace(Point_2 &p){
	// Check whether the point is inside the bounding rect
	if(m_pBoundingRect->bounded_side(p) == CGAL::ON_UNBOUNDED_SIDE){
		return false;
	}

    for(Polygon2_list::iterator pli = m_obstaclePolyList.begin(); pli != m_obstaclePolyList.end(); pli++){
		Polygon_2 &tp = *(pli);
		if(tp.bounded_side(p) == CGAL::ON_BOUNDED_SIDE){
			return false;
		}
	}
	return true;
}


void Roadmap::drawVertexIds(QGraphicsScene& scene){
	// Draw the id text of the vretex
	QFont font;
	QPainterPath path;
	font.setPointSizeF(m_radius/1.5);
	font.setBold(false);
	for(std::set<int>::iterator vit = m_graph.getVertexSet().begin(); vit != m_graph.getVertexSet().end(); vit ++){
		ICPoint_2 p = K_ICK_converter(m_vidPointMap[*vit]);
		QGraphicsSimpleTextItem *ti = scene.addSimpleText(QString::number(m_vidGFMap[*vit]), font);
		ti->setPos(p.x() + m_radius/2, p.y() - m_radius/2);
		ti->setPen(QPen(QColor(Qt::green), 0.03*m_radius, Qt::SolidLine, Qt::RoundCap,Qt::RoundJoin));
		ti->setZValue(2);
	}
}

void Roadmap::buildVisibilityGraph(){
	// To build the visibility graph, we use the package VisiLibity, which is not precise arithematic, 
	// but good for our purpose since we do not use visibility graph as part of our main logic. 

	// Only build once per roadmap
	if(m_pVisibilityGraph != 0) return;

	// Vector to hold all polygon for contructing VisiLibity object
	vector<Polygon> polyVec;

	// We assume that the infated obstacles do not intersect each other but may intersect the boundary.
	// Therefore, we test such intersection and obtain an updated boundary (the inside)
	ECPolygon_2 boundary = convertToExactPolygon(*m_pBoundingRect);
	if(boundary.is_clockwise_oriented()) boundary.reverse_orientation();
	for(Polygon2_list::iterator pli = m_obstaclePolyList.begin(); pli != m_obstaclePolyList.end(); pli++){
		ECPolygon_2 ecPoly = convertToExactPolygon(*pli);
		if(boundaryInterset(boundary, ecPoly)){
			vector<ECPolygon_with_holes_2> outVec;

			// Remove the obstacle "from" the boundary polygon
			CGAL::difference(boundary, ecPoly, back_inserter(outVec)); 

			// Update boundary
			boundary = outVec[0].outer_boundary();
		}
		else{
			// Add all non intersecting obstacles 
			Polygon tempPoly;
			converFromCGALtoVisilibity(ecPoly, tempPoly, false);
			polyVec.push_back(tempPoly);
		}
	}

	// Insert the boundary 
	Polygon boundaryPoly;
	converFromCGALtoVisilibity(boundary, boundaryPoly, true);
	polyVec.insert(polyVec.begin(), boundaryPoly);

	// Build environment and visibility graph	
	m_pEnvironment = new Environment(polyVec);
	m_pVisibilityGraph = new Visibility_Graph(*m_pEnvironment, m_epsilon);
}

double Roadmap::computeShortestPath(Point_2& p1, Point_2& p2, std::list<Point_2> & path){
	ICPoint_2 p1x = K_ICK_converter(p1);
	ICPoint_2 p2x = K_ICK_converter(p2);
	return computeShortestPath(p1x.x(), p1x.y(), p2x.x(), p2x.y(), path);
}

double Roadmap::computeShortestPath(double x1, double y1, double x2, double y2, std::list<Point_2> & path){
	Polyline pl = m_pEnvironment->shortest_path(Point(x1, y1), Point(x2, y2), *m_pVisibilityGraph, m_epsilon);
	for(int i = 0; i < pl.size(); i++){
		path.push_back(Point_2(pl[i].x(), pl[i].y()));
	}
	return getPathLength(path);
}

double Roadmap::computeShortestPath(double x1, double y1, double x2, double y2, vector<pair<double, double> >& path){
	// Compute shortest path
	std::list<Point_2> shortestPath;
	double length = computeShortestPath(x1, y1, x2, y2, shortestPath);

	// Convert the path
	for(std::list<Point_2>::iterator vi = shortestPath.begin(); vi != shortestPath.end(); vi++){
		ICPoint_2 p = K_ICK_converter(*vi);
		path.push_back(pair<double, double>(p.x(), p.y()));
	}
	return length;
}

void Roadmap::addToScene(QGraphicsScene& scene, bool drawEdge, QPen edgePen, bool drawVertex, QPen vertexPen){
	// Paint the edges
	if(drawEdge){
		std::set<std::pair<int, int> > edgeSet = m_graph.getEdgeSet();
		for(std::set<std::pair<int, int> >::iterator eit = edgeSet.begin(); eit != edgeSet.end(); eit++){
			int v1 = (*eit).first;
			int v2 = (*eit).second;
			ICPoint_2 p1 = K_ICK_converter(m_vidPointMap[v1]);
			ICPoint_2 p2 = K_ICK_converter(m_vidPointMap[v2]);
			scene.addLine(p1.x(), p1.y(), p2.x(), p2.y(), edgePen);
			/*}*/

			if(m_boundaryBoundingCycle.hasEdge(v1, v2)){
				scene.addLine(p1.x(), p1.y(), p2.x(), p2.y(), QPen(Qt::blue, 0.05, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
			}
			else{
				for(std::vector<Graph*>::iterator git = m_obsBoundingCycleVec.begin(); git != m_obsBoundingCycleVec.end(); git++){
					if((*git)->hasEdge(v1, v2)){
						scene.addLine(p1.x(), p1.y(), p2.x(), p2.y(), QPen(Qt::yellow, 0.05, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
					}
				}
			}
		}
	}

	// Paint vertices if needed and obtain the fill area
	if(!drawVertex){
		for(std::list<Point_2>::iterator it = m_pointList.begin(); it != m_pointList.end(); it ++){
			ICPoint_2 p = K_ICK_converter(*it);
			scene.addEllipse(p.x() - 0.025, p.y() - 0.025, 0.05, 0.05, vertexPen);
		}
	}

	for(Path_Map::iterator pathIt = m_connectingPathMap.begin(); pathIt != m_connectingPathMap.end(); pathIt ++){
		std::list<Point_2>& shortestPath = pathIt->second;
		std::list<Point_2>::iterator vit = shortestPath.begin();
		if(vit != shortestPath.end()){
			ICPoint_2 p = K_ICK_converter(*vit);
			vit++;
			for(; vit != shortestPath.end(); vit++){
				ICPoint_2 p2 = K_ICK_converter(*vit);
				scene.addLine(p.x(), p.y(), p2.x(), p2.y(), QPen(Qt::magenta, 0.025, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
				p = p2;
			}
		}
	}

	// std::cout << shortestPath.size() << std::endl;

}

void Roadmap::getIntersectingEdges(Polygon_2 & poly, Graph& boundingCycle, bool outerBoundary){
	// Detecting all edges of the lattice graph that intersects the boundary of poly

	// =====================================================================================
	// Get a vertex of the poly and find the hexgaon of the lattice that contains the vetex.
	// For each col and row combo, there are three possible hexgaons the point may fall into

	// Locate the hexgon (relative to some arbitrary base choice)
	ICPoint_2 p0 = K_ICK_converter(poly.vertex(0));
	std::pair<int, int> crp = locateHexagonForPoint(p0);
	int col = crp.first;
	int row = crp.second;

	// Compute the current hexagon
	Point_2 p[6];
	Segment_2 e[6];
	Polygon_2 hex; 
	int pIndex[6];
	populateHexagon(col, row, p, e, hex, pIndex);

	// =====================================================================================
	// With a point of the obstacle and a hexgon containing the point, we iteratively check 
	// intersections between the segments of the obstacle and the hexagon. We assume that 
	// the obstacle would have at least three vertices. Note that a obstacle may be smaller than
	// a single hexagon. 

	std::set<std::pair<int, int> > tempEdgeToRemoveSet;
	int numVertices = poly.size();		// 
	int currentVertexIndex = 0;
	while(currentVertexIndex < numVertices){
		// Get the current segment of polygon to be checked
		Point_2 nextVertex = poly[currentVertexIndex==numVertices-1?0:currentVertexIndex+1];

		// Check whether nextVertex is outside of the current hexagon. If we are still in the 
		// same hexagon, then move to the next obstacle vertex
		if(hex.bounded_side(nextVertex) == CGAL::ON_BOUNDED_SIDE) {
			currentVertexIndex ++;
			continue;
		}

		// If we are here, then we jumped outside of a lattice hexagon. Figure out which edge is 
		// being intersected. We keep doing this until we cover the entire obstacle edge
		Point_2 currentVertex = poly[currentVertexIndex];
		Segment_2 obsEdge(currentVertex, nextVertex);

		std::pair<int, int> lastPair(-10, -10);
		while(true){
			int edgeIndex = 0;
			while(edgeIndex < 6){
				if(intersection(e[edgeIndex], obsEdge)){
					std::pair<int, int> edgePair(pIndex[edgeIndex], pIndex[edgeIndex==5?0:edgeIndex+1]);
					std::pair<int, int> edgePairReverse(pIndex[edgeIndex==5?0:edgeIndex+1], pIndex[edgeIndex]);
					if(lastPair != edgePair && lastPair != edgePairReverse){
						lastPair = edgePair;
						break;
					}
				}
				edgeIndex ++;
			};

			if(edgeIndex < 6){
				// We hit an intersection, mark the edge as to be removed
				tempEdgeToRemoveSet.insert(lastPair);
				m_edgeToBeRemovedSet.insert(lastPair);

				// Figure out the next hexgon to be checked 
				if(col%2 == 0){
					switch(edgeIndex){
						case 0: col--; break;
						case 1: row++; break;
						case 2: col++; break;
						case 3: col++; row--; break;
						case 4: row--; break;
						case 5: col--; row--; break;
						default: break;
					}
				}
				else{
					switch(edgeIndex){
						case 0: col--; row++; break;
						case 1: row++; break;
						case 2: col++; row++; break;
						case 3: col++; break;
						case 4: row--; break;
						case 5: col--; break;
						default: break;
					}
				}
				populateHexagon(col, row, p, e, hex, pIndex);
			}
			else{
				// No more intersections, we are dong with the current obstacle edge
				break;
			}
		}
		currentVertexIndex++;
	}

	// =====================================================================================
	// We now have all the edges of the lattice that lie on the polygonal obstacle boundary. 
	// Next, we locate the bounding cycle in the full lattice graph surrounding the obstacle
	if(tempEdgeToRemoveSet.size() > 0){
		// Iterate through edges to be removed and find one with an end vertex that 
		// has an associated edge in the configuraiton space
		int sIndex, t1Index, t2Index; 
		for(std::set<std::pair<int,int> >::iterator eit = tempEdgeToRemoveSet.begin();
		eit != tempEdgeToRemoveSet.end(); eit++){
			std::pair<int, int> edge = *eit;
			// Try one vertex
			sIndex = edge.first;
			if((t2Index = pointBelongToEdgeOutideObstacle(poly, sIndex, tempEdgeToRemoveSet, outerBoundary))!= -1)
			{
				t1Index = edge.second;
				break;
			}

			// Try the other vertex
			sIndex = edge.second;
			if((t2Index = pointBelongToEdgeOutideObstacle(poly, sIndex, tempEdgeToRemoveSet, outerBoundary))!= -1)
			{
				t1Index = edge.first;
				break;
			}

		}

		// Now edge (sIndex, t1Index) crosses the boundary, and (sIndex, t2Index) it outside the boundary 
		// We simply locate a hexagon containing these two edges to extend the cycle 
		Point_2 p1 = m_vidPointMap[sIndex]; 
		Point_2 p2 = m_vidPointMap[t1Index]; 
		Point_2 p3 = m_vidPointMap[t2Index];
		ICPoint_2 midPoint = K_ICK_converter(Point_2((p2.x() + p3.x())/2, (p2.y() + p3.y())/2));
		crp = locateHexagonForPoint(midPoint);
		populateHexagon(crp.first, crp.second, p, e, hex, pIndex);

		int endIndex = sIndex;		// When we see this index again, we are done

		// Add valid edge to cycle
		addEdgeIfNotThere(sIndex, t2Index, boundingCycle);

		// With the first hexagon, we can keep going along the boundary 
		while(true){
			// Get next candidate index
			int nextIndex = getNextNodeInSequence(sIndex, t2Index, pIndex);
			
			p1 = m_vidPointMap[sIndex]; 
			p2 = m_vidPointMap[nextIndex]; 
			p3 = m_vidPointMap[t2Index];

			// Check whether the edge (t2Index, nextIndex) crosses boundary
			if(edgeInSet(t2Index, nextIndex, tempEdgeToRemoveSet)){
				// In this case, we move to the next hexagon
				crp = getBorderHexagon(crp.first, crp.second, nextIndex, t2Index, pIndex);
				populateHexagon(crp.first, crp.second, p, e, hex, pIndex);
				sIndex = nextIndex;
				// t2Index = sIndex;
			}
			else{
				// Edge is valid, update indices
				sIndex = t2Index;
				t2Index = nextIndex;

				// Add valid edge to cycle
				addEdgeIfNotThere(sIndex, t2Index, boundingCycle);
			}

			if(endIndex == t2Index) break;
		}
	}

	// =====================================================================================
	// Process bounding cycle to remove single degree vertices
	std::set<int> vSet = boundingCycle.getVertexSet();
	std::set<int> vSingleSet;
	// Collect single degree vertices
	for(std::set<int>::iterator vit = vSet.begin(); vit != vSet.end(); vit++){
		if(boundingCycle.getNeighborSet(*vit).size() == 1){
			vSingleSet.insert(*vit);
			// Need to check the neighbor as well
			int prevNbr = *vit;
			int nbr = *(boundingCycle.getNeighborSet(*vit).begin());
			while(boundingCycle.getNeighborSet(nbr).size()==2){
				vSingleSet.insert(nbr);
				std::set<int> nbrSet = boundingCycle.getNeighborSet(nbr);
				nbrSet.erase(prevNbr);

				prevNbr = nbr;
				nbr = *(nbrSet.begin());
			}
		}
	}
	// Delete them
	for(std::set<int>::iterator vit = vSingleSet.begin(); vit != vSingleSet.end(); vit++){
		boundingCycle.removeVertex(*vit);

		// It seems good to remove these single degree vertices 
		if(!outerBoundary)m_vertexToBeRemovedSet.insert(*vit);
	}
}

std::pair<int, int> Roadmap::getBorderHexagon(int col, int row, int i1, int i2, int pIndex[]){
	// First locate the index of i1
	int i1Index = 0;
	for(int i = 0; i < 6; i ++){
		if(pIndex[i] == i1){
			i1Index = i;
			break;
		}
	}

	// Figure out the edge index
	int edgeIndex = 0;
	int nextIndex = (i1Index == 5? 0 : i1Index + 1);
	if(pIndex[nextIndex] == i2){
		// Going clockwise
		edgeIndex = i1Index;
	}
	else{
		// Going counterclockwise
		edgeIndex = (i1Index + 6 - 1)%6;
	}

	// Locate the next hexagon
	if(col%2 == 0){
		switch(edgeIndex){
			case 0: col--; break;
			case 1: row++; break;
			case 2: col++; break;
			case 3: col++; row--; break;
			case 4: row--; break;
			case 5: col--; row--; break;
			default: break;
		}
	}
	else{
		switch(edgeIndex){
			case 0: col--; row++; break;
			case 1: row++; break;
			case 2: col++; row++; break;
			case 3: col++; break;
			case 4: row--; break;
			case 5: col--; break;
			default: break;
		}
	}
	return std::pair<int, int>(col, row);
}

bool Roadmap::edgeInSet(int v1, int v2, std::set<std::pair<int, int> >& edgeSet){
	std::pair<int, int> edge(v1, v2);
	std::pair<int, int> rEdge(v2, v1);
	return (edgeSet.find(edge) != edgeSet.end() || edgeSet.find(rEdge) != edgeSet.end());
}

void Roadmap::addEdgeIfNotThere(int v1, int v2, Graph & graph){
	if(!graph.hasEdge(v1, v2)){
		graph.addEdge(v1, v2);
	}
}

int Roadmap::getNextNodeInSequence(int i1, int i2, int pIndex[]){
	// First locate the index of i1
	int i1Index = 0;
	for(int i = 0; i < 6; i ++){
		if(pIndex[i] == i1){
			i1Index = i;
			break;
		}
	}
	int nextIndex = (i1Index == 5? 0 : i1Index + 1);
	if(pIndex[nextIndex] == i2){
		// Going clockwise
		return pIndex[(i1Index + 2)%6];
	}
	else{
		// Going counterclockwise
		return pIndex[(6 + i1Index - 2)%6];
	}
}

std::pair<int, int> Roadmap::locateHexagonForPoint(ICPoint_2 &p0){
	// Compute the rectangular (col, row) 
	int col = (int)(floor((p0.x() - xs)/(m_edgeLength*1.5)));
	int row = (int)(floor((p0.y() - ys)/(m_edgeLength*sqrt3)));

	// Shift everything to the "origin"
	double dx = p0.x() - (xs + m_edgeLength*1.5*col);
	double dy = p0.y() - (ys + m_edgeLength*sqrt3*row);

	// For each col/row combo, which correspond to a rectangular region, there can be
	// three hexagons corrsponding to a point in that rectangle
	if(col%2 == 0){
		if(dx <= 0.5*m_edgeLength && dy <= sqrt3*m_edgeLength  - sqrt3*dx && dy >= sqrt3*dx){
			// We need to move left
			col --;
			// row ++;
		}
		else if(dy >= m_edgeLength*0.5*sqrt3){
			// Off by one row
			row ++;
		}
		else{
			// Already in the right place, do nothing
		}
	}
	else{
		if(dy >= sqrt3*0.5*m_edgeLength + sqrt3*dx){
			col--;
			row++;
		}
		else if(dy <= sqrt3*0.5*m_edgeLength - sqrt3*dx){
			col--;
		}
		else{
			// Already in the right place, do nothing
		}
	}
	return std::pair<int, int>(col, row);
}

int Roadmap::pointBelongToEdgeOutideObstacle(Polygon_2 & poly, int pIndex, 
	std::set<std::pair<int,int> > &tempEdgeToRemoveSet, bool outerBoundary){
	// Grab the point and check whether it is in the c-space 
	if(pointOutsideObstacle(poly, pIndex, outerBoundary)){
		// The point is not "in" the obstacle, get edges and check that they are not 
		// intersecting with the boundary. Because pIndex is outside obstacle, as long
		// as one edge from pIndex is not intersecting boundary, the edge must be outside 
		// the boundary
		std::set<int>& nbrSet = m_graph.getNeighborSet(pIndex);
		for(std::set<int>::iterator nit = nbrSet.begin(); nit != nbrSet.end(); nit++){
			int et = *nit;
			if(!edgeInSet(pIndex, et, tempEdgeToRemoveSet)){
				return et;
			}
		}
	}
	return -1;
}

bool Roadmap::pointOutsideObstacle(Polygon_2 & poly, int pIndex, bool outerBoundary){
	Point_2 p = m_vidPointMap[pIndex];
	int boundedStatus = poly.bounded_side(p);
	if((boundedStatus == CGAL::ON_BOUNDED_SIDE && outerBoundary) ||
		(boundedStatus == CGAL::ON_UNBOUNDED_SIDE && !outerBoundary)){
		return true;
	}
	return false;
}

void Roadmap::populateHexagon(int col, int row, Point_2* p, Segment_2* e, Polygon_2 &hex, int* pIndex){
	hex.clear();
	if(col%2 == 0){
		pIndex[0] = (col + row*n_w)*2;
		pIndex[1] = (col + row*n_w)*2 + 1;
		pIndex[2] = (col + 1 + row*n_w)*2 + 1;
		pIndex[3] = (col + 1 + row*n_w)*2;
		pIndex[4] = (col + 1 + (row-1)*n_w)*2 + 1;
		pIndex[5] = (col + (row-1)*n_w)*2 + 1;
	}
	else{
		pIndex[0] = (col + row*n_w)*2 + 1;
		pIndex[1] = (col + (row+1)*n_w)*2;
		pIndex[2] = (col + 1 + (row+1)*n_w)*2;
		pIndex[3] = (col + 1 + row*n_w)*2 + 1;
		pIndex[4] = (col + 1 + row*n_w)*2;
		pIndex[5] = (col + row*n_w)*2;
	}
	// Then we get all vertices
	for(int i = 0; i < 6; i ++){
		p[i] = m_vidPointMap[pIndex[i]];
	}

	// Then the edges and the hexgaon
	for(int i = 0; i < 6; i ++){
		e[i] = Segment_2(p[i], p[i == 5? 0 : i+1]);
		hex.push_back(p[i]);
	}
}

void Roadmap::buildHexgaonLattice(){

	// Start the lattice at bottomLeftX, bottomLeftY. 
	// Assume that there are n_w and n_h hexgons inside the rectangle, then we should have
	// (1/2)*sideLength + n_w*sideLength*3/2 <= width and width < (1/2)*sideLength + (n_w + 1)*sideLength*3/2 
	// (sqrt(3)/2)*sideLength + n_h*(sqrt(3)/2)*sideLength <= height and height < (sqrt(3)/2)*sideLength + (n_h + 1)*(sqrt(3)/2)*sideLength
	// From these we can compute n_w and n_h. We add a few to make sure that we cover everything. 

	// Compute lattice nodes
	for(int i = 0; i < n_w; i ++){
		for(int j = 0; j < n_h; j ++){
			// Build one vertical "wave" of vertices
			Point_2 v0(xs + i*m_edgeLength*1.5 + (i%2 == 0? 0 : m_edgeLength*0.5), 
				ys + j*m_edgeLength*sqrt3);
			Point_2 v1(xs + m_edgeLength/2 + i*m_edgeLength*1.5 + (i%2 == 0? 0 : - m_edgeLength*0.5), 
				ys + sqrt3*m_edgeLength/2 + j*m_edgeLength*sqrt3);

			m_vidPointMap[(i + j*n_w)*2] = v0;
			m_pointVidMap[v0] = (i + j*n_w)*2;
			m_vidPointMap[(i + j*n_w)*2 + 1] = v1;
			m_pointVidMap[v1] = (i + j*n_w)*2 + 1;

			m_pointList.push_back(v0);
			m_pointList.push_back(v1);
		}
	}

	// Build adjacency, for each vertex, check whether its three neighbors are present
	for(int id = 0; id < n_w*n_h*2; id ++){
		// Compute w, h
		int w = (id/2)%n_w;
		int h = (id/2)/n_w;
		bool odd = (id%2==1);

		// If odd is true, check the vertex above, which has index (w, h + 1, 0)
		if(odd){
			if(h + 1 < n_h){
				m_graph.addEdge(id, (w + (h + 1)*n_w)*2);
			}
		}
		else{
			m_graph.addEdge(id, (w + (h)*n_w)*2 + 1);
		}
		// If odd is true and w is even, check (w + 1, h, 1)
		if(odd && w%2 == 0){
			if(w + 1 < n_w){
				m_graph.addEdge(id, (w + 1 + h*n_w)*2 + 1);
			}
		}

		// If odd is false and w is odd, check (w + 1, h, 0)
		if(odd == false && w%2 == 1){
			if(w + 1 < n_w){
				m_graph.addEdge(id, (w + 1 + h*n_w)*2);
			}
		}
	}
}

void Roadmap::drawHexagonLattice(QGraphicsScene& scene, bool drawVetex){
	QPen edgePen = QPen(Qt::gray, 0.25, Qt::DashLine, Qt::RoundCap, Qt::RoundJoin);
	QPen vertexPen = QPen(Qt::blue, 0.4, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
	std::set<std::pair<int, int> > edgeSet = m_graph.getEdgeSet();
	std::set<int> vSet;
	for(std::set<std::pair<int, int> >::iterator eit = edgeSet.begin(); eit != edgeSet.end(); eit++){
		int v1 = (*eit).first;
		int v2 = (*eit).second;
		vSet.insert(v1); vSet.insert(v2);
		ICPoint_2 p1 = K_ICK_converter(m_vidPointMap[v1]);
		ICPoint_2 p2 = K_ICK_converter(m_vidPointMap[v2]);
		scene.addLine(p1.x(), p1.y(), p2.x(), p2.y(), edgePen);
	}

	// Paint vertices if needed and obtain the fill area
	if(drawVetex){
		for(std::set<int>::iterator vit = vSet.begin(); vit != vSet.end(); vit ++){
			ICPoint_2 p = K_ICK_converter(m_vidPointMap[*vit]);
			scene.addEllipse(p.x() - 0.25, p.y() - 0.25, 0.5, 0.5, vertexPen);
		}
	}
}

void Roadmap::snapToGraph(vector<pair<double, double> >&coords, vector<int>&snapped){
	set<int> usedVertices;
	// For each member in coords, locate the hexagon it belongs and then locate the closest 
	// graph vertex to the said coordinate. 
	for(int r = 0; r < coords.size(); r ++){
		double x = coords[r].first;
		double y = coords[r].second;
		Point_2 p0 = Point_2(x, y);

		ICPoint_2 p0x = K_ICK_converter(p0);
		std::pair<int, int> crp = locateHexagonForPoint(p0x);
		int col = crp.first;
		int row = crp.second;

		// Compute the current hexagon
		Point_2 p[6];
		Segment_2 e[6];
		Polygon_2 hex; 
		int pIndex[6];
		populateHexagon(col, row, p, e, hex, pIndex);

		// Find the closest vertex that is not occupied and inside the configuration space
		int bestV = -1;
		double minDist = -1;
		for(int vi = 0; vi < 6; vi ++){
			// Only work with vertices inside the configuration space and not used
			if(m_graph.hasVertex(pIndex[vi]) 
			   && usedVertices.find(pIndex[vi]) == usedVertices.end()){
				if(minDist < 0){
					minDist = getDistance(p0, p[vi]);
					bestV = pIndex[vi];
				}
				else{
					double dist = getDistance(p0, p[vi]);
					if(minDist > dist){
						minDist = dist;
						bestV = pIndex[vi];
					}
				}
			}
		}

		// Snap!
		if(bestV != -1){
#ifdef _DEBUG
			cout << "Snapping point (" << x << ", " << y << ") to vertex " << m_vidGFMap[bestV] << endl;
#endif
			snapped.push_back(m_vidGFMap[bestV]);
			usedVertices.insert(bestV);
		}
		else{
			throw "ERROR: Cannot find suitable vertex for snapping";
		}
	}

}

void Roadmap::createRandomStartGoalPairs(int numRobots, double spacing,  vector<pair<double, double> >& starts, 
	vector<pair<double, double> >& goals){

	createRandomCoords(numRobots, spacing, starts);
	createRandomCoords(numRobots, spacing, goals);
}

void Roadmap::createRandomCoords(int numRobots, double spacing, vector<pair<double, double> >& coords){
	// Uniform sample from the free space and pick 5 start/goal pairs
	CGAL::Bbox_2 bbox = m_pBoundingRect->bbox();
	while(true){
		long trialCount = 0;
		for(int i = 0; i < numRobots; i ++){
			while(trialCount < numRobots*30000){
				trialCount ++;
#ifdef _DEBUG
				if(trialCount %1000 == 0) {cout << trialCount << endl;}
#endif
				double x = (rand()%10000)/10000.*(bbox.xmax() - bbox.xmin())  + bbox.xmin();
				double y = (rand()%10000)/10000.*(bbox.ymax() - bbox.ymin())  + bbox.ymin();
#ifdef _DEBUG
				cout << "Randomly created point with x=" << x <<", y=" << y;
#endif
				Point_2 p(x, y);
				// Check that the point is in free configuration space
				if(isPointInCSpace(p)){
					bool good = true;
					// Check that the point has good distance from existing points
					for(int vi = 0; vi < coords.size(); vi ++){
						if(getDistance(p, coords[vi].first, coords[vi].second) < (spacing > 2 ? spacing : 2)*m_radius){
							good = false;
						}
					}

					if(good){
						// No problem! We have a good set
						coords.push_back(pair<double, double>(x, y));
#ifdef _DEBUG
						cout << " - added to set for robot " << i << endl;
#endif
						break;
					}
					else{
#ifdef _DEBUG
						cout << endl; 
#endif
					}

				}
				else{
#ifdef _DEBUG
					cout << endl; 
#endif
				}
			}
			if(trialCount >= numRobots*30000){
				// cout << coords.size() << endl;
				coords.clear();
				break;
			}
		}
		if(coords.size() == numRobots)return;
	}
}


bool Roadmap::solveProblem(vector<pair<int, int> >& sgVec, map<int, vector<int> >& paths,
	string& fileFolder, string& fileNameExtra){
	// Solve the problem by calling external java solver
	
	// =====================================================================================
	// First write the problem 
	string pfString = fileFolder + "\\p" + fileNameExtra + ".txt";
	string sfString = fileFolder + "\\s" + fileNameExtra + ".txt";
	ofstream ps(pfString.c_str(), ios::out);

	// Write the the number of robots
	ps << m_finalGraph.getVertexSet().size() << endl;

	// Write out all the edges
	set<pair<int, int> >& edgeSet = m_finalGraph.getEdgeSet();
	for(set<pair<int, int> >::iterator ei = edgeSet.begin(); ei != edgeSet.end(); ei ++){
		ps << ei->first << ":" << ei->second << " ";
	}
	ps << endl;

	// Write out the starts and then the goals
	for(int r = 0; r < sgVec.size(); r ++){
		ps << sgVec[r].first << " ";
	}
	ps << endl;
	for(int r = 0; r < sgVec.size(); r ++){
		ps << sgVec[r].second << " ";
	}
	ps << endl;
	ps.close();

	// =====================================================================================
	// Make system call 
	string callStr = "java -cp gurobi.jar;mp.jar projects.multipath.general.algorithms.Main";
	if(fileNameExtra.length() > 0) {
		callStr.append(" ").append(pfString).append(" ").append(sfString);
	}
	int i = system(callStr.c_str());

	// =====================================================================================
	// Read in the solution, there should be sgVec.size() robots
	ifstream ss(sfString.c_str(), ios::in);
	if(ss.good()){
		for(int r = 0; r < sgVec.size(); r++){
			string line; 
			getline(ss, line);
			QString s(line.c_str());
			QStringList qsl = s.split(" ");
			paths[r] = vector<int>();
			for(int t = 0; t < qsl.size(); t ++){
				paths[r].push_back(qsl[t].toInt());
			}
			paths[r].pop_back();
		}
		ss.close();
	}
	else{
		ss.close();
		return false;
	}

	// Solve the problem locally
	// ILPSolver solver(&m_finalGraph);
	// solver.solve(sgVec, paths, -1);


	improvePaths(paths);
	return true;
}

void Roadmap::improvePaths(map<int, vector<int> >& paths){
	// Remove obvious ossilation from the paths
	for(int t = 1; t < paths[0].size() - 1; t ++){
		for(int r = 0; r < paths.size(); r ++){
			// Do we have single step ossilation?
			if(paths[r][t - 1] == paths[r][t + 1] && paths[r][t - 1] != paths[r][t]){
				// Check whether any other robots goes to paths[r][t - 1] at t
				bool conflict = false;
				for(int ori = 0; ori < paths.size(); ori ++){
					if(ori == r) continue;
					if(paths[ori][t] == paths[r][t - 1]){
						conflict = true;
					}
				}

				// If not, let the robot stay
				if(!conflict){
					paths[r][t] = paths[r][t - 1];
				}
			}
			// Do a two step look ahead as well
			else if((t < paths[0].size() - 2) && paths[r][t - 1] == paths[r][t + 2] && 
				paths[r][t] == paths[r][t + 1] && paths[r][t - 1] != paths[r][t]){
				// Check whether any other robots goes to paths[r][t - 1] at t
				bool conflict = false;
				for(int ori = 0; ori < paths.size(); ori ++){
					if(ori == r) continue;
					if(paths[ori][t] == paths[r][t - 1] || paths[ori][t + 1] == paths[r][t - 1]){
						conflict = true;
					}
				}

				// If not, let the robot stay
				if(!conflict){
					paths[r][t+1] = paths[r][t] = paths[r][t - 1];
				}
			}
		}
	}
}

pair<double, double> Roadmap::getVertexLocationFromID(int vid){
	// Get vid in the original graph
	vid = m_vidFGMap[vid];

	// Locate the vid
	Point_2& p = m_vidPointMap[vid];
	ICPoint_2 pp = K_ICK_converter(p);

	return pair<double, double>(pp.x(), pp.y());
}