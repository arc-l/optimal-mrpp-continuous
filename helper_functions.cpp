/*
 *	Helper function implementations
 *
 *  Created on: Jan 30, 2015
 *  Author: Jingjin Yu
 */
#include "helper_functions.h"

void populateApproximateDisc(Polygon_2 &poly, Point_2 &center, double radius, int segments){
	for(int i = 0; i < segments; i ++){
		poly.push_back (Point_2 ( center.x() + radius*cos(i*PI*2/segments), center.y() + radius*sin(i*PI*2/segments)));
	}
}

Polygon_2 growPolygonByRadius(Polygon_2 &poly, double radius, int segments){
	// Get the disk for computing Minkowski sum
	Polygon_2 disc;
    Point_2 p(0, 0);
    populateApproximateDisc(disc, p, radius, segments);

	// Do computation
	if(!poly.is_clockwise_oriented()){
		poly.reverse_orientation();
	}
	return minkowski_sum_2 (poly, disc).outer_boundary();
}

double getPathLength(std::list<Point_2> path){
	double length = 0;
	std::list<Point_2>::iterator vit = path.begin();
	CGAL::Cartesian_converter<K,ICK> K_ICK_converter; 
	if(vit != path.end()){
		ICPoint_2 p = K_ICK_converter(*vit);
		vit++;
		for(; vit != path.end(); vit++){
			ICPoint_2 p2 = K_ICK_converter(*vit);
			length += sqrt((p2.x() - p.x())*(p2.x() - p.x()) + (p2.y() - p.y())*(p2.y() - p.y()));
			p = p2;
		}
	}
	return length;
}

double getDistance(Point_2& p1x, Point_2& p2x){
	static CGAL::Cartesian_converter<K,ICK> K_ICK_converter; 
	ICPoint_2 p1 = K_ICK_converter(p1x);
	ICPoint_2 p2 = K_ICK_converter(p2x);
	return sqrt((p2.x() - p1.x())*(p2.x() - p1.x()) + (p2.y() - p1.y())*(p2.y() - p1.y()));
}

double getDistance(Point_2& p, double x, double y){
	static CGAL::Cartesian_converter<K,ICK> K_ICK_converter; 
	ICPoint_2 pp = K_ICK_converter(p);
	return sqrt((pp.x() - x)*(pp.x() - x) + (pp.y() - y)*(pp.y() - y));
}

void getAllPairsShortestPath (Graph* pGraph, map<int, map<int, int> > &dist){
	int V = pGraph->getVertexSet().size();

    /* dist[][] will be the output matrix that will finally have the shortest 
      distances between every pair of vertices */
    int i, j, k;
 
    /* Initialize the solution matrix same as input graph matrix. Or 
       we can say the initial values of shortest distances are based
       on shortest paths considering no intermediate vertex. */
    for (i = 0; i < V; i++)
	{
        for (j = 0; j < V; j++){
            dist[i][j] = (pGraph->hasEdge(i, j) ? 1 : (i == j?0:100000000));
		}
	}
 
    /* Add all vertices one by one to the set of intermediate vertices.
      ---> Before start of a iteration, we have shortest distances between all
      pairs of vertices such that the shortest distances consider only the
      vertices in set {0, 1, 2, .. k-1} as intermediate vertices.
      ----> After the end of a iteration, vertex no. k is added to the set of
      intermediate vertices and the set becomes {0, 1, 2, .. k} */
    for (k = 0; k < V; k++)
    {
        // Pick all vertices as source one by one
        for (i = 0; i < V; i++)
        {
            // Pick all vertices as destination for the
            // above picked source
            for (j = 0; j < V; j++)
            {
                // If vertex k is on the shortest path from
                // i to j, then update the value of dist[i][j]
                if (dist[i][k] + dist[k][j] < dist[i][j])
                    dist[i][j] = dist[i][k] + dist[k][j];
            }
        }
    }
}

void converFromCGALtoVisilibity(ECPolygon_2 &poly, Polygon &visiPoly, bool ccw){
	// Reverse orientation as needed
	if((poly.is_counterclockwise_oriented() && ccw == false) || 
	   (!poly.is_counterclockwise_oriented() && ccw)){
		poly.reverse_orientation();
	}

	// Get the vertices and convert
	static CGAL::Cartesian_converter<EK,ICK> EK_ICK_converter; 
	vector<Point> vp;
	for(ECPolygon_2::Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); vi++){
		ICPoint_2 p = EK_ICK_converter(*vi);
		vp.push_back(Point(p.x(), p.y()));
	}

	visiPoly.set_vertices(vp);
}

bool boundaryInterset(ECPolygon_2 &poly1, ECPolygon_2 &poly2){
	// Check edge by edge
	for(int i = 0; i < poly1.size(); i++){
		ECSegment_2 seg1 = poly1.edge(i);
		for(int j = 0; j < poly2.size(); j ++){
			ECSegment_2 seg2 = poly2.edge(j);
			if(CGAL::do_intersect(seg1, seg2)){return true;}
		}
	}
	return false;
}

ECPoint_2 convertToExactPoint(Point_2 &p){
	static CGAL::Cartesian_converter<K,ICK> K_ICK_converter; 
	ICPoint_2 pp = K_ICK_converter(p);
	return ECPoint_2(pp.x(), pp.y());
}


ECPolygon_2 convertToExactPolygon(Polygon_2 &poly){
	ECPolygon_2 ecPoly;
	for(Polygon_2::Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); vi ++){
		static CGAL::Cartesian_converter<K,ICK> K_ICK_converter; 
		ICPoint_2 pp = K_ICK_converter(*vi);
		ecPoly.push_back(ECPoint_2(pp.x(), pp.y()));
	}
	return ecPoly;
}
