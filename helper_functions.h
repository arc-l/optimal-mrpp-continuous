/*
 *  Some helper functions for geometric computing
 *
 *  Created on: Jan 30, 2015
 *  Author: Jingjin Yu
 */

#ifndef _O_CGAL_HELPER_H_
#define _O_CGAL_HELPER_H_

#include "types.h"
#include "graph.h"
#include "shortest_path/visilibity.hpp"
#include <map>

using namespace std;
using namespace VisiLibity;

// Populate the polygon as an approximate disc with segments nuber of sides
void populateApproximateDisc(Polygon_2 &poly, Point_2 &center, double radius, int segments = 18);

// Compute Minkowski sum of a polygon with a disc 
Polygon_2 growPolygonByRadius(Polygon_2 &poly, double radius, int segments = 18);

// Compute the distance between two points
double getDistance(Point_2& p1, Point_2& p2);
double getDistance(Point_2& p, double x, double y);

// Compute path length
double getPathLength(std::list<Point_2> path);

// All pairs shortest path computation
void getAllPairsShortestPath (Graph* pGraph, map<int, map<int, int> > &dist);

// Convert a CGAL polygon to a VisiLibity polygon
void converFromCGALtoVisilibity(ECPolygon_2 &poly, Polygon &visiPoly, bool ccw = true);

// Does two polygons have boundary edges that intersect?
bool boundaryInterset(ECPolygon_2 &poly1, ECPolygon_2 &poly2);

// Convert to exact contruction
ECPoint_2 convertToExactPoint(Point_2 &p);
ECPolygon_2 convertToExactPolygon(Polygon_2 &poly);

#endif //_O_CGAL_HELPER_H_
