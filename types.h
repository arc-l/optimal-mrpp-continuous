/*
 *	Type definitions
 *
 *  Created on: Jan 30, 2015
 *  Author: Jingjin Yu
 */

#ifndef _O_CGAL_TYPES_H_
#define _O_CGAL_TYPES_H_

#define _CRT_SECURE_NO_WARNINGS 1
#include "number_type.h"

#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polygon_set_2.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/minkowski_sum_2.h>
#include <CGAL/Point_set_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Vector_2.h>
#include <CGAL/General_polygon_with_holes_2.h>

#include <CGAL/Sweep_line_2_algorithms.h>
#include <CGAL/squared_distance_2.h>

#include <boost/foreach.hpp>

// Kernel used for computation
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Segment_2 Segment_2;
typedef std::list<Segment_2> Seg2_list;
typedef K::Line_2 Line_2;
typedef K::Vector_2 Vector_2;

typedef CGAL::Polygon_2<K>  Polygon_2;
typedef std::list<Polygon_2> Polygon2_list;
typedef CGAL::Polygon_with_holes_2<K> Polygon_with_holes_2;
typedef CGAL::Polygon_set_2<K> Polygon_set_2;

typedef CGAL::Point_set_2<K>::Vertex_handle Vertex_handle;

// Kernel used for extracting some intermediate data
typedef CGAL::Exact_predicates_inexact_constructions_kernel ICK;
typedef ICK::Point_2 ICPoint_2;
typedef CGAL::Polygon_2<ICK> ICPolygon_2;


typedef CGAL::Exact_predicates_exact_constructions_kernel EK;
typedef EK::Point_2 ECPoint_2;
typedef EK::Segment_2 ECSegment_2;
typedef CGAL::Polygon_2<EK> ECPolygon_2;
typedef CGAL::Polygon_with_holes_2<EK> ECPolygon_with_holes_2;


#endif //_OPTMLY_CGAL_TYPES_H_
