/*
 *  A GraphicsItem derived class to enable the filling of the inside of a polygonal region. 
 *	Qt GraphicsItem class does not do filling by default (at least as of 4.8.5). 
 *
 *  Created on: Jan 30, 2015
 *  Author: Jingjin Yu
 */

#ifndef _O_ADV_GRAPHICS_ITEM_H_
#define _O_ADV_GRAPHICS_ITEM_H_

#include "types.h"
#include "roadmap.h"

#include <CGAL/Qt/PolygonGraphicsItem.h>
#include <CGAL/Qt/Converter.h>

#include <QColor>
#include <QPen>
#include <QBrush>

template<class T> 
class AdvancedGraphicsItem: public CGAL::Qt::PolygonGraphicsItem<T> {
public:
	bool m_bShowVertices;		// Draw vertices?
	bool m_bShowEdge;			// Draw edges?
	bool m_bFill;				// Fill the polygon?

	QPen m_vertexPen;			// Vertex pen
	QPen m_edgePen;				// Edge pen
	QBrush m_fillBrush;			// Fill brush

	// Need this constructor to take a T* pointer
	AdvancedGraphicsItem(T* p);

	// We need to overwrite paint method to do custom rendering such as filling the area. This
	// is not provided by default by QGraphicsItem in any other way.
	void paint( QPainter * painter, const QStyleOptionGraphicsItem * option, QWidget * widget);
};

template<class T> AdvancedGraphicsItem<T>
::AdvancedGraphicsItem(T* p): CGAL::Qt::PolygonGraphicsItem<T>(p),
		m_bShowVertices(false),m_bShowEdge(true), m_bFill(true),
		m_vertexPen(QPen(Qt::red, 0.025, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin)),
		m_edgePen(QPen(Qt::red, 0.025, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin)),
		m_fillBrush(QColor(Qt::gray)){}


// We need to overwrite paint method to do custom rendering such as filling the area. This
// is not provided by default by QGraphicsItem in any other way.
template<class T>  void AdvancedGraphicsItem<T>
::paint( QPainter * painter, const QStyleOptionGraphicsItem * option, QWidget * widget){

	// Paint the edges
	if(m_bShowEdge){
		painter->setPen(m_edgePen);
		CGAL::Qt::PolygonGraphicsItem<T>::painterostream  = CGAL::Qt::PainterOstream<typename T::Traits>(painter);
		if(CGAL::Qt::PolygonGraphicsItem<T>::drawEdges()) {
			for(typename T::Edge_const_iterator eit = CGAL::Qt::PolygonGraphicsItem<T>::poly->edges_begin();
				eit != CGAL::Qt::PolygonGraphicsItem<T>::poly->edges_end();
				++eit){
					CGAL::Qt::PolygonGraphicsItem<T>::painterostream << *eit;
			}
		}
	}

	// Paint vertices if needed and obtain the fill area
	if(m_bShowVertices || m_bFill){
		CGAL::Qt::Converter<typename T::Traits> convert;
		painter->setPen(m_vertexPen);
		QMatrix matrix = painter->matrix();
		painter->resetMatrix();
		QPainterPath path;
		for(typename T::Vertex_iterator it = CGAL::Qt::PolygonGraphicsItem<T>::poly->vertices_begin();
			it != CGAL::Qt::PolygonGraphicsItem<T>::poly->vertices_end();
			it++){
				QPointF point = matrix.map(convert(*it));
				// Draw vertices here as needed
				if(m_bShowVertices){
					painter->drawPoint(point);
				}
				if(it == CGAL::Qt::PolygonGraphicsItem<T>::poly->vertices_begin()){
					path.moveTo(point.rx(), point.ry());
				}
				else{
					path.lineTo(point.rx(), point.ry());
				}
		}
		path.closeSubpath();
		// Fill the area if needed
		if(m_bFill){
			painter->fillPath(path, QBrush(m_fillBrush));
		}
	}
}

#endif /* _O_ADV_GRAPHICS_ITEM_H_ */
