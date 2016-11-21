/*
 *  A basic adjacency list based simple, undirected graph in which each vertex is indexed
 *	with an integer ID. 
 *
 *  Created on: Jan 30, 2015
 *  Author: Jingjin Yu
 */

#ifndef _O_GRAPH_H_
#define _O_GRAPH_H_

#include "types.h"

#include <QColor>
#include <QPen>
#include <QBrush>
#include <QGraphicsScene>

class Graph{
private:
	std::set<int>						vSet;
	std::map<int, std::set<int> >		adjSetMap;
	std::set<std::pair<int, int> > 		edgeSet;
	std::set<int>						emptySet;

public:
	// =====================================================================================
	//									Core methods
	// =====================================================================================

	Graph(){}

	void clear(){
		vSet.clear();
		adjSetMap.clear();
		edgeSet.clear();
		emptySet.clear();
	}

	bool hasVertex(int v){
		return (vSet.find(v) != vSet.end());
	}

	std::set<int> &getVertexSet(){return vSet;}

	void removeVertex(int v){
		// Delete all edges from the vertex
		std::set<int> nbrSet = getNeighborSet(v);
		for(std::set<int>::iterator vit = nbrSet.begin(); vit != nbrSet.end(); vit++){
			removeEdge(v, *vit);
		}

		// Erase the vertex itself
		vSet.erase(v);
	}

	void addEdge(int v1, int v2){
		// Need to add vertex?
		if(vSet.find(v1) == vSet.end()){vSet.insert(v1);adjSetMap[v1] = std::set<int>();}
		if(vSet.find(v2) == vSet.end()){vSet.insert(v2);adjSetMap[v2] = std::set<int>();}
		// Adding edges as needed
		adjSetMap[v1].insert(v2);
		adjSetMap[v2].insert(v1);
		edgeSet.insert(std::pair<int, int>(v1 < v2? v1 : v2, v1 < v2? v2 : v1));
	}

	std::set<std::pair<int, int> > &getEdgeSet(){return edgeSet;}

	void removeEdge(int v1, int v2){
		// Check that we have the vertices
		if(vSet.find(v1) == vSet.end() || vSet.find(v2) == vSet.end()) return;

		// Remove edges as needed
		adjSetMap[v1].erase(v2);
		adjSetMap[v2].erase(v1);
		edgeSet.erase(std::pair<int, int>(v1 < v2? v1 : v2, v1 < v2? v2 : v1));

		// Remove v1, v2 as needed
		if(adjSetMap[v1].size() == 0){vSet.erase(v1); adjSetMap.erase(v1);}
		if(adjSetMap[v2].size() == 0){vSet.erase(v2); adjSetMap.erase(v2);}
	}

	bool hasEdge(int v1, int v2){
		// Check we have vertex
		if(vSet.find(v1)!=vSet.end()){
			return (adjSetMap[v1].find(v2) != adjSetMap[v1].end());
		}
		return false;
	}

	std::set<int>& getNeighborSet(int v){
		if(vSet.find(v)!=vSet.end()){
			return adjSetMap[v];
		}
		return emptySet;
	}



	// =====================================================================================
	//									Auxiliary methods
	// =====================================================================================

	// If the graph is a cycle, this function returns the vertices of the cycle 
	// ordered sequentially
	void getCycleVertexVector(std::vector<int>& vVec){
		if(vSet.size() == 0) return;
		std::set<int> processedVertexSet;
		// Start from a vertex and retrive other vertices in order. We assume that
		// the cycle is big enough
		int vStart = *(vSet.begin());
		vVec.push_back(vStart);
		processedVertexSet.insert(vStart);
		int vNext = *(getNeighborSet(vStart).begin());
		while(vSet.size() > vVec.size()){
			vVec.push_back(vNext);
			processedVertexSet.insert(vNext);
			std::set<int>::iterator vit = getNeighborSet(vNext).begin();
			vNext = *vit++;
			if(processedVertexSet.find(vNext) != processedVertexSet.end()){
				vNext = *vit;
			}
		}
	}


};

#endif //_O_GRAPH_H_
