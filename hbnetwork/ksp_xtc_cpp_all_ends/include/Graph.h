#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <vector>
#include <stack>

#include "gmxcpp/Utils.h"
// This header file is required for creating a Trajectory object
#include "Trajectory.h"
// This header file is required for some of the command line things used below
#include "CommandLine.h"
// This Topology file is required for creating a Trajectory object
#include "Topology.h"

//using namespace std;

namespace Graph_all
{
	class Graph
	{
	public:
		typedef stack<unsigned> path_t;			// the top is the first node
		typedef vector<path_t>  path_list;		// the paths are sorted by ASC (ascending)
		typedef vector< vector<double> > dag_t; // the two dimension array storing a DAG(: directed acyclic graph)

		/// Input array: there must be only a start node!
		///	if the value of one element of array < 0, indicates there are no arc.
		Graph();
		Graph( dag_t& array );
		Graph( double **array , unsigned N );
		~Graph();

		/// Reset the directed acyclic graph (DAG).
		void Restart( dag_t& array );
		void Restart( double **array , unsigned N );

		/// Find the shortest path from the named start node to the named end node.
		/// Parameter <start> - start node.
		/// Parameter <end> - end node. 
		/// Make sure <end> != <start>.			
		double Dijkstra( unsigned start , unsigned end , path_t& path );

		/// Find the k shortest paths (KSP) from node 0 to N-1.						
		/// Martins' Algorithm (deletion algorithm) Implementation.
		/// Parameter <paths> return all the shortest paths.			
		/// If fails, return 0, or return the real number of all the shortest paths.								
		int MSAforKSP( unsigned start, unsigned end ,unsigned k, path_list& kpaths, int* path_length );

		/// Output the content of "_array" for debug.
		void Output( ostream& out = cerr );
		double dijkstra( unsigned start , unsigned end, int* paths , double* dists );
		double dijkstra( int* paths , unsigned start = 0 , unsigned end = 0 );
		int dijkstra_all_ends(int start, int* end ,int num_ends,  path_list& kpaths, int* path_length);
		vector<int> prims_MST();
	private:
		/// Default is to compute the shortest distance from node 0 to node N-1.

		/// Add a node to graph.
		/// Return the number of new node.
		unsigned addNode( unsigned ni , int preni );	

	private:
		unsigned _N;	 // original size of "_array", it's fixed.
		unsigned _size;	 // size of "_array", because the "_array" maybe be reallocated. 
		dag_t    _array; // dynamic two dimension array
			
	};
}

#endif	// GRAPH_H

