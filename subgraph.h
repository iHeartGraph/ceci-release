/*
 * Copyright 2019 The George Washington University
 * Written by Bibek Bhattarai
 * https://www.seas.gwu.edu/~howie/
 * Contact: bhattarai_b@gwu.edu
 * Please cite the following paper:
 *  
 * Bibek Bhattarai, Hang Liu, H. Howie Huang 2019. CECI: Compact Embedding Cluster Index for Scalabel Subgraph Matching. SIGMOD 2019.
 *
 * This file is part of CECI.
 * CECI is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CECI is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CECI.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef TURBOISO_H
#define TURBOISO_H
#include <iostream>
#include <functional>
#include <vector>
#include <cstring>
#include <cstdint>
#include <map>
#include <algorithm>

#include "graph.h"

// Some flags for checking the status of the nodes
#define VISITED   true
#define UNVISITED false
#define ACCEPTED  true
#define REJECTED  false

#define WORKTAG 1
#define DIETAG 0

//consts
const int MAX_QUERY_NODE = 64;
extern int automorphId[MAX_QUERY_NODE];

using namespace std;
typedef bool node_status_t;

/*typedef enum{
	NODE_UV = 0,
	NODE_VA,
	NODE_VR
} node_status_t;
*/

// Data Structures
struct queryNode{
	public:
		int id;
		int label;
		int vertexId;

		int parent;
		int childCount;
		int ntEdgeCount;
		int childList[MAX_QUERY_NODE];
		int ntEdges[MAX_QUERY_NODE];

		queryNode(){
			childCount = 0;
			ntEdgeCount = 0;
		}
};

//Query graph based variables
//uint8_t  query_label[MAX_QUERY_NODE];
//uint8_t vertexId[MAX_QUERY_NODE];
//bool 	 visited[MAX_QUERY_NODE];
//uint8_t  level[MAX_QUERY_NODE];
//uint8_t	 parent[MAX_QUERY_NODE];



// Thread Status block, this block is used as status in multithreaaded backtracking
class thread_status{
	public:
		uint32_t mybeg, myend;
		bool go;
		int *embedding;
		uint8_t curr_query_node;

	public:
		thread_status()
		{
			mybeg = 0;
			myend = 0;
			go = true;
			embedding = NULL;
			curr_query_node = 0;
		}
		void init(int qsize)
		{
			embedding  = new int[qsize];
			memset(embedding, 0, qsize*sizeof(int));
		}
		
		~thread_status()
		{
	//		delete[] isComplete;
	//		delete[] isStarted;
	//		delete[] embedding;
		}
};

typedef std::vector<int> vecNode;
typedef std::pair<int, std::vector<int> > Data;
//typedef std::vector< Data > vecEdge;

typedef std::vector< Data >::iterator edgeCandIterator;

// Class for comparing two vector entries according to key
class DataCompare
{
	public:
		// Comparison function for sorting
		bool operator()(const Data& lhs, const Data& rhs) const
		{
			return keyLess(lhs.first, rhs.first);
		}

		// COmparison fuunctions for lookup
		bool operator()(const Data& lhs, const Data::first_type& k) const
		{
			return keyLess(lhs.first, k);
		}

		bool operator()(const Data::first_type& k, const Data& rhs) const
		{
			return keyLess(k, rhs.first);
		}
	private:
		bool keyLess(const Data::first_type& k1, const Data::first_type& k2) const
		{
			return k1 < k2;
		}
};

typedef std::pair<double, int> Score;
class ScoreCompareRev
{
	public:
		// comparison function for descending sorting
		bool operator()(const Score& lhs, const Score& rhs) const
		{
			return keyMore(lhs.first, rhs.first);
		}

		// comparison functions for lookup, Although I am not using them here
		bool operator()(const Score& lhs, const Score::first_type& k) const
		{
			return keyMore(lhs.first, k);
		}

		bool operator()(const Score::first_type& k, const Score& rhs) const
		{
			return keyMore(k, rhs.first);
		}
	private:
		bool keyMore(const Score::first_type& k1, const Score::first_type& k2) const
		{
			return k1 > k2;
		}
};

typedef std::pair<int, int> Degree;
class SortbyDegree
{
	public:
		bool operator()(const Degree& lhs, const Degree& rhs) const
		{
			return keyMore(lhs.first, rhs.first);
		}
		bool operator()(const Degree& lhs, const Degree::first_type& k) const
		{
			return keyMore(lhs.first, k);
		}
		bool operator()(const Degree::first_type& k, const Degree& rhs) const
		{
			return keyMore(k, rhs.first);
		}
	private:
		bool keyMore(const Degree::first_type& k1, const Degree::first_type& k2) const
		{
			return k1 > k2;
		}
};

//candidate region representation
struct CRI{
	unsigned long long int *cardinality;
	std::vector<int> node_cands;
	std::vector<std::pair<int, std::vector<int> > >edge_cands;
	std::map<int, std::vector<Data> >nte_cands;
};


class subgraph{
	private:
		// @OpenMP things
		int num_thrds = 56;
		int subregion_count;

		// graph file names
		const char* dataGraphFileName;
		const char* queryGraphFileName;

		graph *dataGraph;
		graph *queryGraph;

		// Querry Graph tree levelwise Representation
		int* level_beg_pos;
		int tree_depth = 0;


		// Statistic counters 
		unsigned long long int sufficientNumberOfEmbeddings;
		unsigned long long int totalNumberOfEmbeddings;
		unsigned long long int totalNumberOfRecursion;

		unsigned long long int sufficientEmbeddingsPerThread;
		unsigned long long int *numberofEmbeddings;
		unsigned long long int *numberofRecursiveCalls;

		//For single threaded solution
		unsigned long long int numberOfEmbeddings_sngl;
		unsigned long long int numberOfRecursiveCalls_sngl;
		
		// Status flags
		node_status_t* visit_status;
		node_status_t* acceptance_status;
		node_status_t* go_status;
		uint8_t* 	   parent_count;

		// Candidate Region Index
		CRI* myCR;
	
		// Query graph spanning Tree
		int startQueryVertex;
		std::vector<queryNode> queryTree;
		int* visitTreeOrder;

		// Query Matching Order/ visit Order
		uint8_t* queryMatchingSequence; // vertexID, ParentID
		uint8_t* queryMatchingSequence_inv;

		int* qGraphtoTreeMap;

		// Thread control/ Status Block
		thread_status *myTCB;
		thread_status myTCB_sngl;

		// to store each embeddings
		std::vector<int*> Embedding_list;

		unsigned long long int threshold;
		unsigned long long int totalwork;

		std::vector<int> preset;
		std::vector<std::vector<int> > work_units;

		std::vector<int> vm;
		std::vector<std::pair<double, int> > workload;
		double* similarity;

	public:
		subgraph(const char* dataGraphFile, const char* queryGraphFile, int num_thrds, bool break_auto, int embedding_count = 1000);
	//	~subgraph();
	

		//utility  functions
		void execute();
		void clean();
		void clearCR();

		// turboIso query processing algorithm
		void genericQueryProc();
		void automorphFreeProc();
		void distributedQueryProc();


		int chooseStartVertex();
		bool NLCFilter(int u, int v);
		bool NTEFilter(int u, int v, int pNode);
		void generateQueryTree();

		void exploreGraph();

		// single threaded implementation
		void exploreGraph_sngl();
		bool exploreCR(int u, vector<int>& vm, int parentDFS);
		void computeMatchingOrder();
		void divideNoTreeEdges(int u, int numberOfNoTreeEdges, float * queryVertexScore);
		void getOrderByBFSScore(int u, float * queryVertexScore, vector<int> & nextVertex);


		//update and reset state
//		void updateState(int u, int v, map<int, int> embedding, map<int, int> inverseEmbedding);
		void subgraphSearch(int thid);
		void subgraphSearch_topk(int thid);

		// single threaded implementation
		void subgraphSearch_sngl(unsigned long long int cardinalityParent);

		bool isJoinable(queryNode& u, int v, int thid);
//		void restoreState(int u, int v, map<int, int> embedding, map<int, int> inverseEmbedding);

		//getter and setter functions
		void showEmbeddings(int thid);
		long totalComputation();
		long totalEmbeddings();
		long totalRecursiveCalls();
		int  numberOfDataGraphs();
		int  numberOfQueryGraphs();
	
		//heuristic based functions 
		void calculateNumThread();
		void workLoadEstimator();
		void calculateSimilarityMatrix();
		std::vector<int> intersection(std::vector<int> &v1, std::vector<int> &v2);

		// Added functions for spanning tree
		int getLeadNode(int index);
		int addToTree(std::vector<int>&, std::vector<float>&, std::vector<int>&, int node, std::vector<int>&);
		int maxEdgeRank(std::vector<float>&, std::vector<float>&);
		void CreateSpanningTree();
		bool allVerticesAdded(std::vector<int>&);
};

#endif
