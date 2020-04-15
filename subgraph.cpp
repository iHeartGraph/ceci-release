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

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <queue>
#include <limits.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <omp.h>
//#include <mpi.h>

#include "subgraph.h"
#include "wtime.h"

bool break_automorph;
int automorphId[MAX_QUERY_NODE];

using namespace std;

subgraph::subgraph(const char* dataGraphFile, const char* queryGraphFile, int count_limit, bool break_auto, int embedding_count){
	dataGraphFileName = dataGraphFile;
	queryGraphFileName = queryGraphFile;
	
	break_automorph = break_auto;
	// Form the data structures
	dataGraph = new graph(dataGraphFileName);
//	dataGraph->test();
	queryGraph = new graph(queryGraphFileName);
//	queryGraph->test();

	// Build inverted label index for graph
	dataGraph->buildLabelVertexList();
	dataGraph->buildVertexLabelVertexList();
	//dataGraph->build_sorted_csr();

	queryGraph->buildLabelVertexList();
	queryGraph->buildVertexLabelVertexList();
	queryGraph->build_sorted_csr();
	
	// Embedding Infomartions
	totalNumberOfRecursion = 0;
	totalNumberOfEmbeddings = 0;
	sufficientNumberOfEmbeddings = embedding_count;
	num_thrds = count_limit;

	numberOfRecursiveCalls_sngl = 0;
	numberOfEmbeddings_sngl = 0;

	// INFO: whether a node is visited or not
	visit_status = new node_status_t[dataGraph->vert_count];
	
	// INFO: whether the visited vertex passed candidate test or not 
	acceptance_status = new node_status_t[dataGraph->vert_count];

	// INFO: count the no of vertex that are adjacent to this one and is candidate of parent
	parent_count = new uint8_t[dataGraph->vert_count];

	// Query Matching Sequence 
	visitTreeOrder				= new int[queryGraph->vert_count];
	queryMatchingSequence 		= new uint8_t[queryGraph->vert_count];
	queryMatchingSequence_inv 	= new uint8_t[queryGraph->vert_count];
	qGraphtoTreeMap 			= new int[queryGraph->vert_count];

	// Counters for search statistics
	numberofEmbeddings 			= new unsigned long long int[num_thrds];
	numberofRecursiveCalls 		= new unsigned long long  int[num_thrds];
	memset(numberofEmbeddings, 0, num_thrds*sizeof(unsigned long long int));
	memset(numberofRecursiveCalls, 0, num_thrds*sizeof(unsigned long long  int));
	
	myTCB = new thread_status[num_thrds]; 
	for(int i = 0; i< num_thrds; ++i){
		myTCB[i].init(queryGraph->vert_count);
	}
	//MPI_Barrier(MPI_COMM_WORLD);
}

/*
void subgraph::execute(){
//	for(int dataGraphIndex = 0; dataGraphIndex < dataGraphVector.size(); ++dataGraphIndex){
//		dataGraph = &dataGraphVector[dataGraphIndex];
//
//		for(int queryGraphIndex = 0; queryGraphIndex < queryGraphVector.size(); ++queryGraphIndex){
//			queryGraph = &queryGraphVector[queryGraphIndex];
//			std::cout << "the number of query graph vertices is" << queryGraph->getNumberOfVertexes() << std::endl;
//			std::cout << "the number of query graph edges is" << queryGraph->getNumberOfEdges() << std::endl;
			genericQueryProc();
//			cout << queryGraphIndex << " : finish one query with embeddings : " << newEmbeddingsForEachQuery << endl;
//		}
//		cout << dataGraphIndex << " finish one data graph ******* " << endl;
//	}
}
*/

/*
// Distributed Query Processing Unit
void subgraph::distributedQueryProc()
{
	//printf("The query tree is generated\n");
	clean();
	double time;

	//Find the node to get started with 
	startQueryVertex = chooseStartVertex();
	if(startQueryVertex == -1){
		return;	
	}

	generateQueryTree();


	myCR = new CRI[queryTree.size()];
	for(int qnode = 0; qnode < queryTree.size(); ++qnode){
		myCR[qnode].node_cands.reserve(10000);
		myCR[qnode].edge_cands.reserve(10000);
		myCR[qnode].cardinality = new unsigned long long int[dataGraph->vert_count];
		std::fill_n(myCR[qnode].cardinality, dataGraph->vert_count, 1);
	}

	vector<int>* startVertexCandidates = &(dataGraph->getLabelVertexList()->find(queryGraph->csr_label[startQueryVertex])->second);
	vector<int>::iterator startVertexCandidateIterator = startVertexCandidates->begin();

	vm.clear();

	int world_size, world_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	for(; startVertexCandidateIterator != startVertexCandidates->end(); ++startVertexCandidateIterator){
		// Degree filter
		if(queryGraph->degree[queryTree[0].vertexId] > dataGraph->degree[*startVertexCandidateIterator]){
			continue;
		}
		// neighborhood label count filter
		if(!NLCFilter(queryTree[0].vertexId, *startVertexCandidateIterator)){
			continue;
		}
		// passes all the filters, add to the list of nodes to be traversed
		vm.push_back(*startVertexCandidateIterator);
	}
	//explore the portion of the graph
	myCR[0].node_cands.clear();
	for(std::vector<int>::iterator cItr = vm.begin(); cItr != vm.end(); ++cItr){
		myCR[0].node_cands.push_back(*cItr);
	}
	
	time = wtime();
	exploreGraph_sngl();

//	if(world_rank == 0){
//		printf("The subregion exploration on Worker %d has completed on %f sec \n", world_rank, wtime()-time);
//	}
	for(int i = 0; i< queryTree.size(); ++i){
		queryMatchingSequence[i] = queryTree[i].id;
		queryMatchingSequence_inv[queryTree[i].id] = i;
	}

//	visit order check
//	printf("The query Matching Sequence is: ");
//	for(int i = 0; i < queryTree.size(); ++i){
//		printf(" %d", queryMatchingSequence[i]);
//	}
//	printf("\n");

	subregion_count = myCR[queryMatchingSequence[0]].node_cands.size();
	int result = 0;
	myTCB_sngl.curr_query_node = 0;
	myTCB_sngl.init(queryTree.size());
	//printf("The total subregion count on worker is %d\n", myCR[0].node_cands.size());
	// Master node, prepares work, distribute it, and collect the result
	if(world_rank == 0){
		time= wtime();
		// single threaded search operation
		std::vector<std::pair<int,int> > works;
		works.reserve(myCR[0].node_cands.size());
		threshold = 0;
		std::vector<int> residue;
		for(std::vector<int>::iterator itr = myCR[0].node_cands.begin(); itr != myCR[0].node_cands.end(); ++itr){
			if(myCR[0].cardinality[*itr] != 0){
				works.push_back(std::pair<int, int>(myCR[0].cardinality[*itr], *itr));
				threshold += myCR[0].cardinality[*itr];
				//printf("%llu\n", myCR[0].cardinality[*itr]);
			}
		}

		std::sort(works.begin(), works.end(),SortbyDegree());
		threshold = threshold/(2*world_size);

		//creating works
	//	time = wtime();
		work_units.reserve(2 * myCR[0].node_cands.size());

		for(int i = 0; i< works.size(); ++i){
			preset.clear();
			if(works[i].first > threshold){
				preset.push_back(works[i].second);
				myTCB_sngl.embedding[queryMatchingSequence[0]] = works[i].second;
				myTCB_sngl.curr_query_node++;
				subgraphSearch_sngl(myCR[0].cardinality[works[i].second]);
				myTCB_sngl.curr_query_node--;
				continue;
			}
			preset.push_back(works[i].second);
			work_units.push_back(preset);
		}
		//printf("number of works %d\t %d\n",works.size(), work_units.size());
		printf("Work prepared in %f seconds\n", wtime()-time);

		std::vector<int> curr_work;
		MPI_Status status;
		// Send the first batch of work
		for(int rank = 1; rank < world_size; ++rank){
			curr_work = work_units[rank-1];
			//printf("work size%d \n", curr_work.size());
			//printf("%d\t%d\n", curr_work[0], curr_work[1]);
			MPI_Send(&curr_work[0], curr_work.size(), MPI_INT, rank, WORKTAG, MPI_COMM_WORLD);
		}

		//while there is still work,
		//send them to slaves one by one
		//count total by receiving result from each slave
		int count = world_size-1;
		curr_work = work_units[count];
		while(count < work_units.size()-1){
			count++;
			MPI_Recv(&result, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			//printf("Received back count of %d\n", result);
			//printf("Sending work of size %d: \t", curr_work.size());
		//	for(int i = 0; i< curr_work.size(); ++i){printf("%d\t", curr_work[i]);}
		//	printf("%llu\n", numberOfEmbeddings_sngl);

			numberOfEmbeddings_sngl += result;

			MPI_Send(&curr_work[0],curr_work.size() , MPI_INT, status.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD);
			curr_work = work_units[count];
		}

		//done working, kill all the slaves
		for(int rank = 1; rank < world_size; ++rank){
			MPI_Send(&count, 1, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
		}
		printf("Final subgraph matching completed at %f Seconds\n", wtime()-time);
		printf("Total Number of Embeddings found is : %llu\n\n",numberOfEmbeddings_sngl);
	}
	else{//slave nodes; receive, do, send, and ask
		std::vector<int> mywork;
		mywork.reserve(10);
		int work_size;
		int thid = 0;
		MPI_Status status;
		for(;;){
			//printf("LADO!!!\n");
			numberOfEmbeddings_sngl = 0;
			//Receive the work
			MPI_Recv(&mywork[0], 10, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			MPI_Get_count(&status, MPI_INT, &work_size);
			if(status.MPI_TAG == DIETAG){
				return;
			}
			// finish the work
			for(int i = 0; i < work_size ; ++i){
		//		printf("work %d is %d \n", i,mywork[i]);
			//	printf("work is %d %d\n", mywork[i], myTCB[thid].curr_query_node);
				myTCB[thid].embedding[queryMatchingSequence[i]] = mywork[i];
				myTCB[thid].curr_query_node++;
			}
			subgraphSearch(thid);
			
			for(int i = 0; i < work_size; ++i){
				myTCB[thid].curr_query_node--;
			}
			result = numberOfEmbeddings_sngl;
			//printf("%d\n", result);

			// send the result back to master and ask for more works
			MPI_Send(&result, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}
	}
}
*/

// Single Machine query procesing engine
void subgraph::genericQueryProc(){
	clean();
	double time;
//	printf("The automorph groups are \n");
//	for(int i  = 0; i < queryGraph->vert_count; ++i){
//		printf("%d ", automorphId[i]);
//	}
//	printf("\n");
	// Find the node to get started with
	startQueryVertex = chooseStartVertex();//implemented
	std::cout << "The starting query Vertex is " << startQueryVertex << endl;
	if(startQueryVertex == -1){
		return;
	}
	
	// rewrite the querygraph as a BFS Tree 
	generateQueryTree(); //implemented

	// allocate the space for candidate region, reserve reasonable space in advance to avoid reallocation over and over again 
	myCR = new CRI[queryTree.size()];
	for(int qnode = 0; qnode < queryTree.size(); ++qnode){
		myCR[qnode].node_cands.reserve(10000);
		myCR[qnode].edge_cands.reserve(10000);
		myCR[qnode].cardinality = new unsigned long long int[dataGraph->vert_count];
		std::fill_n(myCR[qnode].cardinality, dataGraph->vert_count, 1);
	}
	
	// Finding the candidates for start vertex
	vector<int> startVertexCandidates;
	vm.clear();
	// if the label is specific
	if(queryGraph->csr_label[startQueryVertex] != -1){
		startVertexCandidates = dataGraph->getLabelVertexList()->find(queryGraph->csr_label[startQueryVertex])->second;
		vector<int>::iterator startVertexCandidateIterator = startVertexCandidates.begin();
		for(; startVertexCandidateIterator != startVertexCandidates.end(); ++startVertexCandidateIterator){
			// Degree filter
			if(queryGraph->degree[queryTree[0].vertexId] > dataGraph->degree[*startVertexCandidateIterator]){
				continue;
			}
	
			// neighborhood label count filter
			if(!NLCFilter(queryTree[0].vertexId, *startVertexCandidateIterator)){
				continue;
			}
		
			// passes all the filters, add to the list of nodes to be traversed
			vm.push_back(*startVertexCandidateIterator);

		}
	}
	else{
		for(int it = 0; it != dataGraph->vert_count; ++it){
			// Degree filter
			if(queryGraph->degree[queryTree[0].vertexId] > dataGraph->degree[it]){
				continue;
			}
			
			// neighborhood label count filter
			if(!NLCFilter(queryTree[0].vertexId, it)){
				continue;
			}
				
			// passes all the filters, add to the list of nodes to be traversed
			vm.push_back(it);	
		}
	}

	// std::cout << "The vm size " << vm.size() << std::endl;
	//myDistributor = new distributor(vm, dataGraph, queryGraph, startQueryVertex, 64);
	//	workLoadEstimator();
	// Explore the whole graph
	myCR[0].node_cands.clear();
	for(std::vector<int>::iterator cItr = vm.begin(); cItr != vm.end(); ++cItr){

		if(numberOfEmbeddings_sngl >= sufficientNumberOfEmbeddings){
			return;
		}

		//dont need no map, just use this vector
		myCR[0].node_cands.push_back(*cItr);
	}
	time = wtime();
	exploreGraph_sngl();
	//exploreGraph();

//	cout << "The index creation completed at:  " << wtime()-time <<" Seconds" << endl;

	//TODO: implement this function
	//computeMatchingOrder();
	for(int i = 0; i< queryTree.size(); ++i){
		queryMatchingSequence[i] = queryTree[i].id;
		queryMatchingSequence_inv[queryTree[i].id] = i;
	}

	// visit order check
	printf("The query Matching Sequence is: ");
	for(int i = 0; i < queryTree.size(); ++i){
		printf(" %d", queryMatchingSequence[i]);
	}
	printf("\n");

	subregion_count = myCR[queryMatchingSequence[0]].node_cands.size();
//	printf("The total subregion count is %d\n", subregion_count);

	// single threaded search operation
	myTCB_sngl.curr_query_node = 0;
	myTCB_sngl.init(queryTree.size());
	
	/* 
	 * Execution Plan:
	 * put the different embedding/Partial solution for each of the threads
	 * for Each of the subregion, maintain the current status i.e, how many nodes are visited for each query node
	 * record if a subregion is finished
	 * for the first pass, assign equal no of embedding for each node,
	 * if some subregion is finished before completing that no of embedding, go to next subregion
	 * ISSUES:
	 * 		What if some subregion is really huge and all other are very small?
	 * 		What if a thread finishes before collecting the desired number of embeddings?
	 *		What if total no of embedding is less than the desired number?? Dude, exit
	 *		How to find the optimal number of thread necessary for given number of embeddings??
	 *			--> Find the threshold number of embedding for one thread
	 *			--> assign at least that number of thread to the single node
	 */
/*
	// Sort the candidates for first node by cardinality in dcescending order
	// So that biggest works can be started early
	std::vector<std::pair<int,int> > works;
	works.reserve(myCR[0].node_cands.size());
	threshold = 0;
	std::vector<int> residue;
	for(std::vector<int>::iterator itr = myCR[0].node_cands.begin(); itr != myCR[0].node_cands.end(); ++itr){
		if(myCR[0].cardinality[*itr] != 0){
			works.push_back(std::pair<int, int>(myCR[0].cardinality[*itr], *itr));
			threshold += myCR[0].cardinality[*itr];
		//	printf("%llu\n", myCR[0].cardinality[*itr]);
		}
	}
	std::sort(works.begin(), works.end(),SortbyDegree());
	threshold = threshold/(22*num_thrds);

	//creating works
	time = wtime();
	work_units.reserve(2 * myCR[0].node_cands.size());
	for(int i = 0; i< works.size(); ++i){
		preset.clear();
		if(works[i].first > threshold){
			preset.push_back(works[i].second);
			myTCB_sngl.embedding[queryMatchingSequence[0]] = works[i].second;
			myTCB_sngl.curr_query_node++;
			subgraphSearch_sngl(myCR[0].cardinality[works[i].second]);
			myTCB_sngl.curr_query_node--;
			continue;
		}
		preset.push_back(works[i].second);
		work_units.push_back(preset);
	}
*/
	//printf("number of works %d\t %d\n",works.size(), work_units.size());
//	printf("Work prepared in %f seconds\n", wtime()-time);
	time = wtime();
#pragma omp parallel num_threads(num_thrds)
	{
		double xtime = omp_get_wtime();
		//int myCardinality = 0;
	//#pragma omp for schedule(dynamic) nowait
		for(uint32_t cnode = 0; cnode < myCR[0].node_cands.size(); ++cnode)
		{
			int thid = omp_get_thread_num();
			//std::vector<int> mywork = work_units[cnode];
			//for(int i = 0; i < mywork.size(); ++i){
				//printf("work %d is %d \n", i,mywork[i]);
			myTCB[thid].embedding[queryMatchingSequence[0]] = myCR[0].node_cands[cnode];
			myTCB[thid].curr_query_node++;
			//}
			subgraphSearch(thid);
			
			//for(int i = 0; i < mywork.size(); ++i){
				myTCB[thid].curr_query_node--;
			//}
			if(totalEmbeddings() >= 100000){
				break;
			}
		}
	//	printf("%d %llu %f\n", omp_get_thread_num(), omp_get_wtime()-xtime);
	}

	printf("Final subgraph matching completed at %f Seconds\n", wtime()-time);	
	printf("Total Number of Embeddings found is : %llu\n\n",totalEmbeddings());
/*	assert(0);
	//printf("The average is %llu\n", threshold);
	time = wtime();
#pragma omp parallel num_threads(num_thrds)
	{
		double xtime =  omp_get_wtime();
		int myCardinality = 0;
		// get the biggest ones; check if they are bigger than average
		// If yes treat them separately: for all thread, call the backtracking function
		// pass true to break flag; break the program dynamically in next level
		// after it is done, treat the remaining ones
	#pragma omp for schedule(dynamic, 1) nowait
		for(uint32_t cnode = 0; cnode < works.size(); ++cnode)
		{
			if(works[cnode].first > threshold){
				#pragma omp critical
				residue.push_back(works[cnode].second);
				continue;
			}
			int thid = omp_get_thread_num();
		//	int match_data_vert = myCR[0].candi dates[mybeg];
		//	if(myTCB[thid].go == false){
		//		return;
		//	}
			myTCB[thid].embedding[queryMatchingSequence[0]] = works[cnode].second;
			//myCardinality += works[cnode].first;
			myTCB[thid].curr_query_node++;

			subgraphSearch(thid);
			
			myTCB[thid].curr_query_node--;
		}
		//printf("%d %llu %f\n", omp_get_thread_num(), myCardinality, omp_get_wtime()-xtime);
	}
	for(uint32_t cnode = 0; cnode < residue.size(); ++cnode){
		for(int thid = 0; thid < num_thrds; ++thid){
			myTCB[thid].embedding[queryMatchingSequence[0]] = residue[cnode];
			myTCB[thid].curr_query_node++;
		}
		
		myTCB_sngl.embedding[queryMatchingSequence[0]] = residue[cnode];
		myTCB_sngl.curr_query_node++;
		subgraphSearch_sngl(myCR[0].cardinality[residue[cnode]]);
		myTCB_sngl.curr_query_node--;

		for(int thid = 0; thid < num_thrds; ++thid){
			myTCB[thid].curr_query_node--;
		}
	}
	printf("Final subgraph matching completed at %f Seconds\n", wtime()-time);
	
	printf("Total Number of Embeddings found is : %llu\n\n",totalEmbeddings());*/
//	printf("Estimations: \n");
//	for(int x = 0; x < num_thrds; ++x){
//		printf("%d  ", myestimation[x]);
//	}
//	printf("Actual Values: \n");
//	totalRecursiveCalls();
}

void subgraph::clearCR()
{
	for(int queryNode = 0; queryNode < queryTree.size(); ++queryNode){
		myCR[queryNode].node_cands.clear();
		myCR[queryNode].edge_cands.clear();
		myCR[queryNode].nte_cands.clear();
		std::fill_n(myCR[queryNode].cardinality, dataGraph->vert_count, 1);
	}
}

void subgraph::automorphFreeProc(){
	clean();
	double time;
	
	//CreateSpanningTree();
	
	// Find the node to get started with
	startQueryVertex = chooseStartVertex();
	if(startQueryVertex == -1){
		return;
	}
	//std::cout << "The start node is " << startQueryVertex << std::endl;	
	// rewrite the querygraph as a BFS Tree 
	//std::cout <<"Size of the query tree is "<< queryTree.size() << std::endl;
	generateQueryTree(); 

	//std::cout <<"Size of the query tree is "<< queryTree.size() << std::endl;

	// allocate the space for candidate region, reserve reasonable space in advance to avoid reallocation over and over again 
	myCR = new CRI[queryTree.size()];
	for(int qnode = 0; qnode < queryTree.size(); ++qnode){
		myCR[qnode].node_cands.reserve(10000);
		myCR[qnode].edge_cands.reserve(10000);
		myCR[qnode].cardinality = new unsigned long long int[dataGraph->vert_count];
		std::fill_n(myCR[qnode].cardinality, dataGraph->vert_count, 1);
	}
	
	// Finding the candidates for start vertex 
	vector<int>* startVertexCandidates = &(dataGraph->getLabelVertexList()->find(queryGraph->csr_label[startQueryVertex])->second);
	vector<int>::iterator startVertexCandidateIterator = startVertexCandidates->begin();
	
	vm.clear();
	
	for(; startVertexCandidateIterator != startVertexCandidates->end(); ++startVertexCandidateIterator){
		// Degree filter
		if(queryGraph->degree[queryTree[0].vertexId] > dataGraph->degree[*startVertexCandidateIterator]){
			continue;
		}
	
		// neighborhood label count filter
		if(!NLCFilter(queryTree[0].vertexId, *startVertexCandidateIterator)){
			continue;
		}
		
		// passes all the filters, add to the list of nodes to be traversed
		vm.push_back(*startVertexCandidateIterator);
	}


	time = wtime();
	numberOfEmbeddings_sngl = 0;
	for(std::vector<int>::iterator cItr = vm.begin(); cItr != vm.end(); ++cItr){
		if(numberOfEmbeddings_sngl >= sufficientNumberOfEmbeddings){
			break;
		}
		
		clearCR();
		// std::cout << *cItr << std::endl;	
		// time = wtime();
		// for(int j = 0; j < 5; j++){
		myCR[0].node_cands.push_back(*cItr);
		//}
	}
	//std::cout << "Building index" << std::endl;
	exploreGraph_sngl();	
	printf("subgraph exploration  completed at %f Seconds\n", wtime() - time);

	for(int i = 0; i< queryTree.size(); ++i){
		queryMatchingSequence[i] = queryTree[i].id;
		queryMatchingSequence_inv[queryTree[i].id] = i;
	}

	//time = wtime();
		
	int thid = 0;
	if(!myCR[0].node_cands.empty()){
		myTCB[thid].embedding[queryMatchingSequence[0]] = myCR[0].node_cands[0];
		myTCB[thid].curr_query_node++;
		subgraphSearch_topk(thid);
		myTCB[thid].curr_query_node--;
	}
	totalEmbeddings();
	//}
	printf("%llu subgraph matches found at %f Seconds\n", numberOfEmbeddings_sngl, wtime()-time);
/*	uint8_t *grouped = new uint8_t[queryGraph->vert_count];
	std::memset(grouped, 0, queryGraph->vert_count*sizeof(int));

	int porderid = 0;

	for(int i = 0; i < queryGraph->vert_count; ++i){
		if(grouped[i]){continue;}
		automorphId[i] = porderid;
		grouped[i] = 1;
		for(std::vector<int*>::iterator itr = Embedding_list.begin(); itr != Embedding_list.end(); ++itr){
			automorphId[(*itr)[i]] = porderid;
			grouped[(*itr)[i]] = 1;
		}
		porderid++;
	}

	// printing and making sure the grouping is correct
	for(int i  = 0; i < queryGraph->vert_count; ++i){
		printf("%d ", automorphId[i]);
	}
	printf("\n");
*/
}

//Heuristic based workload estimation under each subregion
void
subgraph::workLoadEstimator()
{
	int count = 0;
	workload.reserve(vm.size());
	for(std::vector<int>::iterator cItr = vm.begin(); cItr != vm.end(); ++cItr){
		//printf("Node %d\n", vm[*cItr]);
		double score = 0;
		// for all children of given query node
		// find neighbors having their label
		// add degree of all those neighbors and thats the score
		
		for(int childItr = 0; childItr < queryTree[0].childCount; ++childItr){
			// label of child query vertex
			int label = queryTree[queryTree[0].childList[childItr]].label;

			// get all nodes
			std::map<int, vector<int> >::iterator mItr = dataGraph->adjVertices[*cItr].labelCount.find(label);
			if(mItr == dataGraph->adjVertices[*cItr].labelCount.end()){
				score = 0;
				break;
			}
			score += mItr->second.size();
			for(std::vector<int>::iterator vItr = mItr->second.begin(); vItr != mItr->second.end(); ++vItr){
				score += dataGraph->degree[*vItr]; 
			}
		}
		workload.push_back(std::make_pair(score, *cItr));
		//printf("Score is %f\n", score);
		//count++;

		//assert(count != 100);
	}
	// sort the wokload vector in descending order of score
	std::sort(workload.begin(), workload.end(), ScoreCompareRev());
//	calculateSimilarityMatrix();
//	for(int i = 0; i< workload.size(); ++i){
//		printf("%f %d\n", workload[i].first, workload[i].second);
//		assert(i != 1000);
//	}
}

// This Function calculates the similarity score between all pairs of candidate subregion
// Jaccard coefficient on 1-hop neighborhood is used as the measure of similarity
//
//
// Jaccard Coefficient
//

void subgraph::calculateSimilarityMatrix()
{
	// naive version via looping
	// take elements upto 5% of maximum workload, ignore others as they are non-significant
	// If two subregions are similar enough(50% ??), group them into same embedding-region
	// work for the heuristics to tune the 5% and 50% : need to write this on paper
	double time;
	time = wtime();
	std::map<int, uint8_t> checkedLabel;
	unsigned int significantCount;
//	significantCount = 1000;
	if(workload.size()< 1000)
	{
		significantCount = workload.size();
	}
	else if(workload[1000].first <= 0.05*workload[0].first){
		significantCount = 1000;
	}
	else{
		significantCount = 1001;
		while(workload[significantCount].first > 0.05 * workload[0].first){
			significantCount++;
		}
	}
	
	similarity = new double[significantCount*significantCount];
	uint8_t *state = new uint8_t[dataGraph->vert_count];
	for(int firstnode = 0; firstnode < significantCount; ++firstnode){
		memset(state, 0, dataGraph->vert_count*sizeof(uint8_t));
		uint32_t firstsize = 0;
		checkedLabel.clear();

		for(int childItr = 0; childItr < queryTree[0].childCount; ++childItr){
			int label = queryTree[queryTree[0].childList[childItr]].label;
		
			std::map<int, uint8_t>::iterator lItr = checkedLabel.find(label);
			if(lItr != checkedLabel.end()){	continue;}
			
			checkedLabel.insert(std::pair<int, uint8_t>(label, 1));
			std::map<int, std::vector<int> >::iterator mItr = dataGraph->adjVertices[workload[firstnode].second].labelCount.find(label);
			firstsize += mItr->second.size();

			for(std::vector<int>::iterator vItr1 = mItr->second.begin(); vItr1 != mItr->second.end(); ++vItr1){
				state[*vItr1] = 1;
			}
		}
		for(int secondnode = firstnode+1; secondnode < significantCount; ++secondnode){
			// now check if the nodes in the vectors are already encountered
			// if yes, increase the intersect count, else continue
			uint32_t secondsize = 0;
			uint32_t intersectCount = 0;
			checkedLabel.clear();
			
			for(int childItr = 0; childItr < queryTree[0].childCount; ++childItr){
				int label = queryTree[queryTree[0].childList[childItr]].label;

				std::map<int, uint8_t>::iterator lItr = checkedLabel.find(label);
				if(lItr != checkedLabel.end()){continue;}

				checkedLabel.insert(std::pair<int, uint8_t>(label, 1));
				std::map<int, std::vector<int> >::iterator mItr = dataGraph->adjVertices[workload[secondnode].second].labelCount.find(label);
				secondsize += mItr->second.size();
				for(std::vector<int>::iterator vItr2 = mItr->second.begin(); vItr2 != mItr->second.end(); ++vItr2){
					if(state[*vItr2] == 1){++intersectCount;}
				}
			}
			similarity[firstnode*significantCount+secondnode] = (float)intersectCount/(float)(firstsize+secondsize-intersectCount);
		}
	}
	printf("Similarity matrix calculation completed on %f sec", wtime()-time);
	// test if the similarity matrix makes sense
	for(int i = 0; i< significantCount; ++i){
	//	for(int j = 0; j<= i; ++j){	printf("\t ");}
		for(int j = i+1; j < significantCount; ++j){
			if(similarity[i*significantCount+j] >= 0.4 )
				printf("%f ", similarity[i*significantCount+j]);
		}
	//	printf("\n");
	}
}

/* Compute exploration order for second round of exploration */
/*void subgraph::computeMatchingOrder(){
	queryMatchingSequence.clear();

	// For storing the score corresponding to every vertex
	float* queryVertexScore = new float[necTree.size()];
	for(int queryVertexIndex = 0; queryVertexIndex < necTree.size(); queryVertexIndex++){
		queryVertexScore[queryVertexIndex] = 0;
	}

	for(int queryVertexIndex = 0; queryVertexIndex < necTree.size(); queryVertexIndex++){
		map<int, vector<int> >::iterator CR_key_Iterator = CR_key.find(queryVertexIndex);

		//add the candidate numbers
		if(CR_key_Iterator == CR_key.end()){
			continue;
		}

		vector<int>::iterator parentIterator = CR_key_Iterator->second.begin();
		for(; parentIterator != CR_key_Iterator->second.end(); parentIterator++){
		//	map<std::pair<int, int>, vector<int> >::iterator CR_Iterator = CR.find(std::pair<int, int>(queryVertexIndex, *parentIterator));
			map<int, vector<int> >::iterator CR_Iterator = myCR[queryVertexIndex].cand.find(*parentIterator);
			if(CR_Iterator == myCR[queryVertexIndex].cand.end()){
				continue;
			}
			queryVertexScore[queryVertexIndex] += CR_Iterator->second.size();
		}
	}

	for(int queryVertexIndex = 0; queryVertexIndex < necTree.size(); queryVertexIndex++){
		//divide the non-tree edges
		int numberOfNoTreeEdges;
		if(necTree[queryVertexIndex].parent != -1){
			numberOfNoTreeEdges = queryGraph->degree[necTree[queryVertexIndex].vertexId] - necTree[queryVertexIndex].childList.size() - 1;
		}
		else{
			numberOfNoTreeEdges = queryGraph->degree[necTree[queryVertexIndex].vertexId] - necTree[queryVertexIndex].childList.size();
		}
		divideNoTreeEdges(queryVertexIndex, numberOfNoTreeEdges+1, queryVertexScore);
	}

	vector<int> nextVertex;
	vector<int>::iterator childIterator = necTree[0].childList.begin();
	for(; childIterator != necTree[0].childList.end(); childIterator++){
		nextVertex.push_back(*childIterator);
	}
	getOrderByBFSScore(0, queryVertexScore, nextVertex);

	delete [] queryVertexScore;
}

void subgraph::divideNoTreeEdges(int u, int numberOfNoTreeEdges, float * queryVertexScore){
	queryVertexScore[u] /= numberOfNoTreeEdges;
	for(vector<int>::iterator childIterator = necTree[u].childList.begin(); childIterator != necTree[u].childList.end(); ++childIterator){
		divideNoTreeEdges(*childIterator, numberOfNoTreeEdges, queryVertexScore);
	}
}


//This Function runs iteratively finding minimum score at each level of BFS and digging down that particular vertex to next level
void subgraph::getOrderByBFSScore(int u, float * queryVertexScore, vector<int> & nextVertex){
	queryMatchingSequence.push_back(u);

	if(nextVertex.size() == 0){
		return;
	}

	float minimum = FLT_MAX;
	int minimumVertex = 1;
	int vertexNextId  = -1;

	for(int i = 0; i< nextVertex.size(); ++i){
		if(queryVertexScore[nextVertex[i]] == -1){
			continue;
		}
		else{
			if(queryVertexScore[nextVertex[i]] < minimum){
				minimumVertex = nextVertex[i];
				vertexNextId = i;
				minimum = queryVertexScore[nextVertex[i]];
			}
		}
	}

	queryVertexScore[minimumVertex] = -1;
	if(vertexNextId != -1){
		nextVertex.erase(nextVertex.begin() + vertexNextId);
	}

	for(vector<int>::iterator childIterator = necTree[minimumVertex].childList.begin(); childIterator != necTree[minimumVertex].childList.end(); ++childIterator){
		nextVertex.push_back(*childIterator);
	}
	getOrderByBFSScore(minimumVertex, queryVertexScore, nextVertex);
}

*/
int subgraph::chooseStartVertex(){
	std::map<int, vector<int> >::iterator labelVertexIterator;
	float* rankScore = new float[queryGraph->vert_count];

	for(int queryVertexIterator = 0; queryVertexIterator < queryGraph->vert_count; ++queryVertexIterator){
		if(queryGraph->degree[queryVertexIterator] != 0){
			labelVertexIterator = dataGraph->labelVertexList.find(queryGraph->csr_label[queryVertexIterator]);
			if(queryGraph->csr_label[queryVertexIterator] == -1 || labelVertexIterator != dataGraph->labelVertexList.end()){
				rankScore[queryVertexIterator] = (float)labelVertexIterator->second.size()/ queryGraph->degree[queryVertexIterator];
				// std::cout << "Rank Score " << queryGraph->degree[queryVertexIterator] << std::endl;
			}
			else{
				std::cout << "label not found" << std::endl;
				return -1;
			}
		}
	}

	//minimum number of query vertices is
	int topThree[3] = {-1, -1, -1};

	float minScore = FLT_MAX;
	for(int i = 0; i < queryGraph->vert_count; ++i){
		if(rankScore[i] < minScore){
			minScore = rankScore[i];
			topThree[0] = i;
		}
	}

	minScore = FLT_MAX;
	for(int i = 0; i < queryGraph->vert_count; ++i){
		if(i == topThree[0]){
			continue;
		}
		if(rankScore[i] < minScore){
			minScore = rankScore[i];
			topThree[1] = i;
		}
	}

	minScore = FLT_MAX;
	if(queryGraph->vert_count >= 3){
		for(int i = 0; i < queryGraph->vert_count; ++i){
			if(i == topThree[0] || i == topThree[1]){
				continue;
			}
			if(rankScore[i] < minScore){
				minScore = rankScore[i];
				topThree[2] = i;
			}
		}
	}

	int topThreeCandidateCount[3];
	int us;
	int minTopThree = INT_MAX;

	for(int i = 0; i<3; ++i){
		topThreeCandidateCount[i] = 0;
		if(topThree[i] == -1){
			continue;
		}
		map<int, vector<int> >::iterator candidateListIterator = dataGraph->labelVertexList.find(queryGraph->csr_label[topThree[i]]);
		if(queryGraph->csr_label[topThree[i]] != -1 && candidateListIterator == dataGraph->getLabelVertexList()->end()){
			return -1;
		}

		for(vector<int>::iterator candidateIterator = candidateListIterator->second.begin();
				candidateIterator != candidateListIterator->second.end(); ++candidateIterator){
			if(NLCFilter(topThree[i], *candidateIterator)){
				topThreeCandidateCount[i]++;
			}
		}
		if(minTopThree > topThreeCandidateCount[i]){
			us = topThree[i];
			minTopThree = topThreeCandidateCount[i];
		}
	}
	
//	std::cout << "Three candidates are: " << std::endl;
//	for(int i = 0; i < 3; ++i){
//		std::cout << topThree[i] << std::endl;
//	}
	
	return us;
}

/*
 * Second implementation of spanning Tree (other than BFS)
 * It helps to test two issues; effect of different spanning tree and different visit order
 * The nte-edges now may span multiple level in contrast to 
 */

// check if all the vertices has been added to the spanning tree
bool subgraph::allVerticesAdded(std::vector<int>& vertexSet)
{
	for(int index = 0; index < queryGraph->vert_count; ++index){
		if(vertexSet[index] == 0) return false;
	}
	return true;
}

// get the leading node of given edge
int subgraph::getLeadNode(int index)
{
	for(int i = 0; i < queryGraph->vert_count+1; ++i){
		if(queryGraph->beg_pos[i] > index)
			return (i-1);
	}
	return -1;	
}

// add given edge to the spanning tree
int subgraph::addToTree(std::vector<int>& spanTree, std::vector<float>& edgeRanks, std::vector<int>& addedVertices, int node, std::vector<int>& visitOrder)
{
	if(find(spanTree.begin(), spanTree.end(), node) == spanTree.end()){
		spanTree.push_back(node);
		
		queryTree.push_back(queryNode());
		queryNode &qNode 	= *(queryTree.rbegin());
		qNode.vertexId  	= node;
		qNode.parent    	= -1;
		qNode.id		 	= queryTree.size()-1;
		qNode.label 	  	= queryGraph->csr_label[node];
		qGraphtoTreeMap[node] = qNode.id;
		visitTreeOrder[node] = qNode.id;
		// std::cout << "First Node: ID = " << qNode.id << " Vertex = " << qNode.vertexId << " Parent = "<< "-1" << std::endl;
	}


	// Iterate through nodes neighbors
	for(int neighbor = queryGraph->beg_pos[node]; neighbor < queryGraph->beg_pos[node + 1]; ++neighbor){
		int endNode = queryGraph->csr[neighbor];

		// if the endNode of the current edge is in spanTree
		if(addedVertices[endNode] == 0){
			// Then add to vertex set
			spanTree.push_back(endNode);
			visitOrder.push_back(endNode);
			
			// START
			queryTree.push_back(queryNode());
	
			queryNode &qNode 	= *(queryTree.rbegin());
			qNode.vertexId  	= endNode;
			qNode.parent    	= node;
			qNode.id		 	= queryTree.size()-1;
			qNode.label 	  	= queryGraph->csr_label[endNode];
	
			qGraphtoTreeMap[endNode] = qNode.id;
			visitTreeOrder[endNode] = qNode.id;
		
			// std::cout << "First Node: ID = " << qNode.id << " Vertex = " << qNode.vertexId << " Parent: "<< queryTree[qGraphtoTreeMap[node]].vertexId << std::endl;
	
			if(node != -1){
			//	std::cout << "Parents ID = " << queryTree[qGraphtoTreeMap[node]].id << std::endl;
			//	std::cout << "Child Index = " << queryTree[qGraphtoTreeMap[node]].childCount << std::endl;
				queryTree[qGraphtoTreeMap[node]].childList[queryTree[qGraphtoTreeMap[node]].childCount] = endNode;
				queryTree[qGraphtoTreeMap[node]].childCount++;
			}
	
//			for(int vertItr = 0; vertItr < neighbors.size(); ++vertItr){
//				int cNode = neighbors[vertItr].second;
//				if(!flags[cNode]){
//					flags[cNode] = true;
//					S.push(std::pair<int, int>(cNode, qNode.id));
//				}
//				else if(parentId != visitTreeOrder[cNode]){
//					qNode.ntEdges[qNode.ntEdgeCount] = cNode;
//					qNode.ntEdgeCount++;
//				}
//			}

			//END
			// Update addedVertices
			addedVertices[endNode] = 1;
		}
	}

	for(int index = 0; index < queryGraph->edge_count; ++index){
		int lead = getLeadNode(index);
		int follow = queryGraph->csr[index];

		if (addedVertices[lead] == 1 && addedVertices[follow] == 1){
			edgeRanks[index] = -1;
		}
	}
	return 0;
}

// return the edge with maximum rank and is connected to already selected nodes
int subgraph::maxEdgeRank(std::vector<float>& edgeRanks, std::vector<float>& vertexRanks)
{
	// Track the maximum rank edge
	int maxEdgeRankIndex = -1;
	float maxEdgeRank = 0;

	// Iterate through edgeRanks
	for(int index = 0; index < edgeRanks.size(); ++index){
		// Leading node and following node
		int lead = getLeadNode(index);
		int follow = queryGraph->csr[index];

		// ensure valid vertex ranks
		if(vertexRanks[lead] >= vertexRanks[follow]){
			if(edgeRanks[index] >= maxEdgeRank && edgeRanks[index] != -1){
				maxEdgeRank = edgeRanks[index];
				maxEdgeRankIndex = index;
			}
		}
	}
	if(maxEdgeRankIndex == -1)
		return -1;
	else
		edgeRanks[maxEdgeRankIndex] = -1;

	return maxEdgeRankIndex;
}

// construct a spanning tree corresponding to the query graph 
void subgraph::CreateSpanningTree()
{
	int maxVal = 0;
	int maxNode = -1;

	// Selectivity score cor each node and edges
	// More selective nodes goes early
	std::vector<float> vertexRanks;
	std::vector<float> edgeRanks;

	// Vertices already added to the spanning tree
	std::vector<int> addedVertices;
	addedVertices.reserve(queryGraph->vert_count);
	for(int i = 0; i < queryGraph->vert_count; ++i){
		addedVertices[i] = 0;
	}

	std::vector<int> spanTree;
	std::vector<int> visitOrder;

	// COnstruct the vertex scores
	for(int index = 0; index < queryGraph->vert_count; ++index){
		vertexRanks.push_back((queryGraph->beg_pos[index+1] - queryGraph->beg_pos[index]));
	}

//	std::cout << "Vertex Ranks: "<<std::endl;
//	for(int i = 0; i < vertexRanks.size(); ++i){
//		std::cout << vertexRanks[i] << "\t";
//	}
//	std::cout << "\n";

	// Iterate through edges and compute their rank function
	for(int index = 0; index < queryGraph->edge_count; ++index){
		// leading node 
		int lead = getLeadNode(index);
		float leadRank = (queryGraph->beg_pos[lead+1] - queryGraph->beg_pos[lead]);

		// following node
		int follow = queryGraph->csr[index];
		float followRank = (queryGraph->beg_pos[follow+1] - queryGraph->beg_pos[follow]);
		edgeRanks.push_back(leadRank + followRank);
	}

//	std::cout << "Edge Ranks: "<<std::endl;
//	for(int i = 0; i < edgeRanks.size(); ++i){
//		std::cout << edgeRanks[i] << "\t";
//	}
//	std::cout << "\n";

	// Find max edge rank with valid vertex ranks
	int index = 0;
	while(!allVerticesAdded(addedVertices)){
		index = maxEdgeRank(edgeRanks, vertexRanks);
		//std::cout << "Lead Node: " << getLeadNode(index) << std::endl;

		if(find(visitOrder.begin(), visitOrder.end(), getLeadNode(index)) == visitOrder.end()){
			visitOrder.push_back(getLeadNode(index));
		}
		
		addedVertices[getLeadNode(index)] = 1;
		addToTree(spanTree, edgeRanks, addedVertices, getLeadNode(index), visitOrder);

		//std::cout << "Vertex Set: " << std::endl;
		//for(std::vector<int>::iterator iter = spanTree.begin(); iter != spanTree.end(); ++iter){
		//	std::cout << *iter << std::endl;
		//}

		//std::cout << "Visit Order: " << std::endl;
		//for(std::vector<int>::iterator iter = visitOrder.begin(); iter != visitOrder.end(); ++iter){
		//	std::cout << *iter <<std::endl;
		//}
		
	}

	std::cout <<"The size of query Tree is " << queryTree.size() << std::endl;
	for(int i = 0; i < queryTree.size(); ++i){
		std::cout << "VertexId: " << queryTree[i].vertexId << std::endl;
		std::cout << "parentId: " << queryTree[i].parent << std::endl;
		std::cout << "Childrens: ";
		for(int j = 0; j < queryTree[i].childCount; ++j){
			std::cout << queryTree[i].childList[j] << "	";
		}
		std::cout << "\n";
	}

}

/* 
 * Filter if the node doesnot fulfill neighborhood label Count criterion
 * u - query graph node id, v - data graph node id
 */
bool 
subgraph::NLCFilter(int u, int v)
{
	adjLabelFrequency queryVertex = queryGraph->adjVertices[u];
	adjLabelFrequency dataVertex  = dataGraph->adjVertices[v];

	std::map<int, std::vector<int> >::iterator labelVertexListIterator = queryVertex.labelCount.begin();

	for(; labelVertexListIterator != queryVertex.labelCount.end(); ++labelVertexListIterator){
		// if label is not found, return false
		if(labelVertexListIterator->first == -1){
			continue;
		}
		if(dataVertex.labelCount.find(labelVertexListIterator->first) == dataVertex.labelCount.end()){
			return false;
		}
		// else if the count of some label is more in query than on data graph, return false
		else if(labelVertexListIterator->second.size() > dataVertex.labelCount.find(labelVertexListIterator->first)->second.size()){
			return false;
		}
	}
	//if all label satisfy criterion, return true
	return true;
}

/* 
 * NTE Filter, non-tree edge based filtering, find NTE_List first
 * Filter if the candidate node doesnot have neighbor on candidate of NTE neighbor node
 * NETFilter = Non Tree Egde filter
 * u - nec tree query node id, v - data node id
 */
bool 
subgraph::NTEFilter(int u, int v, int pNode)
{
	bool return_value = true;
	for(int cnode = 0; cnode < queryTree[u].ntEdgeCount; ++cnode)
	{
		//if non tree neighbor is already visited
		//if they have common parent, apply that filter
		//as it will reduce the number of comparison by massive amount
		int target_id = qGraphtoTreeMap[queryTree[u].ntEdges[cnode]];
		if(target_id > u){
			continue;
		}
	
		if(queryTree[u].parent == queryTree[target_id].parent){
			edgeCandIterator mItr = std::lower_bound(myCR[target_id].edge_cands.begin(), myCR[target_id].edge_cands.end(), pNode, DataCompare()); //myCR[target_id].edge_cands.find(pNode);// TODO
			if(mItr != myCR[target_id].edge_cands.end()){
				for(std::vector<int>::iterator vItr = mItr->second.begin(); vItr != mItr->second.end(); ++vItr){
					if(dataGraph->isEdge(v, *vItr)){goto NEXT;}
				}
				return false;
			}
		}
		else{
			for(std::vector<int>::iterator vItr = myCR[target_id].node_cands.begin(); vItr != myCR[target_id].node_cands.end(); ++vItr){
				if(dataGraph->isEdge(v, *vItr)){goto NEXT;}
			}
			return false;
		}
NEXT: return_value &= true;
	}

	return return_value;
}

/*
 * Construct the BFS tree for query Graph
 */
void 
subgraph::generateQueryTree()
{
	//std::cout << "Number of Nodes" << queryGraph->vert_count << std::endl;
	bool *flags = new bool[queryGraph->vert_count];
	memset(flags, false, queryGraph->vert_count * sizeof(bool));
	memset(visitTreeOrder, -1, queryGraph->vert_count * sizeof(uint8_t));

	std::queue<std::pair<int, int> >S;
	int vertexId, parentId;

	S.push(std::make_pair(startQueryVertex, -1));
	flags[startQueryVertex] = true;

	while(!S.empty()){
		vertexId = S.front().first;
		parentId = S.front().second;
		S.pop();
		//printf("The vertex %d \n", vertexId);
		//printf("The queryTree size is %d \n", queryTree.size());

		queryTree.push_back(queryNode());

		queryNode &qNode 	= *(queryTree.rbegin());
		qNode.vertexId  	= vertexId;
		qNode.parent    	= parentId;
		qNode.id		 	= queryTree.size()-1;
		qNode.label 	  	= queryGraph->csr_label[vertexId];
		qNode.childCount = 0;
	
		qGraphtoTreeMap[vertexId] = qNode.id;
		
		visitTreeOrder[vertexId] = qNode.id;

		if(parentId != -1){
			queryTree[parentId].childList[queryTree[parentId].childCount] = qNode.id;
			queryTree[parentId].childCount++;
		}
		
		// push 
		std::vector<std::pair<int, int> > neighbors;
		for(int vertItr = queryGraph->beg_pos[vertexId]; vertItr != queryGraph->beg_pos[vertexId+1]; ++vertItr){
			neighbors.push_back(std::pair<int, int>(queryGraph->degree[queryGraph->csr[vertItr]], queryGraph->csr[vertItr]));
		}
		std::sort(neighbors.begin(), neighbors.end(), SortbyDegree());

	//	for(int vertItr = 0; vertItr < neighbors.size(); ++vertItr){
	//		int cNode = neighbors[vertItr].second;
	//		printf("%d  ",cNode);
	//	}

		for(int vertItr = 0; vertItr < neighbors.size(); ++vertItr){
			int cNode = neighbors[vertItr].second;
			if(!flags[cNode]){
				flags[cNode] = true;
				S.push(std::pair<int, int>(cNode, qNode.id));
			}
			else if(parentId != visitTreeOrder[cNode]){
				qNode.ntEdges[qNode.ntEdgeCount] = cNode;
				qNode.ntEdgeCount++;
			}
		}
		/*
		for(int vertItr = queryGraph->beg_pos[vertexId]; vertItr < queryGraph->beg_pos[vertexId+1]; ++vertItr){
		
			if(!flags[queryGraph->csr[vertItr]]){
				flags[queryGraph->csr[vertItr]] = true;
				S.push(std::pair<int, int>(queryGraph->csr[vertItr], qNode.id));
			}
			else if (parentId != visitTreeOrder[queryGraph->csr[vertItr]]){	
				qNode.ntEdges[qNode.ntEdgeCount] = queryGraph->csr[vertItr];
				qNode.ntEdgeCount++;
			}
		}*/
	}		
	delete [] flags;
}

//  Graph Exploration routine
void 
subgraph::exploreGraph_sngl()
{
	//memset(go_status, true , dataGraph->vert_count * sizeof(node_status_t));
	for(int curr_queryNode = 1; curr_queryNode < queryTree.size(); ++curr_queryNode){
		
		//Query Node Specific Data
		int parent = queryTree[curr_queryNode].parent;
		memset(visit_status, UNVISITED , dataGraph->vert_count * sizeof(node_status_t));
		memset(acceptance_status, REJECTED , dataGraph->vert_count * sizeof(node_status_t));

		// std::cout << curr_queryNode << std::endl;

		//memset(parent_count, 0 , dataGraph->vert_count * sizeof(uint8_t));
		for(int mybeg = 0; mybeg < myCR[parent].node_cands.size(); ++mybeg){
			int pNode = myCR[parent].node_cands[mybeg];
			
			if(pNode == -1){ continue;}
			// do not continue if it has been removed from the candidate set earlier

			vector<int> collector;
			collector.reserve(50);

			// GET the neighbors of candidates of parent vertex
			adjLabelFrequency dataGraphChild = dataGraph->adjVertices[pNode];
			std::vector<int> neighbors;
			neighbors.clear();
			std::map<int, std::vector<int>>::iterator it;
			if(queryTree[curr_queryNode].label == -1){
				for(it = dataGraphChild.labelCount.begin(); it != dataGraphChild.labelCount.end(); ++it){
					neighbors.insert(neighbors.end(), it->second.begin(), it->second.end());
				}	
			}
			else{
				it = dataGraphChild.labelCount.find(queryTree[curr_queryNode].label);
				neighbors.insert(neighbors.begin(), it->second.begin(), it->second.end());
				// insert it->second to neighbors 
			}

			if(neighbors.empty()){continue;}
//			if(childLabelList == dataGraphChild.labelCount.end()){ continue;}
			
			//for(vector<int>::iterator childLabelItem = childLabelList->second.begin(); childLabelItem != childLabelList->second.end(); ++childLabelItem){
			for(vector<int>::iterator childLabelItem = neighbors.begin(); childLabelItem != neighbors.end(); ++childLabelItem){
				// Check if the node has already been listed
			//	if(break_automorph){
			//		if(automorphId[queryTree[parent].vertexId] == automorphId[queryTree[curr_queryNode].vertexId]){
			//			if(((queryTree[parent].vertexId < queryTree[curr_queryNode].vertexId) && (pNode > *childLabelItem))||((queryTree[parent].vertexId > queryTree[curr_queryNode].vertexId) && (pNode < *childLabelItem))){
			//				continue;
			//			}
			//		}
			//	}

				if( visit_status[*childLabelItem] == VISITED){
					if(acceptance_status[*childLabelItem] == REJECTED){
						continue;
					}else{
						goto UPDATE;
					}
				}
				
				visit_status[*childLabelItem] = VISITED;
				acceptance_status[*childLabelItem] = REJECTED;

				if(queryGraph->degree[queryTree[curr_queryNode].vertexId] > dataGraph->degree[*childLabelItem]){
					continue;
				}


				/* Neighborhood Label Count based early filtering */
				if(!NLCFilter(queryTree[curr_queryNode].vertexId, *childLabelItem)){
					continue;
				}

				//myCR[curr_queryNode].childrens[mybeg].push_back(*childLabelItem);
				myCR[curr_queryNode].node_cands.push_back(*childLabelItem);
				acceptance_status[*childLabelItem] = ACCEPTED;
				//go_status[*childLabelItem] = true;

UPDATE:			collector.push_back(*childLabelItem);
				//if(collector.size() == 50){
				//	break;
				//}
				//parent_count[*childLabelItem]++;
		 	}

			//If not empty, simply insert
			if(!collector.empty()){
				//insert(std::pair<int, vector<int> >(pNode, collector));
				myCR[curr_queryNode].edge_cands.push_back(std::pair<int, vector<int> >(pNode, collector));
				collector.clear();
				continue;
			}

			
			// If the collector vector is empty
			// for all visited childs of pNode clear those entries from the candidate
			queryNode myParent = queryTree[parent];

			for(int childIndexLooper = 0; childIndexLooper < myParent.childCount; ++childIndexLooper){
				if(myParent.childList[childIndexLooper] == curr_queryNode){
					break;
				}
				
				// search for the entry corresponding to pNode a.k.a. parent node in the edge list of already visited child nodes
				// Nothing has been done to the parent node yet
				edgeCandIterator map_itr = std::lower_bound(myCR[myParent.childList[childIndexLooper]].edge_cands.begin(), 
											myCR[myParent.childList[childIndexLooper]].edge_cands.end(), pNode, DataCompare());

				if(map_itr != myCR[myParent.childList[childIndexLooper]].edge_cands.end()){
					// Clear the entry found
					myCR[myParent.childList[childIndexLooper]].edge_cands.erase(map_itr);
				}
			}
			//now do the dirty job related to pNode in parent itself
			// 1. set cardinality to 0
			// 2. change element to -1 so that it can be ignored in next sibilings
			myCR[parent].cardinality[pNode] = 0;
			myCR[parent].node_cands[mybeg] = -1;
		}
	}
	// now sort the edge_cands list for all the nodes to allow the lower bound and binary search
	for(int curr_queryNode = 0; curr_queryNode < queryTree.size(); ++curr_queryNode){
		std::sort(myCR[curr_queryNode].edge_cands.begin(), myCR[curr_queryNode].edge_cands.end(), DataCompare());
	}

	// create the nte candidate index in ordered fashion
	//float time1 = wtime();
	for(int curr_queryNode = 0; curr_queryNode < queryTree.size(); ++curr_queryNode){
		for(int cnode = 0; cnode < queryTree[curr_queryNode].ntEdgeCount; ++cnode){
			//if non tree neighbor is already visited, create index for it
			int target_id = qGraphtoTreeMap[queryTree[curr_queryNode].ntEdges[cnode]];
			if(target_id > curr_queryNode){continue;}
			std::vector<Data> xnte_cands;
			xnte_cands.reserve(10000);
			
			memset(visit_status, UNVISITED , dataGraph->vert_count * sizeof(node_status_t));
			memset(acceptance_status, REJECTED , dataGraph->vert_count * sizeof(node_status_t));

			//memset(parent_count, 0 , dataGraph->vert_count * sizeof(uint8_t));
			for(int mybeg = 0; mybeg < myCR[target_id].node_cands.size(); ++mybeg){
				int pNode = myCR[target_id].node_cands[mybeg];
				
				if(pNode == -1){ continue;}
				// do not continue if it has been removed from the candidate set earlier
	
				vector<int> collector;
				collector.reserve(50);
	
				// GET the neighbors of candidates of parent vertex
				adjLabelFrequency dataGraphChild = dataGraph->adjVertices[pNode];
				map<int, vector<int> >::iterator childLabelList = dataGraphChild.labelCount.find(queryTree[curr_queryNode].label);
	
				if(childLabelList == dataGraphChild.labelCount.end()){ continue;}
				
				for(vector<int>::iterator childLabelItem = childLabelList->second.begin(); childLabelItem != childLabelList->second.end(); ++childLabelItem){
			//		if(break_automorph){
			//			if(automorphId[queryTree[target_id].vertexId] == automorphId[queryTree[curr_queryNode].vertexId]){
			//				if(((queryTree[target_id].vertexId < queryTree[curr_queryNode].vertexId) && (pNode > *childLabelItem))||((queryTree[target_id].vertexId > queryTree[curr_queryNode].vertexId) && (pNode < *childLabelItem))){
			//					continue;
			//				}
			//			}
			//		}
					// Check if the node has already been listed
					if( visit_status[*childLabelItem] == VISITED){
						if(acceptance_status[*childLabelItem] == REJECTED){continue;}
						else{goto UPDATE1;}
					}

					visit_status[*childLabelItem] = VISITED;
					acceptance_status[*childLabelItem] = REJECTED;

					if(queryGraph->degree[queryTree[curr_queryNode].vertexId] > dataGraph->degree[*childLabelItem]){
						continue;
					}


					// Neighborhood Label Count based early filtering 
					if(!NLCFilter(queryTree[curr_queryNode].vertexId, *childLabelItem)){
						continue;
					}

					acceptance_status[*childLabelItem] = ACCEPTED;

	UPDATE1:		collector.push_back(*childLabelItem);

				}

				if(!collector.empty()){
					xnte_cands.push_back(std::pair<int, vector<int> >(pNode, collector));
					collector.clear();
					continue;
				}

				// Since no childrens would have been visited yet, no need to 
				// 1. set cardinality to 0
				// 2. change element to -1 so that it can be ignored in next sibilings
				myCR[target_id].cardinality[pNode] = 0;
				myCR[target_id].node_cands[mybeg] = -1;
			}
			std::sort(xnte_cands.begin(), xnte_cands.end(), DataCompare());
			myCR[curr_queryNode].nte_cands.insert(std::pair<int, std::vector<Data> >(target_id, xnte_cands));
			xnte_cands.clear();
		}
	}
	//printf("time for nte-index calculation is %f seconds\n", wtime()-time1);
	//exit(0);
	// code for bottom up refinement and cardinality computation
	// Goals: decay the incomplete search tree to their roots
	// Find cardinality associated with each candidate region, subregion and sub-subregion
	// use the cardinality for obverall work estimation as well as to determine the visit order	
	for(int curr_queryNode = queryTree.size()-1; curr_queryNode > 0; --curr_queryNode){
		// get the parent node
		int parent = queryTree[curr_queryNode].parent;
		for(edgeCandIterator edgItr = myCR[curr_queryNode].edge_cands.begin(); edgItr != myCR[curr_queryNode].edge_cands.end(); ++edgItr){
			int score = 0;
			if(edgItr->second.empty()){continue;}
			for(std::vector<int>::iterator sItr = edgItr->second.begin(); sItr != edgItr->second.end(); ++sItr){
				//int id = edgItr->second[*sItr];
				score += myCR[curr_queryNode].cardinality[*sItr];
				if(myCR[curr_queryNode].cardinality[*sItr] == 0){
					edgItr->second[*sItr] = -1;
				}
			}
			myCR[parent].cardinality[edgItr->first] *= score;
		}
	}
}

void subgraph::exploreGraph()
{
	// This function explores the graph in BFS fashion and lists the candidate for each query node
	// use nec tree for exploration
	for(int curr_queryNode = 1; curr_queryNode < queryTree.size(); ++curr_queryNode){
		// Query Node Specific Data
		int parent = queryTree[curr_queryNode].parent;
		
		memset(visit_status, UNVISITED , dataGraph->vert_count * sizeof(node_status_t));
		memset(acceptance_status, REJECTED , dataGraph->vert_count * sizeof(node_status_t));

		int cand_count = myCR[parent].node_cands.size();
		//myCR[curr_queryNode].childrens = new vector<int>[cand_count];

		// Set Number of threads
		num_thrds = 56;
		if(num_thrds > cand_count){
			num_thrds = cand_count;
		}
		omp_set_num_threads(num_thrds);

#pragma omp parallel num_threads(num_thrds)
		{
			int tid = omp_get_thread_num();
			vector<int> myCands;
			myCands.reserve(10000);

			int step  = myCR[parent].node_cands.size()/num_thrds;
			int mybeg = tid*step;
			int myend = (tid == num_thrds-1)? myCR[parent].node_cands.size():(mybeg + step);

			//if(tid == 0){
			//	cout << "Number of thread launched is = " << omp_get_num_threads() << endl;
			//}

			for(; mybeg < myend; ++mybeg){
				int pNode = myCR[parent].node_cands[mybeg];
				
				std::vector<int> collector;
				collector.reserve(50);
				
				// GET the neighbors of candidates of parent vertex 
				adjLabelFrequency dataGraphChild = dataGraph->adjVertices[pNode];
				map<int, vector<int> >::iterator childLabelList = dataGraphChild.labelCount.find(queryTree[curr_queryNode].label);
	
				// Check the each label filtered candidate 
				if(childLabelList != dataGraphChild.labelCount.end()){
					for(vector<int>::iterator childLabelItem = childLabelList->second.begin(); childLabelItem != childLabelList->second.end(); ++childLabelItem){
						
						// Check if the node has already been listed
						if( visit_status[*childLabelItem] == VISITED){
							if(acceptance_status[*childLabelItem] == REJECTED){
								continue;
							}
							else{
								goto UPDATE;
							}
						}

						visit_status[*childLabelItem] = VISITED;
						acceptance_status[*childLabelItem] = REJECTED;

						// Degree based filtering
						if(queryGraph->degree[queryTree[curr_queryNode].vertexId] > dataGraph->degree[*childLabelItem]){
							continue;
						}
	
						// Neighborhood Label Count based early filtering
						if(!NLCFilter(queryTree[curr_queryNode].vertexId, *childLabelItem)){
							continue;
						}

						// same level and cross level non tree edge based filtering
						if(!NTEFilter(queryTree[curr_queryNode].id, *childLabelItem, pNode)){
							continue;
						}

						acceptance_status[*childLabelItem] = ACCEPTED;
						myCands.push_back(*childLabelItem);

UPDATE:					collector.push_back(*childLabelItem);

						//if(collector.size() == 50){
						//	break;
						//}
					}
				
				
					// If the collector vector is empty
					// for all visited childs of pNode clear those entries from the candidate
					if(!collector.empty()){
					#pragma omp critical
						myCR[curr_queryNode].edge_cands.push_back(pair<int, vector<int> >(pNode, collector));
						
						collector.clear();
						continue;
					}

					queryNode myParent = queryTree[parent];
					for(int childIndexLooper = 0; childIndexLooper < myParent.childCount; ++childIndexLooper){
						if(myParent.childList[childIndexLooper] == curr_queryNode){
							break;
						}
						
						// no need to make atomic, we can tolerate a little bit of error
						edgeCandIterator map_itr =std::lower_bound(myCR[myParent.childList[childIndexLooper]].edge_cands.begin(),
								myCR[myParent.childList[childIndexLooper]].edge_cands.end(), pNode, DataCompare());
						if(map_itr != myCR[myParent.childList[childIndexLooper]].edge_cands.end()){
							#pragma omp critical
								myCR[myParent.childList[childIndexLooper]].edge_cands.erase(map_itr);
						}
					}	
				}
			}
			#pragma omp critical
				myCR[curr_queryNode].node_cands.insert(myCR[curr_queryNode].node_cands.end(), myCands.begin(), myCands.end());
		}
	}
}

//clearCR(necTree[*childIndexLooper].id, vm[vmIndex]);

/*
 * Need to design a cache friendly search function
 * Requirements
 * 	1. cache friendly unlike the DFS based one
 * 	2. Should have limited number of intermediate results
 *  3. 
 */

void 
subgraph::subgraphSearch_sngl(unsigned long long int cardinalityParent)
{
	//gives the id of the query node we are matching
	int matching_node = queryMatchingSequence[myTCB_sngl.curr_query_node];
	
	numberOfRecursiveCalls_sngl++;

	if(myTCB_sngl.curr_query_node == queryGraph->vert_count){
		// printf("One Embedding found\n");
		// showEmbeddings(-1);
		numberOfEmbeddings_sngl++;
		return;
	}
	
	// get the node corresponding to our id
	queryNode nectree_u = queryTree[matching_node];

	// the mapping of the parent in the embedding,
	// it is in position equal to id of parent
	int lookup_value = myTCB_sngl.embedding[nectree_u.parent];

	// find the childs of parent we already matched
	edgeCandIterator candidateListIterator = std::lower_bound(myCR[matching_node].edge_cands.begin(),
			myCR[matching_node].edge_cands.end(), lookup_value, DataCompare());

	// if no children found
	if(candidateListIterator == myCR[matching_node].edge_cands.end()){
	//	printf("List Not Found\n");
		return;
	}
	std::vector<int> inters = candidateListIterator->second;

	//form a vector with intersection of children of each previously visited nodes
	for(int adjItr = 0; adjItr < nectree_u.ntEdgeCount; ++adjItr){
		int nteId = qGraphtoTreeMap[nectree_u.ntEdges[adjItr]];
		if(queryMatchingSequence_inv[nectree_u.id] <  queryMatchingSequence_inv[nteId]){ continue;}
		
		int new_lookup = myTCB_sngl.embedding[nteId];
		// Find the corresponding neighbors of that neighbors
		std::map<int, std::vector<Data> >::iterator mItr = myCR[nectree_u.id].nte_cands.find(nteId);
		if(mItr != myCR[nectree_u.id].nte_cands.end()){
			std::vector<Data>::iterator vItr = std::lower_bound(mItr->second.begin(), mItr->second.end(), new_lookup, DataCompare());
			if(vItr != mItr->second.end()){
				inters = intersection(inters, vItr->second);
			}
		}
	}

	// Go over the whole children list
	//#pragma omp parallel for num_threads(num_thrds) schedule(dynamic)
	for(int j = 0; j < inters.size(); ++j){
		std::vector<int>mypreset;
		mypreset = preset;
		int thid = omp_get_thread_num();
		int candidateIterator = inters[j];
		//check if the given vertex is already matched
		int i = 0;
		for( ; i < myTCB_sngl.curr_query_node; ++i){
			if(myTCB_sngl.embedding[queryMatchingSequence[i]] == candidateIterator){
				break;
			}
		}
		if(i != myTCB_sngl.curr_query_node){
			//printf("Already Used node!!\n");
			continue;
		}

		if(break_automorph){
			i = 0;
			for(; i< myTCB_sngl.curr_query_node; ++i){
				if((automorphId[queryTree[queryMatchingSequence[i]].vertexId] == automorphId[nectree_u.vertexId])/* && (!queryGraph->isEdge(queryTree[queryMatchingSequence[i]].vertexId, nectree_u.vertexId))*/){
					if(((queryTree[queryMatchingSequence[i]].vertexId < nectree_u.vertexId) && (myTCB_sngl.embedding[queryMatchingSequence[i]] > candidateIterator))|| ((queryTree[queryMatchingSequence[i]].vertexId > nectree_u.vertexId) && (myTCB_sngl.embedding[queryMatchingSequence[i]] < candidateIterator))){
						break;
					}
				}
			}
			if(i != myTCB_sngl.curr_query_node){
				continue;
			}
		}
	//	if((myCR[nectree_u.id].cardinality[candidateIterator]*num_thrds*1.25) > totalwork){
	//	#pragma omp critical
	//		residue.push_back(candidateIterator);
	//		continue;
	//	}
		mypreset.push_back(candidateIterator);
//		printf("%d %d\n", myTCB_sngl.embedding[queryMatchingSequence[0]], candidateIterator);
//#pragma omp critical
		work_units.push_back(mypreset);
	//	myTCB[thid].embedding[matching_node] = candidateIterator;
	//	myTCB[thid].curr_query_node++;

	//	subgraphSearch(thid);
	//	myTCB[thid].curr_query_node--;
	}
/*	for(int i = 0; i < residue.size(); ++i){
		for(int thid = 0; thid < num_thrds; ++thid){
			myTCB[thid].embedding[matching_node] = residue[i];
			myTCB[thid].curr_query_node++;
		}
		myTCB_sngl.embedding[matching_node] = residue[i];
		myTCB_sngl.curr_query_node++;
		if(totalwork > 0)
			subgraphSearch_sngl((cardinalityParent*myCR[matching_node].cardinality[residue[i]])/totalwork);
		myTCB_sngl.curr_query_node--;
		for(int thid = 0; thid < num_thrds; ++thid){
			myTCB[thid].curr_query_node--;
		}
	}
	*/
}

//function to generate embeddings from the candidate region data
void 
subgraph::subgraphSearch_topk(int thid)
{
	// gives the id of the query node we are matching
	int matching_node = queryMatchingSequence[myTCB[thid].curr_query_node];
	numberofRecursiveCalls[thid]++;
	
	// One embedding is formed
	if(myTCB[thid].curr_query_node == queryGraph->vert_count){
	//#pragma omp critical
		if(!break_automorph){
			int* newEmbedding = new int[queryGraph->vert_count];
			std::memcpy(newEmbedding, myTCB[thid].embedding, queryGraph->vert_count*sizeof(int));
			Embedding_list.push_back(newEmbedding);
			//showEmbeddings(thid);
		}
		showEmbeddings(thid);
		numberOfEmbeddings_sngl++;//[thid]++;
//		printf("%d \n", numberofEmbeddings[thid]);
		return;
	}

	// get the node corresponding to our id
	queryNode nectree_u = queryTree[matching_node];

	// the mapping of the parent in the current embedding 
	// it is in position equal to id of parent
	int lookup_value = myTCB[thid].embedding[nectree_u.parent];
	
	//	List of candidates corresponding to parent query tree in tree 
	edgeCandIterator candidateListIterator = std::lower_bound(myCR[matching_node].edge_cands.begin(), myCR[matching_node].edge_cands.end(), lookup_value, DataCompare());

	if(candidateListIterator == myCR[matching_node].edge_cands.end()){
		//printf("Candidate List not found");
		return;
	}

	std::vector<int> inters = candidateListIterator->second;

	//form a vector with intersection of children of each previously visited nodes
	for(int adjItr = 0; adjItr < nectree_u.ntEdgeCount; ++adjItr){
		int nteId = qGraphtoTreeMap[nectree_u.ntEdges[adjItr]];
		if(queryMatchingSequence_inv[nectree_u.id] <=  queryMatchingSequence_inv[nteId]){ continue;}
		
		int new_lookup = myTCB[thid].embedding[nteId];
		// Find the corresponding neighbors of that neighbors
		std::map<int, std::vector<Data> >::iterator mItr = myCR[nectree_u.id].nte_cands.find(nteId);
		if(mItr != myCR[nectree_u.id].nte_cands.end()){
			std::vector<Data>::iterator vItr = std::lower_bound(mItr->second.begin(), mItr->second.end(), new_lookup, DataCompare());
			if(vItr != mItr->second.end()){
				std::vector<int> x = vItr->second;
				inters = intersection(x, inters);
			}
		}
	}

	// Go over the whole children list
	for(vector<int>::iterator candidateIterator = inters.begin(); candidateIterator != inters.end(); ++candidateIterator){
		//check if task is completed
		if(myTCB[thid].go == false){
			return;
		}

		if(numberOfEmbeddings_sngl > sufficientNumberOfEmbeddings){
			myTCB[thid].go = false;
			return;
		}
		if(*candidateIterator == -1){
			continue;
		}

		// check if the given vertex is already matched
		int i = 0;
		for( ; i < myTCB[thid].curr_query_node; ++i){
			if(myTCB[thid].embedding[queryMatchingSequence[i]] == *candidateIterator){
				break;
			}
		}
		if(i != myTCB[thid].curr_query_node){ 
			continue;
		}
		if(break_automorph){
			i = 0;
			for(; i< myTCB[thid].curr_query_node; ++i){
				if((automorphId[queryTree[queryMatchingSequence[i]].vertexId] == automorphId[nectree_u.vertexId])/* && (!queryGraph->isEdge(queryTree[queryMatchingSequence[i]].vertexId, nectree_u.vertexId))*/){
					if(((queryTree[queryMatchingSequence[i]].vertexId < nectree_u.vertexId) && (myTCB[thid].embedding[queryMatchingSequence[i]] > *candidateIterator))	|| ((queryTree[queryMatchingSequence[i]].vertexId > nectree_u.vertexId) && (myTCB[thid].embedding[queryMatchingSequence[i]] < *candidateIterator))){
						break;
					}
				}
			}
			if(i != myTCB[thid].curr_query_node){
				continue;
			}
		}
		
		//check if it satisfies the neighborhood requirement: aka non tree edges
		myTCB[thid].embedding[myTCB[thid].curr_query_node] = *candidateIterator;
		myTCB[thid].curr_query_node++;

		subgraphSearch_topk(thid);
		myTCB[thid].curr_query_node--;
	}
}


//function to generate embeddings from the candidate region data
void 
subgraph::subgraphSearch(int thid)
{
	// gives the id of the query node we are matching
	int matching_node = queryMatchingSequence[myTCB[thid].curr_query_node];
	numberofRecursiveCalls[thid]++;
	
	// One embedding is formed
	if(myTCB[thid].curr_query_node == queryGraph->vert_count -1){
		//#pragma omp critical
		// get the node corresponding to our id
		queryNode nectree_u = queryTree[matching_node];

		// the mapping of the parent in the current embedding 
		// it is in position equal to id of parent
		int lookup_value = myTCB[thid].embedding[nectree_u.parent];
	
		//	List of candidates corresponding to parent query tree in tree 
		edgeCandIterator candidateListIterator = std::lower_bound(myCR[matching_node].edge_cands.begin(), myCR[matching_node].edge_cands.end(), lookup_value, DataCompare());

		if(candidateListIterator == myCR[matching_node].edge_cands.end()){
			//printf("Candidate List not found");
			return;
		}

		std::vector<int> inters = candidateListIterator->second;
		//form a vector with intersection of children of each previously visited nodes
		for(int adjItr = 0; adjItr < nectree_u.ntEdgeCount; ++adjItr){
			int nteId = qGraphtoTreeMap[nectree_u.ntEdges[adjItr]];
			if(queryMatchingSequence_inv[nectree_u.id] <=  queryMatchingSequence_inv[nteId]){ continue;}
		
			int new_lookup = myTCB[thid].embedding[nteId];
			// Find the corresponding neighbors of that neighbors
			std::map<int, std::vector<Data> >::iterator mItr = myCR[nectree_u.id].nte_cands.find(nteId);
			if(mItr != myCR[nectree_u.id].nte_cands.end()){
				std::vector<Data>::iterator vItr = std::lower_bound(mItr->second.begin(), mItr->second.end(), new_lookup, DataCompare());
				if(vItr != mItr->second.end()){
					std::vector<int> x = vItr->second;
					inters = intersection(x, inters);
				}
			}
		}
		
		// Go over the whole children list
		for(vector<int>::iterator candidateIterator = inters.begin(); candidateIterator != inters.end(); ++candidateIterator){
			//check if task is completed
			if(myTCB[thid].go == false){
				return;
			}
	
			if(*candidateIterator == -1){
				continue;
			}
	
			// check if the given vertex is already matched
			int i = 0;
			for( ; i < myTCB[thid].curr_query_node; ++i){
				if(myTCB[thid].embedding[queryMatchingSequence[i]] == *candidateIterator){
					break;
				}
			}
			if(i != myTCB[thid].curr_query_node){ 
				continue;
			}
			if(break_automorph){
				i = 0;
				for(; i< myTCB[thid].curr_query_node; ++i){
					if((automorphId[queryTree[queryMatchingSequence[i]].vertexId] == automorphId[nectree_u.vertexId])/* && (!queryGraph->isEdge(queryTree[queryMatchingSequence[i]].vertexId, nectree_u.vertexId))*/){
						if(((queryTree[queryMatchingSequence[i]].vertexId < nectree_u.vertexId) && (myTCB[thid].embedding[queryMatchingSequence[i]] > *candidateIterator))	|| ((queryTree[queryMatchingSequence[i]].vertexId > nectree_u.vertexId) && (myTCB[thid].embedding[queryMatchingSequence[i]] < *candidateIterator))){
					//		printf("%d\t%d\n", i, myTCB[thid].curr_query_node);
							break;
						}
					}
				}
				if(i != myTCB[thid].curr_query_node){
					continue;
				}
			}
			
			// add the node to embedding and call recursively
			myTCB[thid].embedding[myTCB[thid].curr_query_node] = *candidateIterator;
			// Increase the embedding count by 1
			numberofEmbeddings[thid]++;

			//myTCB[thid].curr_query_node++;
	
			//subgraphSearch(thid);
			//myTCB[thid].curr_query_node--;
		}
//		if(!break_automorph){
//			int* newEmbedding = new int[queryGraph->vert_count];
//			std::memcpy(newEmbedding, myTCB[thid].embedding, queryGraph->vert_count*sizeof(int));
//			Embedding_list.push_back(newEmbedding);
//			//showEmbeddings(thid);
//		}
	//	showEmbeddings(thid);
	//	numberofEmbeddings[thid]++;
		return;
	}

	// get the node corresponding to our id
	queryNode nectree_u = queryTree[matching_node];

	// the mapping of the parent in the current embedding 
	// it is in position equal to id of parent
	int lookup_value = myTCB[thid].embedding[nectree_u.parent];
	
	//	List of candidates corresponding to parent query tree in tree 
	edgeCandIterator candidateListIterator = std::lower_bound(myCR[matching_node].edge_cands.begin(), myCR[matching_node].edge_cands.end(), lookup_value, DataCompare());

	if(candidateListIterator == myCR[matching_node].edge_cands.end()){
		//printf("Candidate List not found");
		return;
	}

	std::vector<int> inters = candidateListIterator->second;

	//form a vector with intersection of children of each previously visited nodes
	for(int adjItr = 0; adjItr < nectree_u.ntEdgeCount; ++adjItr){
		int nteId = qGraphtoTreeMap[nectree_u.ntEdges[adjItr]];
		if(queryMatchingSequence_inv[nectree_u.id] <=  queryMatchingSequence_inv[nteId]){ continue;}
		
		int new_lookup = myTCB[thid].embedding[nteId];
		// Find the corresponding neighbors of that neighbors
		std::map<int, std::vector<Data> >::iterator mItr = myCR[nectree_u.id].nte_cands.find(nteId);
		if(mItr != myCR[nectree_u.id].nte_cands.end()){
			std::vector<Data>::iterator vItr = std::lower_bound(mItr->second.begin(), mItr->second.end(), new_lookup, DataCompare());
			if(vItr != mItr->second.end()){
				std::vector<int> x = vItr->second;
				inters = intersection(x, inters);
			}
		}
	}

	// Go over the whole children list
	for(vector<int>::iterator candidateIterator = inters.begin(); candidateIterator != inters.end(); ++candidateIterator){
		//check if task is completed
		if(myTCB[thid].go == false){
			return;
		}

		if(*candidateIterator == -1){
			continue;
		}

		// check if the given vertex is already matched
		int i = 0;
		for( ; i < myTCB[thid].curr_query_node; ++i){
			if(myTCB[thid].embedding[queryMatchingSequence[i]] == *candidateIterator){
				break;
			}
		}
		if(i != myTCB[thid].curr_query_node){ 
			continue;
		}
		if(break_automorph){
			i = 0;
			for(; i< myTCB[thid].curr_query_node; ++i){
				if((automorphId[queryTree[queryMatchingSequence[i]].vertexId] == automorphId[nectree_u.vertexId])/* && (!queryGraph->isEdge(queryTree[queryMatchingSequence[i]].vertexId, nectree_u.vertexId))*/){
					if(((queryTree[queryMatchingSequence[i]].vertexId < nectree_u.vertexId) && (myTCB[thid].embedding[queryMatchingSequence[i]] > *candidateIterator))	|| ((queryTree[queryMatchingSequence[i]].vertexId > nectree_u.vertexId) && (myTCB[thid].embedding[queryMatchingSequence[i]] < *candidateIterator))){
				//		printf("%d\t%d\n", i, myTCB[thid].curr_query_node);
						break;
					}
				}
			}
			if(i != myTCB[thid].curr_query_node){
				continue;
			}
		}
		
		// add the node to embedding and call recursively
		myTCB[thid].embedding[myTCB[thid].curr_query_node] = *candidateIterator;
		myTCB[thid].curr_query_node++;

		subgraphSearch(thid);
		myTCB[thid].curr_query_node--;
	}
}

// INFO: since we have already verified the parent and child relationship during candidate region exploration
// We need to check for presence of non-tree edges only in join phase
// If single threaded, use -1 as thread id
bool subgraph::isJoinable(queryNode& u, int v, int thid)
{
	thread_status temp_tcb;
	if(thid == -1){
		temp_tcb = myTCB_sngl;
	//	printf("isJoinable: Single node version activated\n");
	}
	else{
		temp_tcb = myTCB[thid];
	}

	for(int adjItr = 0; adjItr < u.ntEdgeCount; ++adjItr){
		uint8_t nteId = qGraphtoTreeMap[u.ntEdges[adjItr]];
		if(queryMatchingSequence_inv[u.id] > queryMatchingSequence_inv[nteId]){
//			printf("checking the non tree edge %d --> %d\n", u.id, queryMatchingSequence_inv[nteId]);
			if(!dataGraph->isEdge(v, temp_tcb.embedding[queryMatchingSequence_inv[nteId]])){
				return false;
			}
		}
	}
	return true;
}

/*
 * Clear every bit of information related to the query Graph
 * Run the experiment on single data Graph, so no need to clear the data graph
 *
 */
void subgraph::clean(){
	for(int i = 0; i< queryTree.size(); ++i){
		myCR[i].node_cands.clear();
		myCR[i].edge_cands.clear();
	}
	queryTree.clear();
}

//utility, getter and setter functions
void subgraph::showEmbeddings(int thid){
	int* embedding;
	if(thid == -1){
		embedding = myTCB_sngl.embedding;
	}
	else{
		embedding = myTCB[thid].embedding;
	}

	printf("{ %d: ", numberOfEmbeddings_sngl);
	for(int i = 0; i < queryGraph->vert_count; i++ ){
		printf(" %d->%d, ", qGraphtoTreeMap[i], embedding[i]);
	}
	cout << "}" << endl;
}


long subgraph::totalEmbeddings()
{
//	numberOfEmbeddings_sngl = 0;
	for (int thdid = 0; thdid < num_thrds; ++thdid ){
//		printf("%d ", numberofEmbeddings[thdid]);
		numberOfEmbeddings_sngl += numberofEmbeddings[thdid];
	}
	return numberOfEmbeddings_sngl;
}

long subgraph::totalRecursiveCalls(){
	numberOfRecursiveCalls_sngl = 0;
	for (int thdid = 0; thdid < num_thrds; ++thdid){
		printf("%llu\n ", numberofRecursiveCalls[thdid]);
		numberOfRecursiveCalls_sngl += numberofRecursiveCalls[thdid];
	}
	return numberOfRecursiveCalls_sngl;
}

std::vector<int> subgraph::intersection(std::vector<int> &v1, std::vector<int> &v2)
{
	vector<int> v3;
	sort(v1.begin(), v1.end());
	sort(v2.begin(), v2.end());
	set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(v3));
	return v3;
}

