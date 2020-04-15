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

#include <cstdio>
#include <fstream>
#include <cstring>

// #include <mpi.h>

#include "wtime.h"
#include "subgraph.h"

using namespace std;

// global members
char *dataGraphFileName, *queryGraphFileName;
int omp_thread_count = 0;
char architecture;
int embedding_count = 0;
char operationMode;
double time1;

// Function declarations
void help();

int main(int argc, char* argv[])
{
	//MPI_Init(&argc, &argv);
	int worker_count = 1;
	int worker_rank  = 0;

	//MPI_Comm_size(MPI_COMM_WORLD, &worker_count);
	//MPI_Comm_rank(MPI_COMM_WORLD, &worker_rank);

	if(argc != 3){
		printf("Wrong number of parameters\n");
		help();
		exit(1);
	}

	// Read query and data file names
	dataGraphFileName  	= argv[1];
	queryGraphFileName 	= argv[2];
	graph* datagraph = new graph(dataGraphFileName);
	graph* querygraph = new graph(queryGraphFileName);

	omp_thread_count = 32;
	bool break_automorph = false;
	if(worker_rank == 0) printf("Data graph: %s and Query Graph %s\n", dataGraphFileName, queryGraphFileName);
	subgraph mySubgraphFinder(dataGraphFileName, queryGraphFileName, omp_thread_count, break_automorph);
	mySubgraphFinder.genericQueryProc();

/*
	// The form of operation; list all embedding or specific number of embeddings
	if(atoi(argv[3]) == -1){
		operationMode = 'a';
	}
	else{
		operationMode = 'k';
		embedding_count = atoi(argv[3]);
	}

	// Choose the memory architecture 
	architecture  		= argv[4][0];

	// Number of threads per machine 
	omp_thread_count 	= atoi(argv[5]);

	// First K embeddings
	if(operationMode == 'k'){
		if(worker_rank == 0) printf("Data graph: %s and Query Graph %s\n", dataGraphFileName, queryGraphFileName);
		bool break_automorph = true;
		subgraph mySubgraphFinder(dataGraphFileName, queryGraphFileName, omp_thread_count, break_automorph, embedding_count);
		mySubgraphFinder.automorphFreeProc();
	}
	// All Embeddings
	else if(operationMode == 'a'){
		if(worker_rank == 0) printf("Data graph: %s and Query Graph %s\n", dataGraphFileName, queryGraphFileName);

		// Find Automorphs
		bool break_automorph = false;
		subgraph myAutomorphFinder(queryGraphFileName, queryGraphFileName, omp_thread_count, break_automorph);
		myAutomorphFinder.automorphFreeProc();

		// Subgraph Isomorphism module object
		break_automorph = true;
		subgraph mySubgraphFinder(dataGraphFileName, queryGraphFileName, omp_thread_count, break_automorph);
		if(architecture == 's'){
			printf("Shared memory architecture\n");
			mySubgraphFinder.genericQueryProc();
		}
		else if(architecture == 'd'){
			printf("Distributed memory architecture\n");
			//mySubgraphFinder.distributedQueryProc();
		}
		else{
			printf("The architecture flag is not identified");
		}
	}

	if(worker_rank == 0){printf("**** subgraph Matching completed *****\n");}
	time1 = wtime();
*/	
	//MPI_Finalize();
	return 0;
}

void help()
{
	cout << "Usage:" << "1. application.bin\n"<< "2. dataGraph\n" << "3. queryGraph\n" << endl;
}

