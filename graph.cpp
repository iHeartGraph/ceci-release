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

#include "graph.h"
#include "wtime.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

graph::graph(const char * graphFile)
{
    vert_count = 0;
    vert_max = 0;
    edge_count = 0;
    label_count = 0;
    size_t len = 0;
    ssize_t read;
    char *line = NULL;
  
    /// only load the data file
    /// read the fist time, to get vert_count, edge_count, label_count
    FILE * fp = fopen(graphFile, "r");
    while((read = getline(&line, &len, fp)) != -1)
    {
        if(read > 0 && line[0] == 'v')
        {
            vert_count ++;
            
            index_t v = 0;
            index_t i = 2;
            for( ; i<read && line[i] != ' '; ++i)
            { 
                v = v * 10 + line[i] - '0';
            }
//            label_count ++;
            for( ; line[i] != '\n'; ++i)
            { 
                if(line[i] == ' ')
                {
                    label_count ++;
                }
            }
            if(v > vert_max)
                vert_max = v;
        }
        else
            if(read > 0 && line[0] == 'e')
                edge_count ++;
    }
    fclose(fp);

    if(DEBUG)
    {
        printf("vert_count, %d\nvert_max, %d\nedge_count, %d\nlabel_count, %d\n", vert_count, vert_max, edge_count, label_count);
    }
    
    /// make sure vert_max >= vert_count
    if(vert_max < vert_count)
        vert_max = vert_count;
    beg_pos = new index_t[vert_max + 1];
    csr = new vertex_t[edge_count*2]; // undirected graph
    beg_pos_label = new index_t[vert_max + 1];
    csr_label = new vertex_t[label_count];
    degree = new index_t[vert_max + 1];

    /// initialize
    memset(degree, 0, sizeof(index_t) * (vert_max + 1));

    /// read the second time, write to beg_pos_label, csr_label, degree
    fp = fopen(graphFile, "r");
    read = getline(&line, &len, fp);
    beg_pos_label[0] = 0;
    index_t vert_id = 0;
    if(read > 0 && line[0] == 't')
    {
        ///starts from 'v'
        index_t label_id = 0;
        while((read = getline(&line, &len, fp)) != -1)
        {
            if(line[0] != 'v')
            {
                break;
            }
            ///extract the vertex_id
            index_t v = 0;
            index_t i = 2;
            for( ; i<read && line[i] != ' '; ++i)
            { 
                v = v * 10 + line[i] - '0';
            }
            
            index_t w = 0;
            i++;
			bool sign = false;
			if(line[i] == '-'){ sign = true; i++;}
            //printf("len, %zd\n", read);
            for( ; i < read; ++i)
            {
                if(line[i] == ' ' || line[i] == '\n')
                {
//                  printf("w = %d, ", w);
                    csr_label[label_id] =(sign)? 0-w : w;
//                  printf("%d\n", csr_label[label_id]);
                    label_id ++;
                    w = 0;
                }
                if(line[i] >= '0' && line[i] <= '9')
                {
                    w = w * 10 + line[i] - '0';
                }
            }
//            csr_label[label_id] = w;
//            label_id ++;
            
            vert_id ++;
            beg_pos_label[vert_id] = label_id;
        }
        if(DEBUG)
        {
            printf("vert_count, %d\nlabel_id, %d\n", vert_count, label_id);
//            for(index_t i=0; i<vert_count; ++i)
//            {
//                printf("%d %d", beg_pos_label[i], beg_pos_label[i+1]);
//                for(vertex_t j=beg_pos_label[i]; j<beg_pos_label[i+1]; ++j)
//                {
//                    printf(" %d", csr_label[j]);
//                }
//                printf("\n");
//            }
        
        }
        /// starts from 'e', assume it's unweighted

        do        
        {
            index_t v = 0;
            index_t w = 0;
            index_t i = 2;
            /// the first edge
            for( ; i < read && line[i] != ' '; ++i)
            {
                v = v * 10 + line[i] - '0';
            }
            i++;
            degree[v]++;
            for( ; i < read && line[i] != ' '; ++i)
            {
                w = w * 10 + line[i] - '0';
            }
            degree[w]++;
        }while((read = getline(&line, &len, fp)) != -1);
        
//        if(DEBUG)
//        {
//            vertex_t sum = 0;
//            for(index_t i=0; i<vert_count; ++i)
//            {
//                sum += degree[i];
//                printf("vert_id, %d, degree, %d\n", i, degree[i]);
//            }
//            printf("sum, %d\n", sum);
//        }

    }
    fclose(fp);
    beg_pos[0] = 0;
    for(index_t i=1; i<vert_max+1; ++i)
    {
        beg_pos[i] = beg_pos[i-1] + degree[i-1];
    }
//    if(DEBUG)
//    {
//        for(index_t i=0; i<vert_count+1; ++i)
//        {
//            printf("%d ", beg_pos[i]);
//        }
//        printf("\n");
//    }
    /// read the third time, write to beg_pos, csr
    fp = fopen(graphFile, "r");
    read = getline(&line, &len, fp);
    beg_pos_label[0] = 0;
    vert_id = 0;
    memset(degree, 0, sizeof(index_t) * (vert_max + 1));

    if(read > 0 && line[0] == 't')
    {
        while((read = getline(&line, &len, fp)) != -1)        
        {
            if(read > 0 && line[0] == 'e')
                break;
        }
       
        do
        {
            index_t v = 0;
            index_t w = 0;
            index_t i = 2;
            for( ; i < read && line[i] != ' '; ++i)
            {
                v = v * 10 + line[i] - '0';
            }
            i++;
            for( ; i < read && line[i] != ' '; ++i)
            {
                w = w * 10 + line[i] - '0';
            }
            csr[beg_pos[v] + degree[v]] = w;
            csr[beg_pos[w] + degree[w]] = v;
            degree[v] ++;
            degree[w] ++;
        }while((read = getline(&line, &len, fp)) != -1);
        
//        if(DEBUG)
//        {
//            vertex_t sum = 0;
//            
//            for(index_t i=0; i<vert_count; ++i)
//            {
//                sum += degree[i];
//                printf("d_v, %d:", i);
////                printf("%d %d:", i, beg_pos[i+1] - beg_pos[i]);
//                for(vertex_t j=beg_pos[i]; j<beg_pos[i+1]; ++j)
//                {
//                    printf(" %d", csr[j]);
//                }
//                printf("\n");
////                printf("vert_id, %d, degree, %d\n", i, degree[i]);
//            }
//            printf("sum, %d\n", sum);
//        }
    }
    fclose(fp);
    
}

void graph::build_sorted_csr()
{
	sorted_csr.reserve(vert_count);
	//for all vertices
	
	//create a container to hold the adjacency list
	std::vector<int> current_adjs;
	current_adjs.reserve(vert_count);
	
	for(int vertItr = 0; vertItr < vert_count; ++vertItr){
		current_adjs.clear();
		// All nodes in the neighbor 
		for(int adjItr = beg_pos[vertItr]; adjItr < beg_pos[vertItr+1]; ++adjItr){
			current_adjs.push_back(csr[adjItr]);
		}
		// sort the neighbor list
		std::sort(current_adjs.begin(), current_adjs.end());

		//insert it to index
		sorted_csr.push_back(current_adjs);
		//printf("The size of adjacency list is %d\n", current_adjs.size());
	}
}

//function for generating inverse label indexes
void graph::buildLabelVertexList()
{
	for(int vertItr = 0; vertItr < vert_count; vertItr++){
		for(int labelItr = beg_pos_label[vertItr]; labelItr < beg_pos_label[vertItr+1]; labelItr++ ){
			//std::cout << csr_label[labelItr];
			std::map<int, std::vector<int> >::iterator labelListItr = labelVertexList.find(csr_label[labelItr]);
			if(labelListItr != labelVertexList.end()){
				//std::cout << ": Old" << std::endl;
				labelListItr->second.push_back(vertItr);
			}
			else{
				//std::cout << csr_label[labelItr] <<std::endl;
				std::vector<int> labelVector;
				labelVector.push_back(vertItr);
				labelVertexList.insert(std::pair<int, std::vector<int> >(csr_label[labelItr], labelVector));
			}
		}
	}
}

//function for creating label sorted adjacency list
void graph::buildVertexLabelVertexList()
{
	adjVertices = new adjLabelFrequency[vert_count];
	std::map<int,std::vector<int> >::iterator vertexLabelVertexIterator;
 	// For all vertices
	for(int vertItr = 0; vertItr < vert_count; vertItr++){
		// for all adjacent nodes of given vertex
		for(int adjItr = beg_pos[vertItr]; adjItr < beg_pos[vertItr+1]; adjItr++){
			// std::cout << csr[adjItr] << "\t";
			// for all labels of given adjacent vertex
			for(int lblItr = beg_pos_label[csr[adjItr]]; lblItr < beg_pos_label[csr[adjItr]+1]; lblItr++ ){
			//	std::cout << csr_label[lblItr] << "\t";

				vertexLabelVertexIterator = adjVertices[vertItr].labelCount.find(csr_label[lblItr]);
				if(vertexLabelVertexIterator != adjVertices[vertItr].labelCount.end()){
					//label already there
					vertexLabelVertexIterator->second.push_back(csr[adjItr]);
				}
				else{
					//a new label
					std::vector<int> newVertexLabelList;
                    newVertexLabelList.push_back(csr[adjItr]);
                    adjVertices[vertItr].labelCount.insert(std::pair<int,std::vector<int> >(csr_label[lblItr],newVertexLabelList));
				}
			}
		}
	}
}

// check if two nodes form an edge 
// use binary search instead of linear scan to find the data in the list
bool graph::isEdge(int x, int y)
{
	// for all adjacent nodes of y, find if x exists in the list
	if(binary_search(sorted_csr[x].begin(), sorted_csr[x].end(), y)){
		return true;
	}
	return false;
}

graph::~graph()
{
//    delete[] beg_pos;
//    delete[] csr;
//    delete[] beg_pos_label;
//    delete[] csr_label;
//    delete[] degree;
}

void graph::test()
{
    printf("\nTest result:\n");
    printf("Vert_count, %d\nVert_max, %d\nEdge_count, %d\nLabel_count, %d\n", vert_count, vert_max, edge_count*2, label_count);
}
