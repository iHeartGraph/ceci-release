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

#ifndef __GRAPH_H__
#define __GRAPH_H__
#include <iostream>
#include <vector>
#include <map>

#include "util.h"
		
struct adjLabelFrequency{
	std::map<int, std::vector<int> > labelCount;
};


class graph
{
    public:
        index_t *beg_pos;
        vertex_t *csr;
        index_t *beg_pos_label;
        vertex_t *csr_label;
        index_t vert_count;
        index_t vert_max;
        vertex_t edge_count;
        index_t label_count;
        index_t *degree;
		//The data structures for storing the inverse label index
		std::map<int, std::vector<int> > labelVertexList;
		std::vector<std::vector<int> > sorted_csr;

		adjLabelFrequency* adjVertices;
    public:
        graph(const char *graphFile);
        ~graph();
        void test();

		//function for generating inverse label indexes
		void buildLabelVertexList();
		void buildVertexLabelVertexList();
		void build_sorted_csr();
		bool isEdge(int x, int y);

		std::map<int, std::vector<int> >* getLabelVertexList(){
			return &labelVertexList;
		}

		void testIndex(int label){
			std::map<int, std::vector<int> >::iterator itr = adjVertices[4].labelCount.begin();
			for(; itr != adjVertices[4].labelCount.end(); itr++){
				std::vector<int>::iterator j = itr->second.begin();
				std::cout << itr->first <<" " << itr->second.size() << "\n";
				for(; j != itr->second.end(); ++j){
					std::cout << itr->second[*j] << "\t";
				}
				std::cout << "\n";
			}
		}
};

#endif

