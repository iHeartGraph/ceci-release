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

#ifndef __UTIL_H__
#define __UTIL_H__
#include <stdlib.h>
#include <stdint.h>
#include <sys/stat.h>
#include <algorithm>

#define LOCK(vert, lock) while(!__sync_bool_compare_and_swap(lock+vert,0,-1))
#define UNLOCK(vert, lock) lock[vert]=0

///change to int for SCC
//typedef long index_t;
//typedef long vertex_t;
//typedef double path_t;
//typedef long depth_t;

typedef int index_t;
typedef int vertex_t;
typedef double path_t;
typedef int depth_t;
typedef int color_t;

#define INFTY (float)10000000 
#define NEGATIVE (int)-1
#define ORPHAN	(unsigned char)254
#define UNVIS		(long)-1
#define DEBUG 0 
#define VERBOSE 0 
#define OUTPUT_TIME 1
inline off_t fsize(const char *filename) {
	struct stat st; 
	if (stat(filename, &st) == 0)
		return st.st_size;
	return -1; 
}

#endif
