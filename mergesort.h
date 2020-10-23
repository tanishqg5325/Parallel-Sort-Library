#ifndef MERGE_SORT_H_
#define MERGE_SORT_H_

#include <bits/stdc++.h>
#include <mpi.h>
#include "sort.h"

using namespace std;


void mergeSortPar(int procN[], int numProcs, int maxSz, pSort::dataType *data, int ndata, 
                    int ID, MPI_Datatype& pSortType);

#endif
