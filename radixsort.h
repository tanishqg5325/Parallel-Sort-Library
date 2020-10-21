#ifndef RADIX_SORT_H_
#define RADIX_SORT_H_

#include <bits/stdc++.h>
#include <mpi.h>
#include "sort.h"

using namespace std;


void radixSortPar(int procN[], int numProcs, int maxSz, pSort::dataType *data, int ndata, 
                    int ID, MPI_Datatype& pSortType);

#endif
