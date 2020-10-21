#ifndef MERGE_SORT_H_
#define MERGE_SORT_H_

#include <bits/stdc++.h>
#include <mpi.h>
#include "sort.h"

using namespace std;


void mergeSeq(pSort::dataType *ar1, int n1, pSort::dataType *ar2, int n2, pSort::dataType* ext);
void mergeSortSeq(pSort::dataType *data, pSort::dataType *ext, int l, int r);

void mergePar(vector<int>& ind, int l_sz, int procN[], MPI_Datatype& pSortType, 
                pSort::dataType *data, pSort::dataType *left, pSort::dataType *right, 
                pSort::dataType *ans);

void mergeSortPar(int procN[], int numProcs, int maxSz, pSort::dataType *data, int ndata, 
                    int ID, MPI_Datatype& pSortType);

#endif
