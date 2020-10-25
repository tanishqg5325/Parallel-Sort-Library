#include <bits/stdc++.h>
#include <mpi.h>
#include "mergesort.h"
#include "radixsort.h"
#include "quicksort.h"

using namespace std;


MPI_Datatype pSortType;

void pSort::init() {
    MPI_Init(NULL, NULL);
    // Create MPI DataTypes
    MPI_Datatype charArray;
    MPI_Type_contiguous(LOADSIZE, MPI_CHAR, &charArray);
    MPI_Type_commit(&charArray);
    MPI_Datatype baseTypes[] = {MPI_INT, charArray};
    int blockCount[] = {1, 1};
    MPI_Aint byteOffset[] = {0, 4};
    MPI_Type_create_struct(2, blockCount, byteOffset, baseTypes, &pSortType);
    MPI_Type_commit(&pSortType);

    int ID; MPI_Comm_rank(MPI_COMM_WORLD, &ID);
    srand(ID + time(0));
}


void pSort::close() {
    MPI_Finalize();
}


void pSort::sort(pSort::dataType *data, int ndata, pSort::SortType sorter) {
    int numProcs, ID;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &ID);

    int *procN = new int[numProcs], maxSz;
    assert(MPI_Allgather(&ndata, 1, MPI_INT, procN, 1, MPI_INT, MPI_COMM_WORLD) == MPI_SUCCESS);
    assert(MPI_Allreduce(&ndata, &maxSz, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD) == MPI_SUCCESS);

    switch(sorter) {
        case MERGE:
            mergeSortPar(procN, numProcs, maxSz, data, ndata, ID, pSortType);
            break;
        
        case QUICK:
            quickSortPar(procN, numProcs, maxSz, data, ndata, ID, pSortType);
            break;

        case RADIX:
        case BEST:
            radixSortPar(procN, numProcs, maxSz, data, ndata, ID, pSortType);
    };

    delete[] procN;
}
