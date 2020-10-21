#include "sort.h"

using namespace std;


void pSort::init(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
}


void pSort::close() {
    MPI_Finalize();
}


void pSort::sort(pSort::dataType *data, int ndata, pSort::SortType sorter) {
    int numProcs, ID;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &ID);
    
    // Create MPI DataTypes
    MPI_Datatype charArray;
    MPI_Type_contiguous(LOADSIZE, MPI_CHAR, &charArray);
    MPI_Type_commit(&charArray);
    MPI_Datatype baseTypes[] = {MPI_INT, charArray};
    int blockCount[] = {1, 1};
    MPI_Aint byteOffset[] = {0, 4};
    MPI_Datatype pSortType;
    MPI_Type_create_struct(2, blockCount, byteOffset, baseTypes, &pSortType);
    MPI_Type_commit(&pSortType);

    // srand(ID+1);
    // data = new pSort::dataType[ndata];
    // for(int i=0; i<ndata; ++i) data[i] = {rand() % 100, {'a', 'a', 'a', 'a'}};
    // cout << ID << ": ";
    // for(int i=0;i<ndata;++i) cout<<data[i].key<<" "; cout<<endl;

    int nArray[] = {ndata};
    int *procN= new int[numProcs];
    MPI_Allgather(nArray, 1, MPI_INT, procN, 1, MPI_INT, MPI_COMM_WORLD);

    int maxSz = 0;
    for(int i=0; i<numProcs; ++i) maxSz = max(maxSz, procN[i]);

    mergeSortPar(procN, numProcs, maxSz, data, ndata, ID, pSortType);

    delete[] procN;
}
