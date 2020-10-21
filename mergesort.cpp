#include "mergesort.h"


void mergeSeq(pSort::dataType *ar1, int n1, pSort::dataType *ar2, int n2, pSort::dataType* ext) {
    int i = 0, j = 0, k = 0;
    while(i < n1 && j < n2) {
        if(ar1[i].key <= ar2[j].key) ext[k++] = ar1[i++];
        else ext[k++] = ar2[j++]; 
    }
    while(i < n1) ext[k++] = ar1[i++];
    while(j < n2) ext[k++] = ar2[j++];
    for(int i=0;i<n1;i++) ar1[i] = ext[i];
    for(int i=0;i<n2;i++) ar2[i] = ext[i+n1];
}


void mergeSortSeq(pSort::dataType *data, pSort::dataType *ext, int l, int r) {
    if(l < r) {
        int m = (l + r) >> 1;
        mergeSortSeq(data, ext, l, m);
        mergeSortSeq(data, ext, m+1, r);
        mergeSeq(data+l, m-l+1, data+m+1, r-m, ext);
    }
}


void mergePar(vector<int>& ind, int l_sz, int procN[], MPI_Datatype& pSortType, 
                pSort::dataType *data, pSort::dataType *left, pSort::dataType *right, 
                pSort::dataType *ans) {
    
    MPI_Status status;
    for(int i=0; i<procN[ind[0]]; ++i) left[i] = data[i];
    assert(MPI_Recv(right, procN[ind[l_sz]], pSortType, ind[l_sz], 1, MPI_COMM_WORLD, &status) == MPI_SUCCESS);
    int k = 0, r_sz = ind.size(), l_ind = 0, r_ind = l_sz, i = 0, j = 0;
    while(k < procN[ind[0]]) {
        assert(i < procN[ind[0]]);
        if(r_ind == r_sz)
            data[k++] = left[i++];
        else if(j == procN[ind[r_ind]]) {
            ++r_ind; j = 0;
            assert(MPI_Recv(right, procN[ind[r_ind]], pSortType, ind[r_ind], 1, MPI_COMM_WORLD, &status) == MPI_SUCCESS);
        }
        else if(left[i].key <= right[j].key)
            data[k++] = left[i++];
        else
            data[k++] = right[j++];
    }
    for(int a = 1; a < r_sz; ++a) {
        k = 0;
        while(k < procN[ind[a]]) {
            if(l_ind == l_sz) {
                assert(r_ind < r_sz);
                if(j == procN[ind[r_ind]]) {
                    ++r_ind; j = 0;
                    assert(MPI_Recv(right, procN[ind[r_ind]], pSortType, ind[r_ind], 1, MPI_COMM_WORLD, &status) == MPI_SUCCESS);
                }
                else {
                    ans[k++] = right[j++];
                }
            }
            else if(r_ind == r_sz) {
                assert(l_ind < l_sz);
                if(i == procN[ind[l_ind]]) {
                    ++l_ind; i = 0;
                    assert(MPI_Recv(left, procN[ind[l_ind]], pSortType, ind[l_ind], 1, MPI_COMM_WORLD, &status) == MPI_SUCCESS);
                }
                else {
                    ans[k++] = left[i++];
                }
            }
            else {
                if(i == procN[ind[l_ind]]) {
                    ++l_ind; i = 0;
                    if(l_ind == l_sz) continue;
                    assert(MPI_Recv(left, procN[ind[l_ind]], pSortType, ind[l_ind], 1, MPI_COMM_WORLD, &status) == MPI_SUCCESS);
                }
                else if(j == procN[ind[r_ind]]) {
                    ++r_ind; j = 0;
                    if(r_ind == r_sz) continue;
                    assert(MPI_Recv(right, procN[ind[r_ind]], pSortType, ind[r_ind], 1, MPI_COMM_WORLD, &status) == MPI_SUCCESS);
                }
                else if(left[i].key <= right[j].key)
                    ans[k++] = left[i++];
                else
                    ans[k++] = right[j++];
            }
        }
        MPI_Send(ans, procN[ind[a]], pSortType, ind[a], 1, MPI_COMM_WORLD);
    }
}


void mergeSortPar(int procN[], int numProcs, int maxSz, pSort::dataType *data, int ndata, 
                    int ID, MPI_Datatype& pSortType) {
    
    pSort::dataType *extra = new pSort::dataType[maxSz], *left = NULL, *right = NULL;
    mergeSortSeq(data, extra, 0, ndata-1);
    MPI_Status status;

    // cout << ID << ": ";
    // for(int i=0;i<ndata; ++i) cout << data[i].key << " "; cout << endl;

    int height = log2(numProcs), h = 0;
    if((numProcs & (numProcs - 1)) != 0) height++;

    if(ID % 2 == 0) {
        left = new pSort::dataType[maxSz];
        right = new pSort::dataType[maxSz];
    }

    while(h < height) {
        int par = ID & (~((1 << (h+1)) - 1));
        int sz = min(1 << (h+1), numProcs - par);
        if(sz > (1 << h)) {
            if(par == ID) {
                vector<int> v(sz);
                for(int i=0; i<sz; ++i) v[i] = i+ID;
                mergePar(v, 1<<h, procN, pSortType, data, left, right, extra);
            }
            else {
                for(int i=0; i<procN[ID] ; ++i) extra[i] = data[i];
                assert(MPI_Sendrecv(extra, procN[ID], pSortType, par, 1, data, procN[ID], pSortType, par, 1, MPI_COMM_WORLD, &status) == MPI_SUCCESS);
            }
        }
        h++;
    }

    // cout << ID << ": ";
    // for(int i=0;i<ndata;++i) cout<<data[i].key<<" "; cout<<endl;

    delete[] extra;
    if(ID % 2 == 0) {
        delete[] left;
        delete[] right;
    }
}