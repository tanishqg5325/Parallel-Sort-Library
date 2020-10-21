#include "radixsort.h"


void radixSortPar(int procN[], int numProcs, int maxSz, pSort::dataType *data, int ndata, 
                    int ID, MPI_Datatype& pSortType) {

    #define BASE 16
    #define LOGBASE 4

    int intSz = 8 * sizeof(int);

    int cntNeg = 0, totNeg;
    for(int i=0; i<ndata; ++i) if(data[i].key < 0) cntNeg++;
    assert(MPI_Allreduce(&cntNeg, &totNeg, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD) == MPI_SUCCESS);
    
    vector<pSort::dataType> buckets[BASE];

    int *allBucketSizes = NULL, bucketSizes[BASE], offset[BASE+1];
    pSort::dataType *extra = new pSort::dataType[ndata], *recv_buf = NULL, *ans = NULL;
    if(ID == 0) {
        allBucketSizes = new int[BASE * numProcs];
        recv_buf = new pSort::dataType[maxSz];
        ans = new pSort::dataType[maxSz];
    }
    MPI_Status status;

    for(int i=0; i<intSz; i += LOGBASE) {
        int mask = ((BASE-1) << i);
        // negative integers
        if(i + LOGBASE >= intSz)
            mask &= (~(1 << (intSz - 1))); 
        for(int j=0; j<ndata; ++j) {
            buckets[(abs(data[j].key) & mask) >> i].push_back(data[j]);
        }
        for(int j=0; j<BASE; ++j) bucketSizes[j] = buckets[j].size();
        offset[0] = 0;
        for(int j=1; j<=BASE; ++j) offset[j] = offset[j-1] + bucketSizes[j-1];

        assert(MPI_Gather(bucketSizes, BASE, MPI_INT, allBucketSizes, BASE, MPI_INT, 0, MPI_COMM_WORLD) == MPI_SUCCESS);

        int ind = 0;
        for(int j=0; j<BASE; ++j) {
            for(const pSort::dataType& p : buckets[j])
                extra[ind++] = p;
        }
        assert(ind == ndata);
        int base_ind = 0, proc_ind = 0, data_ind = 0;
        if(ID == 0) {
            for(int j=0; j<numProcs; ++j) {
                ind = 0;
                while(ind < procN[j]) {
                    assert(base_ind < BASE);
                    if(proc_ind == numProcs) {
                        base_ind++;
                        proc_ind = 0;
                    }
                    else if(data_ind == allBucketSizes[proc_ind * BASE + base_ind]) {
                        ++proc_ind; data_ind = 0;
                        if(proc_ind == numProcs) continue;
                        int sz = allBucketSizes[proc_ind * BASE + base_ind];
                        if(sz > 0)
                            assert(MPI_Recv(recv_buf, sz, pSortType, proc_ind, 1, MPI_COMM_WORLD, &status) == MPI_SUCCESS);
                    }
                    else if(proc_ind == 0) {
                        if(j == 0)
                            data[ind++] = extra[offset[base_ind] + data_ind];
                        else
                            ans[ind++] = extra[offset[base_ind] + data_ind];
                        data_ind++;
                    }
                    else {
                        if(j == 0)
                            data[ind++] = recv_buf[data_ind++];
                        else
                            ans[ind++] = recv_buf[data_ind++];
                    }
                }
                if(j > 0)
                    assert(MPI_Send(ans, procN[j], pSortType, j, 1, MPI_COMM_WORLD) == MPI_SUCCESS);
            }
        }
        else {
            MPI_Request req;
            MPI_Irecv(data, ndata, pSortType, 0, 1, MPI_COMM_WORLD, &req);
            for(int j=0; j<BASE; ++j) {
                if(offset[j+1] > offset[j])
                    assert(MPI_Send(extra+offset[j], offset[j+1]-offset[j], pSortType, 0, 1, MPI_COMM_WORLD) == MPI_SUCCESS);
            }
            MPI_Wait(&req, &status);
        }

        for(int j=0; j<BASE; ++j) buckets[j].clear();
    }

    if(totNeg > 0) {
        for(int j=0; j<ndata; ++j) {
            if(data[j].key < 0) buckets[0].push_back(data[j]);
            else buckets[1].push_back(data[j]);
        }
        reverse(buckets[0].begin(), buckets[0].end());

        for(int j=0; j<2; ++j) bucketSizes[j] = buckets[j].size();
        offset[0] = 0;
        for(int j=1; j<=2; ++j) offset[j] = offset[j-1] + bucketSizes[j-1];

        assert(MPI_Gather(bucketSizes, 2, MPI_INT, allBucketSizes, 2, MPI_INT, 0, MPI_COMM_WORLD) == MPI_SUCCESS);

        int ind = 0;
        for(int j=0; j<2; ++j) {
            for(const pSort::dataType& p : buckets[j])
                extra[ind++] = p;
        }
        assert(ind == ndata);

        int base_ind = 0, proc_ind = numProcs-1, data_ind = 0;
        if(ID == 0) {
            int sz = allBucketSizes[2 * proc_ind];
            if(proc_ind > 0 && sz > 0)
                assert(MPI_Recv(recv_buf, sz, pSortType, proc_ind, 1, MPI_COMM_WORLD, &status) == MPI_SUCCESS);
            for(int j=0; j<numProcs; ++j) {
                ind = 0;
                while(ind < procN[j]) {
                    if(proc_ind == -1) {
                        assert(base_ind == 0);
                        base_ind = 1;
                        proc_ind = 0;
                    }
                    else if(data_ind == allBucketSizes[proc_ind * 2 + base_ind]) {
                        if(base_ind == 0) --proc_ind; else ++proc_ind;
                        data_ind = 0;
                        if(proc_ind == -1) continue;
                        int sz = allBucketSizes[proc_ind * 2 + base_ind];
                        if(sz > 0 && proc_ind > 0)
                            assert(MPI_Recv(recv_buf, sz, pSortType, proc_ind, 1, MPI_COMM_WORLD, &status) == MPI_SUCCESS);
                    }
                    else if(proc_ind == 0) {
                        if(j == 0)
                            data[ind++] = extra[offset[base_ind] + data_ind];
                        else
                            ans[ind++] = extra[offset[base_ind] + data_ind];
                        data_ind++;
                    }
                    else {
                        if(j == 0)
                            data[ind++] = recv_buf[data_ind++];
                        else
                            ans[ind++] = recv_buf[data_ind++];
                    }
                }
                if(j > 0)
                    assert(MPI_Send(ans, procN[j], pSortType, j, 1, MPI_COMM_WORLD) == MPI_SUCCESS);
            }
        }
        else {
            MPI_Request req;
            MPI_Irecv(data, ndata, pSortType, 0, 1, MPI_COMM_WORLD, &req);
            for(int j=0; j<2; ++j) {
                if(offset[j+1] > offset[j])
                    assert(MPI_Send(extra+offset[j], offset[j+1]-offset[j], pSortType, 0, 1, MPI_COMM_WORLD) == MPI_SUCCESS);
            }
            MPI_Wait(&req, &status);
        }
    }

    if(ID == 0) {
        delete[] allBucketSizes;
        delete[] recv_buf;
        delete[] ans;
    }

    delete[] extra;
}
