#include "radixsort.h"


void radixSortPar(int procN[], int numProcs, int maxSz, pSort::dataType *data, int ndata, 
                    int ID, MPI_Datatype& pSortType) {

    int BASE = 16, LOGBASE = 4, intSz = 8 * sizeof(int);

    int cntNeg = 0, totNeg;
    for(int i=0; i<ndata; ++i) if(data[i].key < 0) cntNeg++;
    assert(MPI_Allreduce(&cntNeg, &totNeg, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD) == MPI_SUCCESS);
    
    vector<pSort::dataType> buckets[BASE];

    int allBucketSizes[BASE * numProcs], allBucketSizesT[BASE * numProcs], bucketSizes[BASE];
    pSort::dataType *extra = new pSort::dataType[ndata];
    MPI_Status status;

    for(int j=1;j<numProcs;j++) procN[j] += procN[j-1];

    for(int i=0; i<=intSz; i += LOGBASE) {
        int mask = ((BASE-1) << i);
        // negative integers
        if(i == intSz - LOGBASE)
            mask &= (~(1 << (intSz - 1)));
        if(i < intSz) {
            for(int j=0; j<BASE; ++j) bucketSizes[j] = 0;
            for(int j=0; j<ndata; ++j) ++bucketSizes[(data[j].key & mask) >> i];
            for(int j=0; j<BASE; ++j) buckets[j].reserve(bucketSizes[j]);
            for(int j=0; j<ndata; ++j)
                buckets[(data[j].key & mask) >> i].emplace_back(data[j]);
        }
        else {
            if(totNeg == 0) break;
            BASE = 2; bucketSizes[0] = 0;
            for(int j=0; j<ndata; ++j) if(data[j].key < 0) ++bucketSizes[0];
            bucketSizes[1] = ndata - bucketSizes[0];
            for(int j=0; j<BASE; ++j) buckets[j].reserve(bucketSizes[j]);
            for(int j=0; j<ndata; ++j) {
                if(data[j].key < 0) buckets[0].emplace_back(data[j]);
                else buckets[1].emplace_back(data[j]);
            }
        }

        assert(MPI_Allgather(bucketSizes, BASE, MPI_INT, allBucketSizes, BASE, MPI_INT, MPI_COMM_WORLD) == MPI_SUCCESS);
        for(int j=0; j<BASE; ++j)
            for(int k=0; k<numProcs; ++k)
                allBucketSizesT[j*numProcs+k] = allBucketSizes[k*BASE+j];

        for(int j=1; j<BASE*numProcs; ++j) allBucketSizesT[j] += allBucketSizesT[j-1];

        int ind = 0;
        for(int j=0; j<BASE; ++j)
            for(const pSort::dataType& p : buckets[j])
                extra[ind++] = p;
        assert(ind == ndata);

        int l = 0, a = 0, st = 0;
        if(ID > 0) l = procN[ID-1];
        while(allBucketSizesT[a] <= l) a++;
        vector<MPI_Request> req_v;
        while(st < ndata) {
            int prev = 0;
            if(a > 0) prev = allBucketSizesT[a-1];
            if(allBucketSizesT[a] > max(prev, l)) {
                int cnt = min(allBucketSizesT[a]-max(prev, l), ndata-st);
                MPI_Request req;
                MPI_Irecv(data+st, cnt, pSortType, a % numProcs, 1, MPI_COMM_WORLD, &req);
                req_v.push_back(req); st += cnt;
            }
            a++;
        }

        int offset = 0;
        for(int k=0; k<BASE; ++k) {
            if(bucketSizes[k] == 0) continue;
            // send data from [offset[k], offset[k+1])
            int r_send = allBucketSizesT[k*numProcs+ID], l_send = r_send - bucketSizes[k], l, r, s, c;
            for(int j=0; j<numProcs; ++j) {
                r = procN[j]; l = 0;
                if(j > 0) l = procN[j-1];
                if(l == r) continue;
                // [l, r), [l_send, r_send)
                if(l_send >= r) continue;
                if(l >= r_send) break;
                s = offset + max(l - l_send, 0); c = min(r_send, r) - max(l_send, l);
                assert(c > 0); assert(s+c <= ndata);
                MPI_Send(extra+s, c, pSortType, j, 1, MPI_COMM_WORLD);
            }
            offset += bucketSizes[k];
        }

        for(auto& r : req_v) MPI_Wait(&r, &status);
        for(int j=0; j<BASE; ++j) buckets[j].clear();
    }

    delete[] extra;
}
