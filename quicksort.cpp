#include "quicksort.h"


int g_id;

// returns the index
int selectPivot(pSort::dataType *data, int ndata) {
    int cnt = 10, idx;
    if(ndata <= 40) cnt = 4;
    set<pair<int, int>> s;
    for(int i=0;i<cnt;i++) {
        idx = rand() % ndata;
        s.insert({data[idx].key, idx});
    }
    auto it = s.begin();
    for(int i=1;i<cnt/2;i++) it++;
    return (*it).second;
}

int partition(pSort::dataType *data, int l, int r, int p_idx) {
    swap(data[l], data[p_idx]);
    p_idx = l+1;
    for(int i=l+1;i<=r;++i) {
        if((data[i].key < data[l].key) || (data[i].key == data[l].key && (rand() & 1))) {
            swap(data[p_idx], data[i]);
            p_idx++;
        }
    }
    if(p_idx == l+1) return l;
    if(p_idx == r+1) {
        swap(data[l], data[r]);
        return r;
    }
    swap(data[l], data[p_idx-1]);
    return p_idx-1;
}

int partition2(pSort::dataType *data, int l, int r, int pivot) {
    int p_idx = l;
    for(int i=l;i<=r;i++) {
        if((data[i].key < pivot) || (data[i].key == pivot && (rand() & 1))) {
            swap(data[i], data[p_idx]);
            p_idx++;
        }
    }
    return p_idx;
}

void quickSortSeq(pSort::dataType *data, int l, int r) {
    if(l >= r) return;
    bool flag = true;
    for(int i=l;i<r;i++) if(data[i].key > data[i+1].key) {flag = false; break;}
    if(flag) return;
    int p_idx = partition(data, l, r, l + selectPivot(data+l, r-l+1));
    quickSortSeq(data, l, p_idx-1);
    quickSortSeq(data, p_idx+1, r);
}

int quickSortParRec(pSort::dataType **data, pSort::dataType **extra, int ndata, int dataSize, 
        int extraSize, MPI_Comm comm, int commRank, int commSize, MPI_Datatype& pSortType, bool par) {
    
    assert(commSize > 0);
    if(commSize == 1) {
        quickSortSeq(*data, 0, ndata-1);
        if(par) return ndata;
        if(ndata > extraSize) {
            delete[] *extra;
            *extra = new pSort::dataType[ndata];
        }
        for(int i=0;i<ndata;i++) (*extra)[i] = (*data)[i];
        return ndata;
    }

    int pivot = (*data)[selectPivot(*data, ndata)].key;

    int minPivot, maxPivot;
    MPI_Allreduce(&pivot, &minPivot, 1, MPI_INT, MPI_MIN, comm);
    MPI_Allreduce(&pivot, &maxPivot, 1, MPI_INT, MPI_MAX, comm);
    pivot = (0.0 + minPivot + maxPivot) / 2;

    int p_idx = partition2(*data, 0, ndata-1, pivot);

    // for(int i=0;i<ndata;i++) cout<<data[i].key<<" "; cout<<endl;

    int p_idxArray[commSize], ndataArray[commSize];
    MPI_Allgather(&p_idx, 1, MPI_INT, p_idxArray, 1, MPI_INT, comm);
    MPI_Allgather(&ndata, 1, MPI_INT, ndataArray, 1, MPI_INT, comm);
    // p_idx = number of elements in first set. 0 <= p_idx <= ndata

    int lessTot = 0, greatTot, tot = 0, nRecv;
    for(int i=0;i<commSize;++i) {lessTot += p_idxArray[i]; tot += ndataArray[i];}
    greatTot = tot - lessTot;
    assert(lessTot > 0 && greatTot > 0);
    
    int commLsize = 1.0 * commSize * lessTot / tot, commRsize;
    if(commLsize == 0) commLsize = 1;
    else if(commLsize == commSize) --commLsize;
    commRsize = commSize - commLsize;

    vector<MPI_Request> req_v; MPI_Status status;

    if(commRank < commLsize) {
        // left half
        nRecv = lessTot / commLsize; int l = nRecv * commRank;
        if(commRank < lessTot % commLsize) nRecv++;
        l += min(lessTot % commLsize, commRank);

        if(nRecv > extraSize) {
            delete[] *extra; *extra = new pSort::dataType[nRecv];
            extraSize = nRecv;
        }
        int proc_idx = 0, st = 0, sum = 0;
        while(st < nRecv) {
            sum += p_idxArray[proc_idx];
            if(sum > l && p_idxArray[proc_idx] > 0) {
                int cnt = min({nRecv-st, sum-l, p_idxArray[proc_idx]});
                MPI_Request req;
                // cout << g_id << " " << st << " " << cnt << " " << proc_idx << endl;
                MPI_Irecv((*extra)+st, cnt, pSortType, proc_idx, 1, comm, &req);
                st += cnt; req_v.push_back(req);
            }
            proc_idx++;
        }
    }
    else {
        // right half
        commRank -= commLsize;
        nRecv = greatTot / commRsize; int l = nRecv * commRank;
        if(commRank < greatTot % commRsize) nRecv++;
        l += min(greatTot % commRsize, commRank);
        commRank += commLsize;

        if(nRecv > extraSize) {
            delete[] *extra; *extra = new pSort::dataType[nRecv];
            extraSize = nRecv;
        }
        int proc_idx = 0, st = 0, sum = 0;
        while(st < nRecv) {
            int curr = ndataArray[proc_idx] - p_idxArray[proc_idx];
            sum += curr;
            if(sum > l && curr > 0) {
                int cnt = min({nRecv-st, sum-l, curr});
                MPI_Request req;
                // cout << g_id << " " << st << " " << cnt << " " << proc_idx << endl;
                MPI_Irecv((*extra)+st, cnt, pSortType, proc_idx, 1, comm, &req);
                st += cnt; req_v.push_back(req);
            }
            proc_idx++;
        }
    }

    // send to left half
    int l = 0;
    for(int i=0;i<commRank;i++) l += p_idxArray[i];
    int proc_idx = 0, st = 0, sum = 0;
    while(st < p_idx) {
        int curr = lessTot / commLsize;
        if(proc_idx < lessTot % commLsize) curr++;
        sum += curr;
        if(sum > l && curr > 0) {
            int cnt = min({p_idx-st, sum-l, curr});
            // cout << g_id << " " << st << " " << cnt << " " << proc_idx << endl;
            MPI_Send((*data)+st, cnt, pSortType, proc_idx, 1, comm);
            st += cnt;
        }
        proc_idx++;
    }

    // send to right half
    l = 0;
    for(int i=0; i<commRank; ++i) l += ndataArray[i] - p_idxArray[i];
    proc_idx = commLsize; st = p_idx; sum = 0;
    while(st < ndata) {
        int curr = greatTot / commRsize;
        if(proc_idx - commLsize < greatTot % commRsize) curr++;
        sum += curr;
        if(sum > l && curr > 0) {
            int cnt = min({ndata-st, sum-l, curr});
            // cout << g_id << " " << st << " " << cnt << " " << proc_idx << endl;
            MPI_Send((*data)+st, cnt, pSortType, proc_idx, 1, comm);
            st += cnt;
        }
        proc_idx++;
    }

    for(auto& r : req_v) MPI_Wait(&r, &status);

    // distribution complete

    // for(int i=0;i<nRecv;i++) cout<<extra[i].key<<" "; cout<<endl;

    if(commRank < commLsize) {
        MPI_Comm_split(comm, 0, 1, &comm);
        return quickSortParRec(extra, data, nRecv, extraSize, dataSize, comm, commRank, commLsize, pSortType, par^1);
    }
    else {
        MPI_Comm_split(comm, 1, 1, &comm);
        return quickSortParRec(extra, data, nRecv, extraSize, dataSize, comm, commRank-commLsize, commRsize, pSortType, par^1);
    }
}

void quickSortPar(int procN[], int numProcs, int maxSz, pSort::dataType *data, int ndata, 
                    int ID, MPI_Datatype& pSortType) {

    g_id = ID;
    // for(int i=0;i<ndata;i++) cout<<data[i].key<<" "; cout<<endl;
    pSort::dataType *extra = new pSort::dataType[maxSz], *data1 = new pSort::dataType[ndata];
    for(int i=0;i<ndata;i++) data1[i] = data[i];
    int new_ndata = quickSortParRec(&data1, &extra, ndata, ndata, maxSz, MPI_COMM_WORLD, ID, numProcs, pSortType, false);

    // cout << ID << " " << ndata << "\n";
    // for(int i=0;i<new_ndata;i++) cout << extra[i].key << " "; cout << endl;
    
    MPI_Status status; vector<MPI_Request> rec_v;

    int newNarray[numProcs];
    MPI_Allgather(&new_ndata, 1, MPI_INT, newNarray, 1, MPI_INT, MPI_COMM_WORLD);

    int l = 0, proc_idx = 0, st = 0, sum = 0;
    for(int i=0;i<ID;++i) l += procN[i];
    while(st < ndata) {
        int curr = newNarray[proc_idx];
        sum += curr;
        if(sum > l && curr > 0) {
            int cnt = min({ndata-st, sum-l, curr});
            MPI_Request req;
            MPI_Irecv(data+st, cnt, pSortType, proc_idx, 1, MPI_COMM_WORLD, &req);
            st += cnt; rec_v.push_back(req);
        }
        proc_idx++;
    }

    l = 0; proc_idx = 0; st = 0; sum = 0;
    for(int i=0;i<ID;i++) l += newNarray[i];
    while(st < new_ndata) {
        int curr = procN[proc_idx];
        sum += curr;
        if(sum > l && curr > 0) {
            int cnt = min({new_ndata-st, sum-l, curr});
            MPI_Send(extra+st, cnt, pSortType, proc_idx, 1, MPI_COMM_WORLD);
            st += cnt;
        }
        proc_idx++;
    }

    for(auto &r : rec_v) MPI_Wait(&r, &status);

    // if(ID == 2) {for(int i=0;i<ndata;i++) cout<<data[i].key<<" "; cout<<endl;}

    delete[] extra;
    delete[] data1;
}
