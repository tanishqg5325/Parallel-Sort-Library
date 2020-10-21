#ifndef SORT_H_
#define SORT_H_

#include <bits/stdc++.h>
#include <mpi.h>

#define LOADSIZE 4


class pSort {

public:
    typedef enum {
        BEST,
        QUICK,
        MERGE,
        RADIX
    }  SortType;

    typedef struct {
        int key;
        char payload[LOADSIZE];

    } dataType;

   void init(int argc, char* argv[]);
   void close();
   void sort(dataType *data, int ndata, SortType sorter=BEST);
};

#include "mergesort.h"

#endif
