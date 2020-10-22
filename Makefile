all:
	/opt/bin/mpicxx -O3 tester.cpp sort.cpp mergesort.cpp radixsort.cpp quicksort.cpp

clean:
	rm a.out
