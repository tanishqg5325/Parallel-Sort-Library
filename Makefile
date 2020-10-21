all:
	/opt/bin/mpicxx tester.cpp sort.cpp mergesort.cpp radixsort.cpp

clean:
	rm a.out
