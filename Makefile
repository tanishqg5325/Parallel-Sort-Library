EXE := libpsort.a
SRC := sort.cpp mergesort.cpp radixsort.cpp quicksort.cpp
INC := $(wildcard *.h)
OBJ := $(SRC:%.cpp=%.o)

CPPFLAGS := -O3

.PHONY: all clean

all: $(EXE)

$(EXE): $(OBJ)
	ar rvs $@ $^

%.o: %.cpp $(INC)
	/opt/bin/mpicxx -c $< $(CPPFLAGS) -o $@

clean:
	rm *.o
	rm $(EXE)
