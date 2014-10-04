PROG=solve
CXX=g++
CFLAGS=-g
LDFLAGS=

all: $(PROG)

$(PROG): matrix.o box.o main.o
	$(CXX) $(CFLAGS) $^ -o $(PROG) $(LDFLAGS)

test1: matrix.o box.o test1.o
	$(CXX) $(CFLAGS) $^ -o test1 $(LDFLAGS)

gauss: matrix.o util.o gausstest.o
	$(CXX) $(CFLAGS) $^ -o gauss $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CFLAGS) -c $?

clean: 
	rm -f *.o $(PROG)
