PROG=solve
CXX=g++
CFLAGS=-g
LDFLAGS=

$(PROG): matrix.o box.o main.o
	$(CXX) $(CFLAGS) $^ -o $(PROG) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CFLAGS) -c $?

clean: 
	rm -f *.o $(PROG)
