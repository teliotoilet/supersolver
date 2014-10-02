PROG=solve
CXX=g++
CFLAGS=-O2
LDFLAGS=

$(PROG): main.o
	$(CXX) $(CFLAGS) $^ -o $(PROG) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CFLAGS) -c $?

clean: 
	rm *.o $(PROG)
