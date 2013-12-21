CC=g++
CFLAGS=-c -O3
LDFLAGS=-O3 -msse2
SOURCES=ex-stoc.cpp
SOURCES+=sfmt.cpp 
SOURCES+=mersenne.cpp 
SOURCES+=userintf.cpp
SOURCES+=stoc1.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=test

all: $(SOURCES) $(EXECUTABLE)
		
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
	
clean:
	rm -rf *.o $(EXECUTABLE)
