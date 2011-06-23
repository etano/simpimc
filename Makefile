CC=g++
CFLAGS=-I eigen3 -I /opt/local/include -c -O3 -Wall -msse3 -fomit-frame-pointer -fstrict-aliasing -ffast-math -g
LDFLAGS= -I /opt/local/include -O3 -Wall -msse3 -fomit-frame-pointer -fstrict-aliasing -ffast-math -g
LIBFLAGS = -larmadillo -framework Accelerate
SOURCES=src/simpimc.cpp
SOURCES+=src/SimulationClass.cpp
SOURCES+=src/PathClass.cpp
SOURCES+=src/BeadClass.cpp
SOURCES+=src/Moves/MoveClass.cpp
SOURCES+=src/Moves/BisectClass.cpp
SOURCES+=src/Moves/PermBisectClass.cpp
SOURCES+=src/Moves/DisplaceBeadClass.cpp
SOURCES+=src/Moves/DisplaceParticleClass.cpp
SOURCES+=src/Moves/DisplaceAllClass.cpp
SOURCES+=src/Moves/RelabelClass.cpp
SOURCES+=src/Moves/SimplePermClass.cpp
SOURCES+=src/Observables/ObservableClass.cpp
SOURCES+=src/Observables/RClass.cpp
SOURCES+=src/Observables/EnergyClass.cpp
SOURCES+=src/Observables/R2Class.cpp
SOURCES+=src/Fermions.cpp
SOURCES+=src/Action.cpp
SOURCES+=src/RNGClass.cpp
SOURCES+=src/Stats.cpp
SOURCES+=src/rng/sfmt.cpp 
SOURCES+=src/rng/mersenne.cpp 
SOURCES+=src/rng/userintf.cpp
SOURCES+=src/rng/stoc1.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=simpimc

all: $(SOURCES) $(EXECUTABLE)
		
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ $(LIBFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
	
clean:
	rm -rf ./src/rng/*.o ./src/Moves/*.o ./src/Observables/*.o ./src/*.o $(EXECUTABLE)
	
scrubData:
	rm -rf ./data/figures/*.png ./data/traces/*.dat ./data/output/*.out
