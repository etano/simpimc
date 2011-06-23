CC=g++
CFLAGS=-I eigen3 -I /opt/local/include -c -O3 -Wall -msse3 -fomit-frame-pointer -fstrict-aliasing -ffast-math -g
LDFLAGS= -I /opt/local/include -O3 -Wall -msse3 -fomit-frame-pointer -fstrict-aliasing -ffast-math -g
LIBFLAGS = -larmadillo -framework Accelerate
SOURCES=brown-pimc.cpp
SOURCES+=SimulationClass.cpp
SOURCES+=PathClass.cpp
SOURCES+=BeadClass.cpp
SOURCES+=Moves/MoveClass.cpp
SOURCES+=Moves/BisectClass.cpp
SOURCES+=Moves/PermBisectClass.cpp
SOURCES+=Moves/DisplaceBeadClass.cpp
SOURCES+=Moves/DisplaceParticleClass.cpp
SOURCES+=Moves/DisplaceAllClass.cpp
SOURCES+=Moves/RelabelClass.cpp
SOURCES+=Moves/SimplePermClass.cpp
SOURCES+=Observables.cpp
SOURCES+=Observables/ObservableClass.cpp
SOURCES+=Observables/EnergyClass.cpp
SOURCES+=Observables/RClass.cpp
SOURCES+=Observables/R2Class.cpp
SOURCES+=Fermions.cpp
SOURCES+=Action.cpp
SOURCES+=RNGClass.cpp
SOURCES+=Stats.cpp
SOURCES+=rng/sfmt.cpp 
SOURCES+=rng/mersenne.cpp 
SOURCES+=rng/userintf.cpp
SOURCES+=rng/stoc1.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=brown-pimc

all: $(SOURCES) $(EXECUTABLE)
		
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ $(LIBFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
	
clean:
	rm -rf ./rng/*.o ./Moves/*.o ./Observables/*.o *.o $(EXECUTABLE)
	
scrubData:
	rm -rf ./data/figures/*.png ./data/traces/*.dat ./data/output/*.out
