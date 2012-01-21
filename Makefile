CXX      = g++
CXXFLAGS = -O3 -Wall -g -pg -I${HOME}/Projects/boost_1_48_0 -I/usr/bin/include
LDFLAGS  = -larmadillo -lblas -llapack -msse4

TARGET = simpimc
SRCS=src/simpimc.cpp
SRCS+=src/SimulationClass.cpp
SRCS+=src/PathClass.cpp
SRCS+=src/BeadClass.cpp
SRCS+=src/Moves/MoveClass.cpp
SRCS+=src/Moves/BisectClass.cpp
SRCS+=src/Moves/PermBisectClass.cpp
SRCS+=src/Moves/DisplaceBeadClass.cpp
SRCS+=src/Moves/DisplaceParticleClass.cpp
SRCS+=src/Moves/DisplaceAllClass.cpp
SRCS+=src/Moves/RelabelClass.cpp
SRCS+=src/Moves/SimplePermClass.cpp
SRCS+=src/Observables/ObservableClass.cpp
SRCS+=src/Observables/RClass.cpp
SRCS+=src/Observables/EnergyClass.cpp
SRCS+=src/Observables/R2Class.cpp
SRCS+=src/Fermions.cpp
SRCS+=src/Action.cpp
SRCS+=src/RNGClass.cpp
SRCS+=src/Stats.cpp
SRCS+=src/rng/sfmt.cpp 
SRCS+=src/rng/mersenne.cpp 
SRCS+=src/rng/userintf.cpp
SRCS+=src/rng/stoc1.cpp
OBJS   = $(SRCS:.cpp=.o)
DEPS   = $(SRCS:.cpp=.depends)

.PHONY: clean all

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJS) -o $(TARGET)

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.depends: %.cpp
	$(CXX) -M $(CXXFLAGS) $< > $@

clean:
	rm -f $(OBJS) $(DEPS) $(TARGET)

-include $(DEPS)
