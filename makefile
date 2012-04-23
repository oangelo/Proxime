# neste template de Makefile so muda a lista
# dos SOURCES e o nome do EXECUTABLE.

CC=mpic++
CFLAGS=-c -Wall -g 
LDFLAGS=  -lm -lgsl -lgslcblas -lUnitTest++ 
#-lglut -lGLU -lGL
SOURCES=  src/main.cpp src/recurrence_relation.cpp src/bifurcation_diagram.cpp src/Exceptions.cpp src/functions.cpp  src/Numerical_Integration.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=chaos

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC)  $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@


