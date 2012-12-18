# neste template de Makefile so muda a lista
# dos SOURCES e o nome do EXECUTABLE.

CC=mpic++
CFLAGS=-c -Wall -g -std=c++0x  -Weffc++ -Wextra -pedantic
LDFLAGS=  -lm -lgsl -lgslcblas
#-lglut -lGLU -lGL
SOURCES=  src/main.cpp src/recurrence_relation.cpp src/bifurcation_diagram.cpp src/exceptions.cpp src/functions.cpp  src/numerical_integration.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=models

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC)  $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@


