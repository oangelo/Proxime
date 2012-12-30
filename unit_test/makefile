# neste template de Makefile so muda a lista
# dos SOURCES e o nome do EXECUTABLE.

CC=mpic++
CFLAGS=-c -Wall -g -Weffc++ -std=c++11
LDFLAGS=  -lm -lgsl -lgslcblas -lgtest  
#-lglut -lGLU -lGL
SOURCES=  ../src/recurrence_relation.cpp ../src/bifurcation_diagram.cpp ../src/numerical_integration/exceptions.cpp ../src/functions/functions.cpp  ../src/numerical_integration/numerical_integration.cpp main.cpp ../src/functions/double_pendulum.cpp ../src/functions/rossler.cpp ../src/functions/lorenz.cpp 

OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=tests

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC)  $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@


