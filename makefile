# neste template de Makefile so muda a lista
# dos SOURCES e o nome do EXECUTABLE.

CC=g++
CFLAGS=-c -Wall -g -std=c++0x  -Weffc++ -Wextra -pedantic
LDFLAGS=  -lm 
SOURCES=  src/main.cpp src/recurrence_relation.cpp src/bifurcation_diagram.cpp src/numerical_integration/exceptions.cpp src/functions/functions.cpp  src/numerical_integration/numerical_integration.cpp src/functions/double_pendulum.cpp  src/functions/simple_pendulum.cpp src/functions/rossler.cpp src/functions/lorenz.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=models

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC)  $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

install:
	sudo cp $(EXECUTABLE) /usr/local/bin

uninstall:
	sudo rm /usr/local/bin/$(EXECUTABLE)
