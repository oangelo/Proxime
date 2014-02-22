# neste template de Makefile so muda a lista
# dos SOURCES e o nome do EXECUTABLE.

CC=g++
CFLAGS=-c -Wall -g -Weffc++ -std=c++11 -L.. 
LDFLAGS=  -lm -lgtest -lproxime 
SOURCES= main.cpp 

OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=tests

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC)  $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

$(OBJECTS): ../libproxime.so

../libproxime.so:
	@cd .. && make all

