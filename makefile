# Generic Makefile for compiling a simple executable.

CC := g++ 
SRCDIR := src
BUILDDIR := build
LIB_BUILDDIR := lib_build
CFLAGS := -g -Wall -std=c++0x  -Weffc++ -Wextra -pedantic
LDFLAGS=  -lm 
TARGET := proxime 
LIB := libproxime.so 
VPATH = src
SRCEXT := cpp

SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
LIB_OBJECTS := $(patsubst $(SRCDIR)/%,$(LIB_BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
 
all:bin lib

bin:$(TARGET)

$(TARGET):$(OBJECTS)
	@echo " Linking..."; $(CC) $? -o $(TARGET) 
 
$(OBJECTS):$(SOURCES)
	@echo $@
	@mkdir -p  $(shell dirname $@)
	$(CC) $(CFLAGS) $(patsubst $(BUILDDIR)/%.o,$(SRCDIR)/%.cpp,$@) -o $@  -c


lib:$(LIB_OBJECTS)
	$(CC)  -shared -Wl,-soname,$(LIB) -o $(LIB) $? 

$(LIB_OBJECTS):
	@echo $@
	@mkdir -p  $(shell dirname $@)
	$(CC) $(CFLAGS) $(patsubst $(LIB_BUILDDIR)/%.o,$(SRCDIR)/%.cpp,$@) -o $@  -c -fPIC

clean:
	@echo " Cleaning..."; $(RM) -r $(BUILDDIR) $(TARGET) $(LIB_BUILDDIR) $(LIB)
 
.PHONY: clean
.PHONY: lib 
.PHONY: bin 
.PHONY: all 
