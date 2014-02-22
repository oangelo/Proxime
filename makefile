# Generic Makefile for compiling a simple executable.

CC := g++ 
SRCDIR := src
BUILDDIR := build
LIB_BUILDDIR := lib_build
CFLAGS := -g -Wall -std=c++0x  -Weffc++ -Wextra -pedantic
LDFLAGS=  -lm 
TARGET := proxime 
LIB := libproxime.so 
SRCEXT := cpp
INCDIR := /usr/include/proxime
LIBDIR := /usr/lib

SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
LIB_OBJECTS := $(patsubst $(SRCDIR)/%,$(LIB_BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))

HEADERS := $(shell find $(SRCDIR) -type f -name *.h)
DIRS = $(shell find $(SRCDIR) -type d)
INCLUDE = $(patsubst $(SRCDIR)/%,$(INCDIR)/%,$(HEADERS))

vpath %.h $(DIRS)

.PHONY: install 
install: $(INCLUDE)
	@cp $(LIB) $(LIBDIR)

$(INCDIR)/%.h: %.h 
	@mkdir -p  $(shell dirname $@)
	@cp $< $@
 

$(TARGET):$(OBJECTS)
	@echo " Linking..."; $(CC) $? -o $(TARGET) 
 
$(OBJECTS):$(SOURCES)
	@echo $@
	@mkdir -p  $(shell dirname $@)
	$(CC) $(CFLAGS) $(patsubst $(BUILDDIR)/%.o,$(SRCDIR)/%.cpp,$@) -o $@  -c


.PHONY: lib 
lib:$(LIB_OBJECTS)
	$(CC)  -shared -Wl,-soname,$(LIB) -o $(LIB) $? 

$(LIB_OBJECTS):
	@echo $@
	@mkdir -p  $(shell dirname $@)
	$(CC) $(CFLAGS) $(patsubst $(LIB_BUILDDIR)/%.o,$(SRCDIR)/%.cpp,$@) -o $@  -c -fPIC

.PHONY: clean
clean:
	@echo " Cleaning..."; $(RM) -r $(BUILDDIR) $(TARGET) $(LIB_BUILDDIR) $(LIB)

.PHONY: all 
all:bin lib

.PHONY: bin 
bin:$(TARGET)
