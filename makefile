# Generic Makefile for compiling a simple executable.

CC := g++ 
SRCDIR := src
BUILDDIR := build
CFLAGS := -g -Wall -std=c++0x  -Weffc++ -Wextra -pedantic
LDFLAGS=  -lm 
TARGET := proxime 
LIB := libproxime.so 
SRCEXT := cpp
INCDIR := /usr/include/proxime
LIBDIR := /usr/lib

SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
DIRS = $(shell find $(SRCDIR) -type d)

HEADERS := $(shell find $(SRCDIR) -type f -name *.h)
INCLUDE = $(patsubst $(SRCDIR)/%,$(INCDIR)/%,$(HEADERS))

vpath %.h $(DIRS)

.PHONY: install 
install: $(INCLUDE)
	@cp $(LIB) $(LIBDIR)

$(INCDIR)/%.h: %.h 
	@mkdir -p  $(shell dirname $@)
	@cp $< $@
 
.PHONY: lib 
lib:$(LIB)
$(LIB):$(OBJECTS)
	$(CC)  -shared -Wl,-soname,$(LIB) -o $(LIB) $? 

vpath %.$(SRCEXT) $(DIRS)
$(BUILDDIR)/%.o: %.$(SRCEXT) 
	@mkdir -p  $(shell dirname $@)
	$(CC) $(CFLAGS) $< -o $@  -c -fPIC

.PHONY: bin 
bin:$(TARGET)
$(TARGET):$(LIB)
	$(CC) -L. -lproxime -o $(TARGET)

.PHONY: all 
all:bin lib

.PHONY: clean
clean:
	@echo " Cleaning..."; $(RM) -r $(BUILDDIR) $(TARGET) $(LIB)
