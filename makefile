# Generic Makefile for compiling a simple executable.

CC := g++ 
SRCDIR := src
BUILDDIR := build
CFLAGS := -g -Wall -std=c++0x  -Weffc++ -Wextra -pedantic
LDFLAGS=  -lm 
TARGET := models 
VPATH = src
SRCEXT := cpp

SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
 
$(TARGET): $(OBJECTS)
	@echo " Linking..."; $(CC) $? -o $(TARGET) 
 
$(OBJECTS):$(SOURCES)
	@echo $@
	@mkdir -p  $(shell dirname $@)
	@$(teste=oi;mkdir $teste;)
	$(CC) $(CFLAGS) $(patsubst $(BUILDDIR)/%.o,$(SRCDIR)/%.cpp,$@) -o $@  -c

clean:
	@echo " Cleaning..."; $(RM) -r $(BUILDDIR) $(TARGET)
 
-include $(DEPS)
 
.PHONY: clean
