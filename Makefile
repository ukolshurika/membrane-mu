CC=g++

SDIR=src
BDIR=bin

SOURCES=$(addprefix $(SDIR)/,membrane.cpp main.cpp)
HEADERS=$(wildcard $(SDIR)/*.h)
OBJECTS=$(addprefix $(BDIR)/,$(notdir $(SOURCES:.cpp=.o)))

CFLAGS=-Wall -std=c++0x -msse2 -DDEBUG -g
LDFLAGS=-lpthread -lrt
GLFLAGS=-lglut -lGL -lGLU

EXECUTABLE=membrane

$(BDIR)/$(EXECUTABLE): $(OBJECTS)
	$(CC) $^ -o $@ $(LDFLAGS)

.PHONY: clean

clean:
	-rm $(OBJECTS) $(BDIR)/$(EXECUTABLE) $(BDIR)/$(ANIMATE) $(BDIR)/$(PLOT_DATA)
