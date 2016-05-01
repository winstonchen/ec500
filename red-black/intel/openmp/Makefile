TARGET = heatsim
LIBS = -lm -lpthread -fopenmp
CC = icpc
WARNINGS = -Wall -Wextra -Werror
CFLAGS = -g -pedantic -pedantic-errors -fstrict-aliasing -fopenmp $(WARNINGS)

SRCDIR = .
OBJDIR = obj
INCDIR = .

.PHONY: default all clean

default: $(TARGET)
all: default

OBJECTS = $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(wildcard $(SRCDIR)/*.cpp))
HEADERS = $(wildcard $(INCDIR)/*.h)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(HEADERS)
	$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) $(LIBS) -o $@

clean:
	-rm -f *.o


