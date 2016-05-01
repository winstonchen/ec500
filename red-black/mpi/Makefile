.SUFFIXES:
.SUFFIXES: .o .cpp
#============================================================
TARGET	=  heat_mpi

C_SOURCES =  main.cpp util.cpp core.cpp
C_OBJS     =  main.o util.o core.o
MY_INCLUDES =  main.h util.h core.h


CCX = mpic++
CXXFLAGS = -g -Wall
LIBDIRS = -lm

#============================================================
all: $(TARGET)

.o:.cpp	$(MY_INCLUDES)
	$(CCX)  -c  $(CXXFLAGS) $<  

$(TARGET) :   $(C_OBJS)
	$(CCX) $(CXXFLAGS)  $^ $(LIBDIRS)  -o $@

# Implicit rules: $@ = target name, $< = first prerequisite name, $^ = name of all prerequisites
#============================================================

ALL_SOURCES = Makefile $(C_SOURCES) $(MY_INCLUDES)


clean:
	rm -f $(TARGET) $(C_OBJS) core 

tar: $(ALL_SOURCES) $(NOTES)
	tar cvf $(TARGET).tar $(ALL_SOURCES)

$(TARGET).ps: $(ALL SOURCES)
	enscript -pcode.ps $(ALL_SOURCES)


