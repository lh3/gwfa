CC=			gcc
CXX=		g++
CFLAGS=		-g -Wall -O2
CXXFLAGS=	$(CFLAGS)
CPPFLAGS=	-DGWF_DEBUG
INCLUDES=
EDLIB_ROOT=	edlib
PROG=		test-edlib
LIBS=		-lz -lpthread -lm

.SUFFIXES:.c .o
.PHONY:all clean depend

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

test-edlib.o:test-edlib.c
		$(CC) -c $(CFLAGS) -I.. -I$(EDLIB_ROOT)/include $< -o $@

edlib.o:$(EDLIB_ROOT)/src/edlib.cpp
		$(CXX) -c $(CXXFLAGS) -I$(EDLIB_ROOT)/include $< -o $@

test-edlib:test-edlib.o edlib.o
		$(CXX) $(CXXFLAGS) -o $@ $^ -lz

clean:
		rm -fr gmon.out *.o a.out $(PROG) *~ *.a *.dSYM
