CC=			gcc
CXX=		g++
CFLAGS=		-g -Wall -O2
CXXFLAGS=	$(CFLAGS)
CPPFLAGS=	-DGWF_DEBUG
INCLUDES=
PROG=		test-edlib test-wfa
LIBS=		-lz -lpthread -lm

EDLIB_ROOT=edlib # point to edlib source code directory
WFA_ROOT=WFA     # point to compiled WFA source code directory

.SUFFIXES:.c .o
.PHONY:all clean depend

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

test-edlib.o:test-edlib.c
		$(CC) -c $(CFLAGS) -I.. -I$(EDLIB_ROOT)/edlib/include $< -o $@

edlib.o:$(EDLIB_ROOT)/edlib/src/edlib.cpp
		$(CXX) -c $(CXXFLAGS) -I$(EDLIB_ROOT)/edlib/include $< -o $@

test-edlib:test-edlib.o edlib.o
		$(CXX) $(CXXFLAGS) -o $@ $^ -lz

test-wfa.o:test-wfa.c
		$(CC) -c $(CFLAGS) -I.. -I$(WFA_ROOT) $< -o $@

test-wfa:test-wfa.o
		$(CXX) $(CXXFLAGS) -o $@ $< $(WFA_ROOT)/build/libwfa.a -lz

clean:
		rm -fr gmon.out *.o a.out $(PROG) *~ *.a *.dSYM
