CC=			gcc
CXX=		g++
CFLAGS=		-g -Wall -O2
CXXFLAGS=	$(CFLAGS)
CPPFLAGS=	-DGWF_DEBUG
INCLUDES=
PROG=		test-wf-ed test-edlib test-wfa
LIBS=		-lz -lpthread -lm

EDLIB_ROOT=edlib
WFA_ROOT=WFA

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

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

test-wf-ed:wf-ed.c
		$(CC) $(CFLAGS) -I.. -DWF_ED_MAIN -I$(WFA_ROOT) $^ -o $@ -lz

clean:
		rm -fr gmon.out *.o a.out $(PROG) *~ *.a *.dSYM
