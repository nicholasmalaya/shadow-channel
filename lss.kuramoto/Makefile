LBITS := $(shell getconf LONG_BIT)
ifeq ($(LBITS),64)
    # do 64 bit stuff here, like set some CFLAGS
    CFLAGS_FPIC = -fPIC
else
    # do 32 bit stuff here
    CFLAGS_FPIC =
endif

CC = gcc
PY_INC = $(shell python-config --includes)
#CFLAGS = -O0 -g --std=c99 $(PY_INC)
CFLAGS = -O3 $(CFLAGS_FPIC) --std=c99 $(PY_INC)
SWIG = swig -python
LD = ld -shared
LIBS = -lm -lblas

default:	_kuramoto.so

test:	testKura

testPrimal:   testKSprimal.o kuramoto.o
	$(CC) $^ -o $@ $(LIBS)

tar:	testKura.c kuramoto.c kuramoto.i kuramoto.h numpy.i wrapper.py Makefile
	tar -cjvf kuramoto.tar.bz2 $^

testKura:	testKura.o kuramoto.o
	$(CC) $^ -o $@ $(LIBS)

clean:
	rm -f *.so *.o *_wrap.c

_%.so:	%.o %_wrap.o
	$(LD) $^ -o $@ $(LIBS)

%.o:	%.c
	$(CC) $(CFLAGS) -c $< -o $@

%_wrap.c:	%.i
	$(SWIG) $<
