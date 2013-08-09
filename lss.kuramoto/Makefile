CC = gcc
PY_INC = /usr/include/python2.7
#CFLAGS = -O0 -g --std=c99 -fPIC -I$(PY_INC)
CFLAGS = -O3 --std=c99 -fPIC -I$(PY_INC)
SWIG = swig -python
LD = ld -shared
LIBS = -lm -lblas

default:	_kuramoto.so

clean:
	rm -f *.so *.o *_wrap.c

_%.so:	%.o %_wrap.o
	$(LD) $^ -o $@ $(LIBS)

%.o:	%.c
	$(CC) $(CFLAGS) -c $< -o $@

%_wrap.c:	%.i
	$(SWIG) $<
