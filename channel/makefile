PYINC = $(shell python-config --includes) $(shell python3-config --includes)

ifeq ($(LBITS),64)
    # do 64 bit stuff here, like set some CFLAGS
    CFLAGS_FPIC = -fPIC
else
    # do 32 bit stuff here
    CFLAGS_FPIC =
endif

# ---------------------------
# machine directive
# ---------------------------
#
#MACHINE = ices
#MACHINE = tacc
MACHINE = mit

######################################################################## 
ifeq ($(MACHINE),mit)
FFTLIB = -lrfftw -lfftw
HDFLIB = -lhdf5

CC = gcc -finline-functions -I$(FFTINC) -I$(HDFDIR)
CFLAGS = -Wall -O3 $(CFLAGS_FPIC) -fPIC
# CFLAGS = -Wall -g -O0 $(CFLAGS_FPIC) -fPIC
LD     = ld -shared
endif
######################################################################## 
ifeq ($(MACHINE),ices)
FFTDIR =/h2/nick/src/lib
FFTLIB = -L$(FFTDIR)/lib -lrfftw -lfftw
FFTINC = -I$(FFTDIR)/include
HDFDIR =/opt/apps/hdf5/1.8.3
HDFLIB =-L$(HDF5_DIR)/lib -lhdf5
HDFINC =-I$(HDF5_DIR)/include

CC     = gcc -Wall
CFLAGS = -O3 $(CFLAGS_FPIC) -fPIC
LD     = ld -shared
endif
######################################################################## 
ifeq ($(MACHINE),tacc)

FFTDIR =$(TACC_FFTW2_DIR)
FFTLIB =-L$(TACC_FFTW2_LIB) -lrfftw -lfftw 
FFTINC =-I$(TACC_FFTW2_INC)

HDFINC=-I$(TACC_HDF5_INC) 
HDFLIB=-L$(TACC_HDF5_LIB) -lhdf5

CC = gcc -Wall -g -O0 -fPIC

endif
######################################################################## 
SWIG = swig -python

OBJS = LegendreSetup.o \
       arrays.o \
       getMem.o \
       comp_gradient.o \
       comp_hess.o \
       increadjproject0.o \
       increadjproject.o \
       mvOps.o \
       adjproject0.o \
       adjproject.o \
       waveNumsAndCFL.o \
       pass2.o \
       increBoundary.o \
       increBoundary_s.o \
       incre_initAlphaBeta.o \
       pass1.o \
       RsQs.o \
       initAlphaBeta.o \
       initAlphaBeta2.o \
       incre_initAlphaBeta2.o \
       project.o \
       project0.o \
       increproject0.o \
       comp_stat.o \
       spec2phys.o \
       increproject.o \
       write_data2.o \
       restart2.o \
       channel.o \
       channel_wrap.o \
       tangent_forcing0.o \
       tangent_forcing.o

.SECONDARY: $(OBJS)

SOURCESS = $(addsuffix .c,$(basename $(OBJS)))

######################################################################## 
# construct exec
######################################################################## 
LIBS=$(FFTLIB) $(HDFLIB) $(OTHERLIBS) -lc -lm
INCLUDE=$(FFTINC) $(HDFINC) $(OTHERINC)

default: _c_channel.so

tar: $(SOURCES) makefile
	tar -cjvf channel.tar.bz2 $^

# clean directive 'make clean'
clean:
	- /bin/rm -f $(EXEC) *.o *.mod *_wrap.c *.so
	@echo 'files cleaned'

## build library ## 
_c_channel.so:  $(OBJS)
	$(LD) $^ $(INCLUDE) -o $@ $(LIBS) 

%.o: %.c
	@echo building $<
	$(CC) $(CFLAGS) -c $(INCLUDE) $< -o $@ $(PYINC)

%_wrap.c:	%.i
	@echo building $<
	$(SWIG) $<  

