# ---------------------------
# machine directive
# ---------------------------
#
MACHINE = ices
#MACHINE = tacc
#MACHINE = mit

######################################################################## 
ifeq ($(MACHINE),mit)
CC = gcc -Wall -g -O0 -finline-functions -I$(FFTINC) -I$(HDFDIR)
endif
######################################################################## 
ifeq ($(MACHINE),ices)
FFTDIR =/h2/nick/src/lib
FFTLIB = -L$(FFTDIR)/lib -lrfftw -lfftw
FFTINC = -I$(FFTDIR)/include
HDFDIR =/opt/apps/hdf5/1.8.3
HDFLIB =-L$(HDF5_DIR)/lib -lhdf5
HDFINC =-I$(HDF5_DIR)/include

PYINC = $(shell python-config --includes)

ifeq ($(LBITS),64)
    # do 64 bit stuff here, like set some CFLAGS
    CFLAGS_FPIC = -fPIC
else
    # do 32 bit stuff here
    CFLAGS_FPIC =
endif

CC     = gcc
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
       myadjointforce.o \
       increadjproject0.o \
       increadjproject.o \
       mvOps.o \
       adjproject0.o \
       adjproject.o \
       waveNumsAndCFL.o \
       pass2.o \
       myforce.o \
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
       comp_stat2.o \
       comp_stat.o \
       increproject.o \
       restart2.o \
       write_data.o \
       write_data2.o 

.SECONDARY: $(OBJS)

SOURCESS = $(addsuffix .c,$(basename $(OBJS)))

######################################################################## 
# construct exec
######################################################################## 
LIBS=$(FFTLIB) $(HDFLIB) $(OTHERLIBS) -lc -lm
INCLUDE=$(FFTINC) $(HDFINC) $(OTHERINC)

#default: _channel.so

simu2 :  simu2.o $(OBJS)
	$(CC) -o simu2 simu2.o $(OBJS) $(LIBS) -lrfftw -lfftw -lhdf5 -lm

# clean directive 'make clean'
clean:
	- /bin/rm -f $(EXEC) *.o *.mod *_wrap.c *.so
	@echo 'files cleaned'

## build library ## 
_%.so:  $(OBJS)
#	$(CC) $(INCLUDE) $(OBJS) -o $@ $(LIBS) 
	$(LD) $^ $(INCLUDE) -o $@ $(LIBS) 

%.o: %.c
	@echo building $<
	$(CC) $(CFLAGS) -c $(INCLUDE) $< -o $@ $(PYINC)

%_wrap.c:	%.i
	@echo building $<
	$(SWIG) $<  

