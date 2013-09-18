# ---------------------------
# machine directive
# ---------------------------
#
MACHINE = ices
#MACHINE = tacc
#MACHINE = mit

######################################################################## 
ifeq ($(MACHINE),mit)
FFTDIR = /usr
FFTLIB = $(FFTDIR)/lib
FFTINC = $(FFTDIR)/include
HDFDIR= $/opt/apps/hdf5/1.8.3/include
HDFSRC= $/opt/apps/hdf5/1.8.3/lib
CC = gcc -Wall -g -O3 -finline-functions -I$(FFTINC) -I$(HDFDIR)
LIBS =  -L$(FFTLIB) -L$(HDFSRC)
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

CC     = gcc -g -O3 -finline-functions
LD     = ld
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

SOURCES = LegendreSetup.c \
          arrays.c \
          getMem.c \
          comp_gradient.c \
          comp_hess.c \
          myadjointforce.c \
          increadjproject0.c \
          increadjproject.c \
          mvOps.c \
          adjproject0.c \
          adjproject.c \
          waveNumsAndCFL.c \
          pass2.c \
          myforce.c \
          increBoundary.c \
          increBoundary_s.c \
          incre_initAlphaBeta.c \
          pass1.c \
          RsQs.c \
          initAlphaBeta.c \
          initAlphaBeta2.c \
          incre_initAlphaBeta2.c \
          project.c \
          project0.c \
          increproject0.c \
          increproject.c \
          restart.c \
          restart2.c \
          write_data.c \
          write_data2.c \
          init.c   

OBJS = $(addsuffix .o,$(basename $(SOURCES)))

######################################################################## 
# construct exec
######################################################################## 
LIBS=$(FFTLIB) $(HDFLIB) $(OTHERLIBS) -lc -lm
INCLUDE=$(FFTINC) $(HDFINC) $(OTHERINC)

default: _channel.so

tar: $(SOURCES) makefile
	tar -cjvf channel.tar.bz2 $^

# clean directive 'make clean'
clean:
	- /bin/rm -f $(EXEC) *.o *.mod *_wrap.c
	@echo 'files cleaned'

## build library ## 
_%.so:  $(OBJS)
#	$(CC) $(INCLUDE) $(OBJS) -o  $@ $(LIBS) 
	$(LD) $^ $(INCLUDE) -o $@ $(LIBS) 

%.o: %.c
	@echo building $<
	@echo $(PYINC)
	$(CC) -c $(INCLUDE) $< -o $@ $(PYINC)

%_wrap.c:	%.i
	@echo building $<
	$(SWIG) $<  

