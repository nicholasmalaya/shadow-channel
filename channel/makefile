
# machine directive
MACHINE = ices
#MACHINE = tacc
#MACHINE = mit

EXEC = channel_dns

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
FFTDIR = /usr
FFTLIB = $(FFTDIR)/lib
FFTINC = $(FFTDIR)/include
HDFDIR =/opt/apps/hdf5/1.8.3/
HDFSRC =-L$(HDFDIR)/lib
HDFINC =-I$(HDFDIR)/include
CC     = gcc -Wall -g -O3 -finline-functions 

#HDFDIR=$TACC_HDF5_INC
#HDFDIR=/opt/apps/hdf5/1.8.3/include
#HDFSRC=$TACC_HDF5_LIB
#HDFSRC= /opt/apps/hdf5/1.8.3/lib
#GCC = cc -xO5 -I$(FFTINC)
endif
######################################################################## 
ifeq ($(MACHINE),tacc)
FFTDIR = /usr
FFTLIB = $(FFTDIR)/lib
FFTINC = $(FFTDIR)/include
HDFDIR =/opt/apps/hdf5/1.8.3/
HDFSRC =-L$(HDFDIR)/lib
HDFINC =-I$(HDFDIR)/include
CC     = gcc -Wall -g -O3 -finline-functions 

#HDFDIR=$TACC_HDF5_INC
#HDFDIR=/opt/apps/hdf5/1.8.3/include
#HDFSRC=$TACC_HDF5_LIB
#HDFSRC= /opt/apps/hdf5/1.8.3/lib
#GCC = cc -xO5 -I$(FFTINC)
endif
######################################################################## 


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


OBJS = $(addsuffix .o,$(basename $(SOURCES)))

######################################################################## 
# construct exec
######################################################################## 

LIBS=$(FFTLIB) $(HDFLIB) $(OTHERLIBS)
INCLUDE= $(FFTINC) $(HDFINC) $(OTHERINC)

## build the executable ##
$(EXEC): $(OBJS)
	$(COMPILE) $(INCLUDE) $(OBJECTS) -o  $@ $(LDFLAGS) $(LIBS) 

# build cpp object files
%.o: %.cpp
	@echo building $< 
	$(CC) -c $(CFLAGS) $(INCLUDE) $< -o $@

# build c object files
%.o: %.c
	@echo building $<
	$(CC) -c $(CFLAGS) $(INCLUDE) $< -o $@

# clean directive 'make clean'
clean:
	- /bin/rm -f $(EXEC) *.o *.mod 
	@echo 'files cleaned'