FFTDIR =/h2/syang/lib/fftw
FFTLIB = $(FFTDIR)/lib
FFTINC = $(FFTDIR)/include
#HDFDIR=$TACC_HDF5_INC
#HDFDIR=/opt/apps/hdf5/1.8.3/include
#HDFSRC=$TACC_HDF5_LIB
#HDFSRC= /opt/apps/hdf5/1.8.3/lib
HDFDIR= $/opt/apps/hdf5/1.8.3/include
HDFSRC=$/opt/apps/hdf5/1.8.3/lib
GCC = gcc -Wall -g -O3 -finline-functions -I$(FFTINC) -I$(HDFDIR)
#GCC = cc -xO5 -I$(FFTINC)

OBJS1 =  LegendreSetup.o arrays.o getMem.o comp_stat2.o comp_stat.o comp_gradient.o comp_hess.o myadjointforce.o increadjproject0.o increadjproject.o mvOps.o adjproject0.o adjproject.o  waveNumsAndCFL.o pass2.o myforce.o increBoundary.o increBoundary_s.o incre_initAlphaBeta.o
OBJS2 = pass1.o   RsQs.o initAlphaBeta.o initAlphaBeta2.o incre_initAlphaBeta2.o project.o project0.o  increproject0.o increproject.o restart.o restart2.o write_data.o write_data2.o
# compute_err.o compute_err2.o
OBJS = $(OBJS1) $(OBJS2)
#all : simu2 simu

#test7 :  increadjoint_test7.o $(OBJS)
#	gcc -o test7 increadjoint_test7.o $(OBJS) -L$(FFTLIB) -L$(HDFSRC) -lrfftw -lfftw -lhdf5 -lm 
simu2 :  simu2.o $(OBJS)
	gcc -o simu2 simu2.o $(OBJS) -L$(FFTLIB) -L$(HDFSRC) -lrfftw -lfftw -lhdf5 -lm

#simu :  simulation.o $(OBJS)
#	gcc -o simu simulation.o $(OBJS) -L$(FFTLIB) -L$(HDFSRC) -lrfftw -lfftw -lhdf5 -lm

#minChnl :  main.o $(OBJS)
#	gcc -o minChnl main.o $(OBJS) -L$(FFTLIB) -L$(HDFSRC) -lrfftw -lfftw -lhdf5 -lm

#increadjoint_test7.o : increadjoint_test7.c  main.h minChnl.h arrays.h 
#	$(GCC) -c increadjoint_test7.c
simu2.o : simu2.c main.h minChnl.h arrays.h
	$(GCC) -c simu2.c

simulation.o : simulation.c main.h minChnl.h arrays.h
	$(GCC) -c simulation.c

comp_stat.o : comp_stat.c minChnl.h 
	$(GCC) -c comp_stat.c
comp_stat2.o : comp_stat2.c minChnl.h 
	$(GCC) -c comp_stat2.c

comp_gradient.o : comp_gradient.c minChnl.h 
	$(GCC) -c comp_gradient.c
comp_hess.o : comp_hess.c minChnl.h 
	$(GCC) -c comp_hess.c
main.o : main.c main.h minChnl.h arrays.h
	$(GCC) -c main.c
arrays.o : arrays.c arrays.h
	$(GCC) -c arrays.c
LegendreSetup.o : LegendreSetup.c LegendreSetup.h minChnl.h arrays.h RsQs.h
	$(GCC) -c LegendreSetup.c
getMem.o : getMem.c
	$(GCC) -c getMem.c
mvOps.o : mvOps.c minChnl.h
	$(GCC) -c mvOps.c
waveNumsAndCFL.o : waveNumsAndCFL.c minChnl.h arrays.h
	$(GCC) -c waveNumsAndCFL.c
pass1.o : pass1.c minChnl.h 
	$(GCC) -c pass1.c
pass2.o : pass2.c minChnl.h 
	$(GCC) -c pass2.c
increBoundary.o :increBoundary.c minChnl.h
	$(GCC) -c increBoundary.c
project.o : project.c minChnl.h mvOps.h
	$(GCC) -c project.c
increproject0.o : increproject0.c minChnl.h mvOps.h
	$(GCC) -c increproject0.c
increadjproject0.o : increadjproject0.c minChnl.h mvOps.h
	$(GCC) -c increadjproject0.c
increadjproject.o : increadjproject.c minChnl.h mvOps.h
	$(GCC) -c increadjproject.c
project0.o : project0.c minChnl.h mvOps.h
	$(GCC) -c project0.c
adjproject0.o : adjproject0.c minChnl.h mvOps.h
	$(GCC) -c adjproject0.c
adjproject.o : adjproject.c minChnl.h mvOps.h
	$(GCC) -c adjproject.c
increproject.o : increproject.c minChnl.h mvOps.h
	$(GCC) -c increproject.c
RsQs.o : RsQs.c minChnl.h mvOps.h
	$(GCC) -c RsQs.c
initAlphaBeta.o : initAlphaBeta.c minChnl.h mvOps.h
	$(GCC) -c initAlphaBeta.c
initAlphaBeta2.o : initAlphaBeta2.c minChnl.h mvOps.h
	$(GCC) -c initAlphaBeta2.c
incre_initAlphaBeta.o : incre_initAlphaBeta.c minChnl.h mvOps.h
	$(GCC) -c incre_initAlphaBeta.c
incre_initAlphaBeta2.o : incre_initAlphaBeta2.c minChnl.h mvOps.h
	$(GCC) -c incre_initAlphaBeta2.c
restart.o: restart.c minChnl.h 
	$(GCC) -c restart.c
restart2.o: restart2.c minChnl.h 
	$(GCC) -c restart2.c
write_data.o: write_data.c minChnl.h 
	$(GCC) -c write_data.c
write_data2.o: write_data2.c minChnl.h 
	$(GCC) -c write_data2.c

myforce.o: myforce.c  minChnl.h
	$(GCC) -c myforce.c
myadjointforce.o: myadjointforce.c  minChnl.h
	$(GCC) -c myadjointforce.c
increBoundary_s.o :increBoundary_s.c minChnl.h
	$(GCC) -c increBoundary_s.c

clean :
	rm *.o
