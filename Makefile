FC=ifort
#OPTS=-132 -O3 -parallel
OPTS=-132 -CB
#FFLAGS=-I${HOME}/include
FFLAGS=-I/home/sakuma/include
#LDFLAGS=-L${HOME}/lib
LDFLAGS=-L/home/sakuma/lib
LIBS=-lfftw3 -lguide -lpthread
#PROGS= calc_wfn calc_wannier
PROGS= calc_wannier2

OBJS=m_DATA4GW.o m_QG.o m_LMTO.o m_FFT3D.o lin3.o mymath.o myYlm.o wfn2dx.o m_MLWF.o
OBJS2=m_DATA4GW.o m_QG.o m_LMTO.o m_FFT3D.o lin3.o mymath.o myYlm.o wfnrho2.o wfn2dx_2.o m_MLWF.o cross.o

all : ${PROGS}

calc_wfn : ${OBJS} calc_wfn.o
	${FC}  ${OPTS} ${FFLAGS} ${LDFLAGS} ${OBJS} calc_wfn.o -o $@ ${LIBS}

calc_wannier : ${OBJS} calc_wannier.o
	${FC}  ${OPTS} ${FFLAGS} ${LDFLAGS} ${OBJS} calc_wannier.o -o $@ ${LIBS}

calc_wfn2 : ${OBJS2} calc_wfn2.o
	${FC}  ${OPTS} ${FFLAGS} ${LDFLAGS} ${OBJS2} calc_wfn2.o	-o $@ ${LIBS}

calc_wannier2 : ${OBJS2} calc_wannier2.o
	${FC}  ${OPTS} ${FFLAGS} ${LDFLAGS} ${OBJS2} calc_wannier2.o	-o $@ ${LIBS}


clean :
	rm -f *.o *~ ${PROGS} *.mod

.SUFFIXES: .F
.F.o   :
	${FC} ${OPTS} ${FFLAGS} $< -c
