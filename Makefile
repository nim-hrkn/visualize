FC=gfortran 
#OPTS=-132 -O3 -parallel
OPTS=
#FFLAGS=-I${HOME}/include
FFLAGS=  -ffixed-line-length-132 -O3   -fopenmp
#LDFLAGS=-L${HOME}/lib
LDFLAGS=
LIBS=-lfftw3 -lguide -lpthread
#PROGS= calc_wfn calc_wannier
PROGS= calc_wannier6

OBJS=m_DATA4GW.o m_QG.o m_LMTO.o m_FFT3D.o lin3.o mymath.o myYlm.o wfn2dx.o m_MLWF.o
OBJS2=m_DATA4GW.o m_QG.o m_LMTO.o m_FFT3D.o lin3.o mymath.o myYlm.o wfnrho2.o wfn2dx_2.o m_MLWF.o cross.o cubeformat.o xsfformat2.o  keyvalue.o  expand_mesh.o  wfn2dx_abc.o wfnrho_abc.o write_mesh2d.o 

all : ${PROGS}

calc_wfn : ${OBJS} calc_wfn.o
	${FC}  ${OPTS} ${FFLAGS} ${LDFLAGS} ${OBJS} calc_wfn.o -o $@ ${LIBS}

calc_wannier : ${OBJS} calc_wannier.o
	${FC}  ${OPTS} ${FFLAGS} ${LDFLAGS} ${OBJS} calc_wannier.o -o $@ ${LIBS}

calc_wfn2 : ${OBJS2} calc_wfn2.o
	${FC}  ${OPTS} ${FFLAGS} ${LDFLAGS} ${OBJS2} calc_wfn2.o	-o $@ ${LIBS}

calc_wannier2 : ${OBJS2} calc_wannier2.o
	${FC}  ${OPTS} ${FFLAGS} ${LDFLAGS} ${OBJS2} calc_wannier2.o	-o $@ ${LIBS}
calc_wannier3 : ${OBJS2} calc_wannier3.o
	${FC}  ${OPTS} ${FFLAGS} ${LDFLAGS} ${OBJS2} calc_wannier3.o	-o $@ ${LIBS}
calc_wannier4 : ${OBJS2} calc_wannier4.o
	${FC}  ${OPTS} ${FFLAGS} ${LDFLAGS} ${OBJS2} calc_wannier4.o	-o $@ ${LIBS}
calc_wannier5 : ${OBJS2} calc_wannier5.o
	${FC}  ${OPTS} ${FFLAGS} ${LDFLAGS} ${OBJS2} calc_wannier5.o	-o $@ ${LIBS}
calc_wannier6 : ${OBJS2} calc_wannier6.o
	${FC}  ${OPTS} ${FFLAGS} ${LDFLAGS} ${OBJS2} calc_wannier6.o	-o $@ ${LIBS}


calc_wannier3.o: cubeformat.o xsfformat.o 
calc_wannier4.o: cubeformat.o xsfformat.o 
calc_wannier5.o: cubeformat.o xsfformat.o  keyvalue.o 
calc_wannier6.o: cubeformat.o xsfformat.o  keyvalue.o 

clean :
	rm -f *.o *~ ${PROGS} *.mod

.SUFFIXES: .F
.F.o   :
	${FC} ${OPTS} ${FFLAGS} $< -c
