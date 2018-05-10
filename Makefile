CC 	= icc -Kc++ -parallel
CFLAGS	= -parallel -openmp
D_UTIL	= ../util/
WARN    = -Wall -Wextra

O_UTIL	= util.o integration.o spline.o
O_ALL	= $(O_UTIL) table_input.o

L_FFTW	= -openmp -lfftw3_omp -lfftw3

L_FFTWF	= -lfftw3f -lfftw3f_threads -lpthread

L_GSL   = -lgsl -lgslcblas

calc_dens_field: $(O_ALL) calc_dens_field.o header.h calc_dens_field.cpp
	$(CC)  $(CFLAGS) -o calc_dens_field calc_dens_field.o $(O_ALL) $(L_FFTW) -lm

main: main.cpp densityfield.h countlines.h powerspec.h physconst.h growth.h cosmoparams.h
	$(CC)  $(CFLAGS) -o main main.cpp $(L_FFTW) -lm -lgsl -lgslcblas

ftempFFTW: ftempFFTW.cpp densityfield.h countlines.h powerspec.h physconst.h growth.h cosmoparams.h
	$(CC)  $(CFLAGS) -o ftempFFTW ftempFFTW.cpp $(O_ALL) $(L_FFTW) $(L_GSL) -lm
