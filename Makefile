CC 	= icc -Kc++
WARN    = -Wall -Wextra

O_ALL	= integration.o table_input.o

L_FFTW	= -lfftw3

L_FFTWF	= -lfftw3f -lfftw3f_threads -lpthread

L_GSL   = -lgsl -lgslcblas

calc_dens_field: $(O_ALL) calc_dens_field.o header.h calc_dens_field.cpp
	$(CC)  -o calc_dens_field calc_dens_field.o $(O_ALL) $(L_FFTW) -lm

main: main.cpp densityfield.h countlines.h powerspec.h physconst.h growth.h cosmoparams.h
	$(CC)  -o main main.cpp $(L_FFTW) -lm -lgsl -lgslcblas

ftempFFTW: ftempFFTW.cpp densityfield.h countlines.h powerspec.h physconst.h growth.h cosmoparams.h
	$(CC)  -o ftempFFTW ftempFFTW.cpp $(O_ALL) $(L_FFTW) $(L_GSL) -lm
