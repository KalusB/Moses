CC 	= g++ -std=c++98

O_ALL	= integration.o table_input.o

L_FFTW	= -lfftw3

L_FFTWF	= -lfftw3f -lfftw3f_threads -lpthread

L_GSL   = -lgsl -lgslcblas

all: all_cpp all_pyx

all_cpp: calc_dens_field main ftempFFTW

calc_dens_field: $(O_ALL) calc_dens_field.cpp header.h
	$(CC)  -o calc_dens_field calc_dens_field.cpp $(O_ALL) $(L_FFTW) -lm

main: main.cpp densityfield.h countlines.h powerspec.h cosmoparams.h
	$(CC)  -o main main.cpp $(L_FFTW) -lm -lgsl -lgslcblas

ftempFFTW: ftempFFTW.cpp densityfield.h countlines.h powerspec.h cosmoparams.h
	$(CC)  -o ftempFFTW ftempFFTW.cpp $(O_ALL) $(L_FFTW) $(L_GSL) -lm

all_pyx: GalaxyvsContaminant.pyx
	python setup.py build_ext --inplace

clear:
	rm -r calc_dens_field ftempFFTW GalaxyvsContaminant.c *.so *.o main build
