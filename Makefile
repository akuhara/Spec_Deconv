FC       = gfortran 
FCFLAGS  = -Wall -pedantic -std=f95 -fbounds-check -O \
	-Wuninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace \
	-fall-intrinsics
FFTW     = -I/opt/intel/composer_xe_2013_sp1.3.174/mkl/include/fftw -lfftw3
BINDIR   = ./bin
TARGET   = $(BINDIR)/aw_rf
F90_OBJS = aw_rf.o

$(TARGET): $(F90_OBJS)
	$(FC) -o $(TARGET) $(F90_OBJS) $(FCFLAGS) $(FFTW)

$(F90_OBJS): %.o: %.f90
	$(FC) $(FCFLAGS) -c $< $(FFTW)

.PHONY: clean
clean:
	rm -f *.o
	rm -f $(BINDIR)/*
	rm -f *.mod
