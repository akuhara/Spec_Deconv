FC       = gfortran 
FCFLAGS  = -Wall -pedantic -std=f95 -fbounds-check -O \
	-Wuninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace \
	-fall-intrinsics
FFTW     = -I/usr/local/include -lfftw3
BINDIR   = ./bin
TARGET   = $(BINDIR)/spec_deconv
F90_OBJS = spec_deconv.o

$(TARGET): $(F90_OBJS)
	@if [ ! -d $(BINDIR) ]; then mkdir $(BINDIR); fi
	$(FC) -o $(TARGET) $(F90_OBJS) $(FCFLAGS) $(FFTW)

$(F90_OBJS): %.o: %.f90
	$(FC) $(FCFLAGS) -c $< $(FFTW)

.PHONY: clean
clean:
	rm -f *.o
	rm -f $(BINDIR)/*
	rm -f *.mod
