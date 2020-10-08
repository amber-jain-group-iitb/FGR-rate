all: aout

aout: mod_fgr.o fgr.o
	ifort -o aout mod_fgr.o fgr.o -O2 -fopenmp -mkl

%.o: %.f90
	ifort -c -fopenmp $<

clean:
	rm *.o aout

