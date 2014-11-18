EXE = diffmed.exe
FC = /opt/solstudio12.3/bin/f95

#LIB = -framework vecLib
OPT1 = -O3
DIRMED = /opt/med-2.3.6/
DIRHDF5 = /opt/hdf5-1.6.9/
##-fortran_int64/
MEDINCLUDE = -I$(DIRMED)include
#LIBSALOME= $(DIRMED)lib/libmed.a -L$(DIRHDF5)/lib -lhdf5 -lstdc++
LIBSALOME = -L$(DIRMED)lib -lmed -L$(DIRHDF5)lib -lhdf5
LIBLAPACK = -lsunperf
LIBBLAS = 



OBJETS = \
varGlobales.o \
inout.o \
lire_med.o 

all: $(OBJETS) compa.o

	$(FC) $(OPT1) $(OBJETS) compa.o $(LIB) $(LIBSALOME) $(LIBLAPACK) $(LIBBLAS) -o $(EXE)


%.o : %.f90
	$(FC) $(OPT1) $(MEDINCLUDE) -c $< 



compa.o : $(OBJETS)
inout.o : varGlobales.o lire_med.o

clean:
	rm -f *~ *.o *.mod $(EXE)
