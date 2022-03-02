#FCOMPL = af90 $(DEBUG)
FCOMPL = gfortran
#LIBS = 
MAIN=kmc_executable.out

COMPILE_FLAGS=

FFLAGS=-ffree-line-length-4000  #-O1 -O0 -D_DB ${COMPILE_FLAGS} -cpp -Wall -Wextra -fimplicit-none -fcheck=all -fbacktrace -finit-real=nan #-ffpe-trap=invalid,zero,overflow,underflow

# List libraries used by the program here
LDFLAGS=

# The possible libraries that we might want to use
LIBS=

# Suffix-rules:  Begin by throwing away all old suffix-
# rules, and then create new ones for compiling
# *.f90-bfiles.
.SUFFIXES:
.SUFFIXES: .f90 .o

.f90.o:
	$(FCOMPL) -c $(FFLAGS) $(LIBS) $<

obj = common_parameters.o lib.o events.o memory.o TPD_lib.o kmc_lib.o main.o

$(MAIN): $(obj)
	$(FCOMPL) -o $@ $(FFLAGS) $(LDFLAGS) $(obj) $(LIBS)
########################################################################################

clean:
	rm -r *.mod *.o *.out

# include file dependencies

main.o : main.f90 common_parameters.f90 events.f90 TPD_lib.f90 lib.f90 kmc_lib.f90 memory.f90
kmc_lib.o : kmc_lib.f90 lib.f90 common_parameters.f90 events.f90
TPD_lib.o : TPD_lib.f90 common_parameters.f90 events.f90 memory.f90
memory.o : memory.f90 common_parameters.f90 events.f90 lib.f90
events.o : events.f90 lib.f90 common_parameters.f90
lib.o : lib.f90 common_parameters.f90
common_parameters.o : common_parameters.f90
