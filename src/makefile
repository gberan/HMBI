# Set up different object names for serial and parallel versions
SER_OBJS = main_ser.o monomer_ser.o atom_ser.o dimer_ser.o cluster_ser.o cluster_ewald_ser.o params_ser.o \
vector_ser.o matrix_ser.o opt_ser.o multipole_ser.o polarizability_ser.o \
dlf_interface_ser.o 

PAR_OBJS = main_par.o monomer_par.o atom_par.o dimer_par.o cluster_par.o cluster_ewald_ser.o params_par.o \
vector_par.o matrix_par.o opt_par.o multipole_par.o polarizability_par.o \
dlf_interface_par.o 

# Serial version macros
CC = g++
#CC = icpc
FC = gfortran

# Parallel version macros
MPICC = mpiCC
PARALLEL = -DPARALLEL

# Macros that apply to all versions
LDFLAGS =  -L/opt/intel/fce/9.1.052/lib/ -lifcore -limf -L/opt/intel/mkl/9.1.023/lib/em64t -lguide -lmkl -lmkl_lapack -lmkl_em64t -lpthread -lifport
DEBUG = -g
CFLAGS = ${DEBUG} -Wno-deprecated 

# Begin targets

serial: dlfind ${SER_OBJS}
	${CC} ${CFLAGS} ${INCLUDES} ${SER_OBJS} ${LDFLAGS} lib_dlfind.a  -o hmbi.serial

parallel: ${PAR_OBJS}
	${MPICC} ${CFLAGS} ${INCLUDES} ${PAR_OBJS} ${LDFLAGS} lib_dlfind.a -o hmbi.parallel

all: dlfind serial parallel 


install: all links
	mv hmbi.serial ../bin
	mv hmbi.parallel ../bin

links:
	ln -sf ../bin/hmbi.run ../bin/hmbi


clean:
	rm -f *.o ../bin/hmbi.serial ../bin/hmbi.parallel ../bin/hmbi
	cd lib_dlfind && $(MAKE) clean

dlfind: 
	cd lib_dlfind && $(MAKE)


# Notation: $@ = target name, $< = source name
%_ser.o: %.C atom.h cluster.h constants.h dimer.h matrix.h monomer.h multipole.h \
	params.h polarizability.h vector.h dlf_interface.h 
	 $(CC) -c -o $@ $(CFLAGS)  $<

%_par.o: %.C atom.h cluster.h constants.h dimer.h matrix.h monomer.h multipole.h \
	params.h polarizability.h vector.h dlf_interface.h  
	 $(MPICC) -c -o $@ $(CFLAGS) $(PARALLEL) $<

%.o: %.C atom.h cluster.h constants.h dimer.h matrix.h monomer.h multipole.h \
	params.h polarizability.h vector.h dlf_interface.h 
	 $(CC) -c -o $@ $(CFLAGS) $<
