.SUFFIX:
.SUFFIXES: .o .mod .F .f .F90 .f90
F90FLAGS= -g -O2 -xhost -qopenmp -qmkl=parallel -traceback -mcmodel=medium
FC= ifort 
DIR= $1
TAG= sqv123
OBJS=pubdata.o space_info.o correlation.o update.o get_ham.o get_den.o mat_vec.o lanczos.o warm_up.o measure.o quantum.o
TAG = sqv123

all : 	$(TAG)

$(TAG):	$(OBJS) 
	$(FC)  $(OBJS) -qopenmp -qmkl=parallel -o $(TAG)  
.f90.o:
	$(FC) -c $< $(F90FLAGS)
clean: 
	rm -f *.o *.mod core 
