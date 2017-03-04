# For Ubuntu on VMWare
#ifeq ($(HOSTNAME),encyclia.seos.uvic.ca)
FC=ifort
CPP=icpc
LD=$(FC)
COMMFLAGS=-O2 
#COMMFLAGS=-g -I$(NETCDF)/include
FFLAGS=$(COMMFLAGS) 
#else
#########################################################################
# For stanhopea
#NETCDF=/opt/netcdf-3.6.0-p1
#CXXLIB=/usr/lib/gcc/i486-linux-gnu/4.0.2
#FC=g77
#LD=$(CXX)
#COMMFLAGS=-O2 -I$(NETCDF)/include
#COMMFLAGS=-g -I$(NETCDF)/include
##FFLAGS=$(COMMFLAGS) -funix-intrinsics-delete -fbadu77-intrinsics-delete #-fbounds-check
#CXXFLAGS=$(COMMFLAGS) -I../utils
#LDFLAGS=-L$(NETCDF)/lib -L$(CXXLIB) # -static
#LDLIBS=-lnetcdf_c++ -lnetcdf -lstdc++ -lfrtbegin -lg2c
#endif

OBJS=main.o convert_coord_module.o

all: easelatlon

easelatlon: $(OBJS)
	$(FC) -o easelatlon main.o convert_coord_module.o

convert_coord_module.mod: convert_coord_module.o convert_coord_module.f90
	$(FC) -c convert_coord_module.f90

convert_coord_module.o: convert_coord_module.f90
	$(FC) -c convert_coord_module.f90

main.o: convert_coord_module.mod main.f90
	$(FC) -c main.f90

.PHONY: clean

clean:
	rm -f $(OBJS) convert_coord_module.mod easelatlon
