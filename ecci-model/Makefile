FOPTS = -std=f95 -O2 -g -fcheck=all -Wall
FC = gfortran

OBJS = double.o zvode.o  betaEdgeHex_g1_240.o  ODE.o cases.o doHW.o
MODULES = doublep.mod cases.mod odes.mod 
LIBS = -lblas -llinpack
LIBPATH = -L./linpack

.PHONY: all clean 

all: main

main: $(OBJS)
	$(FC) $(FOPTS) $(OBJS) $(LIBPATH) -o $@ $(LIBS)

%.o: %.f90 
	$(FC) -c $<

clean:
	rm -f main $(OBJS) $(MODULES) 
