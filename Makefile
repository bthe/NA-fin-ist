# Compiler used
FC=gfortran

# compiler switches
FFlags=-fno-align-commons -O3 

# program name
PN=naf-ist

all : 
	$(FC) $(FFlags) src/naf-p8r.for -o $(PN)
