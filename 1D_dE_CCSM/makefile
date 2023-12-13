######################################################################
# Makefile for libcrm.so
#
# Also for run_crm which is used for testing
######################################################################

ARCH := $(shell uname)
ifeq ($(ARCH), Linux)
	LIBCRM = libcrm.so
else ifeq ($(ARCH), Darwin)
	LIBCRM = libcrm.dylib
endif

FORTRAN = gfortran

build: $(LIBCRM)
	mkdir ../seaicert/lib
	mv $(LIBCRM) ../seaicert/lib

run_crm: run_crm.f $(LIBCRM)
	@echo $(LIBCRM)
	$(FORTRAN) -o run_crm run_crm.f ./$(LIBCRM)

$(LIBCRM): ccsm3_sir_de.for
	$(FORTRAN) -shared -fPIC -fcray-pointer -fallow-argument-mismatch -o $(LIBCRM) ccsm3_sir_de.for


#modify.o: modify.f
#	$(FORTRAN) -c modify.f
