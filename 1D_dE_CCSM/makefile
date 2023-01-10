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
FLAGS = -shared -fPIC -fcray-pointer

run_crm: run_crm.f $(LIBCRM)
	@echo $(LIBCRM)
	$(FORTRAN) -o run_crm run_crm.f ./$(LIBCRM)

$(LIBCRM): ccsm3_sir_de.for
	$(FORTRAN) $(FLAGS) -o $(LIBCRM) ccsm3_sir_de.for

#modify.o: modify.f
#	$(FORTRAN) -c modify.f
