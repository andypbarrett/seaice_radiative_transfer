######################################################################
# Makefile for libcrm.so
#
# Also for run_crm which is used for testing
######################################################################

FORTRAN = gfortran

run_crm: run_crm.f libcrm.so
	$(FORTRAN) -o run_crm run_crm.f ./libcrm.so

libcrm.so: ccsm3_sir_de.for
	$(FORTRAN) -shared -fPIC -fcray-pointer -o libcrm.so ccsm3_sir_de.for

#modify.o: modify.f
#	$(FORTRAN) -c modify.f
