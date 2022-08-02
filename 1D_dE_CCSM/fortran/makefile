######################################################################
# TEST_PARAMETERS make file
######################################################################

FORTRAN = gfortran

test_parameters: test_parameters.o get_parameters.o print_parameters.o modify.o
	$(FORTRAN) -o test_parameters test_parameters.o get_parameters.o print_parameters.o modify.o

test_parameters.o: test_parameters.f
	$(FORTRAN) -c test_parameters.f

get_parameters.o: get_parameters.f
	$(FORTRAN) -c get_parameters.f

print_parameters.o: print_parameters.f
	$(FORTRAN) -c print_parameters.f

modify.o: modify.f
	$(FORTRAN) -c modify.f
