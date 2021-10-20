#
# Compile and run CCSM3_SIR_DE 
#  23 January 2007  Bruce P. Briegleb
# 
# compile, double precision:
#

gfortran -freal-4-real-8 -fcray-pointer ccsm3_sir_de.for

#
# run 
#

#a.out < ccsm3_sir_de_input.dat > ccsm3_sir_de_output.dat

#
# remove executible
#

#rm a.out

