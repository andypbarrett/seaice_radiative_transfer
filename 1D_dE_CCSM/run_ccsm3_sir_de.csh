#
# Compile and run CCSM3_SIR_DE 
#  23 January 2007  Bruce P. Briegleb
# 
# compile, double precision:
#

pgf77 -r8 ccsm3_sir_de.for

#
# run 
#

a.out < ccsm3_sir_de_input.dat > ccsm3_sir_de_output.dat

#
# remove executible
#

rm a.out

