#/*************************************************************************
#*
#*                     MAKEFILE FOR RNS.C    - Computes the structure of one star                            
#*                                                                         
#*************************************************************************/



#######################################################################
#                                                                     #
# 1. Specify C compiler and ANSI option:                              #
#                                                                     #      
####################################################################### 

#Linux
CC=gcc

#/*************************************************************************
#*                     SPECIFY GRID SIZE
#*************************************************************************/

#STANDARD
LOWSIZE=-DMDIV=65 -DSDIV=129

#HIGH
#SIZE=-DMDIV=101 -DSDIV=201

#VERY HIGH
HIGHSIZE=-DMDIV=151 -DSDIV=301
#HIGHSIZE=-DMDIV=101 -DSDIV=201

#VERY VERY HIGH
#SIZE=-DMDIV=201 -DSDIV=401
#HIGHSIZE=-DMDIV=201 -DSDIV=401


#/*************************************************************************
#*                     COMPILING FLAGS
#*************************************************************************/


# DEBUGGING OPTION
MY_OWN =-g3 -Wall 
#MY_OWN = -g3 -Wall -DDEBUG

#MY_OWN =-O3

#/*************************************************************************
#*                    SOURCE AND OBJECT MACROS
#*************************************************************************/


SOBJ=rns.o findmodel.o equil.o equil_util.o nrutil.o stableorbit.o surface.o


#/*************************************************************************
#*                    MAIN COMPILING INSTRUCTIONS
#*************************************************************************/
rns: $(SOBJ)
	$(CC) $(MY_OWN) -lm  $(HIGHSIZE) -o rns $(SOBJ)

maxmass: $(MOBJ)
	$(CC) $(MY_OWN) -lm  $(HIGHSIZE) -o maxmass $(MOBJ)

rns.o: consts.h struct.h nrutil.h equil.h equil_util.h findmodel.h stableorbit.h rns.c makefile
	$(CC) -c $(MY_OWN) $(CFLAGS) $(COPTFLAGS) $(HIGHSIZE)  rns.c 

maxmass.o: consts.h struct.h nrutil.h equil.h equil_util.h findmodel.h maxmass.c makefile
	$(CC) -c $(MY_OWN) $(CFLAGS) $(COPTFLAGS) $(HIGHSIZE)  maxmass.c 

findmodel.o: consts.h struct.h nrutil.h equil.h equil_util.h surface.h findmodel.h findmodel.c makefile
	$(CC) -c $(MY_OWN) $(COPTFLAGS) $(HIGHSIZE)   findmodel.c

equil.o:equil.h equil_util.h nrutil.h consts.h equil.c makefile
	$(CC) -c $(MY_OWN) $(COPTFLAGS) $(HIGHSIZE)   equil.c

equil_util.o:equil_util.h nrutil.h consts.h equil_util.c makefile
	$(CC) -c $(MY_OWN) $(COPTFLAGS) $(HIGHSIZE)   equil_util.c

nrutil.o:nrutil.h nrutil.c makefile
	$(CC) -c $(COPTFLAGS)   nrutil.c

surface.o: consts.h struct.h nrutil.h equil.h equil_util.h findmodel.h surface.h stableorbit.h surface.c makefile
	$(CC) -c $(COPTFLAGS) $(HIGHSIZE)  surface.c 

stableorbit.o: consts.h struct.h nrutil.h equil.h equil_util.h findmodel.h surface.h stableorbit.h stableorbit.c makefile
	$(CC) -c $(COPTFLAGS) $(HIGHSIZE)  stableorbit.c


