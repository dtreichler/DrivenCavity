#########################################
CFLAGS  =   ifort -r8 -c #-O3
LDFLAGS =   ifort -r8 -o #-O3

CMD  = ./stdsolver

##########################################
OBJ = \
 solver2.o\
 poisson.o\
 exa.o\
 subs.o

##########################################################################

.f.o    :
	$(CFLAGS) $<


$(CMD):		$(OBJ)
		@echo Makefile: ... compiling $@
		$(LDFLAGS) $(@) $(OBJ) 

##########################################################################

solver2.o :  solver2.f dimension.h Makefile
poisson.o : poisson.f dimension.h Makefile
exa.o : exa.f dimension.h Makefile
subs.o : subs.f dimension.h Makefile

