CC=gcc
OMP_PARALLEL=TRUE
LIB_DIR=/cm/shared/apps/gromacs/4.6.7/lib
INC_DIRS=/cm/shared/apps/gromacs/4.6.7/include/gromacs/


ifeq ($(OMP_PARALLEL), TRUE) 
CFLAGS=-O3 -Wall -fopenmp -DOMP_PARALLEL=TRUE
LIBFLAGS=-lgmx -ldl -lm -pthread -lgomp
else
CFLAGS=-O3 -Wall
LIBFLAGS=-lgmx -ldl -lm 
endif

OBJECTS=gmx_efield.o

all: $(OBJECTS)
	$(CC) $(OBJECTS) -L $(LIB_DIR) $(LIBFLAGS) -o g_efield

gmx_efield.o: gmx_efield.c
	$(CC) $(CFLAGS) -c -I $(INC_DIRS) gmx_efield.c

clean:
	rm *.o
