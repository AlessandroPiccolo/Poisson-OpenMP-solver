F90C	  = f90
CC		  = gcc   
FFLAGS    = -O3 -xopenmp
CFLAGS    = -O3 -fopenmp
PWD		  = /home/alpi1674/Assigment3/build

poisson_serial:	poisson_serial.c
	$(CC) $(CFLAGS) -o poisson_serial poisson_serial.c -lm

poisson_parallel:	poisson_parallel.c
	$(CC) $(CFLAGS) -lm -o poisson_parallel poisson_parallel.c

task1:	task1.c
	$(CC) $(CFLAGS) -lm -o $(PWD)/task1 task1.c $(LIBS)
