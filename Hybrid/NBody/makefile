CC=mpicc
LD = mpicc
CFLAGS= -Ofast -Wall -std=c99 -fopenmp
LIBS=-lmpi -lm

RM = /bin/rm -f
OBJS = main.o file_operations.o
EXECUTABLE = galsim

all:$(EXECUTABLE)


$(EXECUTABLE): $(OBJS)
	$(LD) $(OBJS) $(LIBS) -o $(EXECUTABLE) -fopenmp

file_operations.o: file_operations.h file_operations.c
	$(CC) $(CFLAGS) -c file_operations.c -fopenmp

main.o: main.c file_operations.h
	$(CC) $(CFLAGS) -c main.c -fopenmp

clean:
	$(RM) $(EXECUTABLE) $(OBJS)

