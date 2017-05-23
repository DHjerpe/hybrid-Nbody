CC=mpicc
LD = mpicc
CFLAGS= -Ofast -Wall -std=c99 -fopenmp
LIBS=-lmpi -lm

RM = /bin/rm -f
OBJS = main.o
EXECUTABLE = nBody

all:$(EXECUTABLE)


$(EXECUTABLE): $(OBJS)
	$(LD) $(OBJS) $(LIBS) -o $(EXECUTABLE) -fopenmp

main.o: main.c main.h
	$(CC) $(CFLAGS) -c main.c -fopenmp

clean:
	$(RM) $(EXECUTABLE) $(OBJS)

