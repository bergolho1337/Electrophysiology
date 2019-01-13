CC = gcc
CFLAGS = -c -Wall

all: main

main: main.o
	$(CC) -o main main.o -lm

main.o: main.c
	$(CC) $(CFLAGS) main.c -lm

clear:
	rm -rf *.o *.dat main
