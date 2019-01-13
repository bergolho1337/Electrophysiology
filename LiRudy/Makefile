CPP = g++
CPP_FLAGS = -Wall -O2

all: PRd

PRd: PRd.c
	$(CPP) $(CPP_FLAGS) -o $@ $<

clear:
	rm -r PRd

plot:
	python plot.py
