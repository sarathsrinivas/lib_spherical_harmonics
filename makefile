CC=gcc
CFLAGS=-std=c89 -pedantic-errors -Wall
DLIBS=-lm -lgsl -lgslcblas
CSOURCES=$(wildcard *.c)
CHEADERS=$(wildcard *.h)
TEST=main_test.c

check: $(CSOURCES) $(CHEADERS)
	clang -c -std=c89 -Weverything -pedantic-errors $(CSOURCES)
	rm *.o

