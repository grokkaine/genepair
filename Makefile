PROGRAM = genepair
CC = gcc
LIBS = -lm

default: all

all: src/genepair_f.o src/genepair.o
	mkdir bin
	$(CC) src/genepair_f.o src/genepair.o $(LIBS) -o bin/$(PROGRAM)

genepair_f.o: src/genepair_f.c
	$(CC) -c src/genepair_f.c -o src/genepair_f.o

genepair.o: src/genepair.c
	$(CC) -c src/genepair.c -o src/genepair.o

.PHONY: clean

clean:
	rm -rf src/*.o src/*\~ ./*\~ ./bin example/*output.txt

