CC = g++

CFLAGS = -Wall `root-config --cflags`

CLIBS = `root-config --glibs` `root-config --cflags`

all: SabAna

SabAna: Sabat_Analysis_new.C Histo_Collection.h
	$(CC) Sabat_Analysis_new.C -o SabAna.exe $(CLIBS) -std=c++17

clean: 
	rm -rf *o

