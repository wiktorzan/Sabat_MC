CC = g++

CFLAGS = -Wall -Iinclude `root-config --cflags`

CLIBS = `root-config --glibs` `root-config --cflags`

SRCS = src/main.C src/EnergySmearer.C src/AnalysisApp.C src/FileAnalyzer.C
OBJS = $(SRCS:.C=.o)

all: SabAna SabEneFit Analyzer

SabAna: Sabat_Analysis_new.C include/HistCollection.h
	$(CC) Sabat_Analysis_new.C -o SabAna.exe $(CLIBS) -std=c++17

SabEneFit: Sabat_Energy_Fit.C
	$(CC) Sabat_Energy_Fit.C -o SabEneFit.exe $(CLIBS) -std=c++17

Analyzer: ${OBJS} 
	$(CC) -o Analyzer.exe ${OBJS} $(CLIBS) -std=c++17

%.o: %.C
	$(CC) -c -o $@ $< $(CFLAGS)

clean: 
	rm -rf *o

