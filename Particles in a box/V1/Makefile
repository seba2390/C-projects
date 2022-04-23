CFLAGS = -O1 -std=gnu11
CFLAGS += $(shell /opt/homebrew/opt/gsl/bin/gsl-config --cflags)
LDLIBS += $(shell /opt/homebrew/opt/gsl/bin/gsl-config --libs)
LDLIBS += -lm
CC = gcc

.PHONEY: default

default: main run 
.PHONEY: run
run:
	./main

.PHONEY: clean
clean:
	rm main *.o *.txt *.jpg 

main: main.o functions.o

functions.o: functions.c functions.h
	$(CC) $(CFLAGS) -c functions.c $(LDLIBS) 

main.o: main.c functions.h
	$(CC) $(CFLAGS) -c main.c $(LDLIBS) 


