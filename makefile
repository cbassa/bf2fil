# Compiling flags
CFLAGS = -c -g -O3 -I/opt/lofardal/include

# Linking flags
LFLAGS = -L/opt/lofardal/lib -llofardal -lhdf5 -lgsl -lgslcblas

# Compiler
CPP = g++
CC = gcc

bf2fil: bf2fil.o h5read.o filwrite.o filterbank_header.o send_stuff.o swap_bytes.o strings_equal.o error_message.o skz.o
	$(CPP) -o bf2fil bf2fil.o h5read.o filwrite.o filterbank_header.o send_stuff.o swap_bytes.o strings_equal.o error_message.o skz.o $(LFLAGS)

%.o : %.cpp
	$(CPP) $(CFLAGS) $< -o $@

%.o: %.c
	$(CC) -c -g -O3 $< -o $@	

clean:
	rm -f *.o
	rm -f *~
