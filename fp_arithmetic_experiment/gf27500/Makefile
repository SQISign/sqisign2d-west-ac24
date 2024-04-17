CC = clang
CFLAGS = -W -Wextra -Wshadow -Wundef -Wunused -O2 -march=native
LD = clang
LDFLAGS =
TESTLIBS = -lgmp  # Tests use gmp for comparison purposes

all: testgf

clean:
	-rm -f gf27500.o sha3.o testgf.o testgf

testgf: gf27500.o sha3.o testgf.o
	$(LD) $(LDFLAGS) -o testgf gf27500.o sha3.o testgf.o $(TESTLIBS)

gf27500.o: gf27500.c gf27500.h
	$(CC) $(CFLAGS) -c -o gf27500.o gf27500.c

sha3.o: sha3.c sha3.h
	$(CC) $(CFLAGS) -c -o sha3.o sha3.c

testgf.o: testgf.c gf27500.h sha3.h
	$(CC) $(CFLAGS) -c -o testgf.o testgf.c
