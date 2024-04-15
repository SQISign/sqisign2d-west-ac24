CC = clang
CFLAGS = -W -Wextra -Wshadow -Wundef -Wunused -O2 -march=native
LD = clang
LDFLAGS =
TESTLIBS = -lgmp  # Tests use gmp for comparison purposes

all: testgf

clean:
	-rm -f gf65376.o sha3.o testgf.o testgf

testgf: gf65376.o sha3.o testgf.o
	$(LD) $(LDFLAGS) -o testgf gf65376.o sha3.o testgf.o $(TESTLIBS)

gf65376.o: gf65376.c gf65376.h
	$(CC) $(CFLAGS) -c -o gf65376.o gf65376.c

sha3.o: sha3.c sha3.h
	$(CC) $(CFLAGS) -c -o sha3.o sha3.c

testgf.o: testgf.c gf65376.h sha3.h
	$(CC) $(CFLAGS) -c -o testgf.o testgf.c
