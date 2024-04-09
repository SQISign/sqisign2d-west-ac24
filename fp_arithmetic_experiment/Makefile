# ========================================================================

CC ?= clang
FLAGS = -Wall -Wextra -pedantic -std=gnu99 -O3 -march=native -mtune=native
SRC = $(wildcard src/*.c)
OUTPUT = fp_testing

# ------------------------------------------------------------------------

all:
	$(CC) $(SRC) -o $(OUTPUT) $(FLAGS)

clean:
	rm -r $(OUTPUT)

# ========================================================================
