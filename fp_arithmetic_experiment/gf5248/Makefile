# ========================================================================

CC ?= clang
SRC = $(wildcard src/*.c)

# Detect the arch
uname_p := $(shell uname -p) # store the output of the command in a variable

FLAGS = -Wall -Wextra -pedantic -std=gnu99 -O3 -march=native -mtune=native
ifeq ($(strip $(uname_p)),i386)
	FLAGS += -DTARGET_x86
else ifeq ($(strip $(uname_p)),x86_64)
	FLAGS += -DTARGET_x86
endif

# ------------------------------------------------------------------------

all: test_fp test_fp2

test_fp: src/bin/test_fp.c
	$(CC) $(SRC) src/bin/test_fp.c -o test_fp $(FLAGS)

test_fp2: src/bin/test_fp2.c
	$(CC) $(SRC) src/bin/test_fp2.c -o test_fp2 $(FLAGS)


# ========================================================================
