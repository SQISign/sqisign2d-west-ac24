# ========================================================================

CC ?= clang
SRC = $(wildcard src/*.c)
OUTPUT = fp_testing

# Detect the arch
uname_p := $(shell uname -p) # store the output of the command in a variable

FLAGS = -Wall -Wextra -pedantic -std=gnu99 -O3 -march=native -mtune=native
ifeq ($(strip $(uname_p)),i386)
	FLAGS += -DTARGET_x86
endif

# ------------------------------------------------------------------------

all:
	$(CC) $(SRC) -o $(OUTPUT) $(FLAGS)

clean:
	rm -r $(OUTPUT)

# ========================================================================
