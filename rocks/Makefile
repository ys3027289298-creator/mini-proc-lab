CC ?= gcc
CFLAGS ?= -O2 -std=c11 -D_DEFAULT_SOURCE -Wall -Wextra -Wpedantic
LDLIBS ?= -lm -lpng

SRC := src/gfx.c src/stones.c src/analysis.c src/rocks_main.c
HDR := src/gfx.h src/stones.h src/analysis.h

all: rocks

rocks: $(SRC) $(HDR)
	$(CC) $(CFLAGS) $(SRC) -o $@ $(LDLIBS)

test: rocks
	./rocks --check

render: rocks
	./rocks --out artifacts

clean:
	rm -f rocks

.PHONY: all test render clean
