CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
LDLIBS ?= -lm

all: biome

biome: src/biome.c
	$(CC) $(CFLAGS) $< -o $@ $(LDLIBS)

test: biome
	./biome

clean:
	rm -f biome

.PHONY: all test clean
