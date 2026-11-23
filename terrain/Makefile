CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: terrain stitch

terrain: src/terrain.c
	$(CC) $(CFLAGS) $< -o $@ -lm

stitch: src/stitch.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: terrain stitch
	./terrain
	./stitch

clean:
	rm -f terrain stitch
