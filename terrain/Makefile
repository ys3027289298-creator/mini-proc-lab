CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: terrain

terrain: src/terrain.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: terrain
	./terrain

clean:
	rm -f terrain
