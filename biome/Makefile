CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: biome

biome: src/biome.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: biome
	./biome test

gen: biome
	mkdir -p ../artifacts/biome
	./biome gen ../artifacts/biome

clean:
	rm -f biome
