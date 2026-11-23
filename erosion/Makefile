CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: erosion

erosion: src/erosion.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: erosion
	./erosion --n 64 --drops 4000 --iters 3 --seed 7 --out /tmp/erosion_selftest

clean:
	rm -f erosion

.PHONY: all test clean
