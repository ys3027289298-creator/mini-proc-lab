CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: decal

decal: src/decal.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: decal
	./decal
	python3 check_svg.py unfolding.svg

clean:
	rm -f decal unfolding.svg

.PHONY: all test clean
