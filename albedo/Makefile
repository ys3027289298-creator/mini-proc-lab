CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: albedo

albedo: src/albedo.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: albedo
	./albedo test

build: albedo
	./albedo build

clean:
	rm -f albedo out/*.csv out/*.ppm
