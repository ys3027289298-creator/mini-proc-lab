CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
INC := -Iinclude -I../noise/include

all: libterrain.a

src/terrain.o: src/terrain.c include/terrain.h ../noise/include/noise.h
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

libterrain.a: src/terrain.o
	ar rcs $@ $^

tests/test_terrain: tests/test_terrain.c libterrain.a ../noise/libnoise.a
	$(CC) $(CFLAGS) $(INC) $< libterrain.a ../noise/libnoise.a -o $@ -lm

test: tests/test_terrain
	./tests/test_terrain

clean:
	rm -f src/terrain.o libterrain.a tests/test_terrain

.PHONY: all test clean
