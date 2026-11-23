CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
INC := -Iinclude

all: libnoise.a

src/noise.o: src/noise.c include/noise.h
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

libnoise.a: src/noise.o
	ar rcs $@ $^

tests/test_noise: tests/test_noise.c libnoise.a
	$(CC) $(CFLAGS) $(INC) $< libnoise.a -o $@ -lm

test: tests/test_noise
	./tests/test_noise

clean:
	rm -f src/noise.o libnoise.a tests/test_noise

.PHONY: all test clean
