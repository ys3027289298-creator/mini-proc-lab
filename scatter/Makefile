CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
INC := -Iinclude

all: libscatter.a

src/scatter.o: src/scatter.c include/scatter.h
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

libscatter.a: src/scatter.o
	ar rcs $@ $^

tests/test_scatter: tests/test_scatter.c libscatter.a
	$(CC) $(CFLAGS) $(INC) $< libscatter.a -o $@ -lm

tests/test_internal: tests/test_internal.c src/scatter.c include/scatter.h
	$(CC) $(CFLAGS) $(INC) $< -o $@ -lm

test: tests/test_scatter tests/test_internal
	./tests/test_scatter
	./tests/test_internal

clean:
	rm -f src/scatter.o libscatter.a tests/test_scatter tests/test_internal

.PHONY: all test clean
