CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: hydrology

hydrology: src/hydrology.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: hydrology
	mkdir -p out
	./hydrology

clean:
	rm -f hydrology
	rm -rf out
