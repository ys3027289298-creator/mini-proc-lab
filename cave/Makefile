CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: cave

cave: src/cave.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: cave
	./cave

clean:
	rm -f cave
