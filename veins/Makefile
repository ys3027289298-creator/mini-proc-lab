CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: veins

veins: src/veins.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: veins
	./veins

clean:
	rm -f veins
