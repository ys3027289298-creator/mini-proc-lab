CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: puddle

puddle: src/puddle.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: puddle
	./puddle

clean:
	rm -f puddle
