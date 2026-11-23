CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: scatter

scatter: src/scatter.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: scatter
	./scatter

clean:
	rm -f scatter
