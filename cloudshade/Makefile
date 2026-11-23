CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: cloudshade

cloudshade: src/cloudshade.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: cloudshade
	./cloudshade

clean:
	rm -f cloudshade
