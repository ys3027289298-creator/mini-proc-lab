CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: rockfall

rockfall: src/rockfall.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: rockfall
	mkdir -p out
	./rockfall

clean:
	rm -f rockfall
	rm -rf out
