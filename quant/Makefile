CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: quant

quant: src/quant.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: quant
	./quant

clean:
	rm -f quant
