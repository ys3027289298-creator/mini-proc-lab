CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

DATA ?= ../artifacts/data

all: tunnel

tunnel: src/tunnel.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: tunnel
	./tunnel

# regenerate the CSVs consumed by tools/plot_alignment.py
data: tunnel
	mkdir -p $(DATA)
	./tunnel $(DATA)

clean:
	rm -f tunnel

.PHONY: all test data clean
