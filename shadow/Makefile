CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: shadow

shadow: src/shadow.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: shadow
	./shadow

export: shadow
	./shadow export

clean:
	rm -f shadow out/*.bin out/*.csv

.PHONY: all test export clean
