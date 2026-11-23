CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: thermal

thermal: src/thermal.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: thermal
	./thermal

clean:
	rm -f thermal thermal_compare.pgm
