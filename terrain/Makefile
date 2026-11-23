CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: terrain

terrain: src/terrain.c
	$(CC) $(CFLAGS) $< -o $@ -lm

slope_grid: src/slope_grid.c
	$(CC) $(CFLAGS) $< -o $@ -lm

analysis: slope_grid
	./slope_grid

test: terrain
	./terrain

clean:
	rm -f terrain slope_grid
