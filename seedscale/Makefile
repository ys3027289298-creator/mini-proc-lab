CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: hash_test fieldstats

hash_test: src/hash_test.c src/hash.h src/stats.h
	$(CC) $(CFLAGS) src/hash_test.c -o $@ -lm

fieldstats: src/fieldstats.c src/hash.h src/stats.h
	$(CC) $(CFLAGS) src/fieldstats.c -o $@ -lm

out:
	mkdir -p out

test: hash_test
	./hash_test good
	./hash_test bad-linear >/dev/null 2>&1; test $$? -eq 1

stats: all out
	./fieldstats good vn out/good_vn.csv
	./fieldstats good vn out/good_vn_fbm.csv fbm
	./fieldstats good nn out/good_nn.csv
	./fieldstats good sx out/good_sx.csv
	./fieldstats bad-linear vn out/bad_linear_vn.csv

clean:
	rm -f hash_test fieldstats
	rm -rf out
