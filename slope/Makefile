CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
INC := -Iinclude -I../noise/include -I../terrain/include -I../scatter/include

LIBS := libslope.a ../terrain/libterrain.a ../scatter/libscatter.a ../noise/libnoise.a

all: tools/slope_fit_report

src/slope_fit.o: src/slope_fit.c include/slope_fit.h
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

src/slope_stats.o: src/slope_stats.c include/slope_stats.h include/slope_fit.h
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

libslope.a: src/slope_fit.o src/slope_stats.o
	ar rcs $@ $^

tools/slope_fit_report: tools/slope_fit_report.c libslope.a $(LIBS)
	$(CC) $(CFLAGS) $(INC) tools/slope_fit_report.c $(LIBS) -o $@ -lm

tests/test_slope_fit: tests/test_slope_fit.c libslope.a $(LIBS)
	$(CC) $(CFLAGS) $(INC) tests/test_slope_fit.c $(LIBS) -o $@ -lm

test: tests/test_slope_fit
	./tests/test_slope_fit

report: tools/slope_fit_report
	./tools/slope_fit_report

clean:
	rm -f src/*.o libslope.a tools/slope_fit_report tests/test_slope_fit

.PHONY: all test report clean
