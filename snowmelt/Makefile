CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

OUT ?= ../artifacts/snowmelt

all: snowmelt

snowmelt: src/snowmelt.c
	$(CC) $(CFLAGS) $< -o $@ -lm

# 自检 + 非对称统计；同时把渲染网格写到 $(OUT)/
test: snowmelt
	./snowmelt $(OUT)

# 渲染因子场/非对称统计/三条雪线对比图到 $(OUT)/
# 注：f32 网格由本模块 C 程序自产（非外部下载数据），故不加 -I 以便使用用户站点包
plot: test
	python3 plot/plot_snowmelt.py $(OUT) $(OUT)

clean:
	rm -rf snowmelt
