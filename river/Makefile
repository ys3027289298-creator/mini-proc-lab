CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
OUT ?= ../artifacts/river

all: river

river: src/river.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: river
	./river

generate: river
	mkdir -p $(OUT)
	./river generate $(OUT)

# 逆坡判定演示: 在内部格点 (40,30) 强制指向最陡上坡邻居, 审计应判未完成并给出位置
# 可选参数: make fault TIER=1 X=40 Y=30 DIR=-1  (DIR 罗盘: 0E 1SE 2S 3SW 4W 5NW 6N 7NE, -1=自动选最陡上坡)
TIER ?= 1
X ?= 40
Y ?= 30
DIR ?= -1
fault: river
	./river demo-fault $(TIER) $(X) $(Y) $(DIR); test $$? -eq 2

clean:
	rm -f river
