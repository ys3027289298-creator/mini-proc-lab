#ifndef SEEDSCALE_HASH_H
#define SEEDSCALE_HASH_H
#include <stdint.h>

/*
 * 白晶格哈希 H(i,j,s)。
 *
 * 构造：state = i*P1 + j*P2 + s  (mod 2^64)，z = splitmix64(state)。
 *
 * P1,P2,OCT 均为奇数 => 乘常数在模 2^64 下可逆；
 * splitmix64 终混函数是 uint64 上的双射且具备全雪崩；
 * 种子 s 只以加法进入状态字，因此：
 *   - 固定 (i,j) 时，s -> H 是双射（均匀、满置换）；
 *   - 哈希内部不存在任何长度尺度，唯一尺度由采样端的格距引入。
 */

#define SS_P1  0x9E3779B97F4A7C15ULL
#define SS_P2  0xC2B2AE3D27D4EB4FULL
#define SS_OCT 0x165667B19E3779F9ULL

/* SplitMix64 finalizer：uint64 双射，全雪崩。 */
static inline uint64_t ss_mix64(uint64_t z) {
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

/* 合格哈希：白晶格、种子满置换、无内禀尺度。 */
static inline uint64_t ss_hash2(int64_t i, int64_t j, uint64_t s) {
    uint64_t a = (uint64_t)i * SS_P1 + (uint64_t)j * SS_P2;
    return ss_mix64(a + s);
}

/* fBm 八度分离：八度号以奇数常数加法进入，仍保持每八度内的种子双射。 */
static inline uint64_t ss_hash2_oct(int64_t i, int64_t j, uint64_t s, int o) {
    uint64_t a = (uint64_t)i * SS_P1 + (uint64_t)j * SS_P2;
    return ss_mix64(a + s + (uint64_t)o * SS_OCT);
}

/* [0,1) 双精度均匀量，取高 53 位（低位不浪费给尾数）。 */
static inline double ss_u01(uint64_t h) {
    return (double)(h >> 11) * (1.0 / 9007199254740992.0); /* 2^-53 */
}

/* ---- 对照用坏哈希（仅用于检验程序的归因演示） ---- */

/* 坏法：无终混的裸线性组合 => 有代数结构、无雪崩、跨格点强相关；
   且加法型种子项在双线性插值中被整体减掉 => 换种子不换场。 */
static inline uint64_t bad_linear(int64_t i, int64_t j, uint64_t s) {
    return (uint64_t)i * SS_P1 + (uint64_t)j * SS_P2 + s;
}

/* 注：种子->尺度耦合（违反 C2）无法在纯哈希值层面构造——值噪声的尺度由
   插值核与格距决定，只要晶格值保持白，任何坐标重标号都不改几何尺度。
   真正的 C2 违反发生在管线层：种子参与采样坐标/频率。见 fieldstats 的 sx 模式。 */

#endif
