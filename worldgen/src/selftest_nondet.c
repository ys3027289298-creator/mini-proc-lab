/* Test-only negative control, compiled only with -DSP_SELFTEST_NONDET.
 *
 * It reproduces a classic bug: a "global normalization" computed by splitting
 * an array into one stripe per thread, summing each stripe in float, and then
 * adding the partial sums in *thread completion order*.  Both the stripe
 * boundaries and the addition order depend on the thread count and on
 * scheduling, so the result is not reproducible.
 *
 * worldgen-nondet links this in and calls it from the hydro layer so the
 * determinism checks have a real failure to localize.  Never link this into
 * the shipping worldgen binary. */
#include <pthread.h>
#include <stddef.h>

#include "parallel.h"

typedef struct {
    const float *v;
    size_t n;
    int threads;
    float partial[SP_MAX_THREADS];
    int arrived;
    pthread_mutex_t lock;
} SumJob;

typedef struct {
    SumJob *job;
    int tid;
} SumArg;

static void *sum_worker(void *p) {
    SumArg *a = (SumArg *)p;
    float s = 0.0f;
    for (size_t i = (size_t)a->tid; i < a->job->n; i += (size_t)a->job->threads) {
        /* Scale by a non-representable factor so float rounding is visible.
         * (Flow-accumulation values are small integers, which would sum
         * exactly and hide the bug no matter how the stripes are cut.) */
        s += a->job->v[i] * 0.001f + 1e-7f * (float)(i & 7u);
    }
    pthread_mutex_lock(&a->job->lock);
    a->job->partial[a->job->arrived++] = s;
    pthread_mutex_unlock(&a->job->lock);
    return NULL;
}

/* Returns a mean that is only stable when threads == 1. */
float sp_selftest_nondet_mean(const float *v, size_t n, int threads) {
    if (threads < 1) threads = 1;
    if (threads > SP_MAX_THREADS) threads = SP_MAX_THREADS;
    if (threads == 1) {
        float s = 0.0f;
        for (size_t i = 0; i < n; i++) s += v[i];
        return n ? s / (float)n : 0.0f;
    }

    SumJob job;
    job.v = v;
    job.n = n;
    job.threads = threads;
    job.arrived = 0;
    pthread_mutex_init(&job.lock, NULL);

    SumArg args[SP_MAX_THREADS];
    pthread_t tids[SP_MAX_THREADS];
    for (int i = 0; i < threads; i++) {
        args[i].job = &job;
        args[i].tid = i;
        pthread_create(&tids[i], NULL, sum_worker, &args[i]);
    }
    for (int i = 0; i < threads; i++) pthread_join(tids[i], NULL);
    pthread_mutex_destroy(&job.lock);

    float total = 0.0f;
    for (int k = 0; k < job.arrived; k++) total += job.partial[k];
    return n ? total / (float)n : 0.0f;
}
