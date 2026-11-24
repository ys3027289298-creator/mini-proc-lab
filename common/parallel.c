#include "parallel.h"

#include <pthread.h>
#include <stddef.h>

typedef struct {
    SpRowFn fn;
    void *ctx;
    int rows;
    int threads;
} Job;

typedef struct {
    Job *job;
    int tid;
} WorkerArg;

static void *worker(void *p) {
    WorkerArg *a = (WorkerArg *)p;
    for (int y = a->tid; y < a->job->rows; y += a->job->threads) {
        a->job->fn(a->job->ctx, y);
    }
    return NULL;
}

void sp_parallel_rows(int rows, int threads, SpRowFn fn, void *ctx) {
    if (rows <= 0) return;
    if (threads < 1) threads = 1;
    if (threads > SP_MAX_THREADS) threads = SP_MAX_THREADS;
    if (threads > rows) threads = rows;
    if (threads == 1) {
        for (int y = 0; y < rows; y++) fn(ctx, y);
        return;
    }

    Job job = {fn, ctx, rows, threads};
    WorkerArg args[SP_MAX_THREADS];
    pthread_t tids[SP_MAX_THREADS];
    int created[SP_MAX_THREADS] = {0};

    for (int i = 0; i < threads; i++) {
        args[i].job = &job;
        args[i].tid = i;
        if (pthread_create(&tids[i], NULL, worker, &args[i]) != 0) {
            /* Degrade to inline execution; static ownership keeps the output
             * identical to a fully-threaded run. */
            worker(&args[i]);
        } else {
            created[i] = 1;
        }
    }
    for (int i = 0; i < threads; i++) {
        if (created[i]) pthread_join(tids[i], NULL);
    }
}
