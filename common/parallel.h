#ifndef SP_PARALLEL_H
#define SP_PARALLEL_H

#define SP_MAX_THREADS 64

/* Runs fn(ctx, row) exactly once for every row in [0, rows).
 *
 * Rows are statically owned by thread id (row % threads == tid).  Each
 * callback must write only its own row and must read either immutable data
 * or the previous pass's buffer.  Under that contract the output does not
 * depend on the thread count or on scheduling.  There is no cross-row
 * accumulation in this helper by design; a reduction would need an explicit
 * fixed-order tree. */
typedef void (*SpRowFn)(void *ctx, int row);

void sp_parallel_rows(int rows, int threads, SpRowFn fn, void *ctx);

#endif
