#include <stdio.h>
#include <stdlib.h>
#include "queuea.h"

struct queuea
{
    void ** data;
    size_t size;
    int head;
    int tail;
};

queuea * queuea_create(size_t maxarraysize)
{
    queuea * q = malloc(sizeof *q);
    q->data = malloc(maxarraysize * sizeof(void *));
    q->size = maxarraysize;
    q->head = 0;
    q->tail = 0;
}

void queuea_push(queuea * q, void * data)
{
    q->data[q->tail] = data;

    q->tail++;
    q->tail %= q->size;
}

void * queuea_pop(queuea * q)
{
    if (q->head == q->tail)
        return NULL;

    void * d = q->data[q->head];

    q->head++;
    q->head %= q->size;

    return d;
}

void queuea_free(queuea * q)
{
    free(q->data);
    free(q);
}
