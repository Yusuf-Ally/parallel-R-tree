#include <stdio.h>
#include <stdlib.h>
#include "queue.h"

struct queue_item
{
    void * item;
    struct queue_item * next;
};

struct queue
{
    struct queue_item * head;
    struct queue_item * tail;
};


queue * queue_create()
{
    struct queue * q;
    q = malloc(sizeof *q);
    q->head = NULL;
    q->tail = NULL;

    return q;
}

void queue_push(queue * q, void * data)
{
    struct queue_item * item = malloc(sizeof *item);
    item->item = data;
    item->next = NULL;

    if (q->head == NULL)
        q->head = item;
    else
        q->tail->next = item;
    q->tail = item;
}

void * queue_pop(queue * q)
{
    void * ret = NULL;
    if(q->head != NULL)
    {
        struct queue_item * qitem = q->head;
        ret = qitem->item;

        q->head = qitem->next;

        if (q->head == NULL)
            q->tail = NULL;

        free(qitem);
    }

    return ret;
}

void queue_free(queue * q)
{
    while (q->head != NULL)
    {
        struct queue_item * qitem = q->head;
        q->head = q->head->next;
        free(qitem);
    }

    free(q);
}
