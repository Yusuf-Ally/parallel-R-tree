
#ifndef _QUEUE_H_
#define _QUEUE_H_

struct queue;
typedef struct queue queue;

queue * queue_create();
void queue_push(queue * q, void * data);
void * queue_pop(queue * q);
void queue_free(queue * q);

#endif // _QUEUE_

