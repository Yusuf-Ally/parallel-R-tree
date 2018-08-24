
#ifndef _QUEUEA_H_
#define _QUEUEA_H_

struct queuea;
typedef struct queuea queuea;

queuea * queuea_create(size_t maxarraysize);
void queuea_push(queuea * q, void * data);
void * queuea_pop(queuea * q);
void queuea_free(queuea * q);

#endif // _QUEUEA_
