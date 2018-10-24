#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <pthread.h>
#include <semaphore.h>

#ifndef Thread_error_wait
#define Thread_error_wait 5
#endif


// thread
void xperror(int en, const char *msg);

#define xpthread_join(thread, retval) _xpthread_join(thread, retval, __LINE__, __FILE__)
#define xpthread_create(thread, attr, start_routine, arg) _xpthread_create(thread, attr, start_routine, arg, __LINE__, __FILE__)
int _xpthread_create(pthread_t *thread, const pthread_attr_t *attr,
                          void *(*start_routine) (void *), void *arg, int linea, const char *file);
int _xpthread_join(pthread_t thread, void **retval, int linea, const char *file);

// mutex
int _xpthread_mutex_init(pthread_mutex_t *mutex, const pthread_mutexattr_t *attr, int linea, const char *file);
int _xpthread_mutex_destroy(pthread_mutex_t *mutex, int linea, const char *file);
int _xpthread_mutex_lock(pthread_mutex_t *mutex, int linea, const char *file);
int _xpthread_mutex_unlock(pthread_mutex_t *mutex, int linea, const char *file);
#define xpthread_mutex_init(mutex, attr) _xpthread_mutex_init(mutex, attr, __LINE__, __FILE__)
#define xpthread_mutex_destroy(mutex) _xpthread_mutex_destroy(mutex, __LINE__, __FILE__)
#define xpthread_mutex_lock(mutex) _xpthread_mutex_lock(mutex, __LINE__, __FILE__)
#define xpthread_mutex_unlock(mutex) _xpthread_mutex_unlock(mutex, __LINE__, __FILE__)

//semaphores
int _xsem_init(sem_t *sem, int pshared, unsigned int value, int linea, const char *file);
int _xsem_post(sem_t *sem, int linea, const char *file);

#define xsem_init(sem, pshared, value) _xsem_init(sem, pshared, value, __LINE__, __FILE__)
#define xsem_post(sem) _xsem_post(sem, __LINE__, __FILE__)
#define xsem_wait(sem) _xsem_wait(sem, __LINE__, __FILE__)
#define xsem_destroy(sem) _xsem_destroy(sem, __LINE__, __FILE__)

static inline int _xsem_wait(sem_t *sem, int linea, const char *file) {
    const int e = sem_wait(sem);
    if(__builtin_expect(e, 0)) {
        xperror(e, "Error in sem_wait");
        fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),linea,file);
        sleep(Thread_error_wait);  // do not kill immediately other threads
        exit(1);
    }
    return e;
}
int _xsem_destroy(sem_t *sem, int linea, const char *file);
