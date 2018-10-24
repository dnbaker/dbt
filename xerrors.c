#define _GNU_SOURCE 1
#include "xerrors.h"

// xerrors.c
// system calls with output checking



#define Buflen 128

// write error message associated to code en similarly to perror
void xperror(int en, const char *msg) {
  char buf[Buflen];
#if (_POSIX_C_SOURCE >= 200112L || _XOPEN_SOURCE >= 600) && ! _GNU_SOURCE
  int ec = strerror_r(en, buf, Buflen);
  fprintf(stderr, "[Error:%d] %s: %s\n", en, msg, buf);
  /* using posix strerror_r */
#else
  char *errmsg = strerror_r(en, buf, Buflen);
  if(msg!=NULL)
    fprintf(stderr,"%s: %s\n",msg, errmsg);
  else
    fprintf(stderr,"%s\n",errmsg);
#endif
}


// ----- threads

int _xpthread_create(pthread_t *thread, const pthread_attr_t *attr,
                     void *(*start_routine) (void *), void *arg, int linea, const char *file) {
  int e = pthread_create(thread, attr, start_routine, arg);
  if (e!=0) {
    xperror(e, "Error in pthread_create");
    fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),linea,file);
    sleep(Thread_error_wait);  // do not kill immediately other threads
    exit(1);
  }
  return e;
}


int _xpthread_join(pthread_t thread, void **retval, int linea, const char *file) {
  int e = pthread_join(thread, retval);
  if (e!=0) {
    xperror(e, "Error in pthread_join");
    fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),linea,file);
    sleep(Thread_error_wait);  // do not kill immediately other threads
    exit(1);
  }
  return e;
}

// mutex
int _xpthread_mutex_init(pthread_mutex_t *mutex, const pthread_mutexattr_t *attr, int linea, const char *file) {
  int e = pthread_mutex_init(mutex, attr);
  if (e!=0) {
    xperror(e, "Error in pthread_mutex_init");
    fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),linea,file);
    sleep(Thread_error_wait);  // do not kill immediately other threads
    exit(1);
  }
  return e;
}

int _xpthread_mutex_destroy(pthread_mutex_t *mutex, int linea, const char *file) {
  int e = pthread_mutex_destroy(mutex);
  if (e!=0) {
    xperror(e, "Error in pthread_mutex_destroy");
    fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),linea,file);
    sleep(Thread_error_wait);  // do not kill immediately other threads
    exit(1);
  }
  return e;
}

int _xpthread_mutex_lock(pthread_mutex_t *mutex, int linea, const char *file) {
  int e = pthread_mutex_lock(mutex);
  if (e!=0) {
    xperror(e, "Error in pthread_mutex_lock");
    fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),linea,file);
    sleep(Thread_error_wait);  // do not kill immediately other threads
    exit(1);
  }
  return e;
}

int _xpthread_mutex_unlock(pthread_mutex_t *mutex, int linea, const char *file) {
  int e = pthread_mutex_unlock(mutex);
  if (e!=0) {
    xperror(e, "Error in pthread_mutex_unlock");
    fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),linea,file);
    sleep(Thread_error_wait);  // do not kill immediately other threads
    exit(1);
  }
  return e;
}


// semaphores
int _xsem_init(sem_t *sem, int pshared, unsigned int value, int linea, const char *file) {
  int e = sem_init(sem,pshared,value);
  if (e!=0) {
    xperror(e, "Error in sem_init");
    fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),linea,file);
    sleep(Thread_error_wait);  // do not kill immediately other threads
    exit(1);
  }
  return e;
}

int _xsem_post(sem_t *sem, int linea, const char *file) {
  int e = sem_post(sem);
  if (e!=0) {
    xperror(e, "Error in sem_post");
    fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),linea,file);
    sleep(Thread_error_wait);  // do not kill immediately other threads
    exit(1);
  }
  return e;
}



int _xsem_destroy(sem_t *sem, int linea, const char *file) {
  int e = sem_destroy(sem);
  if (e!=0) {
    xperror(e, "Error in sem_destroy");
    fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),linea,file);
    sleep(Thread_error_wait);  // do not kill immediately other threads
    exit(1);
  }
  return e;
}
