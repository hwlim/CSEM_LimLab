#ifndef PTHREADWRAPPER_H_
#define PTHREADWRAPPER_H_

#include<vector>
#include<pthread.h>

struct PThreadWrapper {
  int num_threads_used;

  std::vector<pthread_t> threads;
  pthread_attr_t attr;
  int rc;

  PThreadWrapper() { 
    num_threads_used = 0; 
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  }

  ~PThreadWrapper() {
    pthread_attr_destroy(&attr);
  }

};

#endif
