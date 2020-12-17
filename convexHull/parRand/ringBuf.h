#ifndef RING_BUF_H
#define RING_BUF_H

template<class T>
struct ringBuffer {
  T* A;
  intT n;

  ringBuffer(intT nn): n(nn) {
    A = newA(T, n);
  }

  ~ringBuffer() {
    free(A);
  }

  T& operator[](intT i) {
    return A[i%n];
  }

  T& at(intT i) {
    return A[i%n];
  }

};

#endif
