# -*- coding: utf-8 -*-

cimport cython

#cdef extern from "stdlib.h":
  #void qsort(void *base, int nmemb, int size,
       #int(*compar)(const void *, const void *)) nogil

cdef inline void int_array_init(int *a,int n,int v) nogil:
  """
  initialize integer array a of length n with integer value v
  """
  cdef int i
  for i in xrange(n):
    a[i] = v
  return

cdef inline void float_array_init(float *a,int n,float v) nogil:
  """
  initialize float array a of length n with float value v
  """
  cdef int i
  for i in xrange(n):
    a[i] = v
  return

#cdef inline int compare(const void *aa, const void *bb):
  #"""
  #compare function used with qsort
  #"""

  #cdef int a = (<int *>(aa))[0]
  #cdef int b = (<int *>(bb))[0]

  #if a < b:
    #return -1
  #elif a > b:
    #return 1
  #else:
    #return 0

