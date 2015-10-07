# -*- coding: utf-8 -*-

cimport cython

#cdef extern from "stdlib.h":
  #void qsort(void *base, long nmemb, long size,
       #int(*compar)(const void *, const void *)) nogil

cdef inline void long_array_init(long *a,long n,long v) nogil:
  """
  initialize integer array a of length n with integer value v
  """
  cdef long i
  for i in xrange(n):
    a[i] = v
  return

cdef inline void float_array_init(float *a,long n,float v) nogil:
  """
  initialize float array a of length n with float value v
  """
  cdef long i
  for i in xrange(n):
    a[i] = v
  return

cdef inline void double_array_init(double *a,long n,double v) nogil:
  """
  """
  cdef long i
  for i in xrange(n):
    a[i] = v
  return

#cdef inline long compare(const void *aa, const void *bb):
  #"""
  #compare function used with qsort
  #"""

  #cdef long a = (<long *>(aa))[0]
  #cdef long b = (<long *>(bb))[0]

  #if a < b:
    #return -1
  #elif a > b:
    #return 1
  #else:
    #return 0

