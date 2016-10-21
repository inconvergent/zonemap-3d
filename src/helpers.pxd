# -*- coding: utf-8 -*-

cimport cython

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.nonecheck(False)
cdef inline void long_array_init(long *a,long n,long v) nogil:
  """
  initialize integer array a of length n with integer value v
  """
  cdef long i
  for i in xrange(n):
    a[i] = v
  return

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.nonecheck(False)
cdef inline void float_array_init(float *a,long n,float v) nogil:
  """
  initialize float array a of length n with float value v
  """
  cdef long i
  for i in xrange(n):
    a[i] = v
  return

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.nonecheck(False)
cdef inline void double_array_init(double *a,long n,double v) nogil:
  """
  """
  cdef long i
  for i in xrange(n):
    a[i] = v
  return

