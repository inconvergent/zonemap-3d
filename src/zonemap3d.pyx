# -*- coding: utf-8 -*-
# cython: profile=True

from __future__ import division

cimport cython
from libc.stdlib cimport malloc, free, realloc
#from cython.parallel import parallel, prange
from libc.math cimport sqrt


from helpers cimport int_array_init
from helpers cimport double_array_init

import numpy as np
cimport numpy as np

from time import time


cdef int SIZE = 64


cdef class Zonemap3d:

  def __init__(self, int nz):
    """
    """

    self.vnum = 0

    self.vsize = SIZE

    self.nz = nz

    self.total_zones = nz*nz*nz

    self.greatest_zone_size = SIZE

    self.__init_zones()

    return

  def __cinit__(self, int nz, *arg, **args):

    cdef int total_zones = nz*nz*nz

    self.VZ = <int *>malloc(SIZE*sizeof(int))

    self.ZONES = <sZ **>malloc(total_zones*sizeof(sZ*))

    return

  def __dealloc__(self):

    cdef int i

    for i in xrange(self.total_zones):

      free(self.ZONES[i].ZV)
      free(self.ZONES[i])

    free(self.ZONES)
    free(self.VZ)

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cdef void __assign_xyz_arrays(self, double *X, double *Y, double *Z) :#nogil:

    self.X = X
    self.Y = Y
    self.Z = Z
    return

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cdef void __init_zones(self) :#nogil:
    # somehow this did not work when executed inside cinit

    cdef int i
    cdef sZ *zone

    for i in xrange(self.total_zones):

      zone = <sZ *>malloc(sizeof(sZ))

      zone.i = i
      zone.size = SIZE
      zone.count = 0
      zone.ZV = <int *>malloc(SIZE*sizeof(int))

      self.ZONES[i] = zone

    return

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cdef int __add_vertex(self, int v1) :#nogil:
    """
    """

    cdef int vnum = self.vnum

    cdef double x = self.X[v1]
    cdef double y = self.Y[v1]
    cdef double z = self.Z[v1]

    cdef int z1 = self.__get_z(x,y,z)

    self.__add_v_to_zone(z1, vnum)
    self.VZ[vnum] = z1

    cdef int* new_vz

    if self.vnum>=self.vsize-1:

      new_vz = <int *>realloc(self.VZ, self.vsize*2*sizeof(int))

      if new_vz is not NULL:
        self.VZ = new_vz;
        self.vsize = self.vsize*2
      else:
        ## this is really, really, bad
        return -1

    self.vnum += 1
    return vnum

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cdef int __del_vertex(self, int v1) :#nogil:
    """
    """

    cdef int z1 = self.VZ[v1]

    self.__remove_v_from_zone(z1, v1)
    self.VZ[v1] = -1

    return 1

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cdef int __add_v_to_zone(self, int z1, int v1) :#nogil:

    cdef sZ *zone = self.ZONES[z1]

    zone.ZV[zone.count] = v1
    zone.count += 1

    if zone.count>=zone.size-1:
      return self.__extend_zv_of_zone(zone)

    return 1

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cdef int __extend_zv_of_zone(self, sZ *zone) :#nogil:

    cdef int new_size = zone.size*2
    cdef int* new_zv = <int *>realloc(zone.ZV, new_size*sizeof(int))

    if new_zv is not NULL:
      zone.ZV = new_zv;
      zone.size = new_size
      if new_size>self.greatest_zone_size:
        self.greatest_zone_size = new_size
    else:
      ## this is really, really, bad
      return -1

    return 1

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cdef int __remove_v_from_zone(self, int z1, int v1) :#nogil:

    cdef sZ *zone = self.ZONES[z1]
    cdef int i

    for i in xrange(zone.count):

      if zone.ZV[i] == v1:
        zone.ZV[i] = zone.ZV[zone.count-1]
        zone.count -= 1
        return 1

    return -1

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cdef int __get_z(self, double x, double y, double z) :#nogil:
    """
    """

    cdef int nz = self.nz
    cdef int k = int(z*nz)
    cdef int j = int(y*nz)
    cdef int i = int(x*nz)

    return nz*nz*k + nz*j + i

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cdef int __update_v(self, int v1) :#nogil:

    cdef double x = self.X[v1]
    cdef double y = self.Y[v1]
    cdef double z = self.Z[v1]
    cdef int new_zone = self.__get_z(x, y, z)

    cdef int old_zone = self.VZ[v1]

    if old_zone<0:
      return -1

    if new_zone != old_zone:

      self.__remove_v_from_zone(old_zone, v1)
      self.__add_v_to_zone(new_zone, v1)
      self.VZ[v1] = new_zone

      return 1

    return -1


  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  #@cython.cdivision(True)
  cdef int __sphere_is_free(self, double x, double y, double z, double rad) :#nogil:
    """
    tests if there is another vertex within rad of x,y. rad must be less than
    the width of each zone.
    """

    #TODO: optimize this

    cdef int i
    cdef int j
    cdef sZ *zone
    cdef int zi = self.__get_z(x,y,z)

    cdef int nz = self.nz

    cdef double dx
    cdef double dy
    cdef double dz
    cdef double rad2 = rad*rad

    cdef int zz = int(z*nz)
    cdef int zy = int(y*nz)
    cdef int zx = int(x*nz)

    cdef int a
    cdef int b
    cdef int c

    for a in xrange(max(zx-1,0),min(zx+2,nz)):
      for b in xrange(max(zy-1,0),min(zy+2,nz)):
        for c in xrange(max(zz-1,0),min(zz+2,nz)):

          zone = self.ZONES[c*nz*nz + b*nz + a]

          for j in xrange(zone.count):

            dx = x-self.X[zone.ZV[j]]
            dy = y-self.Y[zone.ZV[j]]
            dz = z-self.Z[zone.ZV[j]]

            if dx*dx+dy*dy+dz*dz<rad2:
              return -1

    return 1

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  #@cython.cdivision(True)
  cdef int __sphere_vertices(
    self,
    double x,
    double y,
    double z,
    double rad,
    int *vertices
  ) :#nogil:

    cdef int i
    cdef int j
    cdef sZ *zone
    cdef int zi = self.__get_z(x,y,z)

    cdef int nz = self.nz

    cdef int num = 0

    cdef double dx
    cdef double dy
    cdef double dz
    cdef double rad2 = rad*rad

    # cdef int zz = int(zi/nz/nz)
    # cdef int zy = int((zi-nz*nz*zz)/nz)
    # cdef int zx = int(zi-nz*nz*zz-nz*zy)

    cdef int zz = int(z*nz)
    cdef int zy = int(y*nz)
    cdef int zx = int(x*nz)

    cdef int a
    cdef int b
    cdef int c

    for a in xrange(max(zx-1,0),min(zx+2,nz)):
      for b in xrange(max(zy-1,0),min(zy+2,nz)):
        for c in xrange(max(zz-1,0),min(zz+2,nz)):

          zone = self.ZONES[c*nz*nz + b*nz + a]

          for j in xrange(zone.count):

            dx = x-self.X[zone.ZV[j]]
            dy = y-self.Y[zone.ZV[j]]
            dz = z-self.Z[zone.ZV[j]]

            if dx*dx+dy*dy+dz*dz<rad2:
              #print(a,b,c,j,zone.ZV[j])
              vertices[num] = zone.ZV[j]
              num += 1

    #print('next')

    return num

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cpdef list _perftest(self, int nmax, int num_points, int num_lookup):

    cdef np.ndarray[double, mode="c",ndim=2] a
    cdef int i
    cdef double t1
    cdef double t2
    cdef list res = []


    cdef double *X = <double *>malloc(nmax*sizeof(double))
    cdef double *Y = <double *>malloc(nmax*sizeof(double))
    cdef double *Z = <double *>malloc(nmax*sizeof(double))
    self.__assign_xyz_arrays(X,Y,Z)


    a = 0.5 + 0.2*(1.0-2.0*np.random.random((num_points,3)))
    t1 = time()
    for i in xrange(num_points):
      X[i] = a[i,0]
      Y[i] = a[i,1]
      Z[i] = a[i,2]
      self.__add_vertex(i)
    t2 = time()
    res.append(('add',t2-t1))


    a = np.random.random((num_lookup,3))
    t1 = time()
    for i in xrange(num_lookup):
      self.__sphere_is_free(a[i,0], a[i,1], a[i,2], 0.03)
    t2 = time()
    res.append(('free',t2-t1))


    a = np.random.random((num_lookup,3))
    t1 = time()
    cdef int asize = self.__get_greatest_zone_size()*27
    cdef int *vertices = <int *>malloc(asize*sizeof(int))
    for i in xrange(num_lookup):
      self.__sphere_vertices(
        a[i,0],
        a[i,1],
        a[i,2],
        0.03,
        vertices
      )
    t2 = time()
    res.append(('sphere',t2-t1))


    t1 = time()
    for i in xrange(num_points):
      self.__del_vertex(i)
    t2 = time()
    res.append(('del',t2-t1))

    free(X)
    free(Y)
    free(Z)
    free(vertices)

    return res

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cpdef int add_vertex(self, int v1):

    return self.__add_vertex(v1)

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cpdef int del_vertex(self, int v1):

    return self.__del_vertex(v1)

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cdef int __get_greatest_zone_size(self) :#nogil:

    return self.greatest_zone_size

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cpdef int update_v(self, int v1):

    return self.__update_v(v1)

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cpdef int sphere_is_free(self, double x, double y, double z, double rad):

    return self.__sphere_is_free(x, y, z, rad)

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cpdef int get_greatest_zone_size(self):

    return self.greatest_zone_size

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cpdef int get_vnum(self):

    return self.vnum

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cpdef list get_zone_info_dicts(self):

    cdef list res = []
    cdef dict d
    cdef int i

    for i in xrange(self.total_zones):

      d = {
        'i': self.ZONES[i].i,
        'size': self.ZONES[i].size,
        'count': self.ZONES[i].count
      }

      res.append(d)

    return res

