# -*- coding: utf-8 -*-
# cython: profile=True

from __future__ import division

cimport cython
from libc.stdlib cimport malloc, free, realloc
from libc.math cimport sqrt

from helpers cimport int_array_init
from helpers cimport float_array_init

cdef int SIZE = 1024


cdef class Zonemap3d:

  def __init__(self, int nz):
    """
    """

    self.vnum = 0

    self.vsize = SIZE

    self.nz = nz

    self.total_zones = (2+nz)*(2+nz)*(2+nz)

    self.greatest_zone_size = SIZE

    self.__init_zones()

    return

  def __cinit__(self, int nz, *arg, **args):

    cdef int total_zones = (2+nz)*(2+nz)*(2+nz)

    self.VZ = <int *>malloc(SIZE*sizeof(int))

    self.ZONES = <sZ **>malloc(total_zones*sizeof(sZ*))

    return

  def __dealloc__(self):

    free(self.VZ)

    # TODO: is this a memory leak?
    free(self.ZONES)

    return

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cdef void __assign_xyz_arrays(self, float *X, float *Y, float *Z) nogil:

    self.X = X
    self.Y = Y
    self.Z = Z

    return

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cdef void __init_zones(self) nogil:
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

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cdef int __add_vertex(self, int v1) nogil:
    """
    """

    cdef int vnum = self.vnum

    cdef float x = self.X[v1]
    cdef float y = self.Y[v1]
    cdef float z = self.Z[v1]

    cdef int z1 = self.__get_z(x,y,z)

    self.__add_v_to_zone(z1, vnum)
    self.VZ[vnum] = z1

    cdef int* new_vz

    if self.vnum>=self.vsize-1:

      new_vz = <int *>realloc(self.VZ, self.vsize*2*sizeof(int))

      if new_vz:
        self.VZ = new_vz;
        self.vsize = self.vsize*2
      else:
        ## this is really, really, bad
        return -1

    self.vnum += 1
    return vnum

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cdef int __del_vertex(self, int v1) nogil:
    """
    """

    cdef int z1 = self.VZ[v1]

    self.__remove_v_from_zone(z1, v1)
    self.VZ[v1] = -1

    return 1

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cdef int __add_v_to_zone(self, int z1, int v1) nogil:

    cdef sZ *zone = self.ZONES[z1]

    zone.ZV[zone.count] = v1
    zone.count += 1

    if zone.count>=zone.size-1:
      return self.__extend_zv_of_zone(zone)

    return 1

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cdef int __extend_zv_of_zone(self, sZ *zone) nogil:

    cdef int new_size = zone.size*2
    cdef int* new_zv = <int *>realloc(zone.ZV, new_size*sizeof(int))

    if new_zv:
      zone.ZV = new_zv;
      zone.size = new_size
      if new_size>self.greatest_zone_size:
        self.greatest_zone_size = new_size
    else:
      ## this is really, really, bad
      return -1

    return 1

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cdef int __remove_v_from_zone(self, int z1, int v1) nogil:

    cdef sZ *zone = self.ZONES[z1]
    cdef int i

    for i in xrange(zone.count):

      if zone.ZV[i] == v1:
        zone.ZV[i] = zone.ZV[zone.count-1]
        zone.count -= 1
        return 1

    return -1

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cdef int __get_z(self, float x, float y, float z) nogil:
    """
    """

    cdef int nz = self.nz
    cdef int nz2 = self.nz+2

    cdef int i = 1 + <int>(x*nz)
    cdef int j = 1 + <int>(y*nz)
    cdef int k = 1 + <int>(z*nz)

    return nz2*nz2*k + nz2*i + j

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cdef int __update_v(self, int v1) nogil:

    cdef float x = self.X[v1]
    cdef float y = self.Y[v1]
    cdef float z = self.Z[v1]
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


  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  @cython.cdivision(True)
  cdef int __sphere_is_free(self, float x, float y, float z, float rad) nogil:
    """
    tests if there is another vertex within rad of x,y. rad must be less than
    the width of each zone.
    """

    cdef int i
    cdef int j
    cdef sZ *zone
    cdef int zi = self.__get_z(x,y,z)

    cdef int nz = self.nz
    cdef int nz22 = (nz+2)*(nz+2)

    cdef float dx
    cdef float dy
    cdef float dz
    cdef float rad2 = rad*rad

    # im sorry ...
    cdef int *neighbors = [
      zi+nz22, zi-1+nz22, zi+1+nz22,
      zi-nz-2+nz22, zi+nz+2+nz22, zi-nz-1+nz22,
      zi-nz-3+nz22, zi+nz+1+nz22, zi+nz+3+nz22,

      zi, zi-1, zi+1,
      zi-nz-2, zi+nz+2, zi-nz-1,
      zi-nz-3, zi+nz+1, zi+nz+3,

      zi-nz22, zi-1-nz22, zi+1-nz22,
      zi-nz-2-nz22, zi+nz+2-nz22, zi-nz-1-nz22,
      zi-nz-3-nz22, zi+nz+1-nz22, zi+nz+3-nz22,
    ]

    for i in xrange(27):

      zone = self.ZONES[neighbors[i]]

      for j in xrange(zone.count):

        dx = x-self.X[zone.ZV[j]]
        dy = y-self.Y[zone.ZV[j]]
        dz = z-self.Z[zone.ZV[j]]

        if dx*dx+dy*dy+dz*dz<rad2:
          return -1

    return 1

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  @cython.cdivision(True)
  cdef int __sphere_vertices(self, float x, float y, float z, float rad, int *vertices) nogil:
    """
    """

    cdef int i
    cdef int j
    cdef sZ *zone
    cdef int zi = self.__get_z(x,y,z)

    cdef int nz = self.nz
    cdef int nz22 = (nz+2)*(nz+2)

    cdef int num = 0

    cdef float dx
    cdef float dy
    cdef float dz
    cdef float rad2 = rad*rad

    cdef int *neighbors = [
      zi+nz22, zi-1+nz22, zi+1+nz22,
      zi-nz-2+nz22, zi+nz+2+nz22, zi-nz-1+nz22,
      zi-nz-3+nz22, zi+nz+1+nz22, zi+nz+3+nz22,

      zi, zi-1, zi+1,
      zi-nz-2, zi+nz+2, zi-nz-1,
      zi-nz-3, zi+nz+1, zi+nz+3,

      zi-nz22, zi-1-nz22, zi+1-nz22,
      zi-nz-2-nz22, zi+nz+2-nz22, zi-nz-1-nz22,
      zi-nz-3-nz22, zi+nz+1-nz22, zi+nz+3-nz22,
    ]

    for i in xrange(27):

      zone = self.ZONES[neighbors[i]]

      for j in xrange(zone.count):

        dx = x-self.X[zone.ZV[j]]
        dy = y-self.Y[zone.ZV[j]]
        dz = z-self.Z[zone.ZV[j]]

        if dx*dx+dy*dy+dz*dz<rad2:

          vertices[num] = zone.ZV[j]
          num += 1

    return num

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cpdef int add_vertex(self, int v1):

    return self.__add_vertex(v1)

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cpdef int del_vertex(self, int v1):

    return self.__del_vertex(v1)

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cdef int __get_greatest_zone_size(self) nogil:

    return self.greatest_zone_size

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cpdef int update_v(self, int v1):

    return self.__update_v(v1)

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cpdef int sphere_is_free(self, float x, float y, float z, float rad):

    return self.__sphere_is_free(x, y, z, rad)

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cpdef int get_greatest_zone_size(self):

    return self.greatest_zone_size

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cpdef int get_vnum(self):

    return self.vnum

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
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

