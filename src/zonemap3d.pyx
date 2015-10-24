# -*- coding: utf-8 -*-
# cython: profile=True

from __future__ import division

cimport cython
from libc.stdlib cimport malloc, free, realloc
#from cython.parallel import parallel, prange
from libc.math cimport sqrt


from helpers cimport long_array_init
from helpers cimport double_array_init

import numpy as np
cimport numpy as np

from time import time


cdef long SIZE = 64


cdef class Zonemap3d:

  def __init__(self, long nz):
    """
    """

    self.vnum = 0

    self.vsize = SIZE

    self.nz = nz

    self.total_zones = nz*nz*nz

    self.greatest_zone_size = SIZE

    self.__init_zones()

    return

  def __cinit__(self, long nz, *arg, **args):

    cdef long total_zones = nz*nz*nz

    self.VZ = <long *>malloc(SIZE*sizeof(long))

    self.ZONES = <sZ **>malloc(total_zones*sizeof(sZ*))

    return

  def __dealloc__(self):

    cdef long i

    for i in xrange(self.total_zones):

      free(self.ZONES[i].ZV)
      free(self.ZONES[i])

    free(self.ZONES)
    free(self.VZ)

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cdef void __assign_xyz_arrays(self, double *X, double *Y, double *Z) nogil:

    self.X = X
    self.Y = Y
    self.Z = Z
    return

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cdef void __init_zones(self) nogil:
    # somehow this did not work when executed inside cinit

    cdef long i
    cdef sZ *zone

    for i in xrange(self.total_zones):

      zone = <sZ *>malloc(sizeof(sZ))

      zone.i = i
      zone.size = SIZE
      zone.count = 0
      zone.ZV = <long *>malloc(SIZE*sizeof(long))

      self.ZONES[i] = zone

    return

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cdef long __add_vertex(self, long v1) nogil:
    """
    """

    cdef long vnum = self.vnum

    cdef double x = self.X[v1]
    cdef double y = self.Y[v1]
    cdef double z = self.Z[v1]

    cdef long z1 = self.__get_z(x,y,z)

    self.__add_v_to_zone(z1, vnum)
    self.VZ[vnum] = z1

    cdef long* new_vz

    if self.vnum>=self.vsize-1:

      new_vz = <long *>realloc(self.VZ, self.vsize*2*sizeof(long))

      if new_vz is not NULL:
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
  cdef long __del_vertex(self, long v1) nogil:
    """
    """

    cdef long z1 = self.VZ[v1]

    self.__remove_v_from_zone(z1, v1)
    self.VZ[v1] = -1

    return 1

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cdef long __add_v_to_zone(self, long z1, long v1) nogil:

    cdef sZ *zone = self.ZONES[z1]

    zone.ZV[zone.count] = v1
    zone.count += 1

    if zone.count>=zone.size-1:
      return self.__extend_zv_of_zone(zone)

    return 1

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cdef long __extend_zv_of_zone(self, sZ *zone) nogil:

    cdef long new_size = zone.size*2
    cdef long* new_zv = <long *>realloc(zone.ZV, new_size*sizeof(long))

    if new_zv is not NULL:
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
  cdef long __remove_v_from_zone(self, long z1, long v1) nogil:

    cdef sZ *zone = self.ZONES[z1]
    cdef long i

    for i in xrange(zone.count):

      if zone.ZV[i] == v1:
        zone.ZV[i] = zone.ZV[zone.count-1]
        zone.count -= 1
        return 1

    return -1

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cdef long __get_z(self, double x, double y, double z) nogil:
    """
    """

    cdef long nz = self.nz
    cdef long k = <long>(z*nz)
    cdef long j = <long>(y*nz)
    cdef long i = <long>(x*nz)

    return nz*nz*k + nz*j + i

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cdef long __update_v(self, long v1) nogil:

    cdef long new_zone = self.__get_z(
      self.X[v1],
      self.Y[v1],
      self.Z[v1]
    )

    cdef long old_zone = self.VZ[v1]

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
  cdef long __sphere_is_free(self, double x, double y, double z, double rad) nogil:
    """
    tests if there is another vertex within* rad of x,y. rad must be less than
    the width of each zone.
    """

    #TODO: optimize this

    cdef long i
    cdef long j
    cdef sZ *zone
    cdef long zi = self.__get_z(x,y,z)

    cdef long nz = self.nz

    cdef double dx
    cdef double dy
    cdef double dz
    cdef double rad2 = rad*rad

    cdef long zz = <long>(z*nz)
    cdef long zy = <long>(y*nz)
    cdef long zx = <long>(x*nz)

    cdef long zvj

    cdef long a
    cdef long b
    cdef long c

    for a in xrange(max(zx-1,0),min(zx+2,nz)):
      for b in xrange(max(zy-1,0),min(zy+2,nz)):
        for c in xrange(max(zz-1,0),min(zz+2,nz)):

          zone = self.ZONES[c*nz*nz + b*nz + a]

          for j in xrange(zone.count):

            zvj = zone.ZV[j]

            dx = x-self.X[zvj]
            dy = y-self.Y[zvj]
            dz = z-self.Z[zvj]

            if dx*dx+dy*dy+dz*dz<rad2:
              return -1

    return 1

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  @cython.cdivision(True)
  cdef long __sphere_is_free_ignore(self, double x, double y, double z, long v, double rad) nogil:
    """
    tests if there is another vertex withlong* rad of x,y. rad must be less than
    the width of each zone.
    """

    #TODO: optimize this

    cdef long i
    cdef long j
    cdef sZ *zone
    cdef long zi = self.__get_z(x,y,z)

    cdef long nz = self.nz

    cdef double dx
    cdef double dy
    cdef double dz
    cdef double rad2 = rad*rad

    cdef long zz = <long>(z*nz)
    cdef long zy = <long>(y*nz)
    cdef long zx = <long>(x*nz)

    cdef long zvj

    cdef long a
    cdef long b
    cdef long c

    for a in xrange(max(zx-1,0),min(zx+2,nz)):
      for b in xrange(max(zy-1,0),min(zy+2,nz)):
        for c in xrange(max(zz-1,0),min(zz+2,nz)):

          zone = self.ZONES[c*nz*nz + b*nz + a]

          for j in xrange(zone.count):

            zvj = zone.ZV[j]

            if zvj == v:
              continue

            dx = x-self.X[zvj]
            dy = y-self.Y[zvj]
            dz = z-self.Z[zvj]

            if dx*dx+dy*dy+dz*dz<rad2:
              return -1

    return 1

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  @cython.cdivision(True)
  cdef long __sphere_vertices(
    self,
    double x,
    double y,
    double z,
    double rad,
    long *vertices
  ) nogil:

    cdef long i
    cdef long j
    cdef sZ *zone
    cdef long zi = self.__get_z(x,y,z)

    cdef long nz = self.nz

    cdef long num = 0

    cdef double dx
    cdef double dy
    cdef double dz
    cdef double rad2 = rad*rad

    # cdef long zz = long(zi/nz/nz)
    # cdef long zy = long((zi-nz*nz*zz)/nz)
    # cdef long zx = long(zi-nz*nz*zz-nz*zy)

    cdef long zz = <long>(z*nz)
    cdef long zy = <long>(y*nz)
    cdef long zx = <long>(x*nz)

    cdef long a
    cdef long b
    cdef long c

    for a in xrange(max(zx-1,0),min(zx+2,nz)):
      for b in xrange(max(zy-1,0),min(zy+2,nz)):
        for c in xrange(max(zz-1,0),min(zz+2,nz)):

          zone = self.ZONES[c*nz*nz + b*nz + a]

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
  @cython.cdivision(True)
  cdef long __sphere_vertices_dst(
    self,
    double x,
    double y,
    double z,
    double rad,
    long *vertices,
    double *dst) nogil:

    cdef long i
    cdef long j
    cdef sZ *zone
    cdef long zi = self.__get_z(x,y,z)

    cdef long nz = self.nz

    cdef long num = 0
    cdef double rad2 = rad*rad
    cdef double nrm2

    cdef double dx
    cdef double dy
    cdef double dz

    cdef long zz = <long>(z*nz)
    cdef long zy = <long>(y*nz)
    cdef long zx = <long>(x*nz)

    cdef long a
    cdef long b
    cdef long c

    for a in xrange(max(zx-1,0),min(zx+2,nz)):
      for b in xrange(max(zy-1,0),min(zy+2,nz)):
        for c in xrange(max(zz-1,0),min(zz+2,nz)):

          zone = self.ZONES[c*nz*nz + b*nz + a]

          for j in xrange(zone.count):

            dx = x-self.X[zone.ZV[j]]
            dy = y-self.Y[zone.ZV[j]]
            dz = z-self.Z[zone.ZV[j]]

            nrm2 = dx*dx+dy*dy+dz*dz

            if nrm2<rad2:

              vertices[num] = zone.ZV[j]
              dst[4*num] = dx
              dst[4*num+1] = dy
              dst[4*num+2] = dz
              dst[4*num+3] = sqrt(nrm2)
              num += 1

    return num

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cpdef list _perftest(self, long nmax, long num_points, long num_lookup):

    cdef np.ndarray[double, mode="c",ndim=2] a
    cdef long i
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
    cdef long asize = self.__get_max_sphere_count()
    cdef long *vertices = <long *>malloc(asize*sizeof(long))
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

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cpdef long add_vertex(self, long v1):

    return self.__add_vertex(v1)

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cpdef long del_vertex(self, long v1):

    return self.__del_vertex(v1)

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cdef long __get_max_sphere_count(self) nogil:

    return self.greatest_zone_size*27

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cdef long __get_greatest_zone_size(self) nogil:

    return self.greatest_zone_size

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cpdef long update_v(self, long v1):

    return self.__update_v(v1)

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cpdef long sphere_is_free(self, double x, double y, double z, double rad):

    return self.__sphere_is_free(x, y, z, rad)

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cpdef long get_max_sphere_count(self):

    return self.__get_max_sphere_count()

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cpdef long get_greatest_zone_size(self):

    return self.__get_greatest_zone_size()

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cpdef long get_vnum(self):

    return self.vnum

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cpdef list get_zone_info_dicts(self):

    cdef list res = []
    cdef dict d
    cdef long i

    for i in xrange(self.total_zones):

      d = {
        'i': self.ZONES[i].i,
        'size': self.ZONES[i].size,
        'count': self.ZONES[i].count
      }

      res.append(d)

    return res

