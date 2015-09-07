#!/usr/bin/python
# -*- coding: utf-8 -*-

from zonemap3d import Zonemap3d


def main():

  nmax = int(1e7)

  num_points = 1000000
  num_lookup = 1000000
  rads = [0.03, 0.02, 0.01]

  for rad in rads:

    nz = int(1.0/rad)
    print('zones\t\t{:d}'.format(nz))
    print('rad\t\t{:2.4f}'.format(rad))
    print(' ')

    ZM = Zonemap3d(nz)
    res = ZM._perftest(nmax, num_points, num_lookup)
    for expl, perf in res:
      print('{:s}\t\t{:3.10f}'.format(expl, perf))
    del(ZM)

    print('\n\n')

  return

if __name__ == '__main__' :

    main()

