__author__ = 'Mathias Zechmeister'
__version__ = '2018-10-01'

import numpy as np
import os
from pause import pause

from ctypes import c_int, c_double
import ctypes

import gls


ptr = np.ctypeslib.ndpointer

#_cgls = np.ctypeslib.load_library(os.path.join(os.path.dirname(__file__), 'gls'), '.')
_cgls = np.ctypeslib.load_library('gls_c', '/home/raid0/zechmeister/python/mzechmeister/GLS/c/')

# void gls(double *t, double *y, double *e_y, long n, double fbeg, long nk, double fstep, double *f, double *p) {
_cgls.gls.argtypes = [ ptr(dtype=np.float),   # t
                       ptr(dtype=np.float),   # y
                       ptr(dtype=np.float),   # e_y
                       c_int,                 # n
                       c_double,              # fbeg
                       c_int,                 # nk
                       c_double,              # fstep
                       ptr(dtype=np.float),   # f
                       ptr(dtype=np.float) ]   # p


'''
# Example
time = np.random.uniform(54000., 56000., 1000)
time = np.random.uniform(0., 560., 73)
flux = 0.15 * np.sin(2. * np.pi * time / 10.)

#    Add some noise
error = 0.5 * np.ones(time.size)
flux += np.random.normal(0, error)

glspy = gls.Gls((time,flux,error))
f = 1*glspy.freq
p = 1*f
gg = time*1
_cgls.gls(time,flux,error*0+1, time.size, f[0], f.size, glspy.fstep, f, p)

from gplot import *
gplot(f,p, ',', glspy.freq, glspy.p)

pause()
'''

class cGls(gls.Gls):
   '''
   GLS periodogram wrapping the c implementation.

   Example
   -------
   >>> time = np.random.uniform(54000., 56000., 1000)
   >>> flux = 0.15 * np.sin(2. * np.pi * time / 10.)

   # Add some noise
   >>> error = 0.5 * np.ones(time.size)
   >>> flux += np.random.normal(0, error)

   >>> cg = cGls((time,flux,error*0+1))
   >>> gplot(cg.freq, cg.p)

   # Comparsion with Gls
   >>> glspy = gls.Gls((time,flux,error))
   >>> gplot+(glspy.freq, glspy.p)

   '''
   def _peakPeriodogram(self):
      pass
   def _calcPeriodogram(self):
      self.p = self.freq * 1.
      print "c version"
      _cgls.gls(self.th, self.y, self.yerr, self.N, self.freq[0], self.nf, self.fstep, self.freq, self.p)
      self._YY = 1.

if __name__ == "__main__":
    import doctest
    from gplot import *
    exec(doctest.script_from_examples(cGls.__doc__))
    from pause import *; pause()

# hh = cGls((time,flux,error))




