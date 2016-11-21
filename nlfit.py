#!/bin/python
# -*- coding: utf-8 -*-

# Authors : J-O Simoneau, M 'Ardi'??

# from pylab import * # Might be required by the loading script
import numpy as np
from numpy import array, sqrt, diag, concatenate, mean, size
from scipy.optimize import leastsq





class nlfit:
    """ Description to come
    """
    def __init__(self, xs, ys, p0, fs, pmasks=None, fullo=1, scales=None,
                   xerrs=None, yerrs=None, verbose=True):
        # TODO : Add the possibility to choose the weigth function (to generate 'self.scales').
        # I want to take advantage of the numpy.array functionnality to its maximum here.
        # NOTE : May not be as general as it should in term of ndarray but for the 2D case it should be fine.

        # Does not work if numpy array, use 'asanyarray' instead
        #if not isinstance(xs,list):
        #    xs, ys, fs = [xs], [ys], [fs]

        # Abscisse
        self.xs = np.asanyarray(xs)
        # xerrors not implemented yet (would require odr)
        if xerrs is None:
            self.xerrs = [[1.]*len(x) for x in self.xs]
        else:
            self.xerrs = np.array(xerrs)

        # Ordonnée
        self.ys = np.asanyarray(ys)
        if yerrs is None:
            self.yerrs = [[1.]*len(y) for y in self.ys]
        else:
            self.yerrs = np.array(yerrs)

        if pmasks is None:
            pmasks = [[1]*len(p0)]*len(xs)

        # Make shure we can use numpy's mask functionality.
        pmasks = np.asanyarray(pmasks).astype(bool)

        if any(len(i) != len(xs) \
                for i in [ys, fs, pmasks, self.xerrs, self.yerrs]):
            raise(AssertionError,"List size don't match")

        self._p0 = p0
        self.para = p0
        self.fs = np.asanyarray(fs)
        self.pmasks = pmasks
        if scales is None:
            self.scales = [mean(abs(y)) for y in ys]
        else:
            self.scales = scales
        self.fullo = fullo
        self.verbose = verbose

    def __call__(self, *args):
        if len(args) == 1:
            return self.__custom_call__(0,*args)
        else:
            return self.__custom_call__(*args)

    def __custom_call__(self, fct_num, x):
        return self.fs[fct_num](x,self.para[self.pmasks[fct_num]])

    def __getitem__(self,i):
        return self.para[i]

    def __len__(self):
        return len(self.para)

    def _ps(self,p):
        """ Post select the parameters for each function call """
        return [p[mask] for mask in self.pmasks]
                                                                         
    def _residuals(self, y, f, x, p, yerr, scale):                       
        return (y - f(x,*p))/yerr/scale                                   
                                                                         
    def _residuals_global(self, p):                                      
        errs = [self._residuals(y,f,x,mp,yerr,scale) \
                for y,f,x,mp,yerr,scale in zip(self.ys, self.fs, self.xs, 
                                          self._ps(p), self.yerrs, self.scales)]
        return concatenate(errs) # Leastsq accepts a list of floating points numbers.

    def _generate_mask(self, curves, independant, shared):
        """ Generates a mask for the case of data with an homogeneous number of shared and 
            non shared parameters given the following parameters : 
                - n : number of curves to fit
                - m : number of independant parameters per function (should be the same for each 
                      function)
                - p : number of shared parameters
            The row are in the following logistic (p,m) shared parameters then independant 
            parameters.
        """
        n = curves
        m = independant
        p = shared
        I = np.identity(n) # n is the number of curves to fit
        v = np.zeros(m)+1. # m is the number of independant parameters per function
        # (should be homogeneous).
        t = np.tensordot(I,v, axes=0)
        mask = np.zeros((n,n*m+p)) # p is the number of shared parameters
        mask[:,:p] = 1. # fix the common parameters
        mask[:,p:] = t.reshape((n,n*m))
        mask = mask.astype(bool)
        self.pmasks = mask

    def leastsq(self, **kwargs):
        self.lsq = leastsq(self._residuals_global, 
                             self.para, full_output=self.fullo, **kwargs)
        if self.lsq[1] is None:
            if self.verbose: print '\n --- FIT DID NOT CONVERGE ---\n'
            self.errs = None
            self.err = None
            self.chi2rs = None
            return False
        else:
            self.para = self.lsq[0]
            self.cv = self.lsq[1]
            self.it = self.lsq[2]['nfev']
            self.computevalues()
            self.errs = array([self.sdcv*sqrt(chi2r) for chi2r in self.chi2rs])
            self.err = self.errs[0]

            if self.verbose:
                print self
            return True

    def computevalues(self):
        self.sdcv = sqrt(diag(self.cv))
        # Matrice de corrélation
        self.corrM = self.cv/self.sdcv/self.sdcv[:,None]
        self.chi2s = [sum(self._residuals(y,f,x,mp,yerr,scale)**2) \
                      for y,f,x,mp,yerr,scale in zip(self.ys, self.fs, self.xs,
                                                       self._ps(self.para), 
                                                       self.yerrs, self.scales)]
        # Chi^2 réduit
        self.chi2rs = [chi2/(len(y)-len(self.para)) \
                        for chi2,y in zip(self.chi2s, self.ys)]

    def __str__(self):
        s = '\n--- FIT ON FUNCTION{} {} ---'+\
            '\n\nFit parameters are\n{}\nFit errors are\n{}\n'+\
            '\nFit covariance\n{}'+\
            '\nFit correlation matrix\n{}\nReduced chi2s are {}\n\n'
        fmt = ['S' if len(self.xs)>1 else '',', '.join([f.__name__ \
                                                          for f in self.fs]), 
                 self.para, self.errs, self.cv, self.corrM, self.chi2rs]
        tmp = fmt[1].rfind(', ')
        if not tmp == -1:
            fmt[1] = fmt[1][:tmp] + ' and ' + fmt[1][tmp+2:]
        return s.format(*fmt)
