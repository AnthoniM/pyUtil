#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import time
import glob
import os
import re
import platform
import itertools
import inspect
import scipy.constants as c
import math
#import copy
from copy import deepcopy
import scipy.special as sp
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.integrate import quad
from functools import wraps
from pyHegel.fitting import fitcurve
import scipy.fftpack
import matplotlib.pyplot as plt

try:
    import util #Not compatible with python3
    def make_cluster_pdf(date):
        """
        Make one file named cluster_'date' from all the files starting with the same 'date'.
        """
        directory = main_dir+'Results/'
        fnames = glob.glob(directory+date+'*.pdf')
        fnames.sort()
        util.merge_pdf(fnames, directory+'cluster_'+date+'.pdf')
except:
    print(('make_cluster_pdf could not be defined because of incompatibility of pyHegel with python3'))

def apply_along_axis_wrapper(func, deepcopy=None):
    """ Appy func to a set of multidimensional arguments."""
    @wraps(func)
    def apply_along_axis_general(*args, **kwargs):
        axis = kwargs.get('axis',func.func_dict.get('axis',-1))
        shape=np.shape(args[0])
        # Bypass x and y are 1D arrays.
        if len(shape)>1:
            if deepcopy:
                cargs = [deepcopy(i) for i in args]
            else:
                cargs = list(args)
            # Swap the array so that the specified axis is last.
            #---------------------------------------------------
            print shape
            swapcount = 0
            if axis!=-1 and len(shape)>1:
                while True:
                    for i in range(len(cargs)):
                        cargs[i] = np.swapaxes(cargs[i], axis, axis+1)
                    swapcount += 1
                    axis += 1
                    if axis==(len(shape)-1) or axis==-1:
                        break
            # Flatten all dimensions but the last to process more easely.
            #---------------------------------------------------
            newshape = cargs[0].shape
            print newshape
            for i in range(len(cargs)):
                cargs[i] = cargs[i].reshape(-1, newshape[-1])
            # Function call.
            #---------------------------------------------------
            cargs = np.rollaxis(np.array(cargs), 1)
            newargs = [func(*arg, **kwargs) for arg in cargs]
            # Verify if the function returned multiple arguments.
            if len(np.shape(newargs))==3:
                for i in range(len(newargs)):
                    newargs[i] = newargs[i].reshape(*(list(newshape[:-1])+[-1]))
                for i in range(swapcount):
                    for j in range(len(newargs)):
                        newargs[j] = np.swapaxes(newargs[j], -1-i, -2-i)
            else:
                newargs = np.array(newargs).reshape(*(list(newshape[:-1])+[-1]))
                for i in range(swapcount):
                    newargs = np.swapaxes(newargs, -1-i, -2-i)
        else:
            newargs = func(*args, **kwargs)
        return np.squeeze(newargs) # Remove single sized dimensions.
    return apply_along_axis_general

def apply_over_array_wrapper(func, deepcopy=None):
    """ Appy func to a set of multidimensional arguments."""
    @wraps(func)
    def apply_over_array_general(*args, **kwargs):
        shape=np.shape(args[0])
        # Bypass x and y are 1D arrays.
        if np.ndim(args[0]):
            if len(args[0])>1 and isinstance(args[0], basestring)==0:
                if deepcopy:
                    cargs = [deepcopy(i) for i in args]
                else:
                    cargs = list(args)
                # Flatten all dimensions but the last to process more easely.
                #---------------------------------------------------
                for i in range(len(cargs)):
                    cargs[i] = (np.asanyarray(cargs[i]).flatten()).tolist()
                cargs = np.rollaxis(np.array(cargs), 1)
                # Function call.
                #---------------------------------------------------
                newargs = [func(*arg, **kwargs) for arg in cargs]
                # Verify if the function returned multiple arguments.
                if np.shape(len(newargs))==2:
                    for i in range(len(newargs)):
                        newargs[i] = newargs[i].reshape(*shape)
                else:
                    newargs = np.array(newargs).reshape(*shape)
            else:
                newargs = func(*args, **kwargs)
        else:
            newargs = func(*args, **kwargs)
        return newargs
    return apply_over_array_general

def phase_corr(v):
    """
    Return the phase in a range of [0,2*pi] in radians.
    """
    return (v+2*np.pi)%(2*np.pi)-np.pi

scale_prefixes = {'y':1e-24, 'z':1e-21, 'a':1e-18, 'f':1e-15, 'p':1e-12, 'n':1e-9, 'u':1e-6, 'm':1e-3, '':1, 'k':1e3, 'M':1e6,\
        'G':1e9, 'T':1e12, 'P':1e15, 'E':1e18, 'Z':1e21}

scale_val_prefixes = dict(zip(scale_prefixes.values(), scale_prefixes.keys()))

def rel_coord(ax, x, y):
    """
    Return coordinates at x and y percent of some 'x' and 'y' axis.
    scale : linear or log
    """
    xmin, xmax, ymin, ymax =ax.axis()
    if ax.get_xscale()=='linear':
        xt = xmin+(xmax-xmin)*x
    elif ax.get_xscale()=='log':
        xt = (xmax/xmin)**x*xmin

    if ax.get_yscale()=='linear':
        yt = ymin+(ymax-ymin)*y
    elif ax.get_yscale()=='log':
        yt = (ymax/ymin)**y*ymin
    return xt, yt

def figtext(ax, x, y, s):
    """ Adds text to subfigure at the relative coordinates (x,y).
    """
    x_, y_ = rel_coord(ax, x, y)
    ax.text(x_, y_, s)

def savefig(fig, name):
    """ Saves figure and crop the output pdf.
    """
    fig.savefig(name)
    os.system('pdfcrop %s %s'%(name,name))

def logarithm(x, base=10):
    """ Compute the logaritm for a given base """
    y = np.log10(x)/np.log10(base)
    return y

def OrderOfMagnitude(x, base=10):
    """ Gives the smallest integer power of a number expressed in a given base. """
    if isinstance(x, list):
        x = np.array(x)
    else:
        if x==0:
            return 0
    x = abs(x)
    y = logarithm(x, base)
    e = np.where(y<0, -np.ceil(abs(y)), y//1)
    return e

def OrderOfAmplitude(x) :
    """ Find the order of magnitude of x and its associated unit prefix. """
    scale_dict = dict(zip(OrderOfMagnitude(scale_prefixes.values()), scale_prefixes.keys()))
    order = np.sort(list(scale_dict.keys()))
    i = np.min(OrderOfMagnitude(x))
    if np.isinf(i):
        return i, ''
    elif np.isnan(i):
        return i, ''
    else:
        i -= i%3
        i = min(i, 21)
        i = max(i, -24) 
        return 10**int(i), scale_dict[i]

def DivideList(x,n) :
    for i in range(len(x)) :
        x[i] = x[i]/n
    return x

def scale(x, Unit, power=False, toint=0, latex=True) :
    """ Find the good scale for a number """
    if x == None or x =='' :
        return None,''
    if x == 0 :
        return x,Unit
    x = np.array(float(x))
    i = x
    if x.size>1:
        i = max(x.flatten())
    n,pre=OrderOfAmplitude(i)
    if latex:
        if pre=='u':
            pre = '$\mu$'
    x=float(x/n)
    if toint:
        x = np.round(x)
    if abs((x-np.round(x)))<1e-10:
        x = int(x)
    if power:
        return x, r'$\times$ 10$^{%s}$'%int(np.log10(scale_prefixes[pre]))+' '+Unit
    else:
        return x, pre+Unit

@apply_over_array_wrapper
def format_scale(x, units='', sep=' ', latex=True, decimal=None, significant=None, power=False):
    if type(x) in [np.string_,str,basestring]:
        return x
    else:
        if decimal==None:
            return sep.join(('{0} {1}'.format(*scale(x, units, power, latex=latex))).split(' '))
        else:
            cmd = '{0:.%sf} {1}'%decimal
            return cmd.format(*scale(x, units, power, latex=latex))

#def format_data(x, dx, units='', decimal=None, significant=None, power=False):
#    ox = OrderOfAmplitude(x)
#    odx = OrderOfAmplitude(dx)
#    cmd = '{0:.%sf} {1}'%decimal+ox/odx
#    s = cmd.format(x/np.min([10**ox,10**odx]))
#
#    s +=cmd.format(dx/np.min([ox,odx]))

def DegtoRad(deg) :
    """ Convert degres to radians """
    rad = deg/360.*2*np.pi
    return rad

def RadtoDeg(rad) :
    """ Convert radians to degres """
    deg = rad/(2*np.pi)*360
    return deg

#
# Instrument/experiment functions
#
def countdown(t):
    """ Wait a certain number of seconds and display the remaining time.
    """
    t = int(t)
    hours = t/3600
    minutes = (t%3600)/60
    seconds = (t%60)
    t1 = time.time()
    ind = 0
    while hours or minutes or seconds:
        ind +=1
        dec = '{0:02}:{1:02}:{2:02}'.format(hours, minutes, seconds)
        #print((dec,))
        print(dec,)
        #print(time.time())
        if not seconds:
            if minutes:
                minutes -= 1
            elif hours:
                hours -= 1
                minutes = 59
            seconds = 59
        else :
            seconds -= 1
        t2 = time.time()
        time.sleep(t1+ind-t2)
        #wait(t1+ind-t2)
        #print(('\r',))
        print('\r',)
    print('        ', '\r',)
    #print(('          ', '\r',))

#
# Conversions of logarithmic units
#
def dB(x, volts=True) :
    if volts:
        return 20*np.log10(abs(x))
    else:
        return 10*np.log10(abs(x))

def dBc2dBm(dBc,carrier_dBm):
        return dBc+carrier_dBm

def dBc2V(dBc, carrier, R=50):
    """
    Transform x in dBc to Volts.

    x : Value in dBc.
    carrier : Power of the initial signal. (1e-3 W for 'dBm' and 'dBc' for any
    other values.).
    R : resistance over which V is measured.
    """

    return np.sqrt(R*carrier)*10**(dBc/20.)

def dBm2W(x):
    """
    Transform x in dBm to Watts.

    x : Value in 'units'.
    """

    return 10**(x/10.)*1e-3

def dBm2V(dBm, R=50):
    """
    Transform x in dBm to Volts.

    x : Value in 'units'.
    R : resistance on which V is measured.
    """

    return dBc2V(dBm, 1e-3, R=R)

def V2dBc(V, carrier, R=50):
    """
    Transform x in Volts_rms to dBc.

    x : Value in 'units'.
    carrier : Power of the initial signal.
    R : resistance on which V is measured.
    """

    return 10*np.log10(V**2/(R*carrier))

def V2dBm(V, R=50):
    """
    Transform x in Volts_rms to dBm.

    x : Value in 'units'.
    R : resistance on which V is measured.
    """

    return V2dBc(V, 1e-3, R=R)

def C2F(x):
    """ Converts Celcius to Fahrenheit. """
    # celcius, fahrenheit
    x1, y1 = 0.,32.
    x2, y2 = 100.,212.
    a = (y2-y1)/(x2-x1)
    b = 32.
    return a*x+b

def F2C(x):
    """ Converts Fahrenheit to Celcius. """
    # celcius, fahrenheit
    x1, y1 = 32.,0.
    x2, y2 = 212.,100.
    a = (y2-y1)/(x2-x1)
    b = y2-a*x2
    return a*x+b

def Gal2L(x):
    """ Converts gallon to liters. """
    return x*3.78541

def L2Gal(x):
    """ Converts gallon to liters. """
    return x/3.78541

#
# File and data manipulation
#

def read_from_file(filename, directory='', convert=True, sep='\t', endline='\n', header=False, multitypes=False, regex=False):
    """ Extracts data from a file.
        filename : name of the file or pattern.
        convert : If True, convert strings to float when possible.
        sep : separator between each column.
        endline : endline character.
        multitypes : enable to have columns with different types.
        directory : path to look for the file. Can be a list of paths.
        regex : use regular expressions (filename should be a pattern).
    """
    multi_sweep = False
    if regex :
        filenames = search_files(filename, directory)
    else :
        if isinstance(directory, list):
            filenames = []
            for d in directory:
                filenames += glob.glob(d+filename)
        else:
            filenames = glob.glob(directory+filename)
    assert(np.size(filenames) != 0)
    a = []
    for fname in filenames:
        f=open(fname,'r')
        lines=f.readlines()
        _header = ''
        v=[]
        for line in lines:
            if line[0]=='#':
                _header += line
                if 'readback numpy shape for line part:' in line:
                    dim = [int(i) for i in line.split(':')[1].split(', ')]
                    multi_sweep = True
            elif line[0]==endline:
                pass
            else:
                lstr=line.strip().split(sep)
                if convert:
                    lfloat=[float_or_string(x) for x in lstr]
                else:
                    lfloat=[x for x in lstr]
                v.append(lfloat)
        if multitypes:
            # TODO : THERE SHOULD BE NO MULTITYPES WHEN SWEEP_MULTI IS ON.
            a.append([[float_or_string(a[i]) for i in range(len(a))] for a in zip(*v)])
        else:
            if multi_sweep:
                # (line, column) to (shape, column) to  (column, shape)
                tmp = np.array(v).reshape(*(dim+[-1]))
                a.append(np.rollaxis(tmp, -1))
            else:
                # (line, column) to (column, line)
                a.append(np.array(v).T)
    a = np.array(a)
    if len(filenames)==1:
        a = a[0]

    if header:
        return a, _header
    else:
        return a

def write_to_file(filename, v, comment):
    """ Write array of data to a file.
        filename : name of the file.
        v : array of data.
        comment : headline comment (must start with '#').
    """
    f=open(filename,'w')
    end = ''
    if comment:
        if comment[-1]!='\n':
            end = '\n'
    f.writelines(comment+end)
    for j in range(len(v[0])):
        for i in range(len(v)):
            term=''
            if i == len(v)-1:
                term='\n'
            f.write('{0}\t{1}'.format(v[i][j],term))
    f.close()

def find_files(directory):
    """
    Get the path of all files contained in the directory and all its
    subdirectories.
    """
    filenames=[]
    if type(directory)==list:
        pass
    else:
        directory = [directory]

    for dir in directory:
        for path, subdirs, files in os.walk(dir):
            for name in files:
                filenames+=[os.path.join(path,name)]
    return filenames

def search_files(file_pattern, directory):
    """ Search files with a certain pattern that are contained in a specified
        directory and all of its subdirectories.
    """
    list_names=find_files(directory)
    matches=[]
    for name in list_names:
        m=re.match(file_pattern,name.split('/')[-1])
        if m:
            matches+=[name]
    return matches

def no_dupl(v):
    """ Remove all redundancy in a set/list/array of values """
    return np.sort(list(set(v)))

def float_or_string(s):
    """ Convert a string to float in two ways.
        The first is a standard conversion (float(s)).
        The second convert a prefix to the right order of magnitude (digit+letter)
        ex. 8m -> 8e-3, 1M -> 1e6.
        If not possible, let the string untouched.
    """
    try:
        return float(s)
    except ValueError:
        pass
    except TypeError:
        return str(s)

    try:
        return float(s[:-1])*scale_prefixes[s[-1]]
    except (KeyError, ValueError):
        return str(s)

def get_data(patterns, dir_path, sort_order=1, exclude_files=[], relevant_groups=None):
    """ Read data from files and put them into a
        sorted multidimensional array.

        Returns : data, list of parameters, list of filenames
    """
    fnames = find_files(dir_path)
    if not isinstance(patterns, list):
        patterns = [patterns]
    for pattern in exclude_files:
        for name in fnames[:]:
            m = re.search(pattern,name)
            if m:
                fnames.remove(name)
    names_lst = []
    groups_lst = []
    data = []
    no_match_flag = 1
    for name in fnames:
        for pattern in patterns:
            m = re.search(pattern,name)
            if m:
                if no_match_flag:
                    no_match_flag = 0
                if len(m.groups())!=0:
                    #if 'off' in m.groups():
                    #   ind=1
                    # TODO deal with file with all parameters at off
                    groups = [float_or_string(x) for x in m.groups()]
                    groups_lst.append(groups[::sort_order])
                names_lst.append(name)
    if no_match_flag:
        raise Exception("No file found !")
    groups = list(map(np.array, zip(*groups_lst)))
    names = np.array(names_lst)
    # Handle silent groups in the ordering
    if relevant_groups is None:
        relevant_groups = np.array([True]*len(groups)).astype(bool)
    else: 
        relevant_groups = np.asanyarray(relevant_groups).astype(bool)
        if (len(relevant_groups)!=len(groups)):
            raise Exception("Relevant groups list not the right lenght !")
    if len(groups)>0:
        grps = np.array(groups)[relevant_groups][::-1].tolist()
        i = np.lexsort([names]+grps)
        names = names[i]
        groups = [g[i] for g in groups]
    data = np.array([list(read_from_file(nm)) for nm in names])
    if len(data.shape)<3:
        print "List of inhomogeneous files lenght"
    return data, groups, names

def axis_slice(x,first,last,step,axis=-1):
    """
    Slice an array x with respect to an given axis between 'start' and 'stop'
    with steps of lenght 'step'.

    Code taken from Christian's _do_stip() function in derivative.py.
    """
    ind = [slice(None)]*x.ndim
    ind[axis] = slice(first,last,step)
    ind = tuple(ind)
    return x[ind]

def flatten(x, axis=-1):
    """ Reshape an array to a 2D array at most.
    """
    tmp = np.rollaxis(x, axis)
    return np.squeeze(tmp.reshape((len(tmp),-1)).T)

def interp(x, xp, yp):
    """ A wrapper around np.interp that 
        make shure that the array are sorted
        in increasing order.
    """
    xp = np.asanyarray(xp)
    yp = np.asanyarray(yp)
    index_array = np.argsort(xp)
    xp = xp[index_array]
    yp = yp[index_array]
    return np.interp(x, xp, yp)


def interpolate(x, xp, yp, axis=-1):
    """ Returns a linear interpolation along a given axis of the 
        data points (xp, yp) along x.
    """
    xp = np.asanyarray(xp)
    if xp.ndim==1:
        _x = [x]
        _xp = [xp]
        _yp = [yp]
    else:
        _x = flatten(x, axis)
        _xp = flatten(xp, axis)
        _yp = flatten(yp, axis)
    return np.array([interp(i,j,k) for i,j,k in zip(_x,_xp,_yp)]).reshape(x.shape)

def nestedfor_comprehension(lst, func, **kwarg):
    """
    Creates a list of size shape(lst) using func like operation on each iteration.
    lst must be the shape of the final desired list.
    """
    index_lst =\
    ['b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
    command = 'new_lst = '+'['*len(lst)+'{0}'.format(func).format(*index_lst[:len(lst)])
    for a in range(len(lst)):
        command += ' for {0} in range({1})]'.format(index_lst[:len(lst)][::-1][a], lst[::-1][a])
    exec(command)
    return new_lst

#
# Mathematical functions
#
def xcothx(x):
    """ Mathematica function xcoth(x). """
    return np.where(np.tanh(x)==0.0, 1.0, x/np.tanh(x))

def _dxcothx(x):
    """ Mathematica function for the first derivative of xcoth(x). """
    if np.sinh(x)==np.inf:
        return 1.0
    elif np.sinh(x)==-np.inf:
        return -1.0
    elif np.sinh(x)==0.0:
        return 0.0
    else:
        return 1/np.sinh(x)*(np.cosh(x)-x/np.sinh(x))
dxcothx = np.vectorize(_dxcothx)

def Zeq_Parallel(*args):
    """ Computes equivalent resistance of parallel impedances. """
    a = 0
    for arg in args:
        a += 1./arg
    return 1./a

def ZC(f, C):
    """ Computes the impedance of a capacitor. """
    return 1./(1j*2*np.pi*f*C)

def ZL(f, L):
    """ Computes the impedance of an inductor. """
    return 1j*2*np.pi*f*L

def gamma(Z2, Z1=50.):
    """ Reflection coefficient of a medium Z1 to a medium Z2.
    """
    return (Z2-Z1)/(Z2+Z1)

def TN_vs_T(T, f, G, Ta):
    x = c.h*f/(2*c.k*T)
    return G*(T*xcothx(x)+Ta)

def Noise(V, T, f, F):
    x0 = np.where(2*c.k*T<1e-308, np.inf, c.h*f/(2*c.k*T))
    x1 = np.where(2*c.k*T<1e-308, np.inf, (c.e*V - c.h*f)/(2*c.k*T))
    x2 = np.where(2*c.k*T<1e-308, np.inf, (c.e*V + c.h*f)/(2*c.k*T))
    return (2*c.k*T)*( xcothx(x0)*(1-F)+F/2.*(xcothx(x1) + xcothx(x2)) )

def T_N(V, T, f, F):
    return Noise(V, T, f, F)/(2*c.k)

def mid_value_array(s, i, j ,k):
        """
        Change the values of the one dimensional array s[i,j,k] by its mid value.
        """
        ss = s[i,j,k]
        ls = len(ss)
        ll = ls/2
        return np.array([ss[ll]]*ls)

def mid_value(s, i, j ,k):
        """
        Returns the mid value of s[i,j,k] being a one dimensional array.
        """
        ss = s[i,j,k]
        ls = len(ss)
        ll = ls/2
        return ss[ll]

def Calib(s, Gain, Ta):
        return (s/Gain-Ta)

def R_from_Noise_plateau(i, f, r, delta=0):
    """
    Extract the resistance from the edges of the noise vaccum noise plateau.
    i : the current at the edge of the plateau.
    f : the frequency at which the noise is measured.
    """
    if delta:
        i, di = i
        f, df = f
        r, dr = r
        a = (c.e*i)/(c.h*f)
        R = 1/(a-1/r)
    if delta:
        dR = R**2*np.sqrt(a**2*( (di/i)**2+(df/f)**2 ) + (1/r)**2*(dr/r)**2)
        return np.array([R, dR])
    else:
        return R

def Z_from_dIdV(V, Vs, Rp, Rin):
    """ Computes the impedance from a dIdV measurement.
        Rp is the current polarisation impedance.
        Rin is the input impedance of the instrument (voltage reading port)

        Supposes that Rp and Rin are real impedances.
    """
    return 1./((Vs/V-1)*1./Rp-1./Rin)

def Lorentzian(x, A=1, B=1, x0=0):
    """ Lorentzian function 
        A : normalization factor (real and positive)
        B : Width at half height (real and positive)
        x0 : offset (real)

        With A=2/B, it is normalized as a probability distribution.
        With this choice, the parameters are independant, thus making
        it easier to fit.
    """
    A = np.abs(A)
    B = np.abs(B)
    x0 = np.real(x0)
    return A/np.pi*(B/2.)**2/((x-x0)**2+(B/2.)**2)

########
# For fit on dS photoassisted (or not) noise as a function of voltage.
########

class PrecisionError(Exception):
    def __init__(self, message, *errors):
        super(PrecisionError, self).__init__(message)
        self.errors = errors

def X(V, f, T, R, F, p=0, units='V2'):
    """ Computes the noise susceptibility
        of a coherent conductor.
    """
    v = c.e*V/(2*c.k*T)
    w = c.h*f/(2*c.k*T)
    assert(units in ['V2','None'])
    if units=='None':
        a = 1.
    elif units=='V2':
        a = c.k*T*R
    return a*F * (xcothx(v+w)+(-1)**p*xcothx(v-w))

def Xphoto(V, Vac, f, f0, T, R, F, p, order, units='V2'):
    """ Computes the noise susceptibility
        of a photoassisted coherent conductor.
    """
    z = c.e*Vac/(c.h*f0)
    func = lambda n: X(V, f+n*f0, T, R, F, p, units=units)
    return np.sum(BesselWeighting(z, func, order, p), axis=0)

def SI(V, Vac, f, f0, T, R, F, p, order, units='V2'):
    """Computes the in phase quadrature of photoassisted noise
       of a coherent conductor 
    """
    return Sphoto(V, Vac, f, f0, T, R, F, order, units=units) + Xphoto(V, Vac, f, f0, T, R, F, p, order, units=units)

def SQ(V, Vac, f, f0, T, R, F, p, order, units='V2'):
    """Computes the out of phase quadrature of photoassisted noise
       of a coherent conductor 
    """
    return Sphoto(V, Vac, f, f0, T, R, F, order, units=units) - Xphoto(V, Vac, f, f0, T, R, F, p, order, units=units)

def S(V, f, T, R, F, units='V2'):
    """ Compute the noise of a coherent conductor.
    """
    T = np.abs(T)
    R = np.abs(R)
    F = np.abs(F)
    v = c.e*V/(2*c.k*T)
    w = c.h*f/(2*c.k*T)
    assert(units in ['V2','None'])
    if units=='None':
        a = 1.
    elif units=='V2':
        a = c.k*T*R
    return a*(F*(xcothx(v-w) + xcothx(v+w))+(1.-F)*2*xcothx(w))

def dS(V, f, T, R, F, variable='None', units='None'):
    """ Compute the derivative of the noise of a coherent conductor.
    """
    # Note : There is a subtlety with a derivative take with respect to the current.
    #        The link between current and voltage is ohm's law dV = Z dI and Z can
    #        be complex. So the passage from one to another can make the derivative
    #        a complex quantity. Is there anything wrong with that ?
    v = c.e*V/(2*c.k*T)
    w = c.h*f/(2*c.k*T)
    # We derive with respect to Voltage V and v = eV/kT
    assert(variable in ['None','V','I'])
    if variable=='None':
        a = 1.
    elif variable=='V':
        a = c.e/(2*c.k*T)
    elif variable=='I':
        a = R*c.e/(2*c.k*T)

    assert(units in ['V2','None'])
    if units=='None':
        b = 1/2.
    elif units=='V2':
        b = c.k*T*R
    return a*b*F*(dxcothx(v-w) + dxcothx(v+w))

def Sphoto(V, Vac, f, f0, T, R, F, order, units='V2'):
    """
    Computes the photoassisted noise (voltage variance) of a tunnel junction.
        V : bias voltage
        Vac : AC excitation
        f : detection frequency
        f0 : excitation frequency
        T : electronic temperature
        R : sample resistance
        F : Fano factor
        order : maximal order of the bessel function
    """
    #v = c.e*V/(2*c.k*T)
    #w = c.h*f/(2*c.k*T)
    #w0 = c.h*f0/(2*c.k*T)
    z = c.e*Vac/(c.h*f0)
    # Voir mémoire G. Gasse.
    # In G. Gasse M. Thesis, g(x)=xcoth(x/2). Here
    # g(x)=xcothx(x). Therefore I omit the factor 1/2
    # over the sum.
    #f = lambda n: c.k*T*R*(F*(xcothx(v-w+n*w0) + xcothx(v+w+n*w0))+(1.-F)*2*xcothx(w))
    func = lambda n: S(V, f+n*f0, T, R, F, units=units)
    return np.sum(BesselWeighting(z, func, order), axis=0)

def dSphoto(V, Vac, f, f0, T, R, F, order, variable='None', units='None'):
    """
    Computes the derivative of the photoassisted noise (voltage variance) of a tunnel junction. 
        V : bias voltage
        Vac : AC excitation
        f : detection frequency
        f0 : excitation frequency
        T : electronic temperature
        order : maximal order of the bessel function
        normalized : normalized to one at +/-infinity when true
    """
    z = c.e*Vac/(c.h*f0)
    func = lambda n: dS(V, f+n*f0, T, R, F, variable=variable, units=units)
    return np.sum(BesselWeighting(z, func, order), axis=0)

def BesselWeighting(z, f, order, p=0, kind='first', threashold=1e-6):
    """ Applies bessel function as weight up to a certain order.
            f : function to weight (should depend only of a integer parameter n)
            order : order of the maximal bessel function
            kind : kind of bessel funciton (either first of second)
    """
    n = np.arange(-order, order+1)
    # Bessel of the first kind
    if kind=='first':
        Jz = lambda n: sp.jv(n,z)
    elif kind=='second':
        raise NotImplementedError
    cond = np.abs(1-np.sum([Jz(i) for i in n]))<threashold
    if cond:
        return np.array([Jz(i)*Jz(i+p)*f(i) for i in n])
    else: 
        message = "The order (%s) is too low compared to 'z = %s'. This might lead to imprecisions higher then %s."%(order, z, threashold)
        raise PrecisionError(message)

########
# Data manipulation
########
def remove_peaks(x, std=2, counts=3):
    """ Replace dispersed points of a set of data with np.nan.

        Parameters :
            std : Number of standard deviation.
            counts : Number of recursion (each recursion splits the set of points in half).
    """
    if counts==0:
        return x
    x1 = remove_peaks(x[:len(x)/2], std, counts-1)
    x2 = remove_peaks(x[len(x)/2:], std, counts-1)
    def crop(x, std):
        a = np.arange(len(x))
        ind_ = [i for i in range(len(a)) if not np.isnan(x[i])]
        xl = np.poly1d(np.polyfit(a[ind_], x[ind_], 1))
        x0 = x - xl(np.arange(len(x)))
        s = np.nanstd(x0)
        x_new = np.where(abs(x0)<std*s, x, np.nan)
        return x_new
    x1_new = crop(x1, std)
    x2_new = crop(x2, std)
    x_new = np.array(x1_new.tolist()+x2_new.tolist())
    return x_new

def logDiff(a, b):
    """ Computes the difference of two numbers in dB units by converting them into linear units.
        Returns a value in dB units.
    """
    if not a>b:
        print("Please put the higher number first : a>b")
        raise
    return 10*np.log10(10**(a/10.)-10**(b/10.))

def lighten(x, y, logxsep, average=False):
    """ Make uniform the number of points of a logarithmic scale with about one point per 'logxsep'.
        Generally used when data is taken in a linear scale. It is used to lighten the plot by
        decreasing the number of points while keeping the density of point to a minimal value.
    """
    logx = np.log10(x)
    dx = np.diff(logx)
    indx = []
    i = 0
    while i<len(dx):
        ind = [i]
        s = dx[i]
        i += 1
        if i==len(dx):
            break
        while s<logxsep:
            s += dx[i]
            ind.append(i)
            i += 1
            if i==len(dx):
                break
        indx.append(ind)

    y_new = []
    x_new = []
    for ind in indx:
        if average:
            y_new.append(np.mean(y[ind[0]:ind[-1]+1]))
            x_new.append(np.mean([x[ind[0]],x[ind[-1]]]))
        else:
            i = int(np.mean([ind[0],ind[-1]]))
            y_new.append(y[i])
            x_new.append(x[i])

    return [np.array(x_new),np.array(y_new)]

def fftFilter(x, y, cutoff=5, axis=-1):
    """ NOTE : INSERT DOCU HERE
    """
    w = np.apply_along_axis(scipy.fftpack.rfft, axis, y)
    func = lambda x : scipy.fftpack.rfftfreq(len(x), x[1]-x[0])
    f = np.apply_along_axis(func, axis, x)
    spectrum = w**2

    cutoff_idx = spectrum < (spectrum.max()/cutoff)
    w2 = w.copy()
    w2[cutoff_idx] = 0

    return np.apply_along_axis(scipy.fftpack.irfft, axis, w2)

@apply_along_axis_wrapper
def peakfilter(x, y, func=None, p0=None, std=5):
    """ Removes peaks from data. It uses func to fit the data, then
        based on the difference between data and fit, it removes
        values over a certain standard deviation.
    """
    p = curve_fit(func, x, y, p0=p0)
    fit = np.array(func(x, *p[0]))
    # Weight the std so that fluctuation can be propotionnal to the slope
    # of the fitted curve.
    # *25 is a choice so that the weight is well balanced.
    weight = np.array(np.abs(np.array([0]+np.diff(fit/np.std(np.abs(fit))*25).tolist())))+1
    y[np.abs((y-fit))/weight>(np.std(y-fit, axis=-1)*std)[...,None]] = None
    return y

@apply_along_axis_wrapper
def fit(x, y, func=None, p0=None, **kwargs):
    params = curve_fit(func, x, y, p0, **kwargs)
    return params, func(x, *params[0])

@apply_along_axis_wrapper
def interpolate1d(x,y, **kwargs):
    return interp1d(x,y, **kwargs)

def expand_dims(x, extra_shape):
    """ Expand dimension of array giving it a new shape
        equal to the old shape plus the extra shape.
    """
    x = np.asanyarray(x)
    init_shape = x.shape
    return np.squeeze(np.array([(x.T).tolist()]*np.prod(extra_shape)).T.reshape(init_shape+extra_shape))

def change_in_array(x):
    """ Return an array with integer elements starting from 0 
        reproducing cyclic nature of the x array.
    """
    val = 0
    c = [val]
    prev = x[0]
    for i in x[1:]:
        if i==prev:
            pass
        else:
            val += 1
            prev = i
        c.append(val)
    return np.array(c)

def dissociate_groups(g):
    c = [change_in_array(i) for i in g]
    mask = [True]*len(g)
    indexes = mask[:]
    while len(indexes):
        ind = indexes.pop(0)
        tmp_indexes = indexes[:]
        for j in tmp_indexes:
            if np.all(c[ind]==c[j]):
                mask[j] = False
                indexes.remove(j)
    return mask

def unique(x):
    """ Find the unique elements of an array
        and keep the order relative to the 
        first occurence of each element.
    """
    x = np.array(x.flatten()[:])
    _,idx = np.unique(x, return_index=True)
    return x[np.sort(idx)]

def gradient(x,y, axis=-1):
    """ Computes the gradient of y relative to x 
        along an given axis.
    """
    x = np.asanyarray(x)
    y = np.asanyarray(y)
    a = np.swapaxes(np.zeros(x.shape), axis, 0)
    dy = np.diff(y, axis=axis)
    dx = np.diff(x, axis=axis)
    p = np.swapaxes(dy/dx, axis, 0) # The first dimension is the one to work with.
    a[0] = p[0]
    a[-1] = p[-1]
    a[1:-1] = (p[:-1]+p[1:])/2.
    # Reshape the array to the shape of x and y.
    return np.swapaxes(a, axis, 0)

def I2V(I, R, zero=0, axis=-1):
    """ Transforms differential resistance measurement (R)
        using a current source to voltage by integrating (RI)
        
        parameters: 
            -zero : specifies the position of the zero voltage point
                 relative to the current (I).
            -axis : axis upon which to integrate
    """
    I = np.asanyarray(I)
    R = np.asanyarray(R)
    assert(I.shape==R.shape)
    if I.ndim>1:
        I = np.rollaxis(I, axis=axis)
        R = np.rollaxis(R, axis=axis)
    shpe = np.shape(I)
    a = np.array([np.trapz(R[:i], I[:i], axis=0) for i in range(1,np.shape(I)[0]+1)])
    # NOTE : could interpolate I to get the zero offset
    if I.ndim>1:
        # Give the initial shape back
        a -= a[I==zero].reshape(shpe[1:])[None]
        return np.transpose(a, range(1,a.ndim)+[0])
    else:
        a -= a[I==zero]
        return a


#def fit(x0, y0, deg=25):
#    _fit = np.poly1d(np.polyfit(x0, y0, deg))
#    return _fit

#def lockin_averaging(x, a, x0, y0, deg=20, type='SIN', axis=-1):
def lockin_averaging(x, a, fit, xlim=None, type='SIN', axis=-1, npts=1001):
    """ Makes an a*sin(w*t) averaging of the fit function evaluated around x.
        The 'axis' argument is not implemented yet.
        Note : The fit function has an unknown validity range, so be carefull when
               interpreting the output for larger values to a compared to the x0 value
               used for the fitting.
    """
    if xlim==None:
        xlim = (np.min(x), np.max(x))
    def basefunc(x, a, fit, lim):
        if (x+a/2.)>=lim[1] or (x-a/2.)<=lim[0]:
            return np.nan
        else:
            if type == 'LIN':
                t = np.linspace(x-a/2., x+a/2., npts)
                return np.average(fit(t))
            if type == 'SIN':
                t = np.linspace(0, 2*np.pi, npts)
                return np.average(fit(x+a/2.*np.sin(2*np.pi*t)))
                #return quad(lambda t: fit(x+a/2.*np.sin(2*np.pi*t)), 0, 1)[0]
    return np.vectorize(lambda _x: basefunc(_x, a, fit, lim=xlim))(x)

def cut3d(X,Y,Z,N, zdir='y'):
    """ Returns an array with the elements of 3d cuts through the 3d graph.
        
        Arguments :
            X, Y, Z : 2d array representing the x, y and z coordinates
            N : Either the number of cuts or an array specifying the value of each cuts.
            zdir : the axis upon which cut ('x', 'y', 'z')
    """
    # Create fig
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection='3d')
    l = ax.contour(X, Y, Z, N, zdir=zdir)
    mask = np.array(['x', 'y', 'z'])!=zdir
    # Close fig
    plt.close(fig)
    # Transform into 2d list to prevent ValueError related to broadcasting the list of array
    # to an single array when shape is not homogeneous.
    cols = np.array([(np.array(map(list, i._segments3d[0]))[:,mask]).T.tolist() for i in l.collections])

    return np.rollaxis(cols, 1), l.levels

def fwhm(x,y):
    """ Estimates the FWHM, and maximum value from raw data """
    mx = np.max(y, axis=-1)
    mn = np.min(y, axis=-1)
    mask = y>(mx+mn)[:,None]/2. # Generate mask
    l = len(mask)
    for i in range(len(mask)):
        indx = np.where(mask[i]==True)[0][1:-1] # Only keep the outmost values
        mask[i,indx] = False                    # and put the rest to false
    xn = x[mask].reshape((l,2))
    yn = y[mask].reshape((l,2))
    return np.diff(xn, axis=-1) # Get a list of full width at half maximum

def heat_water(T, V1, V2, W, Cg=0.4, T0=20.):
    """ T : Target temperature (C)
        V1 : Input volume (L)
        V2 : Volume of water already in the mash (L)
        W : Weight of grains in the mash (kg)
        Cg : Relative heat capacity of grain over water (kg/L)
        T0 : Initial temperature of the mash
    """
    V1 *= 1.
    return T+(V2/V1+W*Cg/V1)*(T-T0)
