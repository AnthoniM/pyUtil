#! /usr/bin/env python
# -*- coding: utf-8 -*-

# author : Anthoni Manseau

import numpy as np
import matplotlib.pyplot as plt
import corefunctions as cf
import sys
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
# implement the default mpl key bindings
from matplotlib.backend_bases import key_press_handler
from copy import copy
if sys.version_info[0] < 3:
    import Tkinter as Tk
else:
    import tkinter as Tk

class SweepData(object):
    """ Reads data from regular expression pattern using groups and reshape the data.

        Parameters:
            pattern : regular expression. For a sweep that generated multiple files and
                      a header file, it should contain the  following string : '.*readval_(\d+).*'.
                      If the pattern does not contain the '.*readval_(\d+).*' string it will be
                      assumed that there is no header file.
            directory : list of paths pointing to specific directories.

        Optional Parameters:
            -verbose : use that option to get additionnal information during the initiation of
                       the class.
            -field_names : list of name given to each field of the regular expression (pattern).
            -relevant_fields : list of bool to indicate which group is relevant in the shaping of data.
            -column_names : list of the name for each column of the files.
            -param_names : list of the name for the sweep parameters.
            -exclude_files : list of pattern of files to exclude from search.

        Output attributes :
            field_names attributes : attributes with the name associated with the field_names list.
                                    Default values are 'self.field1', 'self.field2', ...
            column_names attributes : attributes with the name associated with the column_names list.
                                     Default values are 'self.col1', 'self.col2', ...
            param_names attributes : attributes with the name associated with the param_names list.
                                     Default values are 'self.param1', 'self.param2', ...
    """
    def __init__(self, pattern, directory, **kwargs):
        # Initiate pattern and directory.
        if pattern[-1]=='$':
            # Make shure that the pattern end with the ends
            # of line search charater.
            self._pattern = pattern
        else:
            self._pattern = pattern+'$'
        self._directory = directory
        self._verboseFlag = kwargs.get('verbose', False)
        self._field_names = kwargs.get('field_names', [])
        self._relevant_fields = kwargs.get('relevant_fields', None)
        if self._relevant_fields is not None:
            #raise NotImplementedError, "This functionnality has not been implemented yet !"
            self._relevant_fields = np.asanyarray(self._relevant_fields).astype(bool)
        self._column_names = kwargs.get('column_names', None)
        self._param_names = kwargs.get('param_names', None)
        self._exclude_files = kwargs.get('exclude_files', [])
        self._homogenous = kwargs.get('homogenous', True)

        # Initialize  some gobal parameters
        self._multi_sweep_shape = tuple([])
        self._e = None

        self._verbose('', '-'*60)
        self._verbose('', '-'*60)
        self._verbose('', '-'*60)
        self.__readval = False # If there is no readval field. It is important to take into account all the fields in the pattern.
        if '.*readval_(\d+).*' in self._pattern:
            self._head_file_pattern = self._pattern.replace('.*readval_(\d+).*', '((?!readval).)*')
            if len(cf.search_files(self._pattern, self._directory)):
                self.__readval = True
                self._verbose("Header files.")
                self._load_head()
            else:
                self._verbose("No header file.")
                self._pattern = self._head_file_pattern
        else:
            self._verbose("No header file.")
        self._load()
        # Don't want to call it twice.
        #if self._param_names and not self.__readval:
        self._set_param_names()
        self._set_fields()
        self._set_column_names()
        self._verbose('', '-'*60)
        self._verbose('', '-'*60)
        self._verbose('', '-'*60)

    @property
    def data(self):
        return self._data
    @data.setter
    def data(self, x): 
        self._data = x
        self._set_column_names()
        self._set_param_names()
        self._set_fields()
        self._transformations()

    @property
    def head_data(self):
        return self._head_data
    @head_data.setter
    def head_data(self, x):
        self._head_data = x
        self._set_column_names()
        self._set_param_names()
        self._set_fields()
        self._transformations()

    def _transformations(self):
        pass

    def _verbose(self, x, name=''):
        if self._verboseFlag:
            if name=='':
                print x
            else:
                print name+' : ', x

    def _reshape(self, data, groups, filenames, head=False):
        """ The final shape is the following (c, e, *g, *m, r)
                c : column dimension
                e : extra dimension (when there is(are) unknown(s) dimension(s))
                g : groups shape from the regular expression in the order they appear
                m : multi-sweep shape
                r : rows
                * : means there could be more then one dimension
        """
        # Extract shape of a single file and get rid of single length dimensions
        single_file_shape = np.squeeze(data[0]).shape
        self._verbose(single_file_shape, 'Single file shape')

        # Extract shape of each group
        if head or self.__readval:
            groups = groups[:-1]
        unique_groups = np.array([np.unique(i) for i in groups]) 
        self._verbose(unique_groups, "unique_groups")
        groups_shape = tuple([])
        # Select specific field relevant for the reshaping
        if len(unique_groups):
            if self._relevant_fields is None:
                g = unique_groups
            else:
                g = unique_groups[self._relevant_fields]
            if self.__readval and not head:
                # Get rid the the readval group for reshaping
                g = unique_groups[:-1]
            groups_shape = tuple([len(i) for i in g])
            self._verbose(groups_shape, 'List of groups lenght (l)')
            grps = [np.reshape(g, groups_shape+(-1,)) for g in groups]
        else:
            self._verbose(groups_shape, 'List of groups lenght (l)')
            grps = []

        if head:
            multi_sweep_shape = tuple([])
        else:
            multi_sweep_shape = self._multi_sweep_shape
            self._verbose(multi_sweep_shape, 'Multi-sweep shape')

        # Extract the extra-dimension lenght if any exist
        extra_shape = (filenames.shape[0]*1./np.prod(groups_shape+multi_sweep_shape))
        assert(extra_shape%1==0)
        if extra_shape==1: 
            extra_shape = tuple([])
        else:
            extra_shape = tuple([int(extra_shape)])
            if len(unique_groups):
                
                grps = [np.rollaxis(g, -1) for g in grps]
        self._verbose(extra_shape, 'Extra-dimension lenght')

        # Reshape filenames
        self._verbose('\n# --- Reshaping filenames --- ')
        # If multi_sweep_shape is present it will always come before extra_shape
        new_filename_shape = groups_shape+multi_sweep_shape+extra_shape
        self._verbose(filenames.shape, 'Initial filenames.shape')
        filenames = np.squeeze(filenames.reshape(new_filename_shape))
        if extra_shape:
            # Put the extra shape right after the columns
            filenames = np.rollaxis(filenames, -1)
        self._verbose(filenames.shape, 'Final filenames.shape')

        # Reshape data
        self._verbose('\n# --- Preshaping data --- ')
        self._verbose(data.shape, 'Initial data shape')
        new_data_shape = new_filename_shape+single_file_shape
        # new shape and remove lenght 1 dimensions.
        data = np.squeeze(data.reshape(new_data_shape))
        # Put the columns first
        # shape before rollaxis : (*groups, column, *multi_sweep)
        if extra_shape:
            # Put the extra shape right after the columns
            data = np.rollaxis(data, -(len(single_file_shape)+1))
        data = np.rollaxis(data, -len(single_file_shape))
        self._verbose(data.shape, 'Final data shape')
        self._extra_shape = extra_shape

        # Get the shape of the parameters used in the multi-sweep
        # shape : (column, *groups, *multi_sweep)
        # Get rid of the group dimensions and the column dimension.
        if head:
            self._multi_sweep_shape = data.shape[-len(single_file_shape[1:]):]
            self._verbose(self._multi_sweep_shape, 'Multi-sweep shape')

        # Adds dimension to groups and parameters to be of the same shape as columns.
        if self._homogenous:
            expand_shape = data.shape[1+len(groups_shape+extra_shape):]
            self._verbose(expand_shape, "expand shape")
            groups = [cf.expand_dims(g, expand_shape) for g in grps]

        return data, groups, filenames

    def _load_head(self):
        # Load data by groups.
        self._verbose('-'*30+'Head Pattern'+'-'*30)
        self._verbose(self._head_file_pattern)
        self._verbose('-'*30+'Head'+'-'*30)
        self._head_data, self._head_groups, self._head_fnames = self._reshape(*cf.get_data(self._head_file_pattern, self._directory, exclude_files=self._exclude_files, relevant_groups=self._relevant_fields), head=True)
        #self._head_data, self._head_groups, self._head_fnames = self._reshape(self._head_data, self._head_groups[:-1], self._head_fnames, head=True)
        self._verbose('-'*30+'Done'+'-'*30)

    def _load(self):
        """ Load all files that match 'pattern' contained in 'directory'.
                Class them by group and reshape it with respect to the number
                of distinct elements in each group.
                The number of values for each fixed parameters must be the same.
        """
        # Load data by groups.
        self._verbose('-'*30+'Pattern'+'-'*30)
        self._verbose(self._pattern)
        self._verbose('-'*30+'Data'+'-'*30)
        #self._data, self._groups, self._fnames = cf.get_data(self._pattern, self._directory, exclude_files=self._exclude_files)
        #self._ncol = np.shape(self._data)[1:][0] # Finds the number of columns.
        #self._col_axis = 1-len(self._data.shape) # Finds the axis of for the column dimension.
        #self._data, self._groups, self._fnames = self._reshape(self._data, self._groups, self._fnames)
        self._data, self._groups, self._fnames = self._reshape(*cf.get_data(self._pattern, self._directory, exclude_files=self._exclude_files, relevant_groups=self._relevant_fields))
        self._verbose('', '-'*30+'Done'+'-'*30)

        if self.__readval and self._homogenous:
            expand_shape = self._data.shape[1+len(self._groups)+len(self._extra_shape)+len(self._multi_sweep_shape):]
            tmp = np.array([cf.expand_dims(g, expand_shape) for g in self._head_data])
            self._head_data = tmp

    def _set_fields(self):
        """ Set the field attributes using the list of names 'field_names'.
            If the list is not existing, the default names will be : field1, field2, etc.
        """
        if self._field_names==[]:
            self._field_names = [ 'field%s' %(i+1) for i in range(len(self._groups))]
        assert(len(self._field_names)==len(self._groups))
        for i in range(len(self._field_names)):
            setattr(self, '_'+self._field_names[i]+'_u', cf.unique(self._groups[i]))
            setattr(self, self._field_names[i], self._groups[i])

    def _set_param_names(self):
        """ Set the parameter iterated through the sweep attribute using the list of names 'param_names'.
            If the list is not existing, the default names will be : param1, param2, etc.
        """
        if self.__readval:
            d = self._head_data
        else:
            d = self._data
        if self._param_names==None:
            self._param_names = [ 'param%s' %(i+1) for i in range(d.shape[0])]
        for i in range(len(self._param_names)):
            setattr(self, '_'+self._param_names[i]+'_u', cf.unique(d[i]))
            setattr(self, self._param_names[i], d[i])

    def _set_column_names(self):
        """ Set name for each column and create attribute.
        """
        if self._column_names==None:
            self._column_names = [ 'col%s' %(i+1) for i in range(self._data.shape[0])]
        self._verbose(self._column_names, 'Column_names')
        for i in range(len(self._column_names)):
            setattr(self, self._column_names[i], self._data[i])

class Plot(object):
    def __init__(self):
        pass

    def _select(self, axes, selector, skip=[1,1]):
        swap = []
        for i,s in enumerate(selector):
            if s not in ['x','y','z']:
                swap.append((i,s))
        def select(x):
            # Revert swap order so that when a dimension vanishes the 
            # next swap index is the right one
            for i,s in swap[::-1]:
                x = np.rollaxis(x, axis=i)[s]
            return x
        X, Y, Z = [None]*3
        for i,s in enumerate(selector):
            if s=='x':
                X = select(axes[i])[::skip[0],::skip[1]]
            if s=='y':
                Y = select(axes[i])[::skip[0],::skip[1]]
            if s=='z':
                Z = select(axes[i])[::skip[0],::skip[1]]
        return X,Y,Z

    def _plot(self, x, y, **kwargs):
        figsize = kwargs.get('figsize', [12,9])
        fig = kwargs.get('fig', None)
        ax = kwargs.get('ax', None)
        if fig is None:
            fig, ax = plt.subplots(1,1,figsize=figsize)
        if ax is None:
            ax = fig.axes[0]

        if kwargs.get('reset_colors', False):
            ax.set_prop_cycle(None)

        xunits = kwargs.get('xunits', 1.)
        yunits = kwargs.get('yunits', 1.)

        # Plot
        params = kwargs.get('params', dict())
        lines = ax.plot(xunits*x, yunits*y, **params)
        
        # Labels
        xlabel = kwargs.get('xlabel', '')
        ax.set_xlabel(xlabel) 
        ylabel = kwargs.get('ylabel', '')
        ax.set_ylabel(ylabel) 
        # Scales
        xscales = kwargs.get('xscales', 'linear')
        ax.set_xscale(xscales) 
        yscales = kwargs.get('yscales', 'linear')
        ax.set_yscale(yscales) 
        # Limits
        xlimits = kwargs.get('xlimits', (None,None))
        ax.set_xlim(xlimits) 
        ylimits = kwargs.get('ylimits', (None,None))
        ax.set_ylim(ylimits) 
        # Title 
        suptitle = kwargs.get('suptitle', '')
        fig.suptitle(suptitle)
        # Save figure
        filename = kwargs.get('filename', '')
        # Legend
        bbox_to_anchor = kwargs.get('bbox_to_anchor', (1,0.7))
        legend_labels = kwargs.get('legend_labels', None)
        figtext = kwargs.get('figtext', None)
        if figtext:
            for args in figtext:
                fig.text(*args)
        if legend_labels:
            #            if fig.legends:
            #    fig.legends.pop()
            fig.subplots_adjust(left=0.1, right=0.86)
            fig.legend(lines, legend_labels, bbox_to_anchor=bbox_to_anchor)
        save = kwargs.get('save', False)
        savepath = kwargs.get('savepath',False)
        if save and filename and savepath:
            fig.savefig(savepath+filename)
        return fig

    def _plot_complex(self, x, y, **kwargs):
        figsize = kwargs.get('figsize', [12,9])
        fig = kwargs.get('fig', None)
        if fig is None:
            fig, [ax1,ax2] = plt.subplots(2,1,figsize=figsize)
        else:
            [ax1, ax2] = fig.axes

        if kwargs.get('reset_colors', False):
            ax1.set_prop_cycle(None)
            ax2.set_prop_cycle(None)
        real_imag = kwargs.get('real_imag', False)
        date = kwargs.get('date', '')
        f1 = np.abs
        f2 = lambda x: np.angle(x, 1)
        if real_imag:
            f1 = np.real
            f2 = np.imag
        # Units
        xunits = kwargs.get('xunits', [1.,1.])
        yunits = kwargs.get('yunits', [1.,1.])
        # Plot
        params1, params2 = kwargs.get('params', [dict(),dict()])
        lines = ax1.plot(xunits[0]*x.T, yunits[0]*f1(y.T), **params1)
        ax2.plot(xunits[1]*x.T, yunits[1]*f2(y.T), **params2)
        # Labels
        xlabels = kwargs.get('xlabels', ['',''])
        [ax.set_xlabel(l) for ax,l in zip([ax1,ax2],xlabels)]
        ylabels = kwargs.get('ylabels', ['',''])
        [ax.set_ylabel(l) for ax,l in zip([ax1,ax2],ylabels)]
        # Scales
        xscales = kwargs.get('xscales', ['linear','linear'])
        [ax.set_xscale(l) for ax,l in zip([ax1,ax2],xscales)]
        yscales = kwargs.get('yscales', ['linear','linear'])
        [ax.set_yscale(l) for ax,l in zip([ax1,ax2],yscales)]
        # Limits
        xlimits = kwargs.get('xlimits', [(None,None),(None,None)])
        [ax.set_xlim(l) for ax,l in zip([ax1,ax2],xlimits)]
        ylimits = kwargs.get('ylimits', [(None,None),(None,None)])
        [ax.set_ylim(l) for ax,l in zip([ax1,ax2],ylimits)]
        # Title 
        suptitle = kwargs.get('suptitle', '')
        fig.suptitle(suptitle)
        subtitles = kwargs.get('subtitles', ['',''])
        [ax.set_title(l) for ax,l in zip([ax1,ax2],subtitles)]
        # Figtext
        figtext = kwargs.get('figtext', None)
        if figtext:
            for args in figtext:
                fig.text(*args)
        # Legend
        bbox_to_anchor = kwargs.get('bbox_to_anchor', (1,0.7))
        legend_labels = kwargs.get('legend_labels', None)
        if legend_labels:
            fig.subplots_adjust(left=0.1, right=0.86)
            fig.legend(lines, legend_labels, bbox_to_anchor=bbox_to_anchor)
        # Save figure
        filename = kwargs.get('filename', '')
        fig.text(0.92, 0.01, date)
        save = kwargs.get('save', False)
        savepath = kwargs.get('savepath',False)
        if save and filename and savepath:
            fig.savefig(savepath+filename)
        return fig

    def _plot3D(self, X, Y, Z, **kwargs):
        figsize = kwargs.get('figsize', [12,9])
        fig = plt.figure(figsize=figsize)

        ax = fig.add_subplot(1,1,1, projection='3d')
        rstride = kwargs.get('rstride', 1)
        cstride = kwargs.get('cstride', 1)
        cmap = kwargs.get('cmap', cm.coolwarm)
        linewidth = kwargs.get('linewidth', 0.1)
        xunits = kwargs.get('xunits', 1.)
        yunits = kwargs.get('yunits', 1.)
        zunits = kwargs.get('zunits', 1.)

        # Plot
        ax.plot_surface(xunits*X, yunits*Y, zunits*Z, 
                        rstride=rstride, 
                        cstride=cstride, 
                        cmap=cmap, 
                        linewidth=linewidth)
        offx, _ = ax.get_xlim()
        _, offy = ax.get_ylim()
        offz, _ = ax.get_zlim()
        ax.contour(X, Y, Z, zdir='z', offset=offz, cmap=cm.coolwarm)
        ax.contour(X, Y, Z, zdir='x', offset=offx, cmap=cm.coolwarm)
        ax.contour(X, Y, Z, zdir='y', offset=offy, cmap=cm.coolwarm)

        
        view = kwargs.get('view',(25,-50))
        ax.view_init(*view)
        #fig.text(0.92, 0.01, date)

        # Labels
        xlabel = kwargs.get('xlabel', '')
        ax.set_xlabel(xlabel) 
        ylabel = kwargs.get('ylabel', '')
        ax.set_ylabel(ylabel) 
        zlabel = kwargs.get('zlabel', '')
        ax.set_zlabel(zlabel) 
        # Scales
        xscales = kwargs.get('xscales', 'linear')
        ax.set_xscale(xscales) 
        yscales = kwargs.get('yscales', 'linear')
        ax.set_yscale(yscales) 
        zscales = kwargs.get('zscales', 'linear')
        ax.set_yscale(zscales) 
        # Limits
        xlimits = kwargs.get('xlimits', (None,None))
        ax.set_xlim3d(xlimits) 
        ylimits = kwargs.get('ylimits', (None,None))
        ax.set_ylim3d(ylimits) 
        zlimits = kwargs.get('zlimits', (None,None))
        ax.set_zlim3d(zlimits) 
        # Title 
        suptitle = kwargs.get('suptitle', '')
        fig.suptitle(suptitle)
        # Save figure
        filename = kwargs.get('filename', '')
        save = kwargs.get('save', False)
        savepath = kwargs.get('savepath', False)
        if save and filename and savepath:
            fig.savefig(savepath+filename)
        return fig

class ExploreParameters(object):
    # TODO : It would be nice to ENABLE TO USE A LIST OF FUNCTION that share the same x axis.
    #        Then to solve the set of independant parameters (by an union operation on sets).
    #        For each function create a mask give the right parameters to the function. 
    # 
    #        Also should be able to DISPLAY A LIST OF INITAL DATA corresponding to each functions.
    #
    #        It could also be usefull to CREATE SUBFIGURES for each function in the case of function
    #        whom output are far from one another.
    #
    #        So far I intend to implement the study of a list of functions given a certain list of 
    #        masks. The defauts parameters need to be the lenght of the union of all the function
    #        parameters. 
    #
    #        I should PUT THE PARAMETERES IN A DICTIONARY and then substitute the values in the 
    #        order of apperence of the parameters in each functions 
    #
    #     x  Find the SET OF ALL THE INDEPENDANT parameters 
    #
    #     x  Create a list of the argument names of each function in the right order
    #
    #     =  Pass effectively the RIGHT ARGUMENTS IN THE RIGHT ORDER TO EACH FUNCTION
    #
    #        Deal with the case that THE ARGUMENTS ARE NOT A LIST
    #
    #        RENAME the variable so that they are more meaningfull
    #
    #        Implement response to the ENTER KEY
    #
    #        Implement AUTOSCALING of the y axis

    def __init__(self, x, func, defaults=[], units=[], 
                 figsize=[16,10], 
                 initial_data=None, 
                 real_imag=False, 
                 lasting=0, 
                 scale=[['linear']*2]*2,
                 plot_params=dict(lw=3, alpha=0.7), 
                 masks=None,
                 fig=[None,None,None]):

        # A list of function arguments
        self.func_args = []
        if isinstance(func, list):
            # The set of all independant parameters for all functions
            self.argnames = set()
            for f in func:
                args = f.func_code.co_varnames[1:f.func_code.co_argcount]
                self.func_args.append(args) 
                self.argnames |= set(args)
            self.argnames = list(self.argnames)
        else:
            # Extract only the arguments passed to the function and exclude the first 
            # one (the x axis).
            self.argnames = func.func_code.co_varnames[1:func.func_code.co_argcount] 
            self.func_args = self.argnames

        self.x = x #TODO be more flexible and enable to set the limits
        self.plot_params = plot_params
        # To improve
        self.fig, self.ax1, self.ax2= fig
        self.real_imag = real_imag
        self.lasting = lasting
        self.initial_data = initial_data
        self.figsize = figsize
        self.func = func
        self.scale = scale
        #self.colors = ['b', 'r', 'g', 'y', 'c', 'm']
        self.colors = [(0,0,1), (1,0,0), (0,0.5,0), (0.75,0.75,0), (0,0.75,0.75), (0.75,0,0.75), (0.5,0.5,0.5)]
        self.color_iterator = 0
        self.params = []
        if units:
            self.units = units
        else:
            self.units = ['']*len(self.argnames)
        if defaults:
            self.defaults = defaults
        else:
            self.defaults = [1]*len(self.argnames)

        self.counter = 0
        self.root = Tk.Tk()
        self._makeform()
        self._init_canvas()
        #self.fig.show()

    def _init_canvas(self):
        if not self.fig:
            self.fig, [self.ax1, self.ax2] = plt.subplots(2, 1, figsize=self.figsize)
            self.fig.subplots_adjust(left=0.075, right=0.85)
            self.ax1.set_xscale(self.scale[0][0])
            self.ax1.set_yscale(self.scale[0][1])
            self.ax1.set_xlabel('Frequency (Hz)')
            self.ax2.set_xscale(self.scale[1][0])
            self.ax2.set_yscale(self.scale[1][1])
            self.ax2.set_xlabel('Frequency (Hz)')
        if self.real_imag:
            self.ax1.set_ylabel(r'Real')
            self.ax2.set_ylabel(r'Imaginary')
        else:
            self.ax1.set_ylabel(r'Amplitude')
            self.ax2.set_ylabel(r'Angle ($^\circ$)')
        self.line1 = []
        self.line2 = []
        self.txt = []
        if self.initial_data!=None:
            if self.real_imag:
                self.ax1.plot(self.x, np.real(self.initial_data), 'k--', **self.plot_params)
                self.ax2.plot(self.x, np.imag(self.initial_data), 'k--', **self.plot_params)
            else:
                self.ax1.plot(self.x, np.abs(self.initial_data), 'k----', **self.plot_params)
                self.ax2.plot(self.x, np.angle(self.initial_data, deg=1), 'k--', **self.plot_params)
        self._plot()

        self.root.wm_title("Test of : "+', '.join([f.func_name for f in self.func]))

        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

        toolbar = NavigationToolbar2TkAgg( self.canvas, self.root )
        toolbar.update()
        self.canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
        #self.canvas.mpl_connect('key_press_event', self._on_key_event)

        #ents = self._makeform(self.root, self.argnames, self.defaults)
        button = Tk.Button(master=self.root, text='Quit', command=self._quit)
        button.pack(side=Tk.BOTTOM)
        Update = Tk.Button(master=self.root, text='Update', command=self._update)
        Update.pack(side=Tk.BOTTOM)
        #root=Tk.Tk()

        Tk.mainloop()
        # If you put root.destroy() here, it will cause an error if
        # the window is closed with the window manager.

    def _makeform(self):
        self.entries = {}
        i = 0
        assert(len(self.argnames)==len(self.defaults))
        for field, default in zip(self.argnames, self.defaults):
            if not i%2:
                row = Tk.Frame(self.root)
            lab = Tk.Label(row, width=22, text=field+": ", anchor='w')
            ent = Tk.Entry(row)
            ent.insert(0,str(default))
            row.pack(side=Tk.TOP, fill=Tk.X, padx=5, pady=5)
            lab.pack(side=Tk.LEFT)
            ent.pack(side=Tk.LEFT, expand=Tk.NO, fill=Tk.X)
            self.entries[field] = ent
            i += 1

    def _evaluate_func(self):
        params = {name:float(self.entries[name].get()) for name in self.argnames}
        return [f(self.x, *[params[a] for a in args]) for f,args in zip(self.func, self.func_args)]

    def _text(self):
        self.params += ['\n'.join([name+' : '+cf.format_scale(self.entries[name].get(), unit) for name,unit in zip(self.argnames, self.units)])]
        if self.fig.legends:
            self.fig.legends.pop(0)
        self.fig.legend(self.line1, self.params, bbox_to_anchor=[0.99, 0.99], prop={'size':14})
        #self.txt += [self.fig.text(0.9, 0.6, params)]

    def _plot(self):
        for y in self._evaluate_func():
            if self.real_imag:
                self.line1 += self.ax1.plot(self.x, np.real(y), **self.plot_params)
                self.line2 += self.ax2.plot(self.x, np.imag(y), **self.plot_params)
            else:
                self.line1 += self.ax1.plot(self.x, np.abs(y), **self.plot_params)
                self.line2 += self.ax2.plot(self.x, np.angle(y, deg=1), **self.plot_params)
        self._text()

    def _scaley(self):
        ya_1, ya_2 = self.ax1.get_ylim()
        yp_1, yp_2 = self.ax2.get_ylim()
        a = np.abs(self.y)
        p = np.angle(self.y, deg=1)
        delta_a = np.max(a)-np.min(a)
        delta_p = np.max(p)-np.min(p)
        delta_ya = ya_2-ya_1
        delta_yp = yp_2-yp_1
        if delta_ya>=0.75*delta_a:
            self.ax1.set_ylim(0.9*np.min(a), 1.1*np.max(a))
        if delta_yp>=0.75*delta_p:
            self.ax2.set_ylim(0.9*np.min(p), 1.1*np.max(p))

    def _update(self):
        #self.canvas.figure.clf()
        if self.lasting:
            for i in range(len(self.func)):
                self.line1[-i+1].set_color(self.colors[self.color_iterator%(len(self.colors)-(1+i))])
                self.line2[-i+1].set_color(self.colors[self.color_iterator%(len(self.colors)-(1+i))])
            self.color_iterator+=1
        else:
            self.ax1.set_prop_cycle(None)
            self.ax2.set_prop_cycle(None)
        if self.counter>=self.lasting:
            self.params.pop(0)
            for i in self.func:
                self.line1.pop(0).remove()
                self.line2.pop(0).remove()
            #self.line1 = self.line1[len(self.func):]
            #self.line2 = self.line2[len(self.func):]
        self._plot()
        self.counter += 1
        #self._scaley()
        self.canvas.show()

    def _on_key_event(event):
        print(('you pressed %s'%event.key))
        key_press_handler(event, canvas, toolbar)


    def _quit(self):
        self.root.quit()     # stops mainloop
        self.root.destroy()  # this is necessary on Windows to prevent
                        # Fatal Python Error: PyEval_RestoreThread: NULL tstate
