#!/usr/bin/python
# -*- coding: utf-8 -*-

#Author : Anthoni Manseau

import numpy as np
import corefunctions as cf
from matplotlib.cbook import silent_list
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import sys

if sys.version_info[0] < 3:
    import Tkinter as tk
else:
    import tkinter as tk
    from tkinter import ttk

class Graphic(tk.Frame):
    """ Creates a frame containing the canvas of a matplotlib figure.
    """
    def __init__(self, master, fig, **kwargs):
        """ Initiate frame containing a canvas that plots x vs y """
        tk.Frame.__init__(self, master)
        self._verbose = kwargs.get('verbose', False)
        self.fig = fig
        children = fig.get_children()
        if self.fig.legends!=[]:
            self.legendLabels = [i.get_text() for i in [i for i in children if 'Legend' in str(type(i))][0].get_texts()]
        self.subplots = [i for i in children if 'AxesSubplot' in str(type(i))]
        self.lines = [s.get_lines() for s in self.subplots]
        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand = True)
        toolbar = NavigationToolbar2TkAgg(self.canvas, self)
        toolbar.update()
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand = True)
    def verbose(self, message):
        if self._verbose:
            print(message)
    def update(self):
        self.canvas.draw()
    def axe_lim(self, xmin=None, xmax=None, ymin=None, ymax=None):
        self.ax.set_xlim((xmin, xmax))
        self.ax.set_ylim((ymin, ymax))
    def legend(self, lbl):
        """ Set legend """
        #TODO set bbox_to_anchor option.
        if self.fig.legends!=[]:
            self.fig.legends = []
            self.fig.legend(self.lines[0], self.legendLabels)#, bbox_to_anchor=self.bbox_to_anchor)
    def intensity(self, i, alpha):
        """ Change the intensity of the curve 'i'.
            i : index of the curve.
            alpha : intensity ('alpha' argument for 'plot') of the curve.
        """
        self.verbose('setting alpha : %s'%alpha)
        if isinstance(self.lines, silent_list):
            self.lines[i].set_alpha(alpha)
        else:
            for lines in self.lines:
                lines[i].set_alpha(alpha)

class Legend(tk.Frame):
    def __init__(self, master, subFrame, lbl, **kwargs):
        """ Initial class that creates legend check buttons """
        tk.Frame.__init__(self,  subFrame, relief=tk.GROOVE, bd=1)
        self.master = master
        self.subFrame = subFrame
        self.lbl = lbl
        self._update = kwargs.get('update', self.master.update_curves)
        self._verbose = kwargs.get('verbose', False)
        self.legendCheckButtons = []
        self.legendVars = []

        self.state=[1]*len(self.lbl)
        self.index=range(len(self.lbl))
        self._singlecheck = 0

        def bindingfunction(event):
            canvas.configure(scrollregion=canvas.bbox("all"), width=100, height=900)

        self.frame1 = tk.Frame(self)
        self.frame1.pack()
        self.frame2 = tk.Frame(self)
        canvas=tk.Canvas(self.frame2)
        self.frame=tk.Frame(canvas)
        scrollbar=tk.Scrollbar(self.frame2, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=scrollbar.set)
        scrollbar.pack(side='right', fill='y')
        canvas.pack(side="bottom")
        canvas.create_window((0,0), window=self.frame, anchor='nw')
        self.frame.bind("<Configure>", bindingfunction)
        self.create_widgets(self.lbl)
        self.frame2.pack()

    def verbose(self, message):
        if self._verbose:
            print(message)

    def select(self, i):
        """ Select the ith checkbutton """
        self.legendCheckButtons[i].select()

    def deselect(self, i):
        """ Deselect the ith checkbutton """
        self.legendCheckButtons[i].deselect()

    def read_state(self, i):
        """ Reads the state of the ith checkbutton. """
        self.verbose("read button's state %s : %s"%(i, self.legendVars[i].get()))
        return self.legendVars[i].get()

    def set_state(self, i, newstate):
        """ Sets the state of the ith checkbutton. """
        self.verbose("setting state %s to %s : "%(i,newstate))
        if newstate != self.read_state(i):
            #self.switch_state(i)
            if newstate:
                self.index.append(i)
                self.select(i)
            else:
                self.index.remove(i)
                self.deselect(i)

    def update(self, i):
        """ Updates curves with the checkbuttons. Is called when a checkbutton is used."""
        #self.switch_state(i)
        if self.read_state(i):
            self.index.append(i)
            #self.select(i)
        else:
            self.index.remove(i)
        if self.read_single_mode():
            self.set_single_mode(1)
        self._update()

    def switch_state(self, i):
        """ Switch the state of the ith checkbutton. """
        # NOTE : the ^1 is there because this function is called
        #       whenever the ith checkbutton is clicked. So the
        #       state should be switched. The reason of this is
        #       that the ith variable associated with the check
        #       button does not update itself after a click.
        newstate = self.read_state(i)^1
        print("switching state %s to %s"%(i, newstate))
        print(self.read_state(i))
        if newstate:
            self.index.append(i)
            #self.select(i)
        else:
            self.index.remove(i)
            #self.deselect(i)
        #self.state[i]=newstate

    def checkall(self):
        """ Activates all legend curves """
        for i in range(len(self.state)):
            #if self.state[i]==0:
            if self.read_state(i)==0:
                self.set_state(i, 1)
        self.set_single_mode(0)
        self._update()

    def uncheckall(self):
        """ Desactivate all legend curves """
        for i in self.index[:]:
            self.set_state(i, 0)
        self.set_single_mode(0)
        self._update()

    def set_single_mode(self, i):
        """ Sets in only one curve on. The last of the index list."""
        if i:
            if len(self.index)>0:
                for j in self.index[:-1]:
                    self.set_state(j, 0)
            if not self.read_single_mode():
                self.cb_single.select()
        else:
            if self.read_single_mode():
                self.switch_single_mode()
                self.cb_single.deselect()

    def read_single_mode(self):
        """ Reads the status of the single mode checkbutton.
        """
        return self._singlecheck
        #a = self.subframe1
        #return a.SingleMode.get()

    def switch_single_mode(self):
        """ Switches the status of the single mode checkbutton
        """
        self._singlecheck ^= 1
        return self._singlecheck

    def singlemode(self):
        """ Activate the single curve mode  """
        self.switch_single_mode()
        if self.read_single_mode():
            self.set_single_mode(1)
        self._update()

    def create_widgets(self, lbl):
        """ Create checkbutton for each element of the legend (i.e. each curve of the plot).
            lbl : a list of label for each curve in the legend (should be strings).
        """
        self.cb_check = tk.Button(self.frame1,
                    text='CheckAll',
                    command=self.checkall
                    )
        self.cb_check.pack(anchor=tk.NW)
        self.cb_uncheck = tk.Button(self.frame1,
                    text='UncheckAll',
                    command=self.uncheckall
                    )
        self.cb_uncheck.pack(anchor=tk.NW)
        self.SingleMode = tk.BooleanVar()
        self.cb_single = tk.Checkbutton(self.frame1,
                    text='SingleMode',
                    variable=self.SingleMode,
                    command=self.singlemode
                    )
        self.cb_single.pack(anchor=tk.NW)

        for i in range(len(lbl)) :
            self.legendVars.append(tk.IntVar())
            func = lambda n=i : self.update(n)
            self.legendCheckButtons.append(tk.Checkbutton(self.frame,
                        text=lbl[i],
                        variable=self.legendVars[i],
                        command= func
                        ))
            self.legendCheckButtons[i].pack(anchor=tk.W)
            self.legendCheckButtons[i].select()

class Axis(tk.Frame):
    def __init__(self, root, master, **kwargs):
        """ Initial class that creates legend check buttons """
        self.master = master
        self.root = root
        tk.Frame.__init__(self,  master)
        self.create_widgets()

    def create_widgets(self):
        """ Create checkbutton for each element of the legend (i.e. each curve of the plot).
            lbl : a list of label for each curve in the legend (should be strings).
        """
        self.lbl_xmin = tk.Label(self.master, text="xmin : ")
        self.lbl_xmin.pack(anchor=tk.W)
        self.e_xmin = tk.Entry(self.master)
        self.e_xmin.pack()

        self.lbl_xmax = tk.Label(self.master, text="xmax : ")
        self.lbl_xmax.pack(anchor=tk.W)
        self.e_xmax = tk.Entry(self.master)
        self.e_xmax.pack()

        self.lbl_ymin = tk.Label(self.master, text="ymin : ")
        self.lbl_ymin.pack(anchor=tk.W)
        self.e_ymin = tk.Entry(self.master)
        self.e_ymin.pack()

        self.lbl_ymax = tk.Label(self.master, text="ymax : ")
        self.lbl_ymax.pack(anchor=tk.W)
        self.e_ymax = tk.Entry(self.master)
        self.e_ymax.pack()

        self.button = tk.Button(self.master, text="Enter", command=self.update_lim)
        self.button.pack()

    def update_lim(self):
        """ Update the axis limit. """
        def readout(self):
            """ Converts x from string to float.
                If x is an empty string, return None
            """
            x = self.get()
            self.delete(0, 'end')
            output =  cf.float_or_string(x)
            if isinstance(output, str):
                return None
            else:
                return output
        xmin = readout(self.e_xmin)
        xmax = readout(self.e_xmax)
        ymin = readout(self.e_ymin)
        ymax = readout(self.e_ymax)
        self.root.Graphic.axe_lim(xmin, xmax, ymin, ymax)
        self.refresh()

    def refresh(self):
        """ Refreshes the canvas. """
        self.root.Graphic.update()

class DynamicPlot(tk.Tk, object):
    """ GUI to manipulate display of Graphics, select range of axis and eventually more...
        arguments :
            fig : Figure object.
            lbl : A list a labels for each dimension of the data set.
                  If you visualise the mutlidimensional array as a tree,
                  the lbl array should start with the last dimension of the
                  tree and the last element should be the first dimension
                  of the tree.
            lowAlpha : intensity of the color when a curve is disabled (default=0.05).
            highAlpha : intensity of the color when a curve is enabled (default=1).
    """
    def __init__(self, fig, lbl, **kwargs):
        """ Defines data (x and y), legend label (lbl) and initialise frames """
        self.lowalpha = kwargs.get('lowAlpha', 0.01)
        self.highalpha = kwargs.get('highAlpha', 1)
        self._verbose = kwargs.get('verbose', False)
        self.lbl = lbl
        tk.Tk.__init__(self)
        self.Graphic = Graphic(self, fig,  **kwargs)
        self.Graphic.pack(side='left', expand=True, fill=tk.BOTH)
        self.widget_frame = tk.Frame(self)
        self.widget_frame.pack(expand=True, fill=tk.BOTH)#side="left")
        self.a = []
        assert(type(self.lbl[0])==np.ndarray or type(self.lbl[0])==list)
        if len(self.lbl[0])>1:
            self.shape = [len(l) for l in self.lbl]
            for l in self.lbl:
                self.a.append(Legend(self, self.widget_frame, l, **kwargs))
                self.a[-1].pack(side='right', anchor=tk.NW)
            self.legend_lst = self.make_legend()
        else:
            self.shape = []
            frame = tk.Frame(self)
            self.a.append(Legend(self, frame, self.lbl, **kwargs))
            self.a[-1].pack(anchor=tk.W)
            self.legend_lst = self.lbl
            frame.pack(side=tk.RIGHT, anchor=tk.NW)
        #self.subframe2 = Axis(self, self.widget_frame, **kwargs)
        #self.subframe2.pack(anchor=tk.W)

        self.legend = self.legend_lst
        self.index = range(len(self.legend_lst))

    #NOTE  : Because it seams like the state does not change when check button is
    #        clicked, I will update the state of the checkbutton manually using
    #        the index of the checkbutton that has called the update function.

    def verbose(self, message):
        if self._verbose:
            print(message)

    def make_legend(self):
        """ Generates the legend from the label of each parameter list. """
        lgd = self.lbl[0]
        for i in self.lbl[1:]:
            lgd = np.array([[k+'/'+l for k in lgd] for l in i]).flatten()
        return lgd

    def update_legend(self):
        """ Update the legend. """
        self.legend = [self.legend_lst[i] for i in self.index]
        self.verbose("Legend size : %s"%len(self.legend))
        self.Graphic.legend(self.legend)

    def update_curves(self):
        """ Update the curve selected by the legend checkbuttons.
            'a' is a list of legend frames.
        """
        a = self.a
        index = a[0].index
        # Generate the indexes from each objects in a.
        # NOTE : when len(a)=1 the loop is not called !
        for i in range(1, len(a)):
            index = np.array([[k+np.prod(self.shape[:i])*l for k in index] for l in a[i].index]).flatten()
        # Curves to turn on and off. Update it on the Graphic.
        off =  set(self.index)-set(index)
        on =  set(index)-set(self.index)
        self.verbose('self.index : %s'%self.index)
        self.verbose('index : %s'%index)
        self.verbose('on : %s'%on+' off : %s'%off)
        for i in off:
            self.transparency(i,0)
        for i in on:
            self.transparency(i,1)
        # Update the new indexes.
        self.index = index[:]
        # Update the legend and refresh.
        self.refresh()

    def refresh(self):
        """ Refreshes the canvas. """
        self.update_legend()
        self.Graphic.update()

    def transparency(self, i, state):
        """ Changes the transparency of the ith curve (high of low) """
        a = self.Graphic
        if state:
            a.intensity(i, self.highalpha)
        else:
            a.intensity(i, self.lowalpha)
