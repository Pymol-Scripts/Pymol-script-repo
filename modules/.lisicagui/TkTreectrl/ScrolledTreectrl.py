# -*- coding: utf-8 -*-
'''Treectrl and MultiListbox widgets with scrollbars.'''

try:
    import tkinter
except ImportError:
    import Tkinter as tkinter

from TkTreectrl.Treectrl import Treectrl
from TkTreectrl.MultiListbox import MultiListbox

# see if we have ttk available
try:
    # typically the Python2 variant:
    import ttk
    _have_ttk = True
except ImportError:
    try:
        from tkinter import ttk
        _have_ttk = True
    except ImportError:
        _have_ttk = False

# use themed widgets, if available
if _have_ttk:
    _Frame = ttk.Frame
    _Scrollbar = ttk.Scrollbar
else:
    _Frame = tkinter.Frame
    _Scrollbar = tkinter.Scrollbar

class _UniScrollbar(_Scrollbar):
    '''Internal class.
    A "universal" Scrollbar widget that uses a ttk.Scrollbar if ttk is
    available, else a tkinter.Scrollbar. If a ttk.Scrollbar is used,
    ttk-incompatible configuration options are silently ignored, to keep
    backwards compatibility intact.'''
    def __init__(self, master=None, **kw):
        if _Scrollbar == ttk.Scrollbar:
            # some opts that aren't legal in ttk, ignore them silently:
            for opt in ('activebackground', 'activerelief', 'background',
                        'bg', 'bd', 'borderwidth', 'elementborderwidth',
                        'highlightbackground', 'highlightcolor',
                        'highlightthickness', 'jump', 'relief', 'repeatdelay',
                        'repeatinterval', 'troughcolor', 'width'):
                if opt in kw:
                    del(kw[opt])
        _Scrollbar.__init__(self, master, **kw)

    def cget(self, key):
        if _Scrollbar == ttk.Scrollbar:
            if key in ('activebackground', 'activerelief', 'background',
                       'bg', 'bd', 'borderwidth', 'elementborderwidth',
                       'highlightbackground', 'highlightcolor',
                       'highlightthickness', 'jump', 'relief',
                       'repeatdelay', 'repeatinterval', 'troughcolor',
                       'width'):
                return None
        return _Scrollbar.cget(self, key)
    __getitem__ = cget

    def configure(self, cnf=None, **kw):
        if _Scrollbar == ttk.Scrollbar:
            for opt in ('activebackground', 'activerelief', 'background',
                        'bg', 'bd', 'borderwidth', 'elementborderwidth',
                        'highlightbackground', 'highlightcolor',
                        'highlightthickness', 'jump', 'relief', 'repeatdelay',
                        'repeatinterval', 'troughcolor', 'width'):
                if not cnf is None and opt in cnf:
                    del(cnf[opt])
                if opt in kw:
                    del(kw[opt])
        return _Scrollbar.configure(self, cnf, **kw)
    config = configure

###############################################################################
###############################################################################

class ScrolledWidget(_Frame):
    '''Base class for Tkinter widgets with scrollbars.
    The widget is a standard Tkinter.Frame or a ttk.Frame if ttk is available,
    with an additional configuration option SCROLLMODE which may be one of
    "x", "y", "both" or "auto".
    If SCROLLMODE is one of "x", "y" or "both", one or two static scrollbars
    will be drawn. If SCROLLMODE is set to "auto", two automatic scrollbars
    that appear only if needed will be drawn.
    In order to ensure backwards compatibility, ttk-incompatible configuration
    options passed to the widget or to the scrollbars will be silently ignored
    if ttk widgets are used.
    The Scrollbar widgets can be accessed with the hbar and vbar class
    attributes. Derived classes must override the _setScrolledWidget() method,
    which must return the widget that will be scrolled and should add a class
    attribute that allows to access this widget, so the _setScrolledWidget()
    method for a ScrolledListbox widget might look like:

        def _setScrolledWidget(self):
            self.listbox = Tkinter.Listbox(self)
            return self.listbox

    Note that although it should be possible to create scrolled widget classes
    for virtually any Listbox or Canvas alike Tkinter widget you can *not*
    safely use this class to add automatic scrollbars to a Text or Text alike
    widget. This is because in a scrolled Text widget the value of the
    horizontal scrollbar depends only on the visible part of the Text, not on
    it's whole contents. Thus it may happen that it is the last visible line of
    text that causes the automatic scrollbar to be mapped which then hides this
    last line so it will be unmapped again, but then it is requested again and
    gets mapped and so on forever. There are strategies to avoid this, but
    usually at the cost that there will be situations where the horizontal
    scrollbar remains mapped although it is actually not needed. In order to
    acomplish this with the ScrolledWidget class, at least the _scrollXNow()
    and _scrollBothNow() methods must be overridden with appropriate handlers.
    '''
    def __init__(self, master=None, **kw):
        if 'scrollmode' in kw:
            scrollmode = kw['scrollmode']
            del(kw['scrollmode'])
        else:
            scrollmode = 'auto'
        if not 'width' in kw:
            kw['width'] = 400
        if not 'height' in kw:
            kw['height'] = 300
        if _Frame == ttk.Frame:
            # "bd" is not valid in ttk
            if 'bd' in kw:
                kw['borderwidth'] = kw['bd']
                del(kw['bd'])
            # some more opts that aren't legal in ttk, ignore them silently:
            for opt in ('background', 'bg', 'colormap', 'container',
                        'highlightbackground', 'highlightcolor',
                        'highlightthickness', 'padx', 'pady', 'visual'):
                if opt in kw:
                    del(kw[opt])
        _Frame.__init__(self, master, **kw)
        # call grid_propagate(0) so the widget will not change its size
        # when the scrollbars are mapped or unmapped; that is why we need
        # to specify width and height in any case
        self.grid_propagate(0)
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)

        self._scrollmode = scrollmode

        self._scrolledWidget = self._setScrolledWidget()
        self._scrolledWidget.grid(row=0, column=0, sticky='news')

        self.vbar = _UniScrollbar(self, orient='vertical',
                                            command=self._scrolledWidget.yview)
        self.vbar.grid(row=0, column=1, sticky='ns')
        self.hbar = _UniScrollbar(self, orient='horizontal',
                                            command=self._scrolledWidget.xview)
        self.hbar.grid(row=1, column=0, sticky='ew')

        # Initialise instance variables.
        self._hbarOn = 1
        self._vbarOn = 1
        self._scrollTimer = None
        self._scrollRecurse = 0
        self._hbarNeeded = 0
        self._vbarNeeded = 0

        self._scrollMode(scrollmode)

    def _setScrolledWidget(self):
        '''This method must be overridden in derived classes.
        It must return the widget that should be scrolled and should
        add a reference to the ScrolledWidget object, so it can be accessed
        by the user. For example, to create a scrolled Listbox, do:
        self.listbox = Tkinter.Listbox(self)
        return self.listbox'''
        pass

    def configure(self, cnf=None, **kw):
        # hack: add "scrollmode" to user configurable options
        if not cnf is None and 'scrollmode' in cnf:
            self._scrollMode(cnf['scrollmode'])
            del(cnf['scrollmode'])
        if 'scrollmode' in kw:
            self._scrollMode(kw['scrollmode'])
            del(kw['scrollmode'])
        if _Frame == ttk.Frame:
            # "bd" is not valid in ttk
            if not cnf is None and 'bd' in cnf:
                cnf['borderwidth'] = cnf['bd']
                del(cnf['bd'])
            if 'bd' in kw:
                kw['borderwidth'] = kw['bd']
                del(kw['bd'])
            # some more opts that aren't legal in ttk,
            # ignore them silently:
            for opt in ('bg', 'background', 'highlightbackground',
                        'highlightcolor', 'highlightthickness',
                        'padx', 'pady'):
                if not cnf is None and opt in cnf:
                    del(cnf[opt])
                if opt in kw:
                    del(kw[opt])
        return _Frame.configure(self, cnf, **kw)
    config = configure

    def cget(self, key):
        if key == 'scrollmode':
            return self._scrollmode
        if _Frame == ttk.Frame:
            if key == 'bd':
                key = 'borderwidth'
            elif key in  ('bg', 'background', 'colormap', 'container',
                          'highlightbackground', 'highlightcolor',
                          'highlightthickness', 'padx', 'pady', 'visual'):
                return None
        return _Frame.cget(self, key)
    __getitem__ = cget

    def keys(self):
        keys = _Frame.keys(self) + ['scrollmode']
        keys.sort()
        return keys

    # methods to control the scrollbars;
    # these are mainly stolen from Pmw.ScrolledListbox

    def _scrollMode(self, mode):
        if mode == 'both':
            if not self._hbarOn:
                self._toggleHbar()
            if not self._vbarOn:
                self._toggleVbar()
        elif mode == 'auto':
            if self._hbarNeeded != self._hbarOn:
                self._toggleHbar()
            if self._vbarNeeded != self._vbarOn:
                self._toggleVbar()
        elif mode == 'x':
            if not self._hbarOn:
                self._toggleHbar()
            if self._vbarOn:
                self._toggleVbar()
        elif mode == 'y':
            if not self._vbarOn:
                self._toggleVbar()
            if self._hbarOn:
                self._toggleHbar()
        else:
            raise ValueError('bad value for option scrollmode: ' +\
                    '"%s"; should be "x", "y", "both" or "auto".' % mode )
        self._scrollmode = mode
        self._configureScrollCommands()

    def _configureScrollCommands(self):
        # Clean up previous scroll commands to prevent memory leak.
        tclCommandName = str(self._scrolledWidget.cget('xscrollcommand'))
        if tclCommandName != '':
            self._scrolledWidget.deletecommand(tclCommandName)
        tclCommandName = str(self._scrolledWidget.cget('yscrollcommand'))
        if tclCommandName != '':
            self._scrolledWidget.deletecommand(tclCommandName)
        # If both scrollmodes are not dynamic we can save a lot of
        # time by not having to create an idle job to handle the
        # scroll commands.
        if self._scrollmode == 'auto':
            self._scrolledWidget.configure(
                                    xscrollcommand=self._scrollBothLater,
                                    yscrollcommand=self._scrollBothLater)
        else:
            self._scrolledWidget.configure(xscrollcommand=self._scrollXNow,
                                    yscrollcommand=self._scrollYNow)

    def _scrollXNow(self, first, last):
        first, last = str(first), str(last)
        self.hbar.set(first, last)
        self._hbarNeeded = ((first, last) not in (('0', '1'), ('0.0', '1.0')))
        if self._scrollmode == 'auto':
            if self._hbarNeeded != self._hbarOn:
                self._toggleHbar()

    def _scrollYNow(self, first, last):
        first, last = str(first), str(last)
        self.vbar.set(first, last)
        self._vbarNeeded = ((first, last) not in (('0', '1'), ('0.0', '1.0')))
        if self._scrollmode == 'auto':
            if self._vbarNeeded != self._vbarOn:
                self._toggleVbar()

    def _scrollBothLater(self, first, last):
        # Called by the listbox to set the horizontal or vertical
        # scrollbar when it has scrolled or changed size or contents.
        if self._scrollTimer is None:
            self._scrollTimer = self.after_idle(self._scrollBothNow)

    def _scrollBothNow(self):
        # This performs the function of _scrollXNow and _scrollYNow.
        # If one is changed, the other should be updated to match.
        self._scrollTimer = None
        # Call update_idletasks to make sure that the containing frame
        # has been resized before we attempt to set the scrollbars.
        # Otherwise the scrollbars may be mapped/unmapped continuously.
        self._scrollRecurse = self._scrollRecurse + 1
        self.update_idletasks()
        self._scrollRecurse = self._scrollRecurse - 1
        if self._scrollRecurse != 0:
            return

        xview, yview = self._scrolledWidget.xview(),\
                                    self._scrolledWidget.yview()
        self.hbar.set(*xview)
        self.vbar.set(*yview)
        self._hbarNeeded = (xview != (0.0, 1.0))
        self._vbarNeeded = (yview != (0.0, 1.0))

        # If both horizontal and vertical scrollmodes are dynamic and
        # currently only one scrollbar is mapped and both should be
        # toggled, then unmap the mapped scrollbar.  This prevents a
        # continuous mapping and unmapping of the scrollbars.
        if (self._scrollmode == 'auto' and self._hbarNeeded != self._hbarOn and
                                    self._vbarNeeded != self._vbarOn and
                                    self._vbarOn != self._hbarOn):
            if self._hbarOn:
                self._toggleHbar()
            else:
                self._toggleVbar()
            return

        if self._scrollmode == 'auto':
            if self._hbarNeeded != self._hbarOn:
                self._toggleHbar()
            if self._vbarNeeded != self._vbarOn:
                self._toggleVbar()

    def _toggleHbar(self):
        self._hbarOn = not self._hbarOn
        if self._hbarOn:
            self.hbar.grid(row=1, column=0, sticky='ew')
        else:
            self.hbar.grid_forget()

    def _toggleVbar(self):
        self._vbarOn = not self._vbarOn
        if self._vbarOn:
            self.vbar.grid(row=0, column=1, sticky='ns')
        else:
            self.vbar.grid_forget()

    def destroy(self):
        if self._scrollTimer is not None:
            self.after_cancel(self._scrollTimer)
            self._scrollTimer = None
        _Frame.destroy(self)

###############################################################################
###############################################################################

class ScrolledTreectrl(ScrolledWidget):
    '''Treectrl widget with one or two static or automatic scrollbars.
    Subwidgets are:
        treectrl - TkTreectrl.Treectrl widget
        hbar - horizontal Tkinter.Scrollbar or ttk.Scrollbar
        vbar - vertical Tkinter.Scrollbar or ttk.Scrollbar
    The widget itself is a Tkinter.Frame or ttk.Frame with one additional
    configuration option:
        scrollmode - may be one of "x", "y", "both" or "auto".'''
    def __init__(self, *args, **kw):
        ScrolledWidget.__init__(self, *args, **kw)
    def _setScrolledWidget(self):
        self.treectrl = Treectrl(self)
        return self.treectrl

###############################################################################
###############################################################################

class ScrolledMultiListbox(ScrolledWidget):
    '''MultiListbox widget with one or two static or automatic scrollbars.
    Subwidgets are:
        listbox - TkTreectrl.MultiListbox widget
        hbar - horizontal Tkinter.Scrollbar or ttk.Scrollbar
        vbar - vertical Tkinter.Scrollbar or ttk.Scrollbar
    The widget itself is a Tkinter.Frame or ttk.Frame with one additional
    configuration option:
        scrollmode - may be one of "x", "y", "both" or "auto".'''
    def __init__(self, *args, **kw):
        ScrolledWidget.__init__(self, *args, **kw)
    def _setScrolledWidget(self):
        self.listbox = MultiListbox(self)
        return self.listbox
