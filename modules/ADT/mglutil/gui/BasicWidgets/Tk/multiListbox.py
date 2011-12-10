# $Header: /opt/cvs/python/packages/share1.5/mglutil/gui/BasicWidgets/Tk/multiListbox.py,v 1.5 2008/10/16 22:09:26 vareille Exp $
# $Id: multiListbox.py,v 1.5 2008/10/16 22:09:26 vareille Exp $
# Source: http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/52266
# Modifed by: Sargis Dallakyan sargis@scripps.edu
from Tkinter import *
import Pmw
import tkFont

from mglutil.util.misc import ensureFontCase

class MultiListbox(Pmw.ScrolledFrame):
    """This is a compound widget that gangs multiple Tk Listboxes to a single
       scrollbar to achieve a simple multi-column scrolled listbox. Most of
       the Listbox API is mirrored to make it act like the normal Listbox
       but with multiple values per row.
    """
    def __init__(self, master, lists, **kw):
        Pmw.ScrolledFrame.__init__(self, master, horizflex='expand',vertflex='elastic',
                                    vscrollmode='none', **kw)
        self.lists = []
        self.colmapping = {}
        self.origData = None
        self.rows = None
        self.myFont = tkFont.Font(font = (ensureFontCase('helvetica'), 11, "bold"))
        frame_main = Frame(self.interior())
        frame_main.pack(side=LEFT, expand=YES, fill=BOTH)
        frame_sb = Frame(self.interior())
        frame_sb.pack(side=LEFT, fill=Y)
        sb = Scrollbar(frame_sb, orient=VERTICAL, command=self._scroll)
        m = PanedWindow(frame_main,bd=0,handlepad=0,handlesize=0,sashpad=0)
        m.pack(fill=BOTH, expand=1)
        self.PanedWindow = m
        for l,w in lists:
            frame = Frame(m,bd=0); frame.pack(side=LEFT, expand=YES, fill=BOTH)
            b = Button(frame, text=unicode(l),  relief=RAISED, font=self.myFont,bd=1)
            b.pack(fill = X)
            b.bind('<Button-1>', self._sort)
            self.colmapping[b]=(len(self.lists),1)
            lb = Listbox(frame, width=w, bd=0,
                         exportselection=FALSE,bg='white')
            lb.pack(expand=YES, fill=BOTH)
            self.lists.append(lb)
            lb.bind('<B1-Motion>', lambda e, s=self: s._select(e.y))
            lb.bind('<Button-1>', lambda e, s=self: s._select(e.y))
            lb.bind('<Double-Button-1>', lambda e, s=self: s.output(e.y))
            lb.bind('<Leave>', lambda e: 'break')
            lb.bind('<B2-Motion>', lambda e, s=self: s._b2motion(e.x, e.y))
            lb.bind('<Button-2>', lambda e, s = self: s._button2(e.x, e.y))
#            import pdb
#            pdb.set_trace()
            m.add(frame, width=lb.winfo_reqwidth())    
        sb.pack(expand=YES, fill=Y)
        self.lists[0]['yscrollcommand']=sb.set        
        self.interior().pack(expand=YES, fill=BOTH)
        
    def _sort(self, e):
        # get the listbox to sort by (mapped by the header button)
        b=e.widget
        col, direction = self.colmapping[b]

        # get the entire table data into mem
        tableData = self.get(0,END)
        if self.origData == None:
            import copy
            self.origData = copy.deepcopy(tableData)

        rowcount = len(tableData)

        #remove old sort indicators if it exists
        for btn in self.colmapping.keys():
            lab = btn.cget('text')
            if lab[-1] == u"\u25BC" or lab[-1] ==u"\u25B2":
                btn.config(text=unicode(lab[:-2]))

        btnLabel = b.cget('text')
        #sort data based on direction
        if direction==0:
            tableData = self.origData
        else:
            if direction==1: b.config(text = unicode(btnLabel + " " + u"\u25B2"), font=self.myFont)
            else: b.config(text=unicode(btnLabel + " " + u"\u25BC"),font=self.myFont)
            # sort by col
            def colsort(x, y, mycol=col, direction=direction):
                return direction*cmp(x[mycol], y[mycol])

            tableData.sort(colsort)

        #clear widget
        self.delete(0,END)

        # refill widget
        for row in range(rowcount):
            self.insert(END, tableData[row])

        # toggle direction flag 
        if(direction==1): direction=-1
        else: direction += 1
        self.colmapping[b] = (col, direction)

    def _select(self, y):
         row = self.lists[0].nearest(y)
         self.selection_clear(0, END)
         self.selection_set(row)
         return 'break'

    def _button2(self, x, y):
         for l in self.lists: l.scan_mark(x, y)
         return 'break'

    def _b2motion(self, x, y):
         for l in self.lists: l.scan_dragto(x, y)
         return 'break'

    def _scroll(self, *args):
         for l in self.lists:
             apply(l.yview, args)

    def curselection(self):
         return self.lists[0].curselection()

    def delete(self, first, last=None):
         for l in self.lists:
             l.delete(first, last)

    def get(self, first, last=None):
         result = []
         for l in self.lists:
             result.append(l.get(first,last))
         if last: return apply(map, [None] + result)
         return result
             
    def index(self, index):
         self.lists[0].index(index)

    def insert(self, index, *elements):
         for e in elements:
             i = 0
             for l in self.lists:
                  l.insert(index, e[i])
                  i = i + 1

    def size(self):
         return self.lists[0].size()

    def see(self, index):
         for l in self.lists:
             l.see(index)

    def selection_anchor(self, index):
         for l in self.lists:
             l.selection_anchor(index)

    def selection_clear(self, first, last=None):
         for l in self.lists:
             l.selection_clear(first, last)

    def selection_includes(self, index):
         return self.lists[0].selection_includes(index)

    def selection_set(self, first, last=None):
         for l in self.lists:
             l.selection_set(first, last)   
             
    def output(self, y):
        row = self.lists[0].nearest(y)
        if self.rows:
            print self.rows[row][1]
        
if __name__ == '__main__':
    tk = Tk()
    mlb = MultiListbox(tk, (('Subject', 40), ('Sender', 20), ('Date', 10)))
    mlb.pack(expand=YES,fill=BOTH)
    for i in range(10):
        mlb.insert(END, ('Important Message: %d' % i, 'John Doe', '10/10/%04d' % (1900+i)))
    
    mlb.reposition()

    tk.mainloop()