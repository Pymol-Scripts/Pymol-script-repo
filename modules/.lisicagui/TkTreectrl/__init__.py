# -*- coding: utf-8 -*-
#v.0.7
#added LEFT,RIGHT,HEADER and CONTENT, fixed typo with LEFTMOST
#v.0.8: added __version__

'''Main module of the TkTreectrl package.
Once the TkTreectrl package is installed, it is safe to do :

    from TkTreectrl import *

This will add a number of constants and the following widget classes to the
current namespace:

Treectrl, MultiListbox, ScrolledTreectrl, ScrolledMultiListbox, ScrolledWidget
'''
__version__ = '2.0'

ABOVE = 'above'
ACTIVE = 'active'
ALL = 'all'
ASCII = 'ascii'
AREA = 'area'
BELOW = 'below'
BITMAP  = 'bitmap'
BORDER = 'border'
CANVAS = 'canvas'
COLUMN = 'column'
CONTENT = 'content'
DECREASING = 'decreasing'
DICTIONARY = 'dictionary'
DOUBLE = 'double'
DOT = 'dot'
DYNAMIC = 'dynamic'
ENABLED = 'enabled'
FIRST = 'first'
FIRSTCHILD = 'firstchild'
FOCUS = 'focus'
HEADER = 'header'
IMAGE = 'image'
INCREASING = 'increasing'
INTEGER = 'integer'
ITEM = 'item'
LAST = 'last'
LASTCHILD = 'lastchild'
LEFT = 'left'
LEFTMOST = 'leftmost'
LONG = 'long'
NEXT = 'next'
NEXTSIBLING = 'nextsibling'
NORMAL = 'normal'
OPEN = 'open'
PARENT = 'parent'
PRESSED = 'pressed'
PREV = 'prev'
PREVSIBLING = 'prevsibling'
REAL = 'real'
RECT = 'rect'
RIGHT = 'right'
RIGHTMOST = 'rightmost'
ROOT = 'root'
SELECT = 'select'
SELECTED = 'selected'
STATIC = 'static'
STRING = 'string'
TAIL = 'tail'
TEXT = 'text'
TIME = 'time'
TREE = 'tree'
WINDOW = 'window'

from TkTreectrl.Treectrl import Treectrl
from TkTreectrl.MultiListbox import MultiListbox
from TkTreectrl.ScrolledTreectrl import ScrolledTreectrl
from TkTreectrl.ScrolledTreectrl import ScrolledMultiListbox
# put the ScrolledWidget class in the global namespace so people can easily
# create custom subclasses
from TkTreectrl.ScrolledTreectrl import ScrolledWidget
