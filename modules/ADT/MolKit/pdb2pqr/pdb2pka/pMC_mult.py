# This file was created automatically by SWIG 1.3.29.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

import _pMC_mult
import new
new_instancemethod = new.instancemethod
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'PySwigObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static) or hasattr(self,name):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError,name

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

import types
try:
    _object = types.ObjectType
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0
del types


class PySwigIterator(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, PySwigIterator, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, PySwigIterator, name)
    def __init__(self): raise AttributeError, "No constructor defined"
    __repr__ = _swig_repr
    __swig_destroy__ = _pMC_mult.delete_PySwigIterator
    __del__ = lambda self : None;
    def value(*args): return _pMC_mult.PySwigIterator_value(*args)
    def incr(*args): return _pMC_mult.PySwigIterator_incr(*args)
    def decr(*args): return _pMC_mult.PySwigIterator_decr(*args)
    def distance(*args): return _pMC_mult.PySwigIterator_distance(*args)
    def equal(*args): return _pMC_mult.PySwigIterator_equal(*args)
    def copy(*args): return _pMC_mult.PySwigIterator_copy(*args)
    def next(*args): return _pMC_mult.PySwigIterator_next(*args)
    def previous(*args): return _pMC_mult.PySwigIterator_previous(*args)
    def advance(*args): return _pMC_mult.PySwigIterator_advance(*args)
    def __eq__(*args): return _pMC_mult.PySwigIterator___eq__(*args)
    def __ne__(*args): return _pMC_mult.PySwigIterator___ne__(*args)
    def __iadd__(*args): return _pMC_mult.PySwigIterator___iadd__(*args)
    def __isub__(*args): return _pMC_mult.PySwigIterator___isub__(*args)
    def __add__(*args): return _pMC_mult.PySwigIterator___add__(*args)
    def __sub__(*args): return _pMC_mult.PySwigIterator___sub__(*args)
    def __iter__(self): return self
PySwigIterator_swigregister = _pMC_mult.PySwigIterator_swigregister
PySwigIterator_swigregister(PySwigIterator)

class IntVector(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, IntVector, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, IntVector, name)
    __repr__ = _swig_repr
    def iterator(*args): return _pMC_mult.IntVector_iterator(*args)
    def __iter__(self): return self.iterator()
    def __nonzero__(*args): return _pMC_mult.IntVector___nonzero__(*args)
    def __len__(*args): return _pMC_mult.IntVector___len__(*args)
    def pop(*args): return _pMC_mult.IntVector_pop(*args)
    def __getslice__(*args): return _pMC_mult.IntVector___getslice__(*args)
    def __setslice__(*args): return _pMC_mult.IntVector___setslice__(*args)
    def __delslice__(*args): return _pMC_mult.IntVector___delslice__(*args)
    def __delitem__(*args): return _pMC_mult.IntVector___delitem__(*args)
    def __getitem__(*args): return _pMC_mult.IntVector___getitem__(*args)
    def __setitem__(*args): return _pMC_mult.IntVector___setitem__(*args)
    def append(*args): return _pMC_mult.IntVector_append(*args)
    def empty(*args): return _pMC_mult.IntVector_empty(*args)
    def size(*args): return _pMC_mult.IntVector_size(*args)
    def clear(*args): return _pMC_mult.IntVector_clear(*args)
    def swap(*args): return _pMC_mult.IntVector_swap(*args)
    def get_allocator(*args): return _pMC_mult.IntVector_get_allocator(*args)
    def begin(*args): return _pMC_mult.IntVector_begin(*args)
    def end(*args): return _pMC_mult.IntVector_end(*args)
    def rbegin(*args): return _pMC_mult.IntVector_rbegin(*args)
    def rend(*args): return _pMC_mult.IntVector_rend(*args)
    def pop_back(*args): return _pMC_mult.IntVector_pop_back(*args)
    def erase(*args): return _pMC_mult.IntVector_erase(*args)
    def __init__(self, *args): 
        this = _pMC_mult.new_IntVector(*args)
        try: self.this.append(this)
        except: self.this = this
    def push_back(*args): return _pMC_mult.IntVector_push_back(*args)
    def front(*args): return _pMC_mult.IntVector_front(*args)
    def back(*args): return _pMC_mult.IntVector_back(*args)
    def assign(*args): return _pMC_mult.IntVector_assign(*args)
    def resize(*args): return _pMC_mult.IntVector_resize(*args)
    def insert(*args): return _pMC_mult.IntVector_insert(*args)
    def reserve(*args): return _pMC_mult.IntVector_reserve(*args)
    def capacity(*args): return _pMC_mult.IntVector_capacity(*args)
    __swig_destroy__ = _pMC_mult.delete_IntVector
    __del__ = lambda self : None;
IntVector_swigregister = _pMC_mult.IntVector_swigregister
IntVector_swigregister(IntVector)

class DoubleVector(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, DoubleVector, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, DoubleVector, name)
    __repr__ = _swig_repr
    def iterator(*args): return _pMC_mult.DoubleVector_iterator(*args)
    def __iter__(self): return self.iterator()
    def __nonzero__(*args): return _pMC_mult.DoubleVector___nonzero__(*args)
    def __len__(*args): return _pMC_mult.DoubleVector___len__(*args)
    def pop(*args): return _pMC_mult.DoubleVector_pop(*args)
    def __getslice__(*args): return _pMC_mult.DoubleVector___getslice__(*args)
    def __setslice__(*args): return _pMC_mult.DoubleVector___setslice__(*args)
    def __delslice__(*args): return _pMC_mult.DoubleVector___delslice__(*args)
    def __delitem__(*args): return _pMC_mult.DoubleVector___delitem__(*args)
    def __getitem__(*args): return _pMC_mult.DoubleVector___getitem__(*args)
    def __setitem__(*args): return _pMC_mult.DoubleVector___setitem__(*args)
    def append(*args): return _pMC_mult.DoubleVector_append(*args)
    def empty(*args): return _pMC_mult.DoubleVector_empty(*args)
    def size(*args): return _pMC_mult.DoubleVector_size(*args)
    def clear(*args): return _pMC_mult.DoubleVector_clear(*args)
    def swap(*args): return _pMC_mult.DoubleVector_swap(*args)
    def get_allocator(*args): return _pMC_mult.DoubleVector_get_allocator(*args)
    def begin(*args): return _pMC_mult.DoubleVector_begin(*args)
    def end(*args): return _pMC_mult.DoubleVector_end(*args)
    def rbegin(*args): return _pMC_mult.DoubleVector_rbegin(*args)
    def rend(*args): return _pMC_mult.DoubleVector_rend(*args)
    def pop_back(*args): return _pMC_mult.DoubleVector_pop_back(*args)
    def erase(*args): return _pMC_mult.DoubleVector_erase(*args)
    def __init__(self, *args): 
        this = _pMC_mult.new_DoubleVector(*args)
        try: self.this.append(this)
        except: self.this = this
    def push_back(*args): return _pMC_mult.DoubleVector_push_back(*args)
    def front(*args): return _pMC_mult.DoubleVector_front(*args)
    def back(*args): return _pMC_mult.DoubleVector_back(*args)
    def assign(*args): return _pMC_mult.DoubleVector_assign(*args)
    def resize(*args): return _pMC_mult.DoubleVector_resize(*args)
    def insert(*args): return _pMC_mult.DoubleVector_insert(*args)
    def reserve(*args): return _pMC_mult.DoubleVector_reserve(*args)
    def capacity(*args): return _pMC_mult.DoubleVector_capacity(*args)
    __swig_destroy__ = _pMC_mult.delete_DoubleVector
    __del__ = lambda self : None;
DoubleVector_swigregister = _pMC_mult.DoubleVector_swigregister
DoubleVector_swigregister(DoubleVector)

class FloatVector(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, FloatVector, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, FloatVector, name)
    __repr__ = _swig_repr
    def iterator(*args): return _pMC_mult.FloatVector_iterator(*args)
    def __iter__(self): return self.iterator()
    def __nonzero__(*args): return _pMC_mult.FloatVector___nonzero__(*args)
    def __len__(*args): return _pMC_mult.FloatVector___len__(*args)
    def pop(*args): return _pMC_mult.FloatVector_pop(*args)
    def __getslice__(*args): return _pMC_mult.FloatVector___getslice__(*args)
    def __setslice__(*args): return _pMC_mult.FloatVector___setslice__(*args)
    def __delslice__(*args): return _pMC_mult.FloatVector___delslice__(*args)
    def __delitem__(*args): return _pMC_mult.FloatVector___delitem__(*args)
    def __getitem__(*args): return _pMC_mult.FloatVector___getitem__(*args)
    def __setitem__(*args): return _pMC_mult.FloatVector___setitem__(*args)
    def append(*args): return _pMC_mult.FloatVector_append(*args)
    def empty(*args): return _pMC_mult.FloatVector_empty(*args)
    def size(*args): return _pMC_mult.FloatVector_size(*args)
    def clear(*args): return _pMC_mult.FloatVector_clear(*args)
    def swap(*args): return _pMC_mult.FloatVector_swap(*args)
    def get_allocator(*args): return _pMC_mult.FloatVector_get_allocator(*args)
    def begin(*args): return _pMC_mult.FloatVector_begin(*args)
    def end(*args): return _pMC_mult.FloatVector_end(*args)
    def rbegin(*args): return _pMC_mult.FloatVector_rbegin(*args)
    def rend(*args): return _pMC_mult.FloatVector_rend(*args)
    def pop_back(*args): return _pMC_mult.FloatVector_pop_back(*args)
    def erase(*args): return _pMC_mult.FloatVector_erase(*args)
    def __init__(self, *args): 
        this = _pMC_mult.new_FloatVector(*args)
        try: self.this.append(this)
        except: self.this = this
    def push_back(*args): return _pMC_mult.FloatVector_push_back(*args)
    def front(*args): return _pMC_mult.FloatVector_front(*args)
    def back(*args): return _pMC_mult.FloatVector_back(*args)
    def assign(*args): return _pMC_mult.FloatVector_assign(*args)
    def resize(*args): return _pMC_mult.FloatVector_resize(*args)
    def insert(*args): return _pMC_mult.FloatVector_insert(*args)
    def reserve(*args): return _pMC_mult.FloatVector_reserve(*args)
    def capacity(*args): return _pMC_mult.FloatVector_capacity(*args)
    __swig_destroy__ = _pMC_mult.delete_FloatVector
    __del__ = lambda self : None;
FloatVector_swigregister = _pMC_mult.FloatVector_swigregister
FloatVector_swigregister(FloatVector)

class MC(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, MC, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, MC, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _pMC_mult.new_MC(*args)
        try: self.this.append(this)
        except: self.this = this
    def calc_pKas(*args): return _pMC_mult.MC_calc_pKas(*args)
    def set_MCsteps(*args): return _pMC_mult.MC_set_MCsteps(*args)
    __swig_destroy__ = _pMC_mult.delete_MC
    __del__ = lambda self : None;
MC_swigregister = _pMC_mult.MC_swigregister
MC_swigregister(MC)



