#! /usr/bin/env python
# $Header$
'''Typecodes for dates and times.
'''

from ZSI import _copyright, _floattypes, _inttypes, _get_idstr, EvaluateException
from ZSI.TC import TypeCode, SimpleType
from ZSI.wstools.Namespaces import SCHEMA
import operator, re, time as _time
from time import mktime as _mktime, localtime as _localtime, gmtime as _gmtime
from datetime import tzinfo as _tzinfo, timedelta as _timedelta,\
    datetime as _datetime
from math import modf as _modf

_niltime = [
    0, 0, 0,    # year month day
    0, 0, 0,    # hour minute second
    0, 0, 0     # weekday, julian day, dst flag
]

#### Code added to check current timezone offset
_zero = _timedelta(0)
_dstoffset = _stdoffset = _timedelta(seconds=-_time.timezone)
if _time.daylight: _dstoffset = _timedelta(seconds=-_time.altzone)
_dstdiff = _dstoffset - _stdoffset


class _localtimezone(_tzinfo):
    """ """
    def dst(self, dt):
        """datetime -> DST offset in minutes east of UTC."""
        tt = _localtime(_mktime((dt.year, dt.month, dt.day,
                 dt.hour, dt.minute, dt.second, dt.weekday(), 0, -1)))
        if tt.tm_isdst > 0: return _dstdiff
        return _zero
    
    #def fromutc(...)
    #datetime in UTC -> datetime in local time.

    def tzname(self, dt):
        """datetime -> string name of time zone."""
        tt = _localtime(_mktime((dt.year, dt.month, dt.day,
                 dt.hour, dt.minute, dt.second, dt.weekday(), 0, -1)))
        return _time.tzname[tt.tm_isdst > 0]

    def utcoffset(self, dt):
        """datetime -> minutes east of UTC (negative for west of UTC)."""
        tt = _localtime(_mktime((dt.year, dt.month, dt.day,
                 dt.hour, dt.minute, dt.second, dt.weekday(), 0, -1)))
        if tt.tm_isdst > 0: return _dstoffset
        return _stdoffset

class _fixedoffset(_tzinfo):
    """Fixed offset in minutes east from UTC.
    
    A class building tzinfo objects for fixed-offset time zones.
    Note that _fixedoffset(0, "UTC") is a different way to build a
    UTC tzinfo object.
    """
    #def __init__(self, offset, name):
    def __init__(self, offset):
        self.__offset = _timedelta(minutes=offset)
        #self.__name = name
        
    def dst(self, dt):
        """datetime -> DST offset in minutes east of UTC."""
        return _zero
    
    def tzname(self, dt):
        """datetime -> string name of time zone."""
        #return self.__name
        return "server"
    
    def utcoffset(self, dt):
        """datetime -> minutes east of UTC (negative for west of UTC)."""
        return self.__offset

    
def _dict_to_tuple(d):
    '''Convert a dictionary to a time tuple.  Depends on key values in the
    regexp pattern!
    '''    
    # TODO: Adding a ms field to struct_time tuples is problematic 
    # since they don't have this field.  Should use datetime
    # which has a microseconds field, else no ms..  When mapping struct_time 
    # to gDateTime the last 3 fields are irrelevant, here using dummy values to make
    # everything happy.
    # 

    retval = _niltime[:]
    for k,i in ( ('Y', 0), ('M', 1), ('D', 2), ('h', 3), ('m', 4), ):
        v = d.get(k)
        if v: retval[i] = int(v)
        
    v = d.get('s')
    if v:
        msec,sec = _modf(float(v))
        retval[6],retval[5] = int(round(msec*1000)), int(sec)
            
    v = d.get('tz')
    if v and v != 'Z':
        h,m = map(int, v.split(':'))
        # check for time zone offset, if within the same timezone, 
        # ignore offset specific calculations
        offset=_localtimezone().utcoffset(_datetime.now())
        local_offset_hour = offset.seconds/3600
        local_offset_min = (offset.seconds%3600)%60
        if local_offset_hour > 12: 
            local_offset_hour -= 24
            
        if local_offset_hour != h or local_offset_min != m:                
            if h<0:
                #TODO: why is this set to server
                #foff = _fixedoffset(-((abs(h)*60+m)),"server")
                foff = _fixedoffset(-((abs(h)*60+m)))
            else:
                #TODO: why is this set to server
                #foff = _fixedoffset((abs(h)*60+m),"server")
                foff = _fixedoffset((abs(h)*60+m)) 
                
            dt = _datetime(retval[0],retval[1],retval[2],retval[3],retval[4],
                           retval[5],0,foff)
            
            # update dict with calculated timezone
            localdt=dt.astimezone(_localtimezone())
            retval[0] = localdt.year
            retval[1] = localdt.month
            retval[2] = localdt.day
            retval[3] = localdt.hour
            retval[4] = localdt.minute
            retval[5] = localdt.second
            
    if d.get('neg', 0):
        retval[0:5] = map(operator.__neg__, retval[0:5])
    return tuple(retval)


class Duration(SimpleType):
    '''Time duration.
    '''
    parselist = [ (None,'duration') ]
    lex_pattern = re.compile('^' r'(?P<neg>-?)P' \
                    r'((?P<Y>\d+)Y)?' r'((?P<M>\d+)M)?' r'((?P<D>\d+)D)?' \
                    r'(?P<T>T?)' r'((?P<h>\d+)H)?' r'((?P<m>\d+)M)?' \
                    r'((?P<s>\d*(\.\d+)?)S)?' '$')
    type = (SCHEMA.XSD3, 'duration')


    def text_to_data(self, text, elt, ps):
        '''convert text into typecode specific data.
        '''
        if text is None:
            return None
        m = Duration.lex_pattern.match(text)
        if m is None:
            raise EvaluateException('Illegal duration', ps.Backtrace(elt))
        d = m.groupdict()
        if d['T'] and (d['h'] is None and d['m'] is None and d['s'] is None):
            raise EvaluateException('Duration has T without time')
        try:
            retval = _dict_to_tuple(d)
        except ValueError, e:
            raise EvaluateException(str(e))
    
        if self.pyclass is not None:
            return self.pyclass(retval)
        return retval  
    
    def get_formatted_content(self, pyobj):
        if type(pyobj) in _floattypes or type(pyobj) in _inttypes:
            pyobj = _gmtime(pyobj)
        
        d = {}
        pyobj = tuple(pyobj)
        if 1 in map(lambda x: x < 0, pyobj[0:6]):
            pyobj = map(abs, pyobj)
            neg = '-'
        else:
            neg = ''
            
        val = '%sP%dY%dM%dDT%dH%dM%dS' % \
            ( neg, pyobj[0], pyobj[1], pyobj[2], pyobj[3], pyobj[4], pyobj[5])

        return val
        

class Gregorian(SimpleType):
    '''Gregorian times.
    '''
    lex_pattern = tag = format = None

    def text_to_data(self, text, elt, ps):
        '''convert text into typecode specific data.
        '''
        if text is None:
            return None
        
        m = self.lex_pattern.match(text)
        if not m:
            raise EvaluateException('Bad Gregorian: %s' %text, ps.Backtrace(elt))
        try:
            retval = _dict_to_tuple(m.groupdict())
        except ValueError, e:
            #raise EvaluateException(str(e))
            raise
        
        if self.pyclass is not None:
            return self.pyclass(retval)
        return retval    

    def get_formatted_content(self, pyobj):
        if type(pyobj) in _floattypes or type(pyobj) in _inttypes:
            pyobj = _gmtime(pyobj)
        
        d = {}
        pyobj = tuple(pyobj)
        if 1 in map(lambda x: x < 0, pyobj[0:6]):
            pyobj = map(abs, pyobj)
            d['neg'] = '-'
        else:
            d['neg'] = ''

        ms = pyobj[6]
        if not ms or not hasattr(self, 'format_ms'): 
            d = { 'Y': pyobj[0], 'M': pyobj[1], 'D': pyobj[2],
                'h': pyobj[3], 'm': pyobj[4], 's': pyobj[5], }
            return self.format % d

        if  ms > 999:
            raise ValueError, 'milliseconds must be a integer between 0 and 999'

        d = { 'Y': pyobj[0], 'M': pyobj[1], 'D': pyobj[2],
            'h': pyobj[3], 'm': pyobj[4], 's': pyobj[5], 'ms':ms, }
        return self.format_ms % d
 

class gDateTime(Gregorian):
    '''A date and time.
    '''
    parselist = [ (None,'dateTime') ]
    lex_pattern = re.compile('^' r'(?P<neg>-?)' \
                        '(?P<Y>\d{4,})-' r'(?P<M>\d\d)-' r'(?P<D>\d\d)' 'T' \
                        r'(?P<h>\d\d):' r'(?P<m>\d\d):' r'(?P<s>\d*(\.\d+)?)' \
                        r'(?P<tz>(Z|([-+]\d\d:\d\d))?)' '$')
    tag, format = 'dateTime', '%(Y)04d-%(M)02d-%(D)02dT%(h)02d:%(m)02d:%(s)02dZ'
    format_ms = format[:-1] + '.%(ms)03dZ'
    type = (SCHEMA.XSD3, 'dateTime')

class gDate(Gregorian):
    '''A date.
    '''
    parselist = [ (None,'date') ]
    lex_pattern = re.compile('^' r'(?P<neg>-?)' \
                        '(?P<Y>\d{4,})-' r'(?P<M>\d\d)-' r'(?P<D>\d\d)' \
                        r'(?P<tz>Z|([-+]\d\d:\d\d))?' '$')
    tag, format = 'date', '%(Y)04d-%(M)02d-%(D)02dZ'
    type = (SCHEMA.XSD3, 'date')

class gYearMonth(Gregorian):
    '''A date.
    '''
    parselist = [ (None,'gYearMonth') ]
    lex_pattern = re.compile('^' r'(?P<neg>-?)' \
                        '(?P<Y>\d{4,})-' r'(?P<M>\d\d)' \
                        r'(?P<tz>Z|([-+]\d\d:\d\d))?' '$')
    tag, format = 'gYearMonth', '%(Y)04d-%(M)02dZ'
    type = (SCHEMA.XSD3, 'gYearMonth')

class gYear(Gregorian):
    '''A date.
    '''
    parselist = [ (None,'gYear') ]
    lex_pattern = re.compile('^' r'(?P<neg>-?)' \
                        '(?P<Y>\d{4,})' \
                        r'(?P<tz>Z|([-+]\d\d:\d\d))?' '$')
    tag, format = 'gYear', '%(Y)04dZ'
    type = (SCHEMA.XSD3, 'gYear')

class gMonthDay(Gregorian):
    '''A gMonthDay.
    '''
    parselist = [ (None,'gMonthDay') ]
    lex_pattern = re.compile('^' r'(?P<neg>-?)' \
                        r'--(?P<M>\d\d)-' r'(?P<D>\d\d)' \
                        r'(?P<tz>Z|([-+]\d\d:\d\d))?' '$')
    tag, format = 'gMonthDay', '---%(M)02d-%(D)02dZ'
    type = (SCHEMA.XSD3, 'gMonthDay')


class gDay(Gregorian):
    '''A gDay.
    '''
    parselist = [ (None,'gDay') ]
    lex_pattern = re.compile('^' r'(?P<neg>-?)' \
                        r'---(?P<D>\d\d)' \
                        r'(?P<tz>Z|([-+]\d\d:\d\d))?' '$')
    tag, format = 'gDay', '---%(D)02dZ'
    type = (SCHEMA.XSD3, 'gDay')
    
class gMonth(Gregorian):
    '''A gMonth.
    '''
    parselist = [ (None,'gMonth') ]
    lex_pattern = re.compile('^' r'(?P<neg>-?)' \
                        r'---(?P<M>\d\d)' \
                        r'(?P<tz>Z|([-+]\d\d:\d\d))?' '$')
    tag, format = 'gMonth', '---%(M)02dZ'
    type = (SCHEMA.XSD3, 'gMonth')
    
class gTime(Gregorian):
    '''A time.
    '''
    parselist = [ (None,'time') ]
    lex_pattern = re.compile('^' r'(?P<neg>-?)' \
                        r'(?P<h>\d\d):' r'(?P<m>\d\d):' r'(?P<s>\d*(\.\d+)?)' \
                        r'(?P<tz>Z|([-+]\d\d:\d\d))?' '$')
    tag, format = 'time', '%(h)02d:%(m)02d:%(s)02dZ'
    format_ms = format[:-1] + '.%(ms)03dZ'
    type = (SCHEMA.XSD3, 'time')

if __name__ == '__main__': print _copyright
