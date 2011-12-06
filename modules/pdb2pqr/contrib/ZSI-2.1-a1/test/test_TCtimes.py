#!/usr/bin/env python
import unittest, sys, tests_good, tests_bad, time
from ZSI import *
try:
    import cStringIO as StringIO
except ImportError:
    import StringIO


class TestCase(unittest.TestCase):
    '''Examples from "Definitive XML Schema, Priscilla Walmsley, p237-246
    '''
    def check_dateTime_local_offset(self):
        # UTC with local timezone offset
        # 
        typecode = TC.gDateTime()
        off_hour = time.altzone/60/60
        off_min = time.altzone%60
        stamp_offset = '1968-04-02T13:20:00+%02d:%02d' %(off_hour,off_min)
        data = typecode.text_to_data(stamp_offset, None, None)
        stamp = typecode.get_formatted_content(data)

        correct = "1968-04-01T22:20:00Z"
        self.failUnless(stamp == correct, 
            'dateTime with local offset(%s), expecting "%s" got "%s"' %(
            stamp_offset, correct, stamp))

    def check_valid_dateTime(self):
        typecode = TC.gDateTime()
        for i in ('1968-04-02T13:20:00', '1968-04-02T13:20:15.5', 
            '1968-04-02T13:20:00-05:00', '1968-04-02T13:20:00Z'):
            data = typecode.text_to_data(i, None, None)
            text = typecode.get_formatted_content(data)

    def check_parse_microseconds(self):
        good = (1968, 4, 2, 13, 20, 15, 511, 0, 0)
        typecode = TC.gDateTime()
        data = typecode.text_to_data('1968-04-02T13:20:15.511', None, None)
        self.failUnless(data == good,
            'did not parse something %s, not equal %s' %(data,good))

    def check_serialize_microseconds(self):
        dateTime = '1968-04-02T13:20:15.511Z'
        typecode = TC.gDateTime()
        text = typecode.get_formatted_content((1968, 4, 2, 13, 20, 15, 511, 0, 0))
        self.failUnless(text == dateTime,
            'did not serialze correctly %s, not equal %s' %(text, dateTime))

    def check_serialize_microseconds_1000(self):
        bad = (1968, 4, 2, 13, 20, 15, 1000, 0)
        typecode = TC.gDateTime()
        self.failUnlessRaises(ValueError, typecode.get_formatted_content, bad)

    def check_serialize_microseconds_lessZero(self):
        '''ignore negative microseconds
        '''
        bad = (1968, 4, 2, 13, 20, 15, -1, 0)
        typecode = TC.gDateTime()
        text = typecode.get_formatted_content(bad)
        typecode.get_formatted_content(bad)

    def check_parse_microseconds2(self):
        good = (1968, 4, 2, 13, 20, 15, 500, 0, 0)
        typecode = TC.gDateTime()
        data = typecode.text_to_data('1968-04-02T13:20:15.5Z', None,None)
        self.failUnless(data == good,
            'did not serialze correctly %s, not equal %s' %(data, good))

        #text = typecode.get_formatted_content((1968, 4, 2, 13, 20, 15, 5, 0, 500))
        #self.failUnless(text == dateTime,
        #    'did not serialze correctly %s, not equal %s' %(text, dateTime))

    def check_invalid_dateTime(self):
        typecode = TC.gDateTime()

    def check_valid_time(self):
        typecode = TC.gTime()
        for i in ('13:20:00', '13:20:30.5555', '13:20:00Z'):
            data = typecode.text_to_data(i, None, None)
            text = typecode.get_formatted_content(data)

    def broke_valid_time(self):
        typecode = TC.gTime()
        data = typecode.text_to_data('13:20:00-05:00', None, None)

    def check_invalid_time(self):
        typecode = TC.gTime()
        for i in ('5:20:00', '13:20.5:00',):
            self.failUnlessRaises(Exception, typecode.text_to_data, i, None, None),

    def broke_invalid_time_no_seconds(self):
        typecode = TC.gTime()
        i = '13:20:'
        self.failUnlessRaises(Exception, typecode.text_to_data, i, None, None)

    def broke_invalid_time_bad_timeofday(self):
        typecode = TC.gTime()
        i = '13:65:00'
        self.failUnlessRaises(Exception, typecode.text_to_data, i, None, None)

    def check_valid_date(self):
        typecode = TC.gDate()
        for i in ('1968-04-02', '-0045-01-01', '11968-04-02', '1968-04-02+05:00', '1968-04-02Z'):
            data = typecode.text_to_data(i, None, None)
            text = typecode.get_formatted_content(data)

    def check_invalid_date(self):
        typecode = TC.gDate()
        for i in ('68-04-02', '1968-4-2', '1968/04/02', '04-02-1968',):
            self.failUnlessRaises(Exception, typecode.text_to_data, i, None, None),

    def broke_invalid_date_april31(self):
        # No checks for valid date April 30 days
        typecode = TC.gDate()
        self.failUnlessRaises(Exception, typecode.text_to_data, '1968-04-31', None, None),

#
# Creates permutation of test options: "check", "check_any", etc
#
_SEP = '_'
for t in [i[0].split(_SEP) for i in filter(lambda i: callable(i[1]), TestCase.__dict__.items())]:
    test = ''
    for f in t:
        test += f
        if globals().has_key(test): test += _SEP; continue
        def _closure():
            name = test
            def _makeTestSuite():
                suite = unittest.TestSuite()
                suite.addTest(unittest.makeSuite(TestCase, name))
                return suite
            return _makeTestSuite

        globals()[test] = _closure()
        test += _SEP


makeTestSuite = check
def main():
    unittest.main(defaultTest="makeTestSuite")
if __name__ == "__main__" : main()


