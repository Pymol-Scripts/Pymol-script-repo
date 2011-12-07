#!/usr/bin/env python
############################################################################
# Joshua R. Boverhof, LBNL
# See LBNLCopyright for copyright notice!
###########################################################################
import sys, time
from ZSI.ServiceContainer import AsServer
from DateService_server import simple_Date_Service as _DateService 

class Service(_DateService):
    def soap_getCurrentDate(self, ps):
        request,response = _DateService.soap_getCurrentDate(self, ps)
        response.Today = today = response.new_today()
        request._input
        _SetCurrentDate(today)
        return request,response

    def soap_getDate(self, ps):
        request,response = _DateService.soap_getDate(self, ps)
        response.Day = day = response.new_day()
        _SetDay(day, offset=request.Offset, 
                            date=request.Someday)       
        return request,response


## ADDED WORKER CODE
def _SetCurrentDate(today):
    dt = time.localtime(time.time())
    today.Year = dt[0]
    today.Month = dt[1]
    today.Day = dt[2]
    today.Hour = dt[3]
    today.Minute = dt[4]
    today.Second = dt[5]
    today.Weekday = dt[6]
    today.DayOfYear = dt[7]
    today.Dst = dt[8]

def _SetDay(someDay, offset=None, date=None):
    sec = 3600 * 24  ## seconds/hour * 24h
    providedDate_tuple = (date._year, date._month, date._day,
                          date._hour, date._minute, date._second,
                          date._weekday, date._dayOfYear, date._dst)
    providedDate_sec = time.mktime(providedDate_tuple)
    offset_sec = sec * offset
    newDate_sec = providedDate_sec + offset_sec
    newDate_tuple = time.localtime(newDate_sec)
    if not offset:
        offset = 0

    if not date:
        raise RuntimeError, "Date is required"

    someDay._year = newDate_tuple[0]
    someDay._month = newDate_tuple[1]
    someDay._day = newDate_tuple[2] 
    someDay._hour = newDate_tuple[3]
    someDay._minute = newDate_tuple[4]
    someDay._second = newDate_tuple[5]
    someDay._weekday = newDate_tuple[6]
    someDay._dayOfYear = newDate_tuple[7]
    someDay._dst = newDate_tuple[8]      
    return someDay


if __name__ == "__main__" :
    port = int(sys.argv[1])
    AsServer(port, (Service('test'),))
