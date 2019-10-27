import datetime
import numpy

"""
Written in python3
"""



def UnixTime2PltLabels(unixTime):
    t1970 = datetime.datetime(1970,1,1,0,0,0)
    dtStart = datetime.datetime.utcfromtimestamp(unixTime[0])
    dtStartDay = datetime.datetime(dtStart.year,dtStart.month,dtStart.day,0,0,0)
    tUnixStartTime = (dtStartDay-t1970).total_seconds()

    dtEnd = datetime.datetime.utcfromtimestamp(unixTime[-1])
    dtEndDay = datetime.datetime(dtEnd.year,dtEnd.month,dtEnd.day+1,0,0,0)
    tUnixEndTime = (dtEndDay-t1970).total_seconds()

    TotalTimeDifferenceSeconds = unixTime[-1]-unixTime[0]

    SecondsInDay = 24.*3600.
    SecondsIn2Days = 2.*24.*3600.
    SecondsIn4Days = 4.*24.*3600.

    if (TotalTimeDifferenceSeconds < SecondsInDay):
        dt = 4.*3600. # set to 4 hours
    if (TotalTimeDifferenceSeconds > SecondsInDay) & (TotalTimeDifferenceSeconds < SecondsIn2Days):
        dt = 6.*3600.
    if (TotalTimeDifferenceSeconds > SecondsIn2Days) & (TotalTimeDifferenceSeconds < SecondsIn4Days):
        dt = 12.*3600.
    if (TotalTimeDifferenceSeconds > SecondsIn4Days):
        dt = 24.*3600.
    locs = numpy.arange(tUnixStartTime,tUnixEndTime+dt,dt)
    labels = list()
    for itime in locs:
        tmpdt = datetime.datetime.utcfromtimestamp(itime)
        tmplabel = '%02d/%02d\n%02d:%02d'%(tmpdt.month,tmpdt.day,tmpdt.hour,tmpdt.minute)
        #print(tmplabel)
        labels.append(tmplabel)
    return locs,labels

t1970 = datetime.datetime(1970,1,1,0,0,0)
dtTest = datetime.datetime(2010,1,1,0,0)
tUnixTest = (dtTest-t1970).total_seconds()

# now test first 1 day of data
tUnixTest1 = tUnixTest+numpy.arange(0,22.*3600.,3600.)
print(tUnixTest1)
testLoc,testLabel = UnixTime2PltLabels(tUnixTest1)
print('testLoc',testLoc)
print('testLabel',testLabel)
print('\n\n')

tUnixTest1 = tUnixTest+numpy.arange(0,40.*3600.,3600.)
print(tUnixTest1)
testLoc,testLabel = UnixTime2PltLabels(tUnixTest1)
print('testLoc',testLoc)
print('testLabel',testLabel)
print('\n\n')

tUnixTest1 = tUnixTest+numpy.arange(0,24.*3.*3600.,3600.)
print(tUnixTest1)
testLoc,testLabel = UnixTime2PltLabels(tUnixTest1)
print('testLoc',testLoc)
print('testLabel',testLabel)
print('\n\n')

tUnixTest1 = tUnixTest+numpy.arange(0,24.*5.*3600.,3600.)
print(tUnixTest1)
testLoc,testLabel = UnixTime2PltLabels(tUnixTest1)
print('testLoc',testLoc)
print('testLabel',testLabel)
print('\n\n')
