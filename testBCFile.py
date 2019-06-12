import tables
import sys
sys.path.append('/Users/srkaeppler/Dropbox/research/data/NSF_Dregion_ParticlePrecipitation/Models')

import IRI2016
iri2016 = IRI2016.IRI2016()

import MSIS
msis = MSIS.MSIS()

import numpy
from scipy import interpolate
import time



fname_ac = '/Users/srkaeppler/Dropbox/research/data/CCMC_Conductivites/17march2013/20130316.006_ac_5min-cal.h5'

# this will open and parse data for a straight AC file.
with tables.open_file(fname_ac) as h5:


    ne1=h5.root.FittedParams.Ne.read()#dat1['/FittedParams']['Ne']
    dne1=h5.root.FittedParams.dNe.read()# dne1=dat1['/FittedParams']['dNe']
    time1=h5.root.Time.UnixTime.read()#dat1['/Time']['UnixTime']
    dtime1=h5.root.Time.dtime.read()#dat1['/Time']['dtime']
# (Nrecs1,Nbeams1,Nhts1)=vlos1.shape
    Altitude=h5.root.FittedParams.Altitude.read()#dat1['/FittedParams']['Altitude']

    Babs1=h5.root.Geomag.Babs.read()#dat1['/Geomag']['Babs'];
    if Babs1[0,0]<1.0e-5:
        Babs1=Babs1*1.0e5 # convert to nanoTelsa
# '/FittedParams/IonMass' : [('TITLE','Ion Mass'),('Description', 'Mass of ions used in fitting'),('Unit','amu')],\
# '/FittedParams/Fits' : [('TITLE','Fits'),('Description','Fitted parameters'),('Size','Nrecords x Nbeams x Nranges x Nions+1 x 4 (fraction, temperature, collision frequency, LOS speed)'),('Unit','N/A, Kelvin, s^{-1}, m/s')],\
    mass=h5.root.FittedParams.IonMass.read()#dat1['/FittedParams']['IonMass']
    Fits=h5.root.FittedParams.Fits.read()#dat1['/FittedParams']['Fits'][:,:,:,0:2,0] # O+ fraction
    BeamCodes = h5.root.BeamCodes.read()
    # convert to cm^-3
    nO = h5.root.MSIS.nO.read()/1e6
    nO2 = h5.root.MSIS.nO2.read()/1e6
    nN2 = h5.root.MSIS.nN2.read()/1e6
    Tn = h5.root.MSIS.Tn.read()
    IonMass = h5.root.FittedParams.IonMass.read()
    #nuin=dat1['/FittedParams']['Fits'][:,:,:,0:2,2]

fraction = Fits[:,:,:,0:-1,0] # fraction including electrons as index = -1
nuin = Fits[:,:,:,0:-1,2] # collision frequency, including electrons as index = -1
nuen = Fits[:,:,:,-1,2]
Ti = Fits[:,:,:,0:-1,1]
Te = Fits[:,:,:,-1,1]

fname_bc = '/Users/srkaeppler/Dropbox/research/data/CCMC_Conductivites/17march2013/20130317.004_bc_2min-Ne-cal.h5'

bch5 = tables.open_file(fname_bc)

nebc=bch5.root.NeFromPower.Ne_Mod.read()#dat1['/FittedParams']['Ne']
dnebc=bch5.root.NeFromPower.dNeFrac.read()# dne1=dat1['/FittedParams']['dNe']
timebc=bch5.root.Time.UnixTime.read()#dat1['/Time']['UnixTime']
dtimebc=bch5.root.Time.dtime.read()#dat1['/Time']['dtime']
Altitude=bch5.root.NeFromPower.Altitude.read()
BeamCodes = bch5.root.BeamCodes.read()
glat = bch5.root.Site.Latitude.read()
glon = bch5.root.Site.Longitude.read()
# time by beam by altitude
nObc = numpy.zeros([timebc.shape[0],BeamCodes.shape[0],Altitude.shape[1]])
nO2bc = numpy.zeros([timebc.shape[0],BeamCodes.shape[0],Altitude.shape[1]])
nN2bc = numpy.zeros([timebc.shape[0],BeamCodes.shape[0],Altitude.shape[1]])
Tnbc = numpy.zeros([timebc.shape[0],BeamCodes.shape[0],Altitude.shape[1]])
Tebc = numpy.zeros([timebc.shape[0],BeamCodes.shape[0],Altitude.shape[1]])
Tibc = numpy.zeros([timebc.shape[0],BeamCodes.shape[0],Altitude.shape[1]])

Babs1=bch5.root.Geomag.Babs.read()#dat1['/Geomag']['Babs'];
if Babs1[0,0]<1.0e-5:
    Babs1=Babs1*1.0e5

# I now need the IRI Temperature and MSIS to finish this up
# technically should do this over all beams, but at this point makes sense
# to select the beam
#testDict = iri.IRI2016(tUnix,Lat,Lon,Heibeg,Heiend,step)
timebcMean = numpy.mean(timebc,axis=1)
for itime in range(timebcMean.shape[0]):
    tmpIRI = iri2016.IRI2016(timebcMean[itime],glat,glon,80,140,0.5)
    tmpMSIS = msis.MSIS2(timebcMean[itime],glat,glon,80,140.,0.5,CGSorSI='SI')
    AltGrid = numpy.arange(80.,140.,0.5)*1000.
    for ibeam in range(Altitude.shape[0]):
        f = interpolate.interp1d(AltGrid, tmpMSIS['nO'],bounds_error=False,fill_value=numpy.nan)
        nObc[itime,ibeam,:] = f(Altitude[ibeam,:])
        del f

        f = interpolate.interp1d(AltGrid, tmpMSIS['nO2'],bounds_error=False,fill_value=numpy.nan)
        nO2bc[itime,ibeam,:] = f(Altitude[ibeam,:])

        del f

        f = interpolate.interp1d(AltGrid, tmpMSIS['nN2'],bounds_error=False,fill_value=numpy.nan)
        nN2bc[itime,ibeam,:] = f(Altitude[ibeam,:])

        del f

        f = interpolate.interp1d(AltGrid, tmpMSIS['Tn'],bounds_error=False,fill_value=numpy.nan)
        Tnbc[itime,ibeam,:] = f(Altitude[ibeam,:])
        del f

        f = interpolate.interp1d(AltGrid, tmpIRI['Te'],bounds_error=False,fill_value=numpy.nan)
        Tebc[itime,ibeam,:] = f(Altitude[ibeam,:])
        del f

        f = interpolate.interp1d(AltGrid, tmpIRI['Ti'],bounds_error=False,fill_value=numpy.nan)
        Tibc[itime,ibeam,:] = f(Altitude[ibeam,:])
        del f
        print 'ibeam',ibeam


    print 'itime',itime
