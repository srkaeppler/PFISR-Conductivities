import os
import sys
import math
import scipy
import datetime
import glob
import tables
# import matplotlib
# matplotlib.use('Agg')
# import pylab
import numpy
from scipy.integrate import simps

# some natural constants
v_lightspeed=299792458
v_Boltzmann=1.380658e-23
v_electronmass=9.1093897e-31
v_amu=1.6605402e-27
v_electronradius=2.81794092e-15
v_epsilon0=8.854188E-12
v_elemcharge=1.602177E-19
pi=math.pi

"""
This program is designed to quickly calculate the conductivity and conductance

"""

# read in fitted data


def ProcessFAConductivity(fname_ac):
    # bring this up to speed with python2.7
    # read in the data
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
    #nuin=dat1['/FittedParams']['Fits'][:,:,:,0:2,2]

    fraction = Fits[:,:,:,0:-1,0] # fraction including electrons as index = -1
    nuin = Fits[:,:,:,0:-1,2] # collision frequency, including electrons as index = -1
    nuen = Fits[:,:,:,-1,2]
    PedCond = numpy.zeros(ne1.shape)
    HallCond = numpy.zeros(ne1.shape)
    # define the mobility
    mob=v_elemcharge/(v_amu*nuin); #mob[:,:,:,0]=mob[:,:,:,0]/mass[0]; mob[:,:,:,1]=mob[:,:,:,1]/mass[
    nuinScaler = 1.0
    Babs = numpy.tile(Babs1,ne1.shape[0]).reshape(ne1.shape) # generate a time x beam x altitude array of magnetic field

    for i in range(mass.shape[0]):
        mob[:,:,:,i] = mob[:,:,:,i]/mass[i]
        kappa = mob*1.0
        kappa[:,:,:,i] = kappa[:,:,:,i]*Babs1
        # equation 2.40a in Kelley's textbook
        sp1 = ne1*v_elemcharge*v_elemcharge/(v_amu*mass[i]*nuinScaler*nuin[:,:,:,i]*(1.0+(kappa[:,:,:,i]/nuinScaler)**2.0))*fraction[:,:,:,i]
        # equation 10 from Bostrom 1964 in my notebook
        # should be good > 85 km
        sh1 = ne1*v_elemcharge*fraction[:,:,:,i]/(Babs*(1.0+(kappa[:,:,:,i]/nuinScaler)**2.0)) # fraction in effect weights things
        print sp1
        PedCond = PedCond + sp1
        HallCond = HallCond + sh1


    HallConductance = numpy.zeros(HallCond.shape[0]) # time axis
    PedConductance = numpy.zeros(PedCond.shape[0])
    FA_index = numpy.where((BeamCodes[:,1] == -154.3) & (BeamCodes[:,2] == 77.5))[0][0]
    if FA_index:
        for i in range(HallCond.shape[0]):
            # get out nan values
            qh = numpy.where((~numpy.isnan(HallCond[i,FA_index,:])) & (Altitude[FA_index,:] > 85e3))
            qp = numpy.where((~numpy.isnan(PedCond[i,FA_index,:])) &  (Altitude[FA_index,:] > 85e3))
            # print HallCond[i,1,qh],Altitude[1,qh]
            if (len(qh[0]) > 0):
                print 'Hall', simps(HallCond[i,FA_index,qh],Altitude[FA_index,qh])
                HallConductance[i] = simps(HallCond[i,FA_index,qh],Altitude[FA_index,qh])[0]
            else:
                HallConductance[i] = -1.
            if (len(qp[0]) > 0):
                print 'Ped', simps(PedCond[i,FA_index,qp],Altitude[FA_index,qp])
                PedConductance[i] = simps(PedCond[i,FA_index,qp],Altitude[FA_index,qp])[0]
            else:
                PedConductance[i] = -1.
        HallConductivity = HallCond[:,FA_index,:]
        PedConductivity = PedCond[:,FA_index,:]
        TimeUTHour = numpy.mean(dtime1, axis=1)
        FA_Altitude = Altitude[FA_index,:]
    else:
        HallConductivity = []
        PedConductivity = []
        TimeUTHour = []
        FA_Altitude = []

    return HallConductivity, PedConductivity,PedConductance, HallConductance,TimeUTHour,FA_Altitude

fname_ac=['./20130316.006.done/20130316.006_ac_5min-cal.h5',\
          './20130317.003/20130317.003_ac_5min-cal.h5', \
          './20130317.005.done/20130317.005_ac_5min-cal.h5']
        #   './20170317.005/20170317.005_ac_3min-fitcal.h5',\
        #   './20170317.007/20170317.007_ac_3min-fitcal.h5']
HallConductivity = numpy.zeros([1,25])
HallConductance = numpy.zeros([0])
PedersenConductivity = numpy.zeros([1,25])
PedersenConductance = numpy.zeros([0])
TimeUTHour = numpy.zeros([0])
Altitude = []

for ifile in fname_ac:
    tH1,tP1,tP2,tH2,t,alt = ProcessFAConductivity(ifile)
    HallConductivity = numpy.concatenate((HallConductivity,tH1),axis=0)
    PedersenConductivity = numpy.concatenate((PedersenConductivity,tP1),axis=0)
    PedersenConductance = numpy.concatenate((PedersenConductance,tP2), axis=0)
    HallConductance = numpy.concatenate((HallConductance,tH2), axis=0)
    TimeUTHour = numpy.concatenate((TimeUTHour,t), axis=0)

# remove first element of Hall conductivity
HallConductivity = HallConductivity[1:,:]
PedersenConductivity = PedersenConductivity[1:,:]

print PedersenConductivity.shape

# X = numpy.array([TimeUTHour,PedersenConductance,HallConductance])
# with open('20130317_Conductances.txt', 'wb') as f:
#     f.write('Time (Decimal UT Hours), Pedersen Conductance (mho), Hall Conductance (mho) \n')
#     numpy.savetxt(f,X.T)

# sort by time

# qsort = numpy.argsort(TimeUTHour)
# HallConductance = HallConductance[qsort]
# PedersenConductance = PedersenConductance[qsort]
# HallConductivity = HallConductivity[qsort,:]
# PedersenConductivity = PedersenConductivity[qsort,:]
# TimeUTHour = TimeUTHour[qsort]
# Altitude = alt

# save out to an idl save file


# now calculate the electron part for Pedersen
# sp1e=ne1*v_elemcharge*v_elemcharge/(v_amu*mass[i]*nuinScaler*nuin[:,:,:,i]*(1.0+(kappa[:,:,:,i]/nuinScaler)**2.0))*fraction[:,:,:,i]

# mob=v_elemcharge/(v_amu*nuin); mob[:,:,:,0]=mob[:,:,:,0]/mass[0]; mob[:,:,:,1]=mob[:,:,:,1]/mass[1]
# tmob = mob
# mob = tmob*srkMob # see if this makes a difference in the kappa
# # print 'stop and check the mobility factor of 2'
# # sys.exit()
# kappa=mob*1.0; kappa[:,:,:,0]=kappa[:,:,:,0]*Babs1; kappa[:,:,:,1]=kappa[:,:,:,1]*Babs1
# # mass-weighted params
# nuin1=nuin[:,:,:,0]*fraction[:,:,:,0] + nuin[:,:,:,1]*fraction[:,:,:,1] # mass weighted collision frequency
# mob1=mob[:,:,:,0]*fraction[:,:,:,0]+mob[:,:,:,1]*fraction[:,:,:,1] # mass weighted mobility
# kappa1=kappa[:,:,:,0]*fraction[:,:,:,0]+kappa[:,:,:,1]*fraction[:,:,:,1] # mass weighted kappa
# sp1=ne1*v_elemcharge*v_elemcharge/(v_amu*mass[0]*nuinScaler*nuin[:,:,:,0]*(1.0+(kappa[:,:,:,0]/nuinScaler)**2.0))*fraction[:,:,:,0]
# +ne1*v_elemcharge*v_elemcharge/(v_amu*mass[1]*nuinScaler*nuin[:,:,:,1]*(1.0+(kappa[:,:,:,1]/nuinScaler)**2.0))*fraction[:,:,:,1] # Pedersen conductance
