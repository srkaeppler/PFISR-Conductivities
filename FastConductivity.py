import os
import sys
import math
import scipy
import datetime
import glob
import tables
import matplotlib
#matplotlib.use('Agg')
import pylab as plt
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


def compute_collfreq(nO,nN2,nO2,Tn, Ti=None,Te=1000.0, mj=30.0):

#   nu_in = 0.0
#   nu_in = nu_in + d[1]*scipy.sqrt(0.79/16) # O
#   nu_in = nu_in + d[2]*scipy.sqrt(1.76/28) # N2
#   nu_in = nu_in + d[3]*scipy.sqrt(1.59/32) # O2
#   nu_in = nu_in*2.6e-9
# d[0] - HE NUMBER DENSITY(CM-3)
#   d[1] - O NUMBER DENSITY(CM-3)
#   d[2] - N2 NUMBER DENSITY(CM-3)
#   d[3] - O2 NUMBER DENSITY(CM-3)

    if len(Ti) > 0:
        if Tn.shape[0] == Ti.shape[0]:
            Tr = (Tn+Ti)/2.0
        else:
            raise Exception( 'altitude arrays are not equal for Ti and Tn')
    else:
        Ti=Tn
        Tr = (Tn+Ti)/2.0
    nu_in = 0.0
    if mj==16.0: # O
        #O+ + O
        # n(O) = d[1]
        # Schunk and Nagy Table 4.5
        nu_in = nu_in + nO*3.67e-11*numpy.sqrt(Tr)*(1.0-0.064*numpy.log10(Tr))**2.0
    else:
        nu_in = nu_in + 1.0e6*nO*9.14e-15/numpy.sqrt(mj*(mj+16.0)) # not entirely sure where this comes from
    if mj==28.0: # N2
        # n(N2) = d[2]
        #N2+ + N2
        # Schunk and Nagy Table 4.5
        nu_in = nu_in + nN2*5.14e-11*numpy.sqrt(Tr)*(1.0-0.069*numpy.log10(Tr))**2.0
    else:
        nu_in = nu_in + 1.0e6*nN2*1.80e-14/numpy.sqrt(mj*(mj+28.0))
    if mj==32.0 and Tn>800.0: # O2
        #n(O2) = d[3]
        # O2+ + O2
        # Schunk and Nagy Table 4.5
        nu_in = nu_in + nO2*2.59e-11*numpy.sqrt(Tr)*(1.0-0.073*numpy.log10(Tr))**2.0
    else:
        nu_in = nu_in + 1.0e6*nO2*1.83e-14/numpy.sqrt(mj*(mj+32.0))

    nu_en = 0.0
    nu_en = nu_en +nO*8.2e-10*numpy.sqrt(Te) # O
    nu_en = nu_en + nN2*2.33e-11*(1-1.2e-4*Te)*Te # N2
    nu_en = nu_en + nO2*1.8e-10*(1+3.6e-2*numpy.sqrt(Te))*numpy.sqrt(Te)

    return nu_in, nu_en





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
    PedCond = numpy.zeros(ne1.shape)
    HallCond = numpy.zeros(ne1.shape)

    print nO.shape, nN2.shape, nO2.shape, Tn.shape, Ti[:,:,:,0].shape

    tmpnuin,tmpnuen = compute_collfreq(nO,nN2,nO2,Tn, Ti=Ti[:,:,:,0])
    print 'tmp nu in',tmpnuin.shape
    print 'nu_in',nuin.shape
    # define the mobility
    #mob=v_elemcharge/(v_amu*nuin); #mob[:,:,:,0]=mob[:,:,:,0]/mass[0]; mob[:,:,:,1]=mob[:,:,:,1]/mass[
    mob = v_elemcharge/(v_amu*tmpnuin)
    nuinScaler = 1.0
    Babs = numpy.tile(Babs1,ne1.shape[0]).reshape(ne1.shape) # generate a time x beam x altitude array of magnetic field

    print 'nO.shape', nO.shape
    print 'IonMass', IonMass.shape, IonMass
    # print 'Ti', Ti
    print 'nu in.shape', nuin.shape
    for i in range(mass.shape[0]):
        mob[:,:,:,i] = mob[:,:,:,i]/mass[i]
        kappa = mob*1.0
        kappa[:,:,:,i] = kappa[:,:,:,i]*Babs1
        # equation 2.40a in Kelley's textbook
        sp1 = ne1*v_elemcharge*v_elemcharge/(v_amu*mass[i]*nuinScaler*nuin[:,:,:,i]*(1.0+(kappa[:,:,:,i]/nuinScaler)**2.0))*fraction[:,:,:,i]
        # equation 10 from Bostrom 1964 in my notebook
        # should be good > 85 km
        sh1 = ne1*v_elemcharge*fraction[:,:,:,i]/(Babs*(1.0+(kappa[:,:,:,i]/nuinScaler)**2.0)) # fraction in effect weights things
        # print sp1
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
                # print 'Hall', simps(HallCond[i,FA_index,qh],Altitude[FA_index,qh])
                HallConductance[i] = simps(HallCond[i,FA_index,qh],Altitude[FA_index,qh])[0]
            else:
                HallConductance[i] = -1.
            if (len(qp[0]) > 0):
                # print 'Ped', simps(PedCond[i,FA_index,qp],Altitude[FA_index,qp])
                PedConductance[i] = simps(PedCond[i,FA_index,qp],Altitude[FA_index,qp])[0]
            else:
                PedConductance[i] = -1.
        HallConductivity = HallCond[:,FA_index,:]
        PedConductivity = PedCond[:,FA_index,:]
        TimeUTHour = numpy.mean(dtime1, axis=1)
        UnixTime = numpy.mean(time1,axis=1)
        FA_Altitude = Altitude[FA_index,:]
    else:
        HallConductivity = []
        PedConductivity = []
        TimeUTHour = []
        UnixTime = []
        FA_Altitude = []

    return HallConductivity, PedConductivity,PedConductance, HallConductance,\
           TimeUTHour,FA_Altitude, UnixTime

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
UnixTime = numpy.zeros([0])
Altitude = []

for ifile in fname_ac:
    tH1,tP1,tP2,tH2,t,alt,tunixTime = ProcessFAConductivity(ifile)
    HallConductivity = numpy.concatenate((HallConductivity,tH1),axis=0)
    PedersenConductivity = numpy.concatenate((PedersenConductivity,tP1),axis=0)
    PedersenConductance = numpy.concatenate((PedersenConductance,tP2), axis=0)
    HallConductance = numpy.concatenate((HallConductance,tH2), axis=0)
    TimeUTHour = numpy.concatenate((TimeUTHour,t), axis=0)
    UnixTime = numpy.concatenate((UnixTime,tunixTime), axis=0)

# remove first element of Hall conductivity
HallConductivity = HallConductivity[1:,:]
PedersenConductivity = PedersenConductivity[1:,:]

# sort out the data before printing.
# sort out the times and just make sure everything is ordered by time
#unix time should be increasing
qsort = numpy.argsort(UnixTime)
UnixTime = UnixTime[qsort]
TimeUTHour = TimeUTHour[qsort]
PedersenConductivity = PedersenConductivity[qsort,:]
HallConductivity = HallConductivity[qsort,:]
PedersenConductance = PedersenConductance[qsort]
HallConductance = HallConductance[qsort]

# time gaps program
plt.figure(1)
plt.plot((UnixTime-UnixTime[0])/3600., PedersenConductance, 'b-', label='P')
plt.plot((UnixTime-UnixTime[0])/3600., HallConductance, 'r-', label='H')
plt.legend()
plt.show()


print PedersenConductivity.shape

# plt.figure(1)
# plt.plot(PedersenConductivity)
# plt.show()

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
