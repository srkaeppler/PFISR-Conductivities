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
import pickle

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

Assumptions:
    1. only examine field-aligned radar look direction
    2. The ion neutral collision frequency is weighted by
    concentration, which I believe is from the FLIP model.
    3. I am not entirely sure where one of the ion neutral collision
    frequencies is coming from in the "else" part of the call.
    Contacted Ashton Reimer about that, he didn't know either.
    Mike Nicolls was the original author.
    4. Only looking at conductivities > 85 km altitude


"""


class Conductivity:
    def __init__(self,location='PFISR'):

        self.location = location
        # can configure for other radars too...
        if location == 'PFISR':
            self.azFA = -154.3
            self.elFA = 77.5




    def compute_collfreq(self,nO,nN2,nO2,Tn, Ti=None,Te=1000.0, mj=30.0):
        """
        This was taken from model_utils.py as part of the standard ISR Fitter
        SRK made edits to some of the variables to make it more readable


    #   nu_in = 0.0
    #   nu_in = nu_in + d[1]*scipy.sqrt(0.79/16) # O
    #   nu_in = nu_in + d[2]*scipy.sqrt(1.76/28) # N2
    #   nu_in = nu_in + d[3]*scipy.sqrt(1.59/32) # O2
    #   nu_in = nu_in*2.6e-9
    # d[0] - HE NUMBER DENSITY(CM-3)
    #   d[1] - O NUMBER DENSITY(CM-3)
    #   d[2] - N2 NUMBER DENSITY(CM-3)
    #   d[3] - O2 NUMBER DENSITY(CM-3)
        """
        if Ti is not None:
            if Tn.shape[0] == Ti.shape[0]:
                Tr = (Tn+Ti)/2.0

            else:
                raise Exception( 'altitude arrays are not equal for Ti and Tn')
        else:
            print 'in else'
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

        if mj==32.0:# and Tn>800.0: # O2
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

    def ProcessFAConductivity(self,fname_ac, opt='Original'):
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

        # print nO.shape, nN2.shape, nO2.shape, Tn.shape, Ti[:,:,:,0].shape
        # print 'fraction', fraction
        nuinISR = numpy.zeros(Ti.shape)*numpy.nan
        for i in range(mass.shape[0]):
            print mass[i]
            tmpnuin,tmpnuen = self.compute_collfreq(nO,nN2,nO2,Tn,Ti=Ti[:,:,:,i],mj=mass[i])
            # print 'tmp nu in',tmpnuin.shape
            nuinISR[:,:,:,i] = tmpnuin


        if opt == 'Original':
            mob=v_elemcharge/(v_amu*nuin); # original form in Keely textbook 2009 p. 42
        if opt == 'ISR':
            mob = v_elemcharge/(v_amu*nuinISR)
        nuinScaler = 1.0 # scalar to increase value of ion-neutral collision frequency
        Babs = numpy.tile(Babs1,ne1.shape[0]).reshape(ne1.shape) # generate a time x beam x altitude array of magnetic field


        # loop over the ion type.
        # assume that since nuIn is weighted by fractional amount that
        # performs any weighting in the conductivity equation
        for i in range(mass.shape[0]):
            mob[:,:,:,i] = mob[:,:,:,i]/mass[i] # equivalent to the form on p 42 of Kelley 2009 textbook
            kappa = mob*1.0
            kappa[:,:,:,i] = mob[:,:,:,i]*Babs1 # equivalent to the form on p 42 of Kelley 2009 textbook

            # equation 2.40a in Kelley's textbook
            # note that the weighting is on the outside equivalent to Evans 1977, JGR
            # weighting is fractional concentration, so from 0-1.
            sp1 = fraction[:,:,:,i]*ne1*v_elemcharge*v_elemcharge/(v_amu*mass[i]*nuinScaler*nuin[:,:,:,i]*(1.0+(kappa[:,:,:,i]/nuinScaler)**2.0)) #*fraction[:,:,:,i]

            # equation 10 from Bostrom 1964 in my notebook
            # should be good > 85 km
            # note again the weighting on the outside
            sh1 = fraction[:,:,:,i]*ne1*v_elemcharge/(Babs*(1.0+(kappa[:,:,:,i]/nuinScaler)**2.0)) # fraction in effect weights things
            # print sp1
            PedCond = PedCond + sp1
            HallCond = HallCond + sh1


        HallConductance = numpy.zeros(HallCond.shape[0]) # time axis
        PedConductance = numpy.zeros(PedCond.shape[0])
        NeOut = numpy.zeros(PedCond.shape[0])
        FA_index = numpy.where((BeamCodes[:,1] == self.azFA) & (BeamCodes[:,2] ==self.elFA))[0][0]
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
            NeOut = ne1[:,FA_index,:]
            TimeUTHour = numpy.mean(dtime1, axis=1)
            UnixTime = numpy.mean(time1,axis=1)
            FA_Altitude = Altitude[FA_index,:]
        else:
            HallConductivity = []
            PedConductivity = []
            TimeUTHour = []
            UnixTime = []
            FA_Altitude = []
            NeOut = []

        return HallConductivity, PedConductivity,PedConductance, HallConductance,NeOut,\
               TimeUTHour,FA_Altitude, UnixTime


    def CalculateConductance(self,fname_ac,outLocation):

        # need to figure out how to dynamically allocate this...

        k = True
        for ifile in fname_ac:
            print ifile


            tH1,tP1,tP2,tH2,tNe,t,alt,tunixTime = self.ProcessFAConductivity(ifile)
            if k:
                print tH1.shape
                HallConductivity = numpy.zeros([1,50])
                HallConductance = numpy.zeros([0])
                HallConductanceISR = numpy.zeros([0])
                PedersenConductivity = numpy.zeros([1,50])
                NeOut = numpy.zeros([1,50])
                PedersenConductance = numpy.zeros([0])
                PedersenConductanceISR = numpy.zeros([0])
                TimeUTHour = numpy.zeros([0])
                UnixTime = numpy.zeros([0])
                Altitude = []
                k=False
            tmpHall = numpy.zeros([tH1.shape[0],50])
            tmpHall[:,0:tH1.shape[1]] = tH1
            HallConductivity = numpy.concatenate((HallConductivity,tmpHall),axis=0)
            tmpPed = numpy.zeros([tP1.shape[0],50])
            tmpPed[:,0:tP1.shape[1]] = tP1
            PedersenConductivity = numpy.concatenate((PedersenConductivity,tmpPed),axis=0)
            tmpNe = numpy.zeros([tNe.shape[0],50])
            tmpNe[:,0:tNe.shape[1]] = tNe
            NeOut = numpy.concatenate((NeOut,tmpNe),axis=0)
            PedersenConductance = numpy.concatenate((PedersenConductance,tP2), axis=0)
            HallConductance = numpy.concatenate((HallConductance,tH2), axis=0)
            TimeUTHour = numpy.concatenate((TimeUTHour,t), axis=0)
            UnixTime = numpy.concatenate((UnixTime,tunixTime), axis=0)

            tH1,tP1,tP2,tH2,tNe,t,alt,tunixTime = self.ProcessFAConductivity(ifile,opt='ISR')
            HallConductanceISR = numpy.concatenate((HallConductanceISR,tH2), axis=0)
            PedersenConductanceISR = numpy.concatenate((PedersenConductanceISR,tP2), axis=0)

        # remove first element of Hall conductivity
        HallConductivity = HallConductivity[1:,:]
        PedersenConductivity = PedersenConductivity[1:,:]
        NeOut = NeOut[1:,:]
        # HallConductance = HallConductance[1:]
        # PedersenConductance = PedersenConductance[1:]
        # HallConductanceISR = HallConductanceISR[1:]
        # PedersenConductanceISR = PedersenConductanceISR[1:]

        # sort out the data before printing.
        # sort out the times and just make sure everything is ordered by time
        #unix time should be increasing
        qsort = numpy.argsort(UnixTime)
        UnixTime = UnixTime[qsort]
        TimeUTHour = TimeUTHour[qsort]
        PedersenConductivity = PedersenConductivity[qsort,:]
        HallConductivity = HallConductivity[qsort,:]
        PedersenConductance = PedersenConductance[qsort]
        NeOut = NeOut[qsort,:]
        HallConductance = HallConductance[qsort]
        HallConductanceISR = HallConductanceISR[qsort]
        PedersenConductanceISR = PedersenConductanceISR[qsort]

        outDict = dict()
        outDict['UnixTime'] = UnixTime
        outDict['TimeUTHour'] = TimeUTHour
        outDict['PedersenConductivity'] = PedersenConductivity
        outDict['HallConductivity'] = HallConductivity
        outDict['Ne'] = NeOut
        outDict['HallConductance'] = HallConductanceISR
        outDict['PedersenConductance'] = PedersenConductanceISR



        # do some file name making.
        t0 = datetime.datetime.utcfromtimestamp(int(UnixTime[0]))
        startTimeStr = str(t0.year)+str(t0.month).zfill(2)+str(t0.day).zfill(2)+'T'+str(t0.hour).zfill(2)+str(t0.minute).zfill(2)+'UT'

        t1 = datetime.datetime.utcfromtimestamp(int(UnixTime[-1]))
        endTimeStr = str(t1.year)+str(t1.month).zfill(2)+str(t1.day).zfill(2)+'T'+str(t1.hour).zfill(2)+str(t1.minute).zfill(2)+'UT'

        outFile = self.location+'_Conductance_NoMod_'+startTimeStr+'_'+endTimeStr+'AFOSRYIP.pkl'
        outPng =  self.location+'_Conductance_NoMod_'+startTimeStr+'_'+endTimeStr+'AFOSRYIP.png'
        X = numpy.array([UnixTime,TimeUTHour,PedersenConductance,HallConductance, PedersenConductanceISR, HallConductanceISR])
        with open(outFile,'w') as f:
            pickle.dump(outDict,f)
        # with open(os.path.join(outLocation,outFile), 'wb') as f:
        #     f.write('############################################################################################# \n')
        #     f.write('Author: S.R. Kaeppler \n')
        #     f.write('Email: skaeppl@clemson.edu \n')
        #     f.write('Version number: v0.1\n')
        #     f.write('File Created: %s \n'%datetime.datetime.now())
        #     f.write('Valid Time: %s-%s\n'%(startTimeStr,endTimeStr))
        #
        #     f.write('\n')
        #     f.write('Columns (values of -1 indicate failed altitude integration) \n')
        #     f.write('UnixTime: seconds since 01 January 1970 00:00:00 UT, absolute time reference \n')
        #     f.write('Time: UT Decimal Hours (hours) \n ')
        #     f.write('Pedersen Conductance (mho) using MSIS Tn=Ti \n')
        #     f.write('Hall Conductance (mho) using MSIS Tn=Ti \n')
        #     f.write('Pedersen Conductance (mho) using ISR Ti and MSIS Tn \n')
        #     f.write('Hall Conductance (mho) using ISR Ti and MSIS Tn \n')
        #     f.write('############################################################################################# \n')
        #     numpy.savetxt(f,X.T, fmt='%d,%0.2f, %0.1f, %0.1f,%0.1f,%0.1f')


        # time gaps program
        plt.figure(dpi=300,figsize=(11,8.5))
        tmpT = (UnixTime-UnixTime[0])/3600.
        tmpT = numpy.arange(0,PedersenConductance.shape[0],1)
        plt.plot(tmpT, PedersenConductance, 'b.', label='P Org')
        plt.plot(tmpT, HallConductance, 'r.', label='H Org')
        plt.plot(tmpT, PedersenConductanceISR, 'bx', label='P ISR')
        plt.plot(tmpT, HallConductanceISR, 'rx', label='H ISR')
        plt.xlabel('Time Hours after T0')
        plt.ylabel('Conductance (mho)')
        plt.legend()
        plt.title('%s %s-%s'%(self.location,startTimeStr,endTimeStr))
        plt.tight_layout()
        plt.savefig(os.path.join(outLocation,outPng))
        plt.close()
        # plt.show()




if __name__ == '__main__':
    conduct = Conductivity()


    # 17 March 2013
    # fname_ac=['/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2013/03/PINOT_Daytime31/20130316.006.done/20130316.006_ac_5min-cal.h5',\
    #           '/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2013/03/PINOT_Daytime31/20130317.005.done/20130317.005_ac_5min-cal.h5', \
    #           '/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2013/03/PINOT_Nighttime31/20130317.003/20130317.003_ac_5min-cal.h5',\
    #           '/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2013/03/PINOT_Dregion31/20130317.001.done/20130317.001_ac_5min-cal.h5',\
    #           '/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2013/03/PINOT_Dregion31/20130317.011.done/20130317.011_ac_5min-cal.h5']
    #     #   './20170317.005/20170317.005_ac_3min-fitcal.h5',\
    #     #   './20170317.007/20170317.007_ac_3min-fitcal.h5']
    fname_ac = ['/Users/srkaeppler/research/data/pfisr/conductance_allsky/events/20121106/ISR/20121106.004_ac_3min-cal.h5']


    conduct.CalculateConductance(fname_ac, './')

    # # June 21-24, 2015
    # fname_ac =['/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2015/06/IPY27_Tracking_v03/20150621.002/20150621.002_ac_5min-fitcal.h5',\
    #            '/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2015/06/IPY27_Tracking_v03/20150622.002/20150622.002_ac_5min-fitcal.h5', \
    #            '/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2015/06/IPY27_Tracking_v03/20150623.001/20150623.001_ac_5min-fitcal.h5', \
    #            '/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2015/06/IPY27_Tracking_v03/20150624.003/20150624.003_ac_5min-fitcal.h5', \
    #            '/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2015/06/IPY27_Tracking_v03/20150625.001/20150625.001_ac_5min-fitcal.h5', \
    #            '/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2015/06/WorldDay35/20150622.004/20150622.004_ac_3min-fitcal.h5',\
    #            '/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2015/06/WorldDay35/20150623.002/20150623.002_ac_3min-fitcal.h5',\
    #            '/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2015/06/WorldDay35/20150624.002/20150624.002_ac_3min-fitcal.h5', \
    #            '/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2015/06/WorldDay35/20150624.004/20150624.004_ac_3min-fitcal.h5']
    # conduct.CalculateConductance(fname_ac, './')
    #
    # #Oct 13-15, 2016
    # fname_ac = ['/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2016/10/IPY27_Tracking_v03/20161012.002/20161012.002_ac_5min-fitcal.h5', \
    #             '/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2016/10/IPY27_Tracking_v03/20161012.005/20161012.005_ac_5min-fitcal.h5', \
    #             '/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2016/10/IPY27_Tracking_v03/20161013.002/20161013.002_ac_5min-fitcal.h5', \
    #             '/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2016/10/IPY27_Tracking_v03/20161015.006/20161015.006_ac_5min-fitcal.h5', \
    #             '/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2016/10/WorldDay35/20161013.004/20161013.004_ac_3min-fitcal.h5', \
    #             '/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2016/10/WorldDay35/20161014.002/20161014.002_ac_3min-fitcal.h5', \
    #             '/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2016/10/WorldDay35/20161014.005/20161014.005_ac_3min-fitcal.h5', \
    #             '/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2016/10/WorldDay35/20161015.002/20161015.002_ac_3min-fitcal.h5', \
    #             '/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2016/10/WorldDay35/20161015.005/20161015.005_ac_3min-fitcal.h5',\
    #             '/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2016/10/MSWinds23/20161013.003/20161013.003_ac_5min-fitcal.h5']
    # conduct.CalculateConductance(fname_ac, './')
