from FastConductivityBC import Conductivity
from FastConductivityNe_noTR import Conductivity as Conductivity2
import datetime
import glob

if __name__ == '__main__':
    conduct = Conductivity()
    conductNe_noTR = Conductivity2()


    # 30Jan2016_4Feb2016
    dt1970 = datetime.datetime(1970,01,01,00,00,00)
    dtInitial = datetime.datetime(2010,2,16,00,00,00)
    dtFinal = datetime.datetime(2010,2,17,00,00,00)
    tInitial = (dtInitial-dt1970).total_seconds()
    tFinal = (dtFinal-dt1970).total_seconds()

    timeInterval = [tInitial,tFinal]
    # test case
    # fname_ac = ['/data0/LongTermStorage/research/data/pfisr/conductance_allsky/events/20121106/ISR/20121106.004_ac_3min-cal.h5']

    fname_ac = sorted(glob.glob('/Users/srkaeppler/Dropbox/research/data/CCMC_Conductivites/16Feb2010/20100216.001.done/20100216.001_bc_2min-Ne-cal.h5'))
    conduct.CalculateBCConductance(fname_ac, './ChristineEvents/20100216.001_bc_2min-Ne-cal_el89_az15_04272021/', timeInterval=timeInterval)

    fname_ac = sorted(glob.glob('/Users/srkaeppler/Dropbox/research/data/CCMC_Conductivites/16Feb2010/20100216.001_bc_nenotr_1min.h5'))
    conductNe_noTR.CalculateBCConductance(fname_ac, './ChristineEvents/20100216.001_bc_nenotr_1min_el89_az15_04272021/', timeInterval=timeInterval)

    fname_ac = sorted(glob.glob('/Users/srkaeppler/Dropbox/research/data/CCMC_Conductivites/16Feb2010/20100216.001_bc_nenotr_30sec.h5'))
    conductNe_noTR.CalculateBCConductance(fname_ac, './ChristineEvents/20100216.001_bc_nenotr_30sec_el89_az15_04272021/', timeInterval=timeInterval)


    """
        #
        # # Feb 20, 2014.
        # dt1970 = datetime.datetime(1970,01,01,00,00,00)
        # dtInitial = datetime.datetime(2014,2,20,00,00,00)
        # dtFinal = datetime.datetime(2014,2,21,00,00,00)
        # tInitial = (dtInitial-dt1970).total_seconds()
        # tFinal = (dtFinal-dt1970).total_seconds()
        #
        # timeInterval = [tInitial,tFinal]
        # # test case
        # # fname_ac = ['/data0/LongTermStorage/research/data/pfisr/conductance_allsky/events/20121106/ISR/20121106.004_ac_3min-cal.h5']
        # fname_ac = sorted(glob.glob('/Users/srkaeppler/Dropbox/research/data/CCMC_Conductivites/20Feb2014/*.h5'))
        # conduct.CalculateConductance(fname_ac, './ChristineEvents', timeInterval=timeInterval)
        #
        #
        #
        # # 30Jan2016_4Feb2016
        # dt1970 = datetime.datetime(1970,01,01,00,00,00)
        # dtInitial = datetime.datetime(2016,1,30,00,00,00)
        # dtFinal = datetime.datetime(2016,2,4,00,00,00)
        # tInitial = (dtInitial-dt1970).total_seconds()
        # tFinal = (dtFinal-dt1970).total_seconds()
        #
        # timeInterval = [tInitial,tFinal]
        # # test case
        # # fname_ac = ['/data0/LongTermStorage/research/data/pfisr/conductance_allsky/events/20121106/ISR/20121106.004_ac_3min-cal.h5']
        # fname_ac = sorted(glob.glob('/Users/srkaeppler/Dropbox/research/data/CCMC_Conductivites/30Jan2016_4Feb2016/*.h5'))
        # conduct.CalculateConductance(fname_ac, './ChristineEvents', timeInterval=timeInterval)
    # April 5-6 2010
    dtInitial = datetime.datetime(2010,04,05,00,00,00)
    dtFinal = datetime.datetime(2010,04,07,00,00,00)
    tInitial = (dtInitial-dt1970).total_seconds()
    tFinal = (dtFinal-dt1970).total_seconds()

    timeInterval = [tInitial,tFinal]
    # test case
    fname_ac = ['/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2010/04/Lyons30/20100404.001/20100404.001_ac_3min.h5',
                '/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2010/04/Lyons30/20100405.001/20100405.001_ac_3min.h5',
                '/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2010/04/Lyons30/20100406.001/20100406.001_ac_3min.h5',
                '/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2010/04/IPY17/20100404.002/20100404.002_ac_5min-fitcal.h5',
                '/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2010/04/IPY17/20100405.002/20100405.002_ac_5min-fitcal.h5',
                '/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2010/04/IPY17/20100405.005/20100405.005_ac_5min-fitcal.h5',
                '/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2010/04/IPY17/20100406.002/20100406.002_ac_5min-fitcal.h5',
                '/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2010/04/IPY17/20100407.001/20100407.001_ac_5min-fitcal.h5',
                '/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2010/04/Erickson30/20100405.004/20100405.004_ac_3min-fitcal.h5'
                ]
    conduct.CalculateConductance(fname_ac, './PFISR-AMPERE_Events', timeInterval=timeInterval)



    # August 3-4 2010
    dtInitial = datetime.datetime(2010,8,03,00,00,00)
    dtFinal = datetime.datetime(2010,8,05,00,00,00)
    tInitial = (dtInitial-dt1970).total_seconds()
    tFinal = (dtFinal-dt1970).total_seconds()

    timeInterval = [tInitial,tFinal]
    # test case
    fname_ac = ['/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2010/08/IPY17/20100803.001.done/20100803.001_ac_5min-cal.h5',
                '/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2010/08/IPY17/20100805.002.done/20100805.002_ac_5min-cal.h5',
                '/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2010/08/MSWinds23/20100804.001.done/20100804.001_ac_5min-cal.h5',
                '/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2010/08/MSWinds23/20100805.001.done/20100805.001_ac_5min-cal.h5',
                '/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2010/08/Erickson30/20100804.004.done/20100804.004_ac_5min-cal.h5'
                ]
    conduct.CalculateConductance(fname_ac, './PFISR-AMPERE_Events', timeInterval=timeInterval)


    # Feb 4 2011
    dtInitial = datetime.datetime(2011,2,04,00,00,00)
    dtFinal = datetime.datetime(2011,2,05,00,00,00)
    tInitial = (dtInitial-dt1970).total_seconds()
    tFinal = (dtFinal-dt1970).total_seconds()

    timeInterval = [tInitial,tFinal]
    # test case
    fname_ac = ['/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2011/02/WorldDay30/20110201.002.done/20110201.002_ac_5min-cal.h5'
                ]
    conduct.CalculateConductance(fname_ac, './PFISR-AMPERE_Events', timeInterval=timeInterval)



    # March 10-11 2011
    dtInitial = datetime.datetime(2011,3,10,00,00,00)
    dtFinal = datetime.datetime(2011,3,12,00,00,00)
    tInitial = (dtInitial-dt1970).total_seconds()
    tFinal = (dtFinal-dt1970).total_seconds()

    timeInterval = [tInitial,tFinal]
    # test case
    fname_ac = ['/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2011/03/IPY17/20110306.002.done/20110306.002_ac_5min-cal.h5'
                ]
    conduct.CalculateConductance(fname_ac, './PFISR-AMPERE_Events', timeInterval=timeInterval)


    # May 28-29 2011
    dtInitial = datetime.datetime(2011,5,28,00,00,00)
    dtFinal = datetime.datetime(2011,5,30,00,00,00)
    tInitial = (dtInitial-dt1970).total_seconds()
    tFinal = (dtFinal-dt1970).total_seconds()

    timeInterval = [tInitial,tFinal]
    # test case
    fname_ac = ['/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2011/05/IPY17/20110527.001.done/20110527.001_ac_5min-cal.h5'
                ]
    conduct.CalculateConductance(fname_ac, './PFISR-AMPERE_Events', timeInterval=timeInterval)


    # June 4-5 2011
    dtInitial = datetime.datetime(2011,6,4,00,00,00)
    dtFinal = datetime.datetime(2011,6,6,00,00,00)
    tInitial = (dtInitial-dt1970).total_seconds()
    tFinal = (dtFinal-dt1970).total_seconds()

    timeInterval = [tInitial,tFinal]
    # test case
    fname_ac = ['/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2011/06/IPY17/20110601.002.done/20110601.002_ac_5min-cal.h5'
                ]
    conduct.CalculateConductance(fname_ac, './PFISR-AMPERE_Events', timeInterval=timeInterval)




    # September 17-18 2011
    dtInitial = datetime.datetime(2011,9,17,00,00,00)
    dtFinal = datetime.datetime(2011,9,18,00,00,00)
    tInitial = (dtInitial-dt1970).total_seconds()
    tFinal = (dtFinal-dt1970).total_seconds()

    timeInterval = [tInitial,tFinal]
    # test case
    fname_ac = ['/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2011/09/IPY17/20110913.001.done/20110913.001_ac_5min-cal.h5',
                '/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2011/09/IPY17/20110913.003.done/20110913.003_ac_5min-cal.h5'
                ]
    conduct.CalculateConductance(fname_ac, './PFISR-AMPERE_Events', timeInterval=timeInterval)


    # April 24-26 2012
    dtInitial = datetime.datetime(2012,4,24,00,00,00)
    dtFinal = datetime.datetime(2012,4,27,00,00,00)
    tInitial = (dtInitial-dt1970).total_seconds()
    tFinal = (dtFinal-dt1970).total_seconds()

    timeInterval = [tInitial,tFinal]
    # test case
    fname_ac = ['/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2012/04/IPY17/20120422.002.done/20120422.002_ac_5min-cal.h5',
                '/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2012/04/IPY17/20120427.001.done/20120427.001_ac_5min-cal.h5'
                ]
    conduct.CalculateConductance(fname_ac, './PFISR-AMPERE_Events', timeInterval=timeInterval)



    # June 16-18 2012
    dtInitial = datetime.datetime(2012,6,16,00,00,00)
    dtFinal = datetime.datetime(2012,6,19,00,00,00)
    tInitial = (dtInitial-dt1970).total_seconds()
    tFinal = (dtFinal-dt1970).total_seconds()

    timeInterval = [tInitial,tFinal]
    # test case
    fname_ac = ['/media/srk/KaepplerAMISRProcessed/AMISR_PROCESSED/processed_data/PFISR/2012/06/IPY17/20120615.001.done/20120615.001_ac_5min-cal.h5'
                ]
    conduct.CalculateConductance(fname_ac, './PFISR-AMPERE_Events', timeInterval=timeInterval)
    """
