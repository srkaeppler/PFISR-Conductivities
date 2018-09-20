from FastConductivity import Conductivity
import datetime

if __name__ == '__main__':
    conduct = Conductivity()

    dt1970 = datetime.datetime(1970,01,01,00,00,00)
    dtInitial = datetime.datetime(2012,11,06,10,00,00)
    dtFinal = datetime.datetime(2012,11,06,14,00,00)
    tInitial = (dtInitial-dt1970).total_seconds()
    tFinal = (dtFinal-dt1970).total_seconds()

    timeInterval = [tInitial,tFinal]
    # test case
    fname_ac = ['/data0/LongTermStorage/research/data/pfisr/conductance_allsky/events/20121106/ISR/20121106.004_ac_3min-cal.h5']
    conduct.CalculateConductance(fname_ac, './PFISR-AMPERE_Events', timeInterval=timeInterval)


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


    
