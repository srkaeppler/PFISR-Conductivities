from FastConductivity import Conductivity
import datetime
import glob

if __name__ == '__main__':
    conduct = Conductivity()

    # Feb 20, 2014.
    dt1970 = datetime.datetime(1970,01,01,00,00,00)
    dtInitial = datetime.datetime(2017,9,5,00,00,00)
    dtFinal = datetime.datetime(2017,9,13,00,00,00)
    tInitial = (dtInitial-dt1970).total_seconds()
    tFinal = (dtFinal-dt1970).total_seconds()

    timeInterval = [tInitial,tFinal]
    # test case
    # fname_ac = ['/data0/LongTermStorage/research/data/pfisr/conductance_allsky/events/20121106/ISR/20121106.004_ac_3min-cal.h5']
    fname_ac = sorted(glob.glob('/Users/srkaeppler/Dropbox/research/data/CCMC_Conductivites/092017_Lay/*.h5'))
    conduct.CalculateConductance(fname_ac, './092017_Lay', timeInterval=timeInterval)
