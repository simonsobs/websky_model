from __future__ import print_function
import numpy as np
import sys


class WebSky:
    """class for websky maps and dark matter halo catalogues

    Parameters
    ----------

    directory_path : str
        directory of websky halo catalogues and maps. Default to nersc
    halo_catalogue : str
        directory of websky halo catalogue. Default to nersc convention
    websky_cosmo : dictionary
        useful cosmological parameters of the websky run
    verbose : bool
        prints information on halo catalogues and maps
    """

    def __init__(self, 
                 directory_path='/global/project/projectdirs/sobs/v4_sims/mbs/websky/v0/',
                 websky_version = 'v0',
                 halo_catalogue = 'halos.pksc',
                 kappa_map_name = 'kappa.fits',
                 comptony_map_name = 'compton-y.fits',
                 websky_cosmo = {'Omega_M': 0.31, 'Omega_B': 0.049, 'Omega_L': 0.69, 
                                 'h': 0.68, 'sigma_8': 0.81, 'n_s':0.965},
                 verbose = True
    ):


        self.directory_path = directory_path 
        self.websky_version = websky_version 
        self.halo_catalogue = halo_catalogue
        self.websky_cosmo   = websky_cosmo
        self.verbose        = verbose

    def load_halo_catalogue(self, mmin=1e0, mmax=1e30, rmin=0., rmax=14.e3):
        """load in peak patch dark matter halo catalogue

        Returns
        -------

        halodata : np.array((10, Nhalo)
            numpy array of halo information, 10 floats per halo
            x [Mpc], y [Mpc], z [Mpc], vx [km/s], vy [km/s], vz [km/s], 
            M [M_sun (M_200,m)], x_lag [Mpc], y_lag [Mpc], z_lag [Mpc]
        """

        rho_mean = 2.775e11 * self.websky_cosmo['Omega_M'] * self.websky_cosmo['h']**2
        halo_catalogue_file = open(self.directory_path+self.halo_catalogue,"rb")
        
        # load catalogue header
        Nhalo            = np.fromfile(halo_catalogue_file, dtype=np.int32, count=1)[0]
        RTHMAXin         = np.fromfile(halo_catalogue_file, dtype=np.float32, count=1)
        redshiftbox      = np.fromfile(halo_catalogue_file, dtype=np.float32, count=1)
        if self.verbose: print("\nNhalo = ", Nhalo)

        nfloats_perhalo = 10
        npkdata         = nfloats_perhalo*Nhalo

        # load catalogue data
        halodata        = np.fromfile(halo_catalogue_file, dtype=np.float32, count=npkdata)
        halodata        = np.reshape(halodata,(Nhalo,nfloats_perhalo))

        # change from R_th to halo mass (M_200,M)
        halodata[:,6]   = 4.0/3*np.pi * halodata[:,6]**3 * rho_mean        
        
        # cut mass range
        dm = [ (halodata[:,6] > mmin) & (halodata[:,6] < mmax)]
        halodata = halodata[dm]

        # cut distance range
        rpp = halodata[:,0]**2 + halodata[:,1]**2 + halodata[:,2]**2
        dm = [ (rpp > rmin**2) & (rpp < rmax**2)]
        halodata = halodata[dm]

        if self.verbose: print("halo catalogue loaded: min, max mass /1e12 = ", np.min(halodata[:,6])/1e12, np.max(halodata[:,6])/1e12 )
        
        return halodata


    def cib_map_file_name(self, freq='545'):
        """get file name of cib map, given a frequency

        Parameters
        ----------

        freq : str or int
            frequency of desired map in GHz

        Returns
        -------

        cib_file_name : str
            name of cib file at given frequency
        """

        cib_file_name = 'cib-'+str(freq)+'GHZ_v'+str(websky_version)+'.fits'

        return self.directory_path+self.websky_version+'/'+cib_file_name

    def kappa_map_file_name(self):
        """get file name of kappa map

        Returns
        -------

        kappa_file_name : str
            name of kappa map file 
        """

        return self.directory_path+self.websky_version+'/'+self.kappa_map_name


    def comptony_map_file_name(self):
        """get file name of compton-y map

        Returns
        -------

        comptony_file_name : str
            name of compton-y map file 
        """
        
        return self.directory_path+self.websky_version+'/'+self.comptony_map_name


