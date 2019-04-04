from __future__ import print_function
import numpy as np
import healpy as hp
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
                 directory_path='/global/project/projectdirs/sobs/v4_sims/mbs/websky/data/',
                 websky_version = 'v0.3',
                 halo_catalogue = 'halos.pksc',
                 kappa_map_name = 'kap.fits',
                 comptony_map_name = 'tsz.fits',
                 ksz_map_name = 'ksz.fits',
                 websky_cosmo = {'Omega_M': 0.31, 'Omega_B': 0.049, 'Omega_L': 0.69, 
                                 'h': 0.68, 'sigma_8': 0.81, 'n_s':0.965},
                 verbose = True
    ):


        self.directory_path = directory_path 
        self.websky_version = websky_version 
        self.halo_catalogue = halo_catalogue
        self.websky_cosmo   = websky_cosmo
        self.verbose        = verbose
        self.kappa_map_name = kappa_map_name
        self.comptony_map_name = comptony_map_name
        self.ksz_map_name = ksz_map_name

    def import_astropy(self):
        """load in astropy. Only used if comoving distance to redshift calculations required
        """

        from astropy.cosmology import FlatLambdaCDM, z_at_value
        import astropy.units as u

        self.z_at_value = z_at_value
        self.u = u
        self.astropy_cosmo = FlatLambdaCDM(H0=self.websky_cosmo['h']*100, Om0=self.websky_cosmo['Omega_M'])

    def load_halo_catalogue(self, mmin=0., mmax=np.inf, zmin=0., zmax=np.inf, rmin=0., rmax=np.inf, practical=True):
        """load in peak patch dark matter halo catalogue

        Requires astropy if using distance to redshift calculations, or redshift cuts

        Returns
        -------

        if practical==True: only generally useful information (including redshifts)
        halodata : np.array((Nhalo, 8))
            numpy array of halo information, 10 floats per halo
            x [Mpc], y [Mpc], z [Mpc], vx [km/s], vy [km/s], vz [km/s], 
            M [M_sun (M_200,m)], redshift (v_pec not included)

        if practical==False: returns everything in halo catalogue
        halodata : np.array((Nhalo, 10+))
            numpy array of halo information, 10 floats per halo
            x [Mpc], y [Mpc], z [Mpc], vx [km/s], vy [km/s], vz [km/s], 
            M [M_sun (M_200,m)], x_lag [Mpc], y_lag [Mpc], z_lag [Mpc], ...



        """

        halo_catalogue_file = open(self.directory_path+self.websky_version+'/'+self.halo_catalogue,"rb")
        
        # load catalogue header
        Nhalo            = np.fromfile(halo_catalogue_file, dtype=np.int32, count=1)[0]
        RTHMAXin         = np.fromfile(halo_catalogue_file, dtype=np.float32, count=1)
        redshiftbox      = np.fromfile(halo_catalogue_file, dtype=np.float32, count=1)
        if self.verbose: print("\nNumber of Halos in full catalogue %d \n " % Nhalo)

        nfloats_perhalo = 10
        npkdata         = nfloats_perhalo*Nhalo

        # load catalogue data
        halodata        = np.fromfile(halo_catalogue_file, dtype=np.float32, count=npkdata)
        halodata        = np.reshape(halodata,(Nhalo,nfloats_perhalo))

        # change from R_th to halo mass (M_200,M)
        rho_mean = 2.775e11 * self.websky_cosmo['Omega_M'] * self.websky_cosmo['h']**2
        halodata[:,6] = 4.0/3*np.pi * halodata[:,6]**3 * rho_mean        
        
        # cut mass range
        if mmin > 0 or mmax < np.inf:
            dm = (halodata[:,6] > mmin) & (halodata[:,6] < mmax) 
            halodata = halodata[dm]

        # cut redshift range
        if zmin > 0 or zmax < np.inf:
            self.import_astropy()
            rofzmin = self.astropy_cosmo.comoving_distance(zmin).value
            rofzmax = self.astropy_cosmo.comoving_distance(zmax).value

            rpp =  np.sqrt( np.sum(halodata[:,:3]**2, axis=1))

            dm = (rpp > rofzmin) & (rpp < rofzmax) 
            halodata = halodata[dm]

        # cut distance range
        if rmin > 0 or rmax < np.inf:
            rpp =  np.sqrt(np.sum(halodata[:,:3]**2, axis=1))

            dm = (rpp > rmin) & (rpp < rmax) 
            halodata = halodata[dm]


        # get halo redshifts and crop all non practical information
        if practical:
            self.import_astropy()

            # set up comoving distance to redshift interpolation table
            rpp =  np.sqrt( np.sum(halodata[:,:3]**2, axis=1))

            zminh = self.z_at_value(self.astropy_cosmo.comoving_distance, rpp.min()*self.u.Mpc)
            zmaxh = self.z_at_value(self.astropy_cosmo.comoving_distance, rpp.max()*self.u.Mpc)
            zgrid = np.linspace(zminh, zmaxh, 10000)
            dgrid = self.astropy_cosmo.comoving_distance(zgrid).value
            
            # first replace 7th column with redshift
            halodata[:,7] = np.interp(rpp, dgrid, zgrid)
            
            # crop un-practical halo information
            halodata = halodata[:,:8]

        Nhalo = halodata.shape[0] 
        Nfloat_perhalo = halodata.shape[1]

        # write out halo catalogue information
        if self.verbose: 
            print("Halo catalogue after cuts: np.array((Nhalo=%d, floats_per_halo=%d)), containing:\n" % (Nhalo, Nfloat_perhalo))
            if practical:
                print("0:x [Mpc], 1:y [Mpc], 2:z [Mpc], 3:vx [km/s], 4:vy [km/s], 5:vz [km/s],\n"+ 
                      "6:M [M_sun (M_200,m)], 7:redshift(chi_halo) \n")
            else:
                print("0:x [Mpc], 1:y [Mpc], 2:z [Mpc], 3:vx [km/s], 4:vy [km/s], 5:vz [km/s],\n"+ 
                      "6:M [M_sun (M_200,m)], 7:x_lag [Mpc], 8:y_lag [Mpc], 9:z_lag [Mpc]\n")

        return halodata


    def catalogue_to_map(self, halodata, nside=512, weight=1):
        """project halos into healpix map of size nside

        Parameters
        ----------

        halodata : array
            halo catalogue of size (Nhalo, 10)
        nside: int
            healpix nside of map created
        weight: str
            weighting to give halos when adding to map
            possibilities: 1=number density, or array of size Nhalo for custom (e.g. mass) 

        Returns
        -------

        map : np.array( (12 * nside**2))
            healpy map of halos, with "flux" proportional to weight

        """
        # create empty healpy map
        map = np.zeros(hp.nside2npix(nside))
        
        # get pixel id from halo x,y,z vector
        pix = hp.vec2pix(nside, halodata[:,0], halodata[:,1], halodata[:,2])

        # add flux to map 
        np.add.at(map, pix, weight)

        return map

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

        cib_file_name = 'cib_'+str(freq).zfill(4)+ '.fits'

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


    def ksz_map_file_name(self):
        """get file name of ksz map

        Returns
        -------

        ksz_file_name : str
            name of ksz map file 
        """
        
        return self.directory_path+self.websky_version+'/'+self.ksz_map_name

