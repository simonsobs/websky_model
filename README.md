# Websky Models

websky_models is a python library for dealing with websky maps and halo catalogues

Some functionalities require healpy or astropy. astropy is only used if comoving distance to redshift calculations required for halo catalogue (if zmin/zmax are used for halo catalogue, or practical=True). Should work without astropy installed if not using these.

## Usage

```python
import numpy as np
import websky_model as wm

# load model
wmodel = wm.WebSky(directory_path='/global/project/projectdirs/sobs/v4_sims/mbs/websky/data/',websky_version = 'v0.1', verbose=True)

# read in halo catalogue
hcat  = wmodel.load_halo_catalogue(mmin=0, mmax=np.inf, zmin=0, zmax=np.inf, rmin=0., rmax=np.inf, practical=True)

# project halos onto healpix map. 
# each halo is given a "flux" of weight=constant, or weight=array(len(Nhalo))
hpmap = wmodel.catalogue_to_map(hcat, nside=512, weight=1)

# get cib map filename
fname_cib = wmodel.cib_map_file_name(freq='545')

# get kappa map filename
fname_kappa = wmodel.kappa_map_file_name()

# get compton-y map filename
fname_comptony = wmodel.comptony_map_file_name()
```

