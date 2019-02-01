# Websky Models

websky_models is a Python library for dealing with websky maps and halo catalogues


## Usage

```python
import websky_model as wm

# load model
wmodel = wm.WebSky(directory_path='/global/project/projectdirs/sobs/v4_sims/mbs/websky/data/',websky_version = 'v0')

# read in halo catalogue
hcat  = wmodel.load_halo_catalogue(mmin=1e0, mmax=1e30, rmin=0., rmax=14.e3)

# get cib map filename
fname_cib = wmodel.cib_map_file_name(freq='545')

# get kappa map filename
fname_kappa = wmodel.kappa_map_file_name()

# get compton-y map filename
fname_comptony = wmodel.comptony_map_file_name()
```

