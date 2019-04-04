from websky_model import cib
import healpy as hp
import matplotlib.pyplot as plt

dirpath = "/global/project/projectdirs/sobs/v4_sims/mbs/websky/data/"

wcib = cib.WebSkyCIB(dirpath,verbose=True)
cib_150 = wcib.signal(150.)

hp.mollview(cib_150[0])
plt.savefig("cib.png")
