from so_pysm_models.interpolating import InterpolatingComponent
from . import websky

class WebSkyCIB(InterpolatingComponent):
    """PySM component interpolating between precomputed maps"""
    def __init__(self, path, input_units="MJysr", target_nside=4096, interpolation_kind="linear",
                 pixel_indices=None,
                 mpi_comm=None, verbose=False):
        InterpolatingComponent.__init__(self, path, input_units, target_nside, interpolation_kind,
                                    False, pixel_indices,
                                    mpi_comm, verbose)

    def get_filenames(self,path):
        # Override this to implement name convention
        ws = websky.WebSky(path)
        ifreqs = [27,39,93,144,225,280]
        freqs = []
        for i in range(len(ifreqs)):
            freqs.append(ifreqs[i]-1)
            freqs.append(ifreqs[i])
            freqs.append(ifreqs[i]+1)
        freqs = freqs + [100,217,353,545,857]
        fnames = {}
        for freq in freqs:
            fnames[freq] = ws.cib_map_file_name(str(freq))
        return fnames
