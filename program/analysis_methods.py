import numpy as np
from astropy.io import fits


def _Phase_Calculator(self, time_array, MJD_ref):
  
    MJD_times = MJD_ref + (time_array / 86400.0)

    
    dt_sec = (MJD_times - self.T0) * 86400.0

    P0_sec = self.orbital_period      
    Pdot   = self.orbital_pdot        

    if Pdot == 0.0:
        return (dt_sec / P0_sec) % 1.0

  
    cycles = (dt_sec / P0_sec) - 0.5 * (Pdot / (P0_sec * P0_sec)) * (dt_sec * dt_sec)
    return cycles % 1.0


def _Phase_Binner(self, phase, rate, error):
    
    bins = np.linspace(0.0, 1.0, self.nbins + 1)
    centres = 0.5 * (bins[:-1] + bins[1:])
    digitized = np.digitize(phase, bins) - 1

    binned_rate = np.full(self.nbins, np.nan)
    binned_err  = np.full(self.nbins, np.nan)

    for i in range(self.nbins):
        mask = digitized == i
        if np.any(mask):
            w = 1.0 / error[mask]**2
            binned_rate[i] = np.average(rate[mask], weights=w)
            binned_err[i]  = np.sqrt(1.0 / np.sum(w))

    return centres, binned_rate, binned_err


def _Epoch_Binning(self):
    self.obsid_epoch = {}

    self._Refresh_ObsID()

    for obsid in self.ObsID_array:
        paths = self._ObsID_Paths(obsid)
        fits_file = self._Light_Curve_File(obsid)

        if not fits_file.exists():
            continue

        with fits.open(fits_file) as hdul:
            header = hdul[1].header

            if "MJDREFI" in header and "MJDREFF" in header:
                mjdref = header["MJDREFI"] + header["MJDREFF"]
                tstart = header["TSTART"]
                mjd = mjdref + tstart / 86400.0

                self.obsid_epoch[obsid] = mjd


def _Flare_Filtering(self, rate, median, std): #Simple flare filter model
    three_sigma = 3 * std
    flares = (rate - median) > three_sigma
    below = (median - rate) > three_sigma
    normal = ~(flares | below)

    return flares, below, normal


def _Comparison_Model(self, phase, rate, err): #Simple comparison model, can be updated with a proper published model


    m = np.isfinite(phase) & np.isfinite(rate) & np.isfinite(err) & (err > 0)

    x = phase[m]
    y = rate[m]
    s = err[m]

    if len(y) < 4:
        return None

    ang = 2*np.pi*x

    X = np.column_stack([
        np.ones(len(x)),
        np.sin(ang),
        np.cos(ang)
    ])

    Xw = X / s[:,None]
    yw = y / s

    beta, _, _, _ = np.linalg.lstsq(Xw, yw, rcond=None)

    C, a, b = beta

    A = np.sqrt(a*a + b*b)

    delta = np.arctan2(b,a)
    phi0 = (-delta/(2*np.pi)) % 1.0

    yfit = C + a*np.sin(ang) + b*np.cos(ang)

    chi2 = np.sum(((y - yfit)/s)**2)
    dof = len(y) - 3

    return C, A, phi0, chi2, dof


def _Light_Curve_File(self, obsid):
    paths = self._ObsID_Paths(obsid)

    if self.use_nicerl3:

        lc_bary = paths["event_cl"] / f"ni{obsid}mpu7_sr_bary.lc"
        lc_raw  = paths["event_cl"] / f"ni{obsid}mpu7_sr.lc"
        return lc_bary if lc_bary.exists() else lc_raw
    else:

        return paths["xti_dir"] / f"{obsid}_bary.fits"


def _Load_Background_File(self, obsid):
    file = self.images / "Background_Files" / f"netlc_{obsid}.npz"
    paths = self._ObsID_Paths(obsid)

    
    if (paths["event_cl"] / "nibackgen3C50.SKIPPED_218").exists() or \
       (paths["event_cl"] / "nibackgen3C50.SKIPPED_255_GTI").exists():
        raise RuntimeError(f"{obsid} skipped: no valid background/GTI.")

    if not file.exists():
        if self.auto_resolve:
            self.Background_Subtraction()

        if not file.exists():
            raise RuntimeError(f"{obsid}: background-subtracted lightcurve not available.")

    d = np.load(file)
    return d["t"], d["rate"], d["error"]
