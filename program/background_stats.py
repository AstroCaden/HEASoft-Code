import numpy as np
from astropy.io import fits
import os
import sys
import subprocess


def Background_Subtraction(self):

    if self.use_nicerl3:
        self.logger.info("L3 is active so Background_Subtraction is skipped")
        return

    self._Refresh_ObsID()
    self.images = self.star_folder / "Outputs"
    self.images.mkdir(parents=True, exist_ok=True)

    bkg_dir = self.images / "Background_Files"
    bkg_dir.mkdir(parents=True, exist_ok=True)

    for obsid in self.ObsID_array:

        light_curve_file = self._Light_Curve_File(obsid)
        paths = self._ObsID_Paths(obsid)

        background_file = paths["event_cl"] / "nibackgen3C50_bkg.pi"

        if not light_curve_file.exists():
            if self.auto_resolve:
                self.Barycorr_Curve(self.lower_energy_limit, self.upper_energy_limit)
            else:
                sys.exit(0)

        if not background_file.exists():

            evt_file = paths["event_cl"] / f"ni{obsid}_0mpu7_cl_bary.evt"

            if not evt_file.exists():
                if self.auto_resolve:
                    self.Barycorr()
                else:
                    sys.exit(0)

            cmd_nibackgen = self.base_dir / "CMD_files" / "cmd_nibackgen3C50.sh"

            result = subprocess.run(
                [
                    cmd_nibackgen.as_posix(),
                    self.obsids_dir.as_posix(),
                    str(obsid),
                ],
                cwd=paths["event_cl"],
                capture_output=True,
                text=True,
                check=False,
            )

            if result.returncode == 0:
                self.logger.info(f"Background spectrum generated for {obsid}")

            elif result.returncode == 218:
                self._Failed_ObsID(obsid, "nibackgen3C50: no good exposure after cuts (218)", where="background")
                continue

            elif result.returncode == 255 and ("GTI re-definition FAILED" in (result.stdout + result.stderr)):
                self._Failed_ObsID(obsid, "nibackgen3C50: GTI re-definition failed / zero exposure (255)", where="background")
                continue

            elif result.returncode != 0:
                reason = (result.stderr or result.stdout or "").strip()
                self._Failed_ObsID(obsid, f"nibackgen3C50 failed ({result.returncode})\n{reason}", where="background")
                continue

        output_file = bkg_dir / f"netlc_{obsid}.npz"

        with fits.open(light_curve_file) as hdul:
            sdat = hdul[1].data
            time_col = "BARYTIME" if "BARYTIME" in sdat.columns.names else "TIME"

            st = np.asarray(sdat[time_col], dtype=float)
            sr = np.asarray(sdat["RATE"], dtype=float)
            se = np.asarray(sdat["ERROR"], dtype=float)

        with fits.open(background_file) as hdul:
            bdat = hdul[1].data

            if "RATE" in bdat.columns.names:
                bkg = np.asarray(bdat["RATE"], dtype=float)
            else:
                bkg = np.asarray(bdat["COUNTS"], dtype=float)

            mean_bkg_rate = np.mean(bkg)
            br_i = np.full_like(st, mean_bkg_rate)
            be_i = np.sqrt(br_i)

        net_r = sr - br_i
        net_e = np.sqrt(se**2 + be_i**2)

        np.savez(output_file, t=st, rate=net_r, error=net_e)

        
       
def Statistics(self):
    self.images = self.star_folder / "Outputs"
    self._Refresh_ObsID()
    self.images.mkdir(parents=True, exist_ok=True)

    kept_obsids = []
    self.medians_bary = []
    self.errors_bary = []
    self.std_bary = []

    for obsid in self.ObsID_array:

        if self.use_nicerl3:
            fits_file = self._Light_Curve_File(obsid)

            if not fits_file.exists():
                if self.auto_resolve:
                    self.Dataset_Processing_l3()
                    self.Barycorr()
                    fits_file = self._Light_Curve_File(obsid)
                else:
                    sys.exit(0)

            if not fits_file.exists():
                self.logger.warning(f"{obsid}: missing L3 light curve even after auto_resolve; skipping.")
                continue

            with fits.open(fits_file) as plot:
                data = plot[1].data
                if len(data) == 0:
                    self.logger.warning(f"{obsid}: empty light curve; skipping.")
                    continue
                rate = np.asarray(data["RATE"], dtype=float)
                error = np.asarray(data["ERROR"], dtype=float)

        else:
            try:
                t, rate, error = self._Load_Background_File(obsid)
            except RuntimeError as e:
                self.logger.warning(f"Skipping {obsid} in Statistics(): {e}")
                continue

        kept_obsids.append(obsid)
        self.medians_bary.append(np.median(rate))
        self.errors_bary.append(np.median(error))
        self.std_bary.append(np.std(rate))

    output_txt = self.images / "Light_Curve_Uncertainties.txt"
    np.savetxt(
        output_txt,
        np.column_stack([kept_obsids, self.medians_bary, self.errors_bary, self.std_bary]),
        fmt="%s",
        header="ObsID Median Error Std",
    )
