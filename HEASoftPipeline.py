import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
from datetime import datetime
from datetime import timedelta
import glob
from astroquery.heasarc import Heasarc
from astropy.time import Time
from astropy.table import Table
import subprocess
from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry
from astropy.coordinates import SkyCoord
import logging
from pathlib import Path
import json
import sys
import gzip, shutil
import re
from concurrent.futures import ThreadPoolExecutor

class Pipeline:

    def __init__(
        self,
        star_name,
        orbital_period,
        working_dir,
        auto_resolve,
        reprocess,
        number_of_datasets,
        obsid_to_exclude,
        obsid_specifically_chosen,
        lower_energy_limit,
        upper_energy_limit,
        nbins,
        epoch_width,
        flare_filter,
        T0,
        use_nicerl3,
        concurrent_downloads,
    ):





        
        self.star_name = star_name
        self.orbital_period = orbital_period
        self.working_dir = working_dir
        self.auto_resolve = auto_resolve
        self.reprocess = reprocess
        self.number_of_datasets = number_of_datasets
        self.ObsID_excluded = obsid_to_exclude
        self.obsid_specifically_chosen = obsid_specifically_chosen
        self.lower_energy_limit = lower_energy_limit
        self.upper_energy_limit = upper_energy_limit
        self.nbins = nbins
        self.epoch_width = epoch_width
        self.flare_filter = flare_filter
        self.T0 = T0
        self.use_nicerl3 = use_nicerl3
        self.concurrent_downloads = concurrent_downloads
        
        
        if working_dir == "working":
            self.base_dir = Path.cwd()
        else:
            path = Path(working_dir)
            if path.exists():
                
                self.base_dir = Path(working_dir).resolve()
            else:
                print("The directory you entered doesn't exist")
                return
                
        
        

        self.obsid_epoch = {}
        safe_folder = re.sub(r"[^A-Za-z0-9_.-]+", "_", self.star_name).strip("_")
        self.star_folder = self.base_dir / safe_folder
        self.star_folder.mkdir(exist_ok=True)
        print(f"Inside {self.star_folder}")

        heasarc = Heasarc()
        self.heasarc = heasarc
        self.pos = SkyCoord.from_name(f"{self.star_name}")
        self.RA = self.pos.ra.deg
        self.Dec = self.pos.dec.deg

        print(f"RA:{self.RA}, Dec {self.Dec}")
        
        self.ObsID_current = [e.name for e in self.star_folder.iterdir() if e.is_dir() and e.name.isdigit()]
        self.orbital_period = orbital_period * 86400
 
        self._Refresh_ObsID()

        self.images = self.star_folder / "Outputs"
        if self.images.exists():
            shutil.rmtree(self.images)
        self.images.mkdir(parents=True, exist_ok=True)




        self.logger = logging.getLogger(self.star_name)
        self.logger.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

        
        ch = logging.StreamHandler()
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)

        
        log_file = self.star_folder / f"{self.star_name}_pipeline.log"

        if log_file.exists():
            log_file.unlink()


        
        fh = logging.FileHandler(log_file)
        fh.setFormatter(formatter)
        self.logger.addHandler(fh)
        self.logger.info(f"Pipeline initialised for {self.star_name} at {self.star_folder}")
        
        self._Reprocessing()



    def _ObsID_Paths(self, obsid):
        obs_dir = self.star_folder / obsid
        xti_dir = obs_dir / "xti"
        event_cl = xti_dir / "event_cl"
        auxil = obs_dir / "auxil"

        return {
            "obs_dir": obs_dir,
            "xti_dir": xti_dir,
            "event_cl": event_cl,
            "auxil": auxil
        }
                

        

    def _Refresh_ObsID(self):
        if not self.star_folder.exists():
            self.ObsID_array = []
            self.ObsID_current = []
            return

        self.ObsID_array = sorted([
            d.name
            for d in self.star_folder.iterdir()
            if d.is_dir() and d.name[:9].isdigit()  
        ])
        self.ObsID_current = self.ObsID_array.copy()



    def _Reprocessing(self):

        if self.reprocess == True:
            self.logger.info("Beginning reprocessing")

            for obsid in self.ObsID_current:
                paths = self._ObsID_Paths(obsid)

                event_dir = paths["event_cl"]

                delete_file = event_dir / f"ni{obsid}_0mpu7_cl.evt"
                if delete_file.exists():
                    delete_file.unlink()
                    self.logger.info(f"Removed {delete_file}")

                delete_file = event_dir / f"ni{obsid}_0mpu7_cl_bary.evt"
                if delete_file.exists():
                    delete_file.unlink()
                    self.logger.info(f"Removed {delete_file}")

                delete_file = paths["xti_dir"] / f"{obsid}_bary.fits"
                if delete_file.exists():
                    delete_file.unlink()
                    self.logger.info(f"Removed {delete_file}")

            self.logger.info("Reprocessing completed")
                    


    def _Phase_Calculator(self, time_array, MJD_ref):

        MJD_times = MJD_ref + (time_array / 86400.0)
        P_orb_days = self.orbital_period / 86400.0
        phase = ((MJD_times - self.T0) / P_orb_days) % 1.0

        return phase



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
        file = self.images / f"netlc_{obsid}.npz"
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


    def _Single_Downloader(self, url: str):
        self.logger.info(f"Downloading {url}")
        subprocess.run(
            [os.path.expanduser("~/download_wget.pl"), url],
            cwd=self.star_folder
        )



    def Dataset_Download(self, ObsID_wanted, ObsID_excluded, ObsID_Chosen):

        heasarc = self.heasarc
        table = heasarc.query_region(self.pos, catalog="nicermastr", radius="3 arcmin")

        target = self.star_name.lower().replace(" ", "").replace("_", "")
        namecol = np.array(table["name"]).astype(str)
        mask = np.array([target in n.lower().replace(" ", "").replace("_", "") for n in namecol])
        table = table[mask]

        today_mjd = Time.now().mjd
        table = table[(table["time"] > 0) & (table["public_date"] <= today_mjd)]
        table.sort("time")

        obsids_time = np.array(table["obsid"]).astype(str)
        _, first_idx = np.unique(obsids_time, return_index=True)
        unique_obsids = obsids_time[np.sort(first_idx)]

        if ObsID_excluded is not None and len(ObsID_excluded) > 0:
            if ObsID_excluded == "current":
                ObsID_excluded = self.ObsID_current
            ObsID_excluded = set(map(str, ObsID_excluded))

            before = list(unique_obsids)
            unique_obsids = np.array([o for o in unique_obsids if o not in ObsID_excluded], dtype=str)
            removed = sorted(set(before) - set(unique_obsids))

            self.logger.info(f"Excluded {len(removed)} ObsIDs by request.")

        if ObsID_Chosen is not None and len(ObsID_Chosen) > 0:
            ObsID_Chosen = set(map(str, ObsID_Chosen))

            unique_obsids = np.array([o for o in unique_obsids if o in ObsID_Chosen], dtype=str)
            self.logger.info(f"Whitelisted {len(unique_obsids)} ObsIDs by request.")

        if len(unique_obsids) == 0:
            self.logger.info("No ObsIDs left after filtering; nothing to download.")
            return

        n = min(int(ObsID_wanted), len(unique_obsids))
        idx = np.linspace(0, len(unique_obsids) - 1, n, dtype=int)
        selected_obsids = unique_obsids[idx]
        ObsID_run = len(selected_obsids)
        self.logger.info(f"Will download {len(selected_obsids)} ObsIDs: {selected_obsids}")

        subset_table = table[np.isin(np.array(table["obsid"]).astype(str), selected_obsids)]
        links = heasarc.locate_data(subset_table)

        urls = []
        seen = set()
        for row in links:
            url = row["access_url"]
            if url in seen:
                continue
            seen.add(url)
            urls.append(url)
            if len(urls) >= ObsID_run:
                break

        if not urls:
            self.logger.info("No download URLs found; nothing to download.")
            return

        max_workers = self.concurrent_downloads
        self.logger.info(f"Downloading {len(urls)} ObsIDs with {max_workers} workers")

        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            pool.map(self._Single_Downloader, urls)

        self._Refresh_ObsID()



      

    def Dataset_Processing_l2(self):

        
        self._Refresh_ObsID()
        if len(self.ObsID_current) == 0:
            self.logger.info("No datasets found")
            sys.exit(0)
        
        for obsid in self.ObsID_current:
            paths = self._ObsID_Paths(obsid)
            event_dir = paths["event_cl"]

            confirmation_file = event_dir / f"ni{obsid}_0mpu7_cl.evt"


            if not confirmation_file.exists() or confirmation_file.stat().st_size == 0:
                subprocess.run(
                    ["nicerl2", f"indir={paths['obs_dir']}", "clobber=yes"],
                    check=True
                )
                self.logger.info(f"{obsid} has been L2 processed")
            else:
                self.logger.info(f"{obsid} is already L2 processed")


        if self.use_nicerl3 == True:
            self.Dataset_Processing_l3()
    

    def Dataset_Processing_l3(self): #Use L2 for main process but L3 is standard so can compare with L2

        filter_min = self.lower_energy_limit
        filter_max = self.upper_energy_limit

        self._Refresh_ObsID()
        if len(self.ObsID_current) == 0:
            if self.auto_resolve == True:
                self.Dataset_Processing_l2()
                self._Refresh_ObsID()
            else:
                sys.exit(0)

        for obsid in self.ObsID_current:
            paths = self._ObsID_Paths(obsid)
            event_dir = paths["event_cl"]

            confirmation_file = event_dir / f"ni{obsid}mpu7_sr.lc"

            if not confirmation_file.exists() or confirmation_file.stat().st_size == 0:

                
                subprocess.run(
                    ["nicerl3-lc", str(paths["xti_dir"].parent), f"pirange={filter_min}-{filter_max}", "timebin=30", "clobber=YES"],
                    check=True
                )


                
                self.logger.info(f"{obsid} has been L3 processed")
            else:
                self.logger.info(f"{obsid} is already L3 processed")
  

               
    def Barycorr(self):
        self._Refresh_ObsID()

        for obsid in self.ObsID_current:
            paths = self._ObsID_Paths(obsid)
            orbitf = paths["auxil"] / f"ni{obsid}.orb.gz"

            if self.use_nicerl3:
              
                infile  = paths["event_cl"] / f"ni{obsid}mpu7_sr.lc"
                outfile = paths["event_cl"] / f"ni{obsid}mpu7_sr_bary.lc"
                barytime = "no"

                if not infile.exists():
                    if self.auto_resolve:
                        self.logger.info(f"{infile} not found, running Dataset_Processing_l3()")
                        self.Dataset_Processing_l3()
                    else:
                        self.logger.info(f"{infile} not found, terminating")
                        sys.exit(0)

            else:
               
                infile  = paths["event_cl"] / f"ni{obsid}_0mpu7_cl.evt"
                outfile = paths["event_cl"] / f"ni{obsid}_0mpu7_cl_bary.evt"
                barytime = "yes"

                if not infile.exists():
                    if self.auto_resolve:
                        self.logger.info(f"{infile} not found, running Dataset_Processing_l2()")
                        self.Dataset_Processing_l2()
                    else:
                        self.logger.info(f"{infile} not found, terminating")
                        sys.exit(0)

            if outfile.exists():
                self.logger.info(f"{outfile.name} already exists")
                continue

            cmd = [
                "barycorr",
                f"infile={infile}",
                f"outfile={outfile}",
                f"orbitfiles={orbitf}",
                f"ra={self.RA}",
                f"dec={self.Dec}",
                "refframe=ICRS",
                "ephem=JPLEPH.430", 
                f"barytime={barytime}",
                "clobber=yes",
            ]

            env = os.environ.copy()
            env["HEADASNOQUERY"] = "1"
            env["HEADASYES"] = "1"

            subprocess.run(cmd, env=env, check=True)
            self.logger.info(f"{outfile.name} generated")


    def Barycorr_Curve(self, filter_min, filter_max):

        if self.use_nicerl3:
            self.logger.info("L3 is active so Barycorr_Curve is skipped")
            return
        
        self._Refresh_ObsID()

        cmd_bary_curve_xco = self.base_dir / "CMD_files" / "cmd_bary_curve.xco"

        for obsid in self.ObsID_array:
            paths = self._ObsID_Paths(obsid)

             
            event_file = paths["event_cl"] / f"ni{obsid}_0mpu7_cl_bary.evt"
            output_curve = paths["xti_dir"] / f"{obsid}_bary.fits"

            if output_curve.exists():
                self.logger.info(f"{output_curve.name} already exists")
                continue

            with open(cmd_bary_curve_xco, "r") as f:
                xco_cmds = f.read()

            xco_cmds = xco_cmds.format(
                event_dir=paths["event_cl"].as_posix(),
                event_file=f"ni{obsid}_0mpu7_cl_bary.evt",
                filter_min=filter_min,
                filter_max=filter_max,
                output_curve=output_curve.as_posix()
            )

            if not event_file.exists():
                if self.auto_resolve:
                    self.logger.info(f"{event_file} not found, running Barycorr")
                    self.Barycorr()
                    if not event_file.exists():
                        self.logger.info(f"{event_file} still missing after Barycorr, skipping {obsid}")
                        continue
                else:
                    self.logger.info(f"{event_file} not found, terminating")
                    sys.exit(0)
                
            subprocess.run(f"xselect << EOF\n{xco_cmds}\nEOF",
                           shell=True, executable="/bin/bash", check=True)

            self.logger.info(f"{output_curve.name} generated")



    def Background_Subtraction(self):

        if self.use_nicerl3:
            self.logger.info("L3 is active so Background_Subtraction is skipped")
            return

        self._Refresh_ObsID()
        self.images = self.star_folder / "Outputs"
        self.images.mkdir(parents=True, exist_ok=True)

        for obsid in self.ObsID_array:

            light_curve_file = self._Light_Curve_File(obsid)
            paths = self._ObsID_Paths(obsid)

            background_file = paths["event_cl"] / "nibackgen3C50_bkg.pi"

            if not light_curve_file.exists():
                if self.auto_resolve:
                    self.Barycorr_Curve(
                        self.lower_energy_limit,
                        self.upper_energy_limit
                    )
                else:
                    sys.exit(0)

            if not background_file.exists():

                evt_file = paths["event_cl"] / f"ni{obsid}_0mpu7_cl_bary.evt"

                if not evt_file.exists():
                    if self.auto_resolve:
                        self.Barycorr()
                    else:
                        sys.exit(0)

                result = subprocess.run(
                    [
                        "nibackgen3C50",
                        f"rootdir={self.star_folder}",
                        f"obsid={obsid}",
                        f"bkgidxdir={os.environ['CALDB']}/data/nicer/xti/bcf/bkg",
                        f"bkglibdir={os.environ['CALDB']}/data/nicer/xti/bcf/bkg",
                        "gainepoch=2020",
                        "clobber=YES",
                    ],
                    cwd=paths["event_cl"],
                    capture_output=True,
                    text=True,
                    check=False,   
                )

                if result.returncode == 0:
                    self.logger.info(f"Background spectrum generated for {obsid}")
                elif result.returncode == 218:
                    self.logger.warning(f"{obsid}: nibackgen3C50 returned 218 (no good exposure after cuts). Skipping.")
                    bad_dir = paths["obs_dir"]
                    new_dir = bad_dir.parent / f"#{obsid}"
                    try:
                        bad_dir.rename(new_dir)
                        self.logger.warning(f"{obsid}: renamed folder to {new_dir.name} so it will be skipped next runs.")
                    except Exception as e:
                        self.logger.warning(f"{obsid}: failed to rename folder ({e}); leaving as-is.")
                    continue

                elif result.returncode == 255 and ("GTI re-definition FAILED" in (result.stdout + result.stderr)):
                    self.logger.warning(f"{obsid}: nibackgen3C50 GTI re-definition failed (255). Skipping (zero exposure).")
                    bad_dir = paths["obs_dir"]
                    new_dir = bad_dir.parent / f"#{obsid}"
                    try:
                        bad_dir.rename(new_dir)
                        self.logger.warning(f"{obsid}: renamed folder to {new_dir.name} so it will be skipped next runs.")
                    except Exception as e:
                        self.logger.warning(f"{obsid}: failed to rename folder ({e}); leaving as-is.")
                    continue

                else:
                    self.logger.error(f"{obsid}: nibackgen3C50 failed with code {result.returncode}")
                    self.logger.error(result.stdout)
                    self.logger.error(result.stderr)
                    result.check_returncode()  

            output_file = self.images / f"netlc_{obsid}.npz"

            
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

            np.savez(
                output_file,
                t=st,
                rate=net_r,
                error=net_e
            )



            
           
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
        
                
    def Time_Plot(self):

        self._Refresh_ObsID()
        self.images = self.star_folder / "Outputs"
      
        stats_file = self.images / "Light_Curve_Uncertainties.txt"
        
        if not stats_file.exists():
            if self.auto_resolve == True:
                self.logger.info(f"{stats_file} not found, running Statistics")
                self.Statistics()
            else:
                self.logger.info(f"{stats_file} not found, terminating")
                sys.exit(0)

        stats = np.loadtxt(stats_file, dtype=str, skiprows=1, ndmin=2)

        self.medians_bary = stats[:, 1].astype(float)
        self.errors_bary  = stats[:, 2].astype(float)
        self.std_bary     = stats[:, 3].astype(float)

        for i, obsid in enumerate(self.ObsID_array):

            paths = self._ObsID_Paths(obsid)
            fits_file = self._Light_Curve_File(obsid)

            if not fits_file.exists():
                self.logger.info(f"{fits_file.name} missing, skipping")
                continue

            with fits.open(fits_file) as plot:
                data = plot[1].data
                header = plot[1].header

                if len(data) == 0:
                    self.logger.info(f"{obsid} is empty, skipping")
                    continue

                three_sigma = 3 * self.std_bary[i]
                median_val  = self.medians_bary[i]

                try:
                    t, rate, error = self._Load_Background_File(obsid)
                except RuntimeError as e:
                    self.logger.warning(f"Skipping {obsid}: {e}")
                    continue

                three_sigma = 3 * self.std_bary[i]
                median_val  = self.medians_bary[i]

                flares_filter  = (rate - median_val) > three_sigma
                below_threshold = (median_val - rate) > three_sigma
                normal = ~(flares_filter | below_threshold)

                t0 = header.get("TSTART", float(t[0]))
                t_rel = t - t0


                plt.scatter(t_rel[flares_filter],   rate[flares_filter],   color="red")
                plt.scatter(t_rel[below_threshold], rate[below_threshold], color="gray")
                plt.scatter(t_rel[normal],          rate[normal],          color="blue")

                plt.errorbar(
                    t_rel,
                    rate,
                    yerr=error,
                    fmt=".",
                    color="gray",
                    capsize=3,
                    capthick=1,
                    alpha=0.5
                )

                plt.xlabel("Time since start (s)")
                plt.ylabel("Count rate (/s)")
                plt.title(f"{self.star_name} {obsid}")



                output_png = self.images / f"barycentered_{self.star_name}_{obsid}.png"

                plt.savefig(output_png)
                plt.clf()
                plt.close()




    def Unbinned_Plot(self):

        self._Refresh_ObsID()
        self.images = self.star_folder / "Outputs"
       
        stats_file = self.images / "Light_Curve_Uncertainties.txt"
        
        if not stats_file.exists():
            if self.auto_resolve == True:
                self.logger.info(f"{stats_file} not found, running Statistics")
                self.Statistics()
            else:
                self.logger.info(f"{stats_file} not found, terminating")
                sys.exit(0)

        stats = np.loadtxt(stats_file, dtype=str, skiprows=1, ndmin=2)

        self.medians_bary = stats[:, 1].astype(float)
        self.errors_bary  = stats[:, 2].astype(float)
        self.std_bary     = stats[:, 3].astype(float)

        for i, obsid in enumerate(self.ObsID_array):

            paths = self._ObsID_Paths(obsid)
            fits_file = self._Light_Curve_File(obsid)

            if not fits_file.exists():
                self.logger.info(f"{fits_file.name} missing, skipping")
                continue

            with fits.open(fits_file) as plot:
                data = plot[1].data

                if len(data) == 0:
                    self.logger.info(f"{obsid} is empty, skipping")
                    continue

                time_array, rate, error = self._Load_Background_File(obsid)

                three_sigma = 3 * self.std_bary[i]
                median_val  = self.medians_bary[i]

                flares_filter  = (rate - median_val) > three_sigma
                below_threshold = (median_val - rate) > three_sigma
                normal = ~(flares_filter | below_threshold)

                header = plot[1].header
                mjdref = header["MJDREFI"] + header["MJDREFF"]

        

                phase = self._Phase_Calculator(time_array, mjdref)

                phase_plot = np.concatenate([phase, phase + 1.0])
                rate_plot  = np.concatenate([rate, rate])
                error_plot = np.concatenate([error, error])

                flares_plot = np.concatenate([flares_filter, flares_filter])
                below_plot  = np.concatenate([below_threshold, below_threshold])
                normal_plot = np.concatenate([normal, normal])

                plt.scatter(phase_plot[flares_plot], rate_plot[flares_plot], color="red")
                plt.scatter(phase_plot[below_plot],  rate_plot[below_plot],  color="gray")
                plt.scatter(phase_plot[normal_plot], rate_plot[normal_plot], color="blue")

                plt.errorbar(
                    phase_plot,
                    rate_plot,
                    yerr=error_plot,
                    fmt=".",
                    color="gray",
                    capsize=3,
                    capthick=1,
                    alpha=0.5
                )

                plt.xlabel("Orbital phase")
                plt.ylabel("Count rate (/s)")
                plt.title(f"{self.star_name} {obsid} phase folded (unbinned)")
                plt.xlim(0, 2)

                output_png = self.images / f"phase_folded_unbinned_{self.star_name}_{obsid}.png"
                plt.savefig(output_png)
                plt.clf()
                plt.close()
                        

    def Binned_Plot(self):
        self._Refresh_ObsID()
        self.images = self.star_folder / "Outputs"
       

        for obsid in self.ObsID_array:

            paths = self._ObsID_Paths(obsid)
            fits_file = self._Light_Curve_File(obsid)

            if not fits_file.exists():
                continue

            with fits.open(fits_file) as plot:
                data = plot[1].data

                if len(data) == 0:
                    continue

                header = plot[1].header
                mjdref = header["MJDREFI"] + header["MJDREFF"]

                time_array, rate, error = self._Load_Background_File(obsid)

                phase = self._Phase_Calculator(time_array, mjdref)

                if self.flare_filter:
                    median = np.median(rate)
                    std = np.std(rate)
                    flares, below, normal = self._Flare_Filtering(rate, median, std)
                    phase = phase[normal]
                    rate = rate[normal]
                    error = error[normal]

                centers, binned_rate, binned_err = self._Phase_Binner(phase, rate, error)
                fit = self._Comparison_Model(centers, binned_rate, binned_err)
                if fit is not None:
                    C, A, phi0, chi2, dof = fit
                    self.logger.info(f"{obsid} chi2/dof = {chi2/dof:.3f} (dof={dof})")

       
                plt.errorbar(
                    centers,
                    binned_rate,
                    yerr=binned_err,
                    fmt="o",
                    capsize=3
                )

                plt.xlabel("Orbital phase")
                plt.ylabel("Mean count rate (/s)")
                plt.title(f"{self.star_name} {obsid} phase folded (binned)")
                plt.xlim(0, 1)

                output_png = self.images / f"phase_folded_binned_{self.star_name}_{obsid}.png"
                plt.savefig(output_png)
                plt.clf()
                plt.close()




    def Stacked_Plot(self):
        self._Refresh_ObsID()
        self.images = self.star_folder / "Outputs"
        

        self._Epoch_Binning()
        epoch_groups = {}
        for obsid in self.ObsID_array:
            mjd = self.obsid_epoch.get(obsid)
            if mjd is None:
                continue
            epoch_start = mjd - (mjd % self.epoch_width)
            epoch_groups.setdefault(epoch_start, []).append(obsid)

        for epoch_start in sorted(epoch_groups.keys()):
            obsids = epoch_groups[epoch_start]
            all_phases, all_rates, all_errors, all_obsids = [], [], [], []

            for obsid in obsids:
                paths = self._ObsID_Paths(obsid)
                fits_file = self._Light_Curve_File(obsid)
                if not fits_file.exists():
                    continue
                with fits.open(fits_file) as plot:
                    data = plot[1].data
                    if len(data) == 0:
                        continue

                    header = plot[1].header
                    mjdref = header["MJDREFI"] + header["MJDREFF"]

            
                    
                    time_array, rates, errors = self._Load_Background_File(obsid)

                    phases = self._Phase_Calculator(time_array, mjdref)

                    if self.flare_filter:
                        median = np.median(rates)
                        std = np.std(rates)
                        flares, below, normal = self._Flare_Filtering(rates, median, std)
                        phases = phases[normal]
                        rates = rates[normal]
                        errors = errors[normal]

                    mask_nonzero = errors > 0
                    phases = phases[mask_nonzero]
                    rates = rates[mask_nonzero]
                    errors = errors[mask_nonzero]

                    if len(phases) == 0:
                        continue

                    all_phases.append(phases)
                    all_rates.append(rates)
                    all_errors.append(errors)
                    all_obsids.append(np.full(len(phases), obsid))

            if len(all_phases) == 0:
                self.logger.info(f"No valid data for epoch starting at MJD {epoch_start:.2f}, skipping")
                continue

            all_phases = np.concatenate(all_phases)
            all_rates = np.concatenate(all_rates)
            all_errors = np.concatenate(all_errors)
            all_obsids = np.concatenate(all_obsids)

            output_txt = self.images / f"stacked_epoch_{self.star_name}_{int(epoch_start)}.txt"
            np.savetxt(
                output_txt,
                np.column_stack([all_obsids, all_phases, all_rates, all_errors]),
                fmt="%s",
                header="ObsID Phase Rate Error"
            )

            centers, binned_rate, binned_err = self._Phase_Binner(all_phases, all_rates, all_errors)

            fit = self._Comparison_Model(centers, binned_rate, binned_err)
            if fit is not None:
                C, A, phi0, chi2, dof = fit
                self.logger.info(f"Epoch {int(epoch_start)} chi2/dof = {chi2/dof:.3f} (dof={dof})")

            

            m = np.isfinite(binned_rate) & np.isfinite(binned_err)
            plt.errorbar(centers[m], binned_rate[m], yerr=binned_err[m], fmt="o", capsize=3)
            plt.xlabel("Orbital phase")
            plt.ylabel("Mean count rate (/s)")
            plt.title(f"{self.star_name} stacked phase-folded light curve (Epoch start MJD {epoch_start:.2f})")
            plt.xlim(0, 1)

            output_png = self.images / f"stacked_epoch_{self.star_name}_{int(epoch_start)}.png"
            plt.savefig(output_png)
            plt.clf()
            plt.close()
            self.logger.info(f"Stacked phase-folded curve saved as {output_png} and data saved as {output_txt}")


##########################################################################################################################



if __name__ == "__main__":
   
    config_file = sys.argv[1] if len(sys.argv) > 1 else "config.json"

    with open(config_file, "r") as f:
        config = json.load(f)

    pipeline = Pipeline(
        star_name=config["star_name"],
        orbital_period=config["orbital_period"],
        working_dir=config["working_dir"],
        auto_resolve=config["auto_resolve"],
        reprocess=config["reprocess"],
        number_of_datasets=config["number_of_datasets"],
        obsid_to_exclude=config["ObsID_to_exclude"],
        obsid_specifically_chosen=config["ObsID_specifically_chosen"],
        lower_energy_limit=config["lower_energy_limit"],
        upper_energy_limit=config["upper_energy_limit"],
        nbins=config["bins"],
        epoch_width=config["epoch_range"],
        flare_filter=config["flare_filtering"],
        T0=config["reference_epoch"],
        use_nicerl3=config["run_nicer_l3"],
        concurrent_downloads=config["concurrent_downloads"],
        
        
    )

    steps = config.get("run_steps", {})

    if steps.get("download", False):
        pipeline.Dataset_Download(
            ObsID_wanted=config["number_of_datasets"],
            ObsID_excluded=config["ObsID_to_exclude"],
            ObsID_Chosen=config["ObsID_specifically_chosen"]
            
        )

    if steps.get("process_l2", False):
        pipeline.Dataset_Processing_l2()



    if steps.get("barycorr", False):
        pipeline.Barycorr()

    if steps.get("barycurve", False):
        pipeline.Barycorr_Curve(
        filter_min=config["lower_energy_limit"],
        filter_max=config["upper_energy_limit"]
        )

    if steps.get("remove_background", False):
        pipeline.Background_Subtraction()


    if steps.get("statistics", False):
        pipeline.Statistics()

    if steps.get("time_plot", False):
        pipeline.Time_Plot()
        
    if steps.get("unbinned_plot", False):
        pipeline.Unbinned_Plot()
        
    if steps.get("binned_plot", False):
        pipeline.Binned_Plot()

    if steps.get("stacked_plot", False):
        pipeline.Stacked_Plot()

    pipeline.logger.info(f"Pipeline complete")

