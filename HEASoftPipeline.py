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
            
        self.star_folder = self.base_dir / self.star_name
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




    
    def Dataset_Download(self, ObsID_wanted, ObsID_excluded, ObsID_Chosen):


        heasarc = self.heasarc
        table = heasarc.query_region(self.pos, catalog="nicermastr", radius="3 arcmin")




        target = self.star_name.lower().replace(" ", "").replace("_", "")
        namecol = np.array(table["name"]).astype(str)

        mask = np.array([
            target in n.lower().replace(" ", "").replace("_", "")
            for n in namecol
        ])
        table = table[mask]





       
        today_mjd = Time.now().mjd
        table = table[(table["time"] > 0) & (table["public_date"] <= today_mjd)]

        table.sort("time")

                

        ObsID_array = np.array(table["obsid"])

        ObsID_date = []
        for mjd in table["public_date"]:
            t = Time(mjd, format="mjd")
            ObsID_date.append(f"{t.datetime.year}_{t.datetime.month:02d}")
        ObsID_date = np.array(ObsID_date)

        Observation_date = []
        for mjd in table["processing_date"]:
            t = Time(mjd, format="mjd")
            Observation_date.append(f"{t.datetime.year}_{t.datetime.month:02d}")
        Observation_array = np.array(Observation_date)

        
        Epoch_array = np.array([self.obsid_epoch[str(o)] for o in self.ObsID_current])



        ObsID_run = min(ObsID_wanted, len(table))
        

        unique_obsids = np.unique(table["obsid"])

        if ObsID_excluded is not None and len(ObsID_excluded) > 0:
            unique_obsids = [obs for obs in unique_obsids if obs not in ObsID_excluded]

        if ObsID_Chosen is not None and len(ObsID_Chosen) > 0:
            unique_obsids = [obs for obs in unique_obsids if obs in ObsID_Chosen]
            
        
        selected_obsids = unique_obsids[:ObsID_wanted]

   
        subset_table = table[np.isin(table["obsid"], selected_obsids)]
        links = heasarc.locate_data(subset_table)


        downloaded = 0
        seen = set()


        for row in links:

            if downloaded >= ObsID_run:
                break

            url = row["access_url"]

            if url in seen:
                continue
            seen.add(url)

            self.logger.info(f"Downloading {url}")
            subprocess.run([os.path.expanduser("~/download_wget.pl"), url], check=True, cwd=self.star_folder)
            downloaded += 1


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
    

    def Dataset_Processing_l3(self): #Use L2 for main process but L3 is standard so can compare L2 results to

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
                barytime = "yes"

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
       
    def Statistics(self): #Statistics is preliminary and not fully complete for proper use
        self.images = self.star_folder / "Outputs"
       
        self._Refresh_ObsID()
        self.images.mkdir(parents=True, exist_ok=True)
        self.medians_bary = []
        self.errors_bary = []
        self.std_bary = []

        for obsid in self.ObsID_array:

            paths = self._ObsID_Paths(obsid)
            fits_file = self._Light_Curve_File(obsid)
            
            if not fits_file.exists():
                if self.auto_resolve:
                    if self.use_nicerl3:
                        self.logger.info(f"{fits_file} not found, running Dataset_Processing_l3 + Barycorr")
                        self.Dataset_Processing_l3()
                        self.Barycorr()
                    else:
                        self.logger.info(f"{fits_file} not found, running Barycorr_Curve")
                        self.Barycorr_Curve(self.lower_energy_limit, self.upper_energy_limit)
                    fits_file = self._Light_Curve_File(obsid) 
                else:
                    self.logger.info(f"{fits_file} not found, terminating")
                    sys.exit(0)

            with fits.open(fits_file) as plot:
                data = plot[1].data

                temp_median = np.median(data["RATE"])
                temp_error  = np.median(data["ERROR"])
                std_temp    = np.std(data["RATE"])

            self.medians_bary.append(temp_median)
            self.errors_bary.append(temp_error)
            self.std_bary.append(std_temp)

        self.logger.info(f"Medians: {self.medians_bary}")
        self.logger.info(f"Errors: {self.errors_bary}")
        self.logger.info(f"Standard Deviation: {self.std_bary}")

        output_txt = self.images / "Light_Curve_Uncertainties.txt"

        np.savetxt(
            output_txt,
            np.column_stack([
                self.ObsID_array,
                self.medians_bary,
                self.errors_bary,
                self.std_bary
            ]),
            fmt="%s",
            header="ObsID Median Error Std"
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

                flares_filter  = (data["RATE"] - median_val) > three_sigma
                below_threshold = (median_val - data["RATE"]) > three_sigma
                normal = ~(flares_filter | below_threshold)

                
                time_col = "BARYTIME" if "BARYTIME" in data.columns.names else "TIME"
                t = np.asarray(data[time_col], dtype=float)

                t0 = header.get("TSTART", float(t[0]))
                t_rel = t - t0

                plt.scatter(t_rel[flares_filter],   data["RATE"][flares_filter],   color="red")
                plt.scatter(t_rel[below_threshold], data["RATE"][below_threshold], color="gray")
                plt.scatter(t_rel[normal],          data["RATE"][normal],          color="blue")

                plt.errorbar(
                    t_rel,
                    data["RATE"],
                    yerr=data["ERROR"],
                    fmt=".",
                    color="gray",
                    capsize=3,
                    capthick=1,
                    alpha=0.5
                )

                plt.xlabel("Time since start (s)")
                plt.ylabel("Count rate (/s)")
                plt.title(f"{self.star_name} {obsid} ({time_col})")

                xmax = min(self.orbital_period, float(np.nanmax(t_rel)))
                plt.xlim(0, xmax)

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

                three_sigma = 3 * self.std_bary[i]
                median_val  = self.medians_bary[i]

                flares_filter  = (data["RATE"] - median_val) > three_sigma
                below_threshold = (median_val - data["RATE"]) > three_sigma
                normal = ~(flares_filter | below_threshold)

                header = plot[1].header
                mjdref = header["MJDREFI"] + header["MJDREFF"]

                time_col = "BARYTIME" if "BARYTIME" in data.columns.names else "TIME"
                time_array = data[time_col]

                phase = self._Phase_Calculator(time_array, mjdref)

                phase_plot = np.concatenate([phase, phase + 1.0])
                rate_plot  = np.concatenate([data["RATE"], data["RATE"]])
                error_plot = np.concatenate([data["ERROR"], data["ERROR"]])

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

                time_col = "BARYTIME" if "BARYTIME" in data.columns.names else "TIME"
                time_array = data[time_col]

                phase = self._Phase_Calculator(time_array, mjdref)



                
                rate = data["RATE"]
                error = data["ERROR"]

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

                    time_col = "BARYTIME" if "BARYTIME" in data.columns.names else "TIME"
                    time_array = data[time_col]

                    phases = self._Phase_Calculator(time_array, mjdref)


                    
                    rates = data["RATE"]
                    errors = data["ERROR"]

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

            

            plt.errorbar(centers, binned_rate, yerr=binned_err, fmt='o', capsize=3)
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
        use_nicerl3=config["run_nicer_l3"]
        
        
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

