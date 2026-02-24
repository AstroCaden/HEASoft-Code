import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
from datetime import datetime
from datetime import timedelta
from astroquery.heasarc import Heasarc
from astropy.time import Time
from astropy.table import Table
import subprocess
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
        orbital_pdot,
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
        self.orbital_pdot = orbital_pdot
        
        
        if working_dir == "working":
            self.base_dir = Path.cwd()
        else:
            path = Path(working_dir)
            if path.exists():
                
                self.base_dir = Path(working_dir).resolve()
            else:
                raise RuntimeError("Working directory does not exist")
                
        
        

        safe_folder = re.sub(r"[^A-Za-z0-9_.-]+", "_", self.star_name).strip("_")

        datasets_dir = self.base_dir / "Datasets"
        datasets_dir.mkdir(exist_ok=True)

        self.star_folder = datasets_dir / safe_folder
        self.star_folder.mkdir(exist_ok=True)
        self.obsids_dir = self.star_folder / "ObsIDs"
        self.obsids_dir.mkdir(exist_ok=True)
        print(f"Inside {self.star_folder}")

        heasarc = Heasarc()
        self.heasarc = heasarc
        self.pos = SkyCoord.from_name(f"{self.star_name}")
        self.RA = self.pos.ra.deg
        self.Dec = self.pos.dec.deg

        print(f"RA:{self.RA}, Dec {self.Dec}")
        

        self.orbital_period = orbital_period * 86400
 
        

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
        self.failed = {}

        f = self.star_folder / "failed_obsids.json"
        if f.exists():
            try:
                self.failed = json.loads(f.read_text())
            except Exception:
                self.failed = {}


        self._Refresh_ObsID()
        self._Reprocessing()



    def _ObsID_Paths(self, obsid):
        obs_dir = self.obsids_dir / obsid
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
        if not self.obsids_dir.exists():
            self.ObsID_array = []
            self.ObsID_current = []
            return

        failed = set(getattr(self, "failed", {}).keys())

        self.ObsID_array = sorted([
            d.name
            for d in self.obsids_dir.iterdir()
            if d.is_dir()
            and d.name[:9].isdigit()
            and d.name not in failed
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
                    


    def _Failed_ObsID(self, obsid: str, reason: str, where: str = ""):

        msg = f"{obsid} FAILED"
        if where:
            msg += f" in {where}"
        msg += f": {reason}"

        self.logger.warning(msg)

        self.failed[str(obsid)] = msg

     
        f = self.star_folder / "failed_obsids.json"
        try:
            f.write_text(json.dumps(self.failed, indent=2))
        except Exception as e:
            self.logger.warning(f"Could not save failed_obsids.json: {e}")


from program.download import Dataset_Download, _Single_Downloader
from program.processing import Dataset_Processing_l2, Dataset_Processing_l3, Barycorr, Barycorr_Curve
from program.background_stats import Background_Subtraction, Statistics
from program.plotting import Time_Plot, Unbinned_Plot, Binned_Plot, Stacked_Plot
from program.analysis_methods import _Phase_Calculator, _Phase_Binner, _Epoch_Binning, _Flare_Filtering, _Comparison_Model, _Light_Curve_File, _Load_Background_File

Pipeline.Dataset_Download = Dataset_Download
Pipeline._Single_Downloader = _Single_Downloader

Pipeline.Dataset_Processing_l2 = Dataset_Processing_l2
Pipeline.Dataset_Processing_l3 = Dataset_Processing_l3
Pipeline.Barycorr = Barycorr
Pipeline.Barycorr_Curve = Barycorr_Curve

Pipeline.Background_Subtraction = Background_Subtraction
Pipeline.Statistics = Statistics

Pipeline.Time_Plot = Time_Plot
Pipeline.Unbinned_Plot = Unbinned_Plot
Pipeline.Binned_Plot = Binned_Plot
Pipeline.Stacked_Plot = Stacked_Plot

Pipeline._Phase_Calculator = _Phase_Calculator
Pipeline._Phase_Binner = _Phase_Binner
Pipeline._Epoch_Binning = _Epoch_Binning
Pipeline._Flare_Filtering = _Flare_Filtering
Pipeline._Comparison_Model = _Comparison_Model
Pipeline._Light_Curve_File = _Light_Curve_File
Pipeline._Load_Background_File = _Load_Background_File
