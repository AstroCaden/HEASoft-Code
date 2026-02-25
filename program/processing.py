import subprocess
import os
import sys
import numpy as np
from astropy.io import fits


def Dataset_Processing_l2(self):
    self._Refresh_ObsID()
    if len(self.ObsID_current) == 0:
        self.logger.info("No datasets found")
        sys.exit(0)

    for obsid in self.ObsID_current:
        paths = self._ObsID_Paths(obsid)

        hk_dir = paths["xti_dir"] / "hk"
        hk_files = list(hk_dir.glob(f"ni{obsid}_?mpu*.hk")) + list(hk_dir.glob(f"ni{obsid}_?mpu*.hk.gz"))
        if len(hk_files) == 0:
            self._Failed_ObsID(obsid, f"missing hk files in {hk_dir}", where="nicerl2")
            continue

        confirmation_file = paths["event_cl"] / f"ni{obsid}_0mpu7_cl.evt"
        if confirmation_file.exists() and confirmation_file.stat().st_size > 0:
            self.logger.info(f"{obsid} is already L2 processed")
            continue

        try:
            cmd_nicerl2 = self.base_dir / "CMD_files" / "cmd_nicer_l2.sh"
            subprocess.run([cmd_nicerl2.as_posix(), paths["obs_dir"].as_posix()], check=True)
            self.logger.info(f"{obsid} has been L2 processed")
        except subprocess.CalledProcessError as e:
            self._Failed_ObsID(obsid, f"nicerl2 failed (code {e.returncode})", where="nicerl2")
            continue

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

            
            cmd_nicerl3 = self.base_dir / "CMD_files" / "cmd_nicer_l3.sh"
            subprocess.run(
                [cmd_nicerl3.as_posix(), str(paths["xti_dir"].parent), f"{filter_min}-{filter_max}"],
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

        cmd_barycorr = self.base_dir / "CMD_files" / "cmd_barycorr.sh"
        log_dir = self.star_folder / "Outputs" / "Barycorr_Logs"
        log_dir.mkdir(parents=True, exist_ok=True)
        log_file = log_dir / f"barycorr_{obsid}.log"
        self.logger.info(f"[{datetime.now()strftime('%H:%M:%S')}] - Barycorr processing {obsid}")
        with open(log_file, "w") as f:
            p = subprocess.run(
                [
                    cmd_barycorr.as_posix(),
                    infile.as_posix(),
                    outfile.as_posix(),
                    orbitf.as_posix(),
                    str(self.RA),
                    str(self.Dec),
                    barytime,
                ],
                stdout=f,
                stderr=subprocess.STDOUT,
                text=True
            )

        self.logger.info(f"[{datetime.now()strftime('%H:%M:%S')}] - Barycorr complete for {obsid}")
        if p.returncode != 0:
            reason = (p.stderr or p.stdout or "unknown error").strip()
            self._Failed_ObsID(obsid, reason, where="barycorr")
            continue

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
