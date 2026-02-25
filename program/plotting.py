import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys


plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
})


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
    if stats.size == 0 or stats.shape[1] < 4:
        self.logger.info(f"{stats_file} is empty, unable to run Time_Plot")
        return

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


            plt.scatter(t_rel[flares_filter],   rate[flares_filter],   s=14, marker="o", linewidths=0.5, edgecolors="k", color="red",  alpha=0.9)
            plt.scatter(t_rel[below_threshold], rate[below_threshold], s=14, marker="o", linewidths=0.5, edgecolors="k", color="gray", alpha=0.6)
            plt.scatter(t_rel[normal],          rate[normal],          s=14, marker="o", linewidths=0.5, edgecolors="k", color="black", alpha=0.8)

            plt.errorbar(
                t_rel,
                rate,
                yerr=error,
                fmt="none",
                ecolor="0.4",
                elinewidth=0.8,
                capsize=0,
                alpha=0.35
            )

            plt.xlabel("Time since start (s)")
            plt.ylabel("Count rate (/s)")
            plt.title(f"{self.star_name} {obsid}")



            output_png = self.images / "Time_Plots" / f"barycentered_{self.star_name}_{obsid}.png"
            plt.grid(True, which="major", alpha=0.2)
            plt.grid(True, which="minor", alpha=0.1)
            plt.tight_layout()
            output_png.parent.mkdir(parents=True, exist_ok=True)
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
    if stats.size == 0 or stats.shape[1] < 4:
        self.logger.info(f"{stats_file} is empty, unable to run Unbinned_Plot")
        return

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

    

            header = plot[1].header
            mjdref = header["MJDREFI"] + header["MJDREFF"]
            phase = self._Phase_Calculator(time_array, mjdref)



            t = time_array
            if np.nanmin(t) < 1e6:
                t = t + header["TSTART"]

            phase = self._Phase_Calculator(t, mjdref)

            
            flares_filter  = (rate - median_val) > three_sigma
            below_threshold = (median_val - rate) > three_sigma
            normal = ~(flares_filter | below_threshold)

            

            phase_plot = np.concatenate([phase, phase + 1.0])
            rate_plot  = np.concatenate([rate, rate])
            error_plot = np.concatenate([error, error])

            flares_plot = np.concatenate([flares_filter, flares_filter])
            below_plot  = np.concatenate([below_threshold, below_threshold])
            normal_plot = np.concatenate([normal, normal])

            plt.scatter(phase_plot[flares_plot], rate_plot[flares_plot], s=14, marker="o", linewidths=0.5, edgecolors="k", color="red",  alpha=0.9)
            plt.scatter(phase_plot[below_plot],  rate_plot[below_plot],  s=14, marker="o", linewidths=0.5, edgecolors="k", color="gray", alpha=0.6)
            plt.scatter(phase_plot[normal_plot], rate_plot[normal_plot], s=14, marker="o", linewidths=0.5, edgecolors="k", color="black", alpha=0.8)

            plt.errorbar(
                phase_plot,
                rate_plot,
                yerr=error_plot,
                fmt="none",
                ecolor="0.4",
                elinewidth=0.8,
                capsize=0,
                alpha=0.35
            )

            plt.xlabel("Orbital phase")
            plt.ylabel("Count rate (/s)")
            plt.title(f"{self.star_name} {obsid} phase folded (unbinned)")
            plt.xlim(0, 2)

            output_png = self.images / "Unbinned_Plots" / f"phase_folded_unbinned_{self.star_name}_{obsid}.png"
            plt.grid(True, which="major", alpha=0.2)
            plt.grid(True, which="minor", alpha=0.1)
            plt.tight_layout()
            output_png.parent.mkdir(parents=True, exist_ok=True)
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

            t = time_array
            if np.nanmin(t) < 1e6:
                t = t + header["TSTART"]

            phase = self._Phase_Calculator(t, mjdref)

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
                markersize=4,
                markerfacecolor="white",
                markeredgecolor="black",
                markeredgewidth=0.8,
                ecolor="0.35",
                elinewidth=0.9,
                capsize=0,
                alpha=0.95
            )
            
            plt.xlabel("Orbital phase")
            plt.ylabel("Mean count rate (/s)")
            plt.title(f"{self.star_name} {obsid} phase folded (binned)")
            plt.xlim(0, 1)

            output_png = self.images / "Binned_Plots" / f"phase_folded_binned_{self.star_name}_{obsid}.png"
            plt.grid(True, which="major", alpha=0.2)
            plt.grid(True, which="minor", alpha=0.1)
            plt.tight_layout()
            output_png.parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(output_png)
            plt.clf()
            plt.close()


def Stacked_Plot(self):
    self._Refresh_ObsID()
    self.images = self.star_folder / "Outputs"
    (self.images / "Stacked_Plots").mkdir(parents=True, exist_ok=True)

    self._Epoch_Binning()
    overlay_centers = []
    overlay_rates = []
    overlay_errors = []
    overlay_labels = []
    epoch_groups = {}
    
    for obsid in self.ObsID_array:
        mjd = self.obsid_epoch.get(obsid)
        if mjd is None:
            continue
        epoch_id = int(mjd // self.epoch_width)
        epoch_groups.setdefault(epoch_id, []).append(obsid)


    overlay_txt = self.images / "Stacked_Plots" / f"epoch_overlay_{self.star_name}.txt"


    with open(overlay_txt, "w") as f:
        f.write("# Epoch overlay file\n")
        f.write("# Each epoch block: PhaseCenter  Rate  Error\n")
        f.write("#\n")
     

    for epoch_id in sorted(epoch_groups.keys()):
        obsids = epoch_groups[epoch_id]
        epoch_start = epoch_id * self.epoch_width
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

                t = time_array
                if np.nanmin(t) < 1e6:
                    t = t + header["TSTART"]

                phases = self._Phase_Calculator(t, mjdref)

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

        output_txt = self.images / "Stacked_Plots" /f"stacked_epoch_{self.star_name}_{epoch_id}.txt"
        np.savetxt(
            output_txt,
            np.column_stack([all_obsids, all_phases, all_rates, all_errors]),
            fmt="%s",
            header="ObsID Phase Rate Error"
        )

        centers, binned_rate, binned_err = self._Phase_Binner(all_phases, all_rates, all_errors)




        with open(overlay_txt, "a") as f:
            f.write(f"# Epoch {epoch_id}  MJD {epoch_start:.0f}-{epoch_start+self.epoch_width:.0f}\n")
            f.write("# PhaseCenter  Rate  Error\n")
            m = np.isfinite(binned_rate) & np.isfinite(binned_err)
            np.savetxt(
                f,
                np.column_stack([centers[m], binned_rate[m], binned_err[m]]),
                fmt="%.8f %.8f %.8f"
            )
            f.write("\n")


        
        overlay_centers.append(centers)
        overlay_rates.append(binned_rate)
        overlay_errors.append(binned_err)
        overlay_labels.append(f"MJD {epoch_start:.0f}-{epoch_start+self.epoch_width:.0f}")

        fit = self._Comparison_Model(centers, binned_rate, binned_err)
        if fit is not None:
            C, A, phi0, chi2, dof = fit
            self.logger.info(f"Epoch {epoch_id} chi2/dof = {chi2/dof:.3f} (dof={dof})")

        m = np.isfinite(binned_rate) & np.isfinite(binned_err)
        
        plt.errorbar(
            centers[m],
            binned_rate[m],
            yerr=binned_err[m],
            fmt="o",
            markersize=4,
            markerfacecolor="white",
            markeredgecolor="black",
            markeredgewidth=0.8,
            ecolor="0.35",
            elinewidth=0.9,
            capsize=0,
            alpha=0.95
        )
        plt.xlabel("Orbital phase")
        plt.ylabel("Mean count rate (/s)")
        plt.title(f"{self.star_name} stacked phase-folded light curve (Epoch start MJD {epoch_start:.2f})")
        plt.xlim(0, 1)

        output_png = self.images / "Stacked_Plots" / f"stacked_epoch_{self.star_name}_{epoch_id}.png"
        plt.grid(True, which="major", alpha=0.2)
        plt.grid(True, which="minor", alpha=0.1)
        plt.tight_layout()
        output_png.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_png)
        plt.clf()
        plt.close()
        self.logger.info(f"Stacked phase-folded curve saved as {output_png} and data saved as {output_txt}")

    if len(overlay_centers) > 0:

        plt.figure()

        for centers, rate, err, label in zip(
            overlay_centers,
            overlay_rates,
            overlay_errors,
            overlay_labels
        ):

            m = np.isfinite(rate) & np.isfinite(err)

            plt.errorbar(
                centers[m],
                rate[m],
                yerr=err[m],
                fmt="o",
                markersize=4,
                markerfacecolor="white",
                markeredgewidth=0.8,
                capsize=0,
                elinewidth=0.9,
                alpha=0.95,
                label=label
            )

        plt.xlabel("Orbital phase")
        plt.ylabel("Mean count rate (/s)")
        plt.title(f"{self.star_name} epoch overlay")
        plt.xlim(0, 1)

        plt.grid(True, which="major", alpha=0.2)
        plt.grid(True, which="minor", alpha=0.1)

        plt.legend(fontsize=8)
        plt.tight_layout()

        overlay_png = self.images / "Stacked_Plots" / f"epoch_overlay_{self.star_name}.png"

        overlay_png.parent.mkdir(parents=True, exist_ok=True)

        plt.savefig(overlay_png)
        plt.clf()
        plt.close()

        self.logger.info(f"Epoch overlay saved as {overlay_png}")
