import json
import sys
from program.pipeline_core import Pipeline


if __name__ == "__main__":
   
    config_file = sys.argv[1] if len(sys.argv) > 1 else "config.json"

    with open(config_file, "r") as f:
        config = json.load(f)

    pipeline = Pipeline(
        star_name=config["star_name"],
        orbital_period=config["orbital_period"],
        orbital_pdot=config["orbital_pdot"],
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
