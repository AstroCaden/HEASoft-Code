import numpy as np
import os
import subprocess
from astropy.time import Time
import sys
from concurrent.futures import ThreadPoolExecutor

def _Single_Downloader(self, url: str):
    self.logger.info(f"Downloading {url}")
    subprocess.run(
        [os.path.expanduser("~/download_wget.pl"), url],
        cwd=self.obsids_dir
    )

def Dataset_Download(self, ObsID_wanted, ObsID_excluded, ObsID_Chosen):

    self._Refresh_ObsID()

    heasarc = self.heasarc
    table = heasarc.query_region(self.pos, catalog="nicermastr", radius="3 arcmin")

    target = self.star_name.lower().replace(" ", "").replace("_", "")
    namecol = np.array(table["name"]).astype(str)
    mask = np.array([target in n.lower().replace(" ", "").replace("_", "") for n in namecol])
    table = table[mask]

    today_mjd = Time.now().mjd
    table = table[(table["time"] > 0) & (table["public_date"] <= today_mjd)]
    table.sort("time")


    if "exposure" in table.colnames:
        exp_col = "exposure"
    else:
        self.logger.info("Cannot find exposure column")
        exp_col = None

    if exp_col is not None:
        before = len(table)
        table = table[table[exp_col] >= 5.0]
        removed = before - len(table)
        if removed > 0:
            self.logger.info(f"Excluded {removed} ObsIDs with exposure < 5 seconds.")

    obsids_time = np.array(table["obsid"]).astype(str)
    _, first_idx = np.unique(obsids_time, return_index=True)
    unique_obsids = obsids_time[np.sort(first_idx)]

    if ObsID_excluded is not None and len(ObsID_excluded) > 0:
        if ObsID_excluded == "current":
            ObsID_excluded = [d.name for d in self.obsids_dir.iterdir() if d.is_dir() and d.name.isdigit()]
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

    already_have = set(d.name for d in self.obsids_dir.iterdir() if d.is_dir() and d.name.isdigit())
    target_total = int(ObsID_wanted)
    current_total = len(already_have)
    need = max(0, target_total - current_total)

    if need == 0:
        self.logger.info(f"Already have {current_total} ObsIDs (target {target_total}); nothing to download.")
        return

    candidates = np.array([o for o in unique_obsids if o not in already_have], dtype=str)

    if len(candidates) == 0:
        self.logger.info("No new ObsIDs available to download.")
        return

    n = min(need, len(candidates))
    idx = np.linspace(0, len(candidates) - 1, n, dtype=int)
    selected_obsids = candidates[idx]

    self.logger.info(
        f"Have {current_total}, target {target_total} => need {need}. "
        f"Will download {len(selected_obsids)} ObsIDs: {selected_obsids}"
    )

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
        if len(urls) >= len(selected_obsids):
            break

    if not urls:
        self.logger.info("No download URLs found; nothing to download.")
        return

    max_workers = self.concurrent_downloads
    self.logger.info(f"Downloading {len(urls)} ObsIDs with {max_workers} workers")

    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        pool.map(self._Single_Downloader, urls)

    self._Refresh_ObsID()
    if os.environ.get("PIPELINE_AFTER_DOWNLOAD") != "1":
        os.environ["PIPELINE_AFTER_DOWNLOAD"] = "1"
        os.execv(sys.executable, [sys.executable, "-m", "program.main"] + sys.argv[1:])
