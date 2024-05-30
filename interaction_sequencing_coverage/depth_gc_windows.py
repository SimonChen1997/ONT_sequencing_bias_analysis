import os
import pandas as pd

## store the chromosome name in a list
chrom = list(range(1, 31))

## function to calculate the sequencing coverage of each region of sequencing kits
def sequencing_depth_windows(window_bed, depth_file):
    depth_in_windows = pd.DataFrame(columns=["sequencing_kit", "chromosome", "sample_id",
                                             "total_depth", "start", "end", "gc"])

    for i in range(len(window_bed)):
        chrom = window_bed.loc[i, "chromosome"]
        start_position = float(window_bed.loc[i, "start"])
        end_position = float(window_bed.loc[i, "end"])
        gc = float(window_bed.loc[i, "gc_percentage"])

        region_data = depth_file[(depth_file["chromosome"] == chrom) &
                                 (depth_file["position"] >= start_position) &
                                 (depth_file["position"] < end_position)]

        region_data_summary = region_data.groupby(["sequencing_kit", "chromosome", "sample_id"])[
            "depth"].sum().reset_index(name="total_depth")

        region_data_summary["start"] = start_position
        region_data_summary["end"] = end_position
        region_data_summary["gc"] = gc

        depth_in_windows = pd.concat([depth_in_windows, region_data_summary])

    return pd.DataFrame(depth_in_windows)

for i in chrom:
    data_folder = f"/path/TSU_30_sequencing_depth/chr_{i}"
    os.chdir(data_folder)  # Change the current working directory

    LSK_RBK_depth = pd.read_csv(f"lsk_rbk_sequencing_depth_chr{i}.tsv", sep="\t", header=None,
                                names=["chromosome", "position", "depth", "sample_id", "sequencing_kit"])
    chrom_window = pd.read_csv(f"reference_10kb_GC_chr{i}.tsv", sep="\t", header=None,
                               names=["chromosome", "start", "end", "gc_percentage"])

    summary_depth = sequencing_depth_windows(chrom_window, LSK_RBK_depth)
    print(f"lsk_rbk_sequencing_depth_chr{i}.tsv is done")

    output_file = f"/path/chr_{i}_depth_gc_window_python.tsv"
    summary_depth.to_csv(output_file, sep="\t", index=False)



    