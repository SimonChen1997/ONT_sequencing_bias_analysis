import os
import pandas as pd

## store the chromosome name in a list
chrom = list(range(1,31))

## function to calculate the insertion number of each region of sequencing kits
def insertion_site_v2(window_bed, insertion_file):
    insertion_in_windows = pd.DataFrame(columns=["sequencing_kit", "chromosome", "sample_id",
                                                 "insertion_number", "start", "end", "gc"])

    ## get the values the row with index i
    for i in range(len(window_bed)):
        chromo = window_bed.loc[i, "chromosome"]
        start = float(window_bed.loc[i, "start"])
        end = float(window_bed.loc[i, "end"])
        gc = float(window_bed.loc[i, "gc_percentage"])

        ## target the defined region
        region_data = insertion_file[(insertion_file["chromosome"] == chromo) &
                                     (insertion_file["insertion"] >= start) &
                                     (insertion_file["insertion"] < end)]

        ## perform the calculation
        region_summary = region_data.groupby(["sequencing_kit", "chromosome", "sample_id"])[
            "insertion"].count().reset_index(name="insertion_number")

        region_summary["start"] = start
        region_summary["end"] = end
        region_summary["gc"] = gc
        insertion_in_windows = pd.concat([insertion_in_windows, region_summary])

    return insertion_in_windows


for i in chrom:
    data_folder = f"/path/TSU_30_insertion_site/chr_{i}"
    os.chdir(data_folder)  # Change the current working directory

    ## read necessary files
    insertion_file = pd.read_csv(f"lsk_rbk_insertion_site_chr{i}.tsv", sep="\t", header=None,
                                names=["chromosome", "insertion", "sample_id", "sequencing_kit"])
    window_bed = pd.read_csv(f"reference_10kb_GC_chr{i}.tsv", sep="\t", header=None,
                               names=["chromosome", "start", "end", "gc_percentage"])

    insertion_position = insertion_site_v2(window_bed, insertion_file)
    print(f"lsk_rbk_insertion_site_chr{i}.tsv is done")

    output_file = f"/path/chr_{i}_insertion_gc_window_python.tsv"
    insertion_position.to_csv(output_file, sep="\t", index=False)