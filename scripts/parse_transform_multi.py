import os
import pandas as pd
from upsetplot import from_memberships
import itertools
import argparse
import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
from venn import venn
from upsetplot import UpSet
import warnings

random.seed(1234)
warnings.simplefilter(action='ignore', category=FutureWarning)






def calculate_median_length(directory, primeSS):
    """
    Calculate the median number of rows in files of specified folder, for specific 5' or 3'.
    """
    prime = '.ss' + primeSS
    line_counts = []

    for filename in os.listdir(directory):
        if filename.endswith(prime):
            file_path = os.path.join(directory, filename)
            try:
                with open(file_path, 'r') as f:
                    line_count = sum(1 for _ in f) - 1 
                    if line_count > 0: 
                        line_counts.append(line_count)
            except Exception as e:
                print(f"Error reading {filename}: {e}")

    if not line_counts:
        print(f"Warning: No files matching {prime} found in the directory or all files have zero lines. Skipping...")
        return None

    # Calculate and return the median. Used for automatically finding a cutoff. (If not sure, initiate with this to get an understanding of the data)
    median_length = int(np.median(line_counts))
    print(f"Calculated median length for {prime}: {median_length}")
    print(f"Samples with the most splice sites: {sorted(line_counts, reverse=True)[:5]}")
    return median_length






def get_label_SS(outdir, directory, sub_df, length, primeSS='3', label=str(), cutoff=10):
    """
    Finds common splice sites between the samples of a given condition.
    """
    common_splice_sites = None
    count = 0
    prime = '.ss' + primeSS

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    subtype_files = sub_df['file_name'].tolist()
    processed_count = 0

    files_with_line_counts = []
    for filename in os.listdir(directory):
        if filename.endswith(prime) and filename.split('.')[0] in subtype_files:
            file_path = os.path.join(directory, filename)
            try:
                with open(file_path, 'r') as f:
                    line_count = sum(1 for _ in f) - 1 
                    files_with_line_counts.append((filename, line_count))
            except Exception as e:
                print(f"Error reading {filename}: {e}")

    # Sort files by line count in descending order (number of splice sites)
    files_with_line_counts.sort(key=lambda x: x[1], reverse=True)
    print(f"Processing files for label '{label}' sorted by splice site count: {[f[0] for f in files_with_line_counts]}")

    for filename, _ in files_with_line_counts:
        file_path = os.path.join(directory, filename)
        processed_count += 1

        try:
            data = pd.read_csv(file_path, sep='\t')
            splice_sites = data.apply(lambda row: f"{row['chrStart']}:{row['chrEnd']}", axis=1).tolist()

            if len(splice_sites) < length:
                count += 1
                print(f"Skipping {filename} due to insufficient splice sites ({len(splice_sites)} < {length}).")
                continue

            if common_splice_sites is None:
                common_splice_sites = set(splice_sites)
                print(f"Initialized common splice sites with {len(common_splice_sites)} sites from {filename}.")
            else:
                before_update = len(common_splice_sites)
                previous_splice_sites = common_splice_sites.copy()  # Save the current list before intersection
                common_splice_sites.intersection_update(splice_sites)
                after_update = len(common_splice_sites)

                if after_update == 0:
                    print(f"Intersection with {filename} resulted in 0 splice sites. Keeping the previous list.")
                    common_splice_sites = previous_splice_sites  # Revert to the previous list if 0 intersection -- failsafe
                else:
                    print(f"Updated common splice sites with {filename}: {before_update} -> {after_update} sites.")
        except Exception as e:
            print(f"Error processing {filename}: {e}")

    print(f"Skipped {count} samples for label '{label}' due to insufficient splice sites.")
    print(f"Final common splice sites for label '{label}': {len(common_splice_sites) if common_splice_sites else 0}")

    if common_splice_sites:
        file_path = os.path.join(outdir, f'common_{label}_SS{primeSS}_L{length}.list')
        with open(file_path, 'w') as file:
            file.write(f'Common splice sites between samples for {label}\n')
            for item in common_splice_sites:
                file.write(f"{item}\n")
        print(f"Common splice sites for {label} saved at {file_path}")

    return common_splice_sites, count, processed_count






def create_SE_df(outdir, common_splice_sites, directory, sub_df, length, primeSS='3', label=str(), subtype_column=None):
    """
    Create a DataFrame for downstream analysis, finding common splice sites between conditions.
    """
    if not common_splice_sites:
        print(f"Warning: No common splice sites found for {label}. Skipping...")
        return None

    prime = '.ss' + primeSS
    data_df = pd.DataFrame(columns=list(common_splice_sites))

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    subtype_files = sub_df['file_name'].tolist()

    for filename in os.listdir(directory):
        if filename.endswith(prime) and filename.split('.')[0] in subtype_files:
            sample_name = filename.split('.')[0]
            file_path = os.path.join(directory, filename)

            try:
                data = pd.read_csv(file_path, sep='\t')
                splice_sites = data.apply(lambda row: f"{row['chrStart']}:{row['chrEnd']}", axis=1)

                if len(splice_sites) < length:
                    continue

                filtered_data = data[splice_sites.isin(common_splice_sites)].copy()
                temp_df = pd.DataFrame([filtered_data['splice_eff'].values],
                                       index=[sample_name],
                                       columns=splice_sites[splice_sites.isin(common_splice_sites)].values).reindex(columns=common_splice_sites)
                if not temp_df.empty and not temp_df.isna().all(axis=None):
                    if subtype_column and subtype_column in sub_df.columns:
                        temp_df[subtype_column] = sub_df.loc[sub_df['file_name'] == sample_name, subtype_column].values[0] if not sub_df.loc[sub_df['file_name'] == sample_name, subtype_column].empty else None
                    data_df = pd.concat([data_df, temp_df], ignore_index=False)
            except Exception as e:
                print(f"Error processing {filename}: {e}")

    if not data_df.empty:
        data_df = data_df.dropna(how='all', axis=1)
        data_df['label'] = label
        output_file = os.path.join(outdir, f'{label}_SE_SS{primeSS}_L{length}.tsv')
        data_df.to_csv(output_file, sep='\t', index=True)
        print(f"Saved file at {output_file}\nShape: {data_df.shape}")
    else:
        print(f"Warning: DataFrame for label '{label}' is empty. Skipping...")
    return data_df






def prepare_upset_data(upset_data):
    """
    Converts a dictionary of sets into a format suitable for UpSet plots.
    """
    if not upset_data:
        print("Warning: No data provided for UpSet plot. Skipping...")
        return None

    memberships = []

    all_splice_sites = sorted(set.union(*upset_data.values()))
    for splice_site in all_splice_sites:
        memberships.append({
            "Splice Site": splice_site,
            "Membership": [subtype for subtype, sites in upset_data.items() if splice_site in sites]
        })

    upset_series = from_memberships(
        [entry["Membership"] for entry in memberships],
        data=[entry["Splice Site"] for entry in memberships]
    )
    return upset_series






def get_label_SS_single(outdir, directory, length, primeSS, label, drop_threshold, individual_drop, cutoff, drop_cut):
    """
    Subfunction for single mode to find common splice sites with additional thresholds.
    """
    prime = '.ss' + primeSS
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    sample_dict = {}
    for f in os.listdir(directory):
        if f.endswith(prime):
            path = os.path.join(directory, f)
            d = pd.read_csv(path, sep='\t')
            splice_sites = d.apply(lambda r: f"{r['chrStart']}:{r['chrEnd']}", axis=1).tolist()
            if len(splice_sites) >= length:
                sample_dict[f.split('.')[0]] = set(splice_sites)

    if not sample_dict:
        print(f"[get_label_SS_single] {label}: No samples found with >= {length} sites.")
        return None

    largest_sample = max(sample_dict, key=lambda x: len(sample_dict[x]))
    common_splice_sites = set(sample_dict[largest_sample])
    allowed_drop = len(common_splice_sites) - drop_threshold
    print(f"[get_label_SS_single] {label}: largest_sample={largest_sample}, base_sites={len(common_splice_sites)}")
    print(f"[get_label_SS_single] {label}: allowed_drop={allowed_drop}, individual_drop={individual_drop}, drop_cut={drop_cut}")

    del sample_dict[largest_sample]
    sorted_samples = sorted(sample_dict.keys(), key=lambda x: len(sample_dict[x]), reverse=True)
    print(f"[get_label_SS_single] {label}: will process {len(sorted_samples)} more samples, largest->smallest.")

    for smp in sorted_samples:
        old_size = len(common_splice_sites)
        proposed = common_splice_sites.intersection(sample_dict[smp])
        new_size = len(proposed)

        # Check if the intersection violates the drop_cut threshold
        if new_size < drop_cut:
            print(f"[get_label_SS_single] {label}: skipping {smp}, would drop below drop_cut={drop_cut}.")
            continue

        # Apply allowed_drop and individual_drop logic
        if old_size > allowed_drop:
            if new_size < allowed_drop:
                print(f"[get_label_SS_single] {label}: skipping {smp}, would drop below allowed_drop={allowed_drop}.")
            else:
                common_splice_sites = proposed
                print(f"[get_label_SS_single] {label}: intersected with {smp}, from {old_size} -> {new_size}")
        else:
            diff = old_size - new_size
            if diff <= individual_drop:
                common_splice_sites = proposed
                print(f"[get_label_SS_single] {label}: intersected with {smp} using individual_drop, {old_size} -> {new_size}")
            else:
                print(f"[get_label_SS_single] {label}: skipped {smp}, diff {diff} > individual_drop={individual_drop}")

    # Final check for drop_cut after all intersections
    if len(common_splice_sites) < drop_cut:
        print(f"[get_label_SS_single] {label}: final intersection drops below drop_cut={drop_cut}. Returning None.")
        return None

    list_path = os.path.join(outdir, f"common_{label}_SS{primeSS}_L{length}_ad{drop_threshold}id{individual_drop}dc{drop_cut}c{cutoff}.list")
    with open(list_path, 'w') as file:
        for site in common_splice_sites:
            file.write(f"{site}\n")

    print(f"[get_label_SS_single] {label}: final intersection has {len(common_splice_sites)} sites, saved at {list_path}")
    return common_splice_sites






def create_SE_df_single(outdir, common_splice_sites, directory, length, primeSS, label, drop_threshold, individual_drop, cutoff, subtype_column=None):
    """
    Subfunction for single mode to create a DataFrame for downstream analysis.
    """
    prime = '.ss' + primeSS
    cols = list(common_splice_sites) if common_splice_sites else []
    df = pd.DataFrame(columns=cols)

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    included_samples = 0
    for f in os.listdir(directory):
        if f.endswith(prime):
            sample_id = f.split('.')[0]
            if sample_id not in valid_file_names:
                continue
            else:
                path = os.path.join(directory, f)
                d = pd.read_csv(path, sep='\t')
                s = d.apply(lambda r: f"{r['chrStart']}:{r['chrEnd']}", axis=1)
                if len(s) < length:
                    continue

            intersect_sites = set(s).intersection(common_splice_sites)
            if len(intersect_sites) < len(common_splice_sites):
                continue

            flt = d[s.isin(common_splice_sites)].copy()
            c = s[s.isin(common_splice_sites)]
            tmp = pd.DataFrame([flt['splice_eff'].values], index=[f.split('.')[0]], columns=c.values)
            tmp = tmp.reindex(columns=cols)

            if tmp.isna().any().any():
                continue

            # Append subtype column if provided
            if subtype_column:
                subtype_value = split_dfs[label].loc[split_dfs[label]['file_name'] == f.split('.')[0], subtype_column]
                tmp[subtype_column] = subtype_value.values[0] if not subtype_value.empty else None

            df = pd.concat([df, tmp])
            included_samples += 1

    df['label'] = label
    tsv_path = os.path.join(outdir, f"{label}_SE_SS{primeSS}_L{length}_ad{drop_threshold}id{individual_drop}dc{drop_cut}c{cutoff}.tsv")
    df.to_csv(tsv_path, sep='\t', index=True)
    print(f"[create_SE_df_single] {label}: shape={df.shape}, included_samples={included_samples}, columns={len(cols)}")
    print(f"[create_SE_df_single] {label}: saved at {tsv_path}")

    return df






def load_splice_site_list(filepath):
    """
    Load a list of splice sites from a file (one per line).
    """
    with open(filepath, 'r') as f:
        sites = [line.strip() for line in f if line.strip()]
    return set(sites)






def process_single_selected_mode(outdir, directory, splice_site_list, ss, length, label):
    """
    Process the single_selected mode: include only files that contain ALL splice sites in the provided list.
    """
    prime = '.ss' + ss
    included_samples = 0
    checked_samples = 0
    included_files = []
    excluded_files = []
    cols = list(splice_site_list)
    df = pd.DataFrame(columns=cols)

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for f in os.listdir(directory):
        if f.endswith(prime):
            path = os.path.join(directory, f)
            try:
                d = pd.read_csv(path, sep='\t')
                s = d.apply(lambda r: f"{r['chrStart']}:{r['chrEnd']}", axis=1)
                checked_samples += 1
                if not splice_site_list.issubset(set(s)):
                    excluded_files.append(f)
                    print(f"[single_selected] Excluding {f}: does not contain all {len(splice_site_list)} provided splice sites.")
                    continue
                flt = d[s.isin(splice_site_list)].copy()
                c = s[s.isin(splice_site_list)]
                tmp = pd.DataFrame([flt['splice_eff'].values], index=[f.split('.')[0]], columns=c.values)
                tmp = tmp.reindex(columns=cols)
                if tmp.isna().any().any():
                    excluded_files.append(f)
                    print(f"[single_selected] Excluding {f}: contains NaN values for provided splice sites.")
                    continue
                df = pd.concat([df, tmp])
                included_samples += 1
                included_files.append(f)
                print(f"[single_selected] Including {f}: contains all {len(splice_site_list)} provided splice sites (may have more).")
            except Exception as e:
                print(f"[single_selected] Error processing {f}: {e}")

    print(f"[single_selected] Checked {checked_samples} files. Included: {included_samples}, Excluded: {len(excluded_files)}")
    if included_files:
        print(f"[single_selected] Included files: {included_files}")
    if excluded_files:
        print(f"[single_selected] Excluded files: {excluded_files}")

    # Save .tsv
    df['label'] = label
    tsv_path = os.path.join(outdir, f"{label}_SE_SS{ss}_L{length}_single_selected.tsv")
    df.to_csv(tsv_path, sep='\t', index=True)
    print(f"[single_selected] {label}: shape={df.shape}, included_samples={included_samples}, columns={len(cols)}")
    print(f"[single_selected] {label}: saved at {tsv_path}")

    # Save .list
    list_path = os.path.join(outdir, f"common_{label}_SS{ss}_L{length}_single_selected.list")
    with open(list_path, 'w') as file:
        for site in splice_site_list:
            file.write(f"{site}\n")
    print(f"[single_selected] {label}: saved list at {list_path}")

    return df






def process_modes(directory, outdir, split_dfs, modes, ss_tup, cutoff, subtype_column=None, condition_identifier=None, drop_threshold=10, individual_drop=2, drop_cut=5, single_selected_file=None):
    """
    Process the modes ('all', 'binary', 'single', 'single_selected').
    """
    summary_info = []

    for mode in modes:
        print(f"Processing mode: {mode}")
        for ss, length in ss_tup:
            print(f"Processing splice site type: SS{ss}, Length: {length}")
            if int(length) == 0:
                median_length = calculate_median_length(directory, ss)
                if (median_length is None):
                    continue
                length_to_use = median_length
            else:
                length_to_use = int(length)

            mode_outdir = os.path.join(outdir, f'SS{ss}/common_{mode}')
            os.makedirs(mode_outdir, exist_ok=True)

            dynamic_prefix = os.path.basename(outdir).split('_')[0]  # Keep dynamic_prefix

            if mode == "single_selected":
                if not single_selected_file:
                    print("Error: --single_selected_file must be provided for single_selected mode.")
                    continue
                splice_site_list = load_splice_site_list(single_selected_file)
                label = "single_selected"
                print(f"Processing single_selected mode with {len(splice_site_list)} splice sites.")
                df = process_single_selected_mode(
                    mode_outdir, directory, splice_site_list, ss, length_to_use, label
                )
                summary_info.append({
                    "Mode": mode,
                    "Condition": label,
                    "Splice_Sites": len(splice_site_list),
                    "File_Length": length_to_use,
                    "Included_Samples": df.shape[0] if df is not None else 0
                })
            elif mode == "single":
                print("Processing single mode...")
                if condition_identifier:
                    # Filter the DataFrame for the specific condition
                    if condition_identifier not in concat_df[condition_column].unique():
                        print(f"Error: Condition '{condition_identifier}' not found in the condition column.")
                        continue
                    single_df = concat_df[concat_df[condition_column] == condition_identifier]
                    print(f"Filtered DataFrame for condition '{condition_identifier}': {single_df.shape}")
                else:
                    # Combine all conditions into a single condition
                    single_df = concat_df
                    print(f"Using entire manifest for single mode: {single_df.shape}")

                # Ensure only rows matching the condition_identifier are included
                valid_file_names = set(single_df['file_name'])
                print(f"Valid file names for condition '{condition_identifier}': {len(valid_file_names)}")

                filtered_files = filter_directory_files(directory, valid_file_names, ss)
                print(f"Filtered files in directory for condition '{condition_identifier}': {len(filtered_files)}")

                # Debug: Check for files in the directory that are not in the manifest
                directory_files = [f.split('.')[0] for f in os.listdir(directory) if f.endswith('.ss' + ss)]
                extra_files = set(directory_files) - valid_file_names

                print(f"Processing single condition: {condition_identifier}")
                print(f"Working on {ss} with readSum > {cutoff}, minimum splice site entries={length_to_use}")

                common_splice_sites = get_label_SS_single(
                    mode_outdir, directory, length=length_to_use, primeSS=ss, label=condition_identifier,
                    drop_threshold=drop_threshold, individual_drop=individual_drop, cutoff=cutoff, drop_cut=drop_cut
                )

                if common_splice_sites:
                    print(f"[process_modes] {condition_identifier}: Found {len(common_splice_sites)} common splice sites.")
                    create_SE_df_single(
                        mode_outdir, common_splice_sites, directory, length_to_use, primeSS=ss, label=condition_identifier,
                        drop_threshold=drop_threshold, individual_drop=individual_drop, cutoff=cutoff, subtype_column=subtype_column
                    )
                else:
                    print(f"[process_modes] {condition_identifier}: No common splice sites found.")

                summary_info.append({
                    "Mode": mode,
                    "Condition": condition_identifier,
                    "Splice_Sites": len(common_splice_sites) if common_splice_sites else 0,
                    "File_Length": length_to_use
                })
            elif mode == "all":
                common_across_all = None
                concatenated_df = pd.DataFrame()
                for condition, condition_df in split_dfs.items():
                    if len(split_dfs.items()) == 2:
                        print("Warning: Only two conditions found. Use binary mode. Skipping...")
                        return
                    print(f"Processing condition: {condition}")
                    common_splice_sites, skipped_count, processed_count = get_label_SS(
                        mode_outdir, directory, condition_df, length=length_to_use, primeSS=ss, label=condition
                    )
                    if common_across_all is None:
                        common_across_all = common_splice_sites
                        print(f"Initialized common splice sites across all conditions with {len(common_across_all) if common_across_all else 0} sites.")
                    else:
                        before_update = len(common_across_all)
                        common_across_all.intersection_update(common_splice_sites)
                        after_update = len(common_across_all)
                        print(f"Updated common splice sites across all conditions: {before_update} -> {after_update} sites.")

                    subtype_df = create_SE_df(
                        mode_outdir, common_across_all, directory, condition_df, length=length_to_use, primeSS=ss, label=condition, subtype_column=subtype_column
                    )
                    if subtype_df is not None:
                        concatenated_df = pd.concat([concatenated_df, subtype_df], axis=0)

                    summary_info.append({
                        "Mode": mode,
                        "Subtype/Combination": condition,
                        "Files_Processed": processed_count,
                        "Files_Skipped": skipped_count,
                        "Final_Samples": processed_count - skipped_count,
                        "Splice_Sites": len(common_splice_sites) if common_splice_sites else 0,
                        "File_Length": length_to_use
                    })

                # Save common splice sites across all conditions, keeping a clean non-NA dataframe (this is our final DF)
                if not concatenated_df.empty:
                    concatenated_df = concatenated_df.dropna(how='all', axis=1)
                    concat_output_file = os.path.join(mode_outdir, f'{condition_column}_{dynamic_prefix}_{mode}_conditions_concat_SE_SS{ss}_L{length_to_use}.tsv')
                    concatenated_df.to_csv(concat_output_file, sep='\t', index=True)
                    print(f"Saved concatenated file for mode 'all' at {concat_output_file}\nShape: {concatenated_df.shape}")

                if common_across_all:
                    common_file = os.path.join(mode_outdir, f'{condition_column}_{dynamic_prefix}common_across_all_conditions_SS{ss}_L{length_to_use}.list')
                    with open(common_file, 'w') as file:
                        file.write("Common splice sites across all conditions:\n")
                        for item in common_across_all:
                            file.write(f"{item}\n")
                    print(f"Saved common splice sites across all conditions at {common_file}")

                # Venn diagram
                venn_data = {condition: get_label_SS(mode_outdir, directory, condition_df, length=length_to_use, primeSS=ss, label=condition)[0]
                             for condition, condition_df in split_dfs.items()}
                venn_dict = {key: value for key, value in venn_data.items() if value is not None and len(value) > 0}  # Ensure non-empty sets

                if len(venn_dict) > 1:
                    plt.figure(figsize=(12, 10))
                    venn(venn_dict)
                    venn_file = os.path.join(mode_outdir, f'{condition_column}_{dynamic_prefix}_Venn_all_SE_SS{ss}_L{length_to_use}.pdf')
                    plt.title(f"Venn Diagram for All Conditions (SS={ss}, L={length_to_use})")
                    plt.savefig(venn_file)
                    print(f"Venn Diagram for all conditions saved at {venn_file}")
                    plt.close()
                else:
                    print("Warning: Venn diagram skipped due to insufficient or empty sets.")

                # UpSet plot
                upset_data = {}
                for condition, condition_df in split_dfs.items():
                    common_splice_sites, _, _ = get_label_SS(
                        mode_outdir, directory, condition_df, length=length_to_use, primeSS=ss, label=condition
                    )
                    upset_data[condition] = common_splice_sites

                print("Subtypes and their splice sites:")
                for subtype, sites in upset_data.items():
                    print(f"{subtype}: {len(sites)} splice sites")

                upset_series = prepare_upset_data(upset_data)

                plt.figure(figsize=(20, 20))
                upset = UpSet(
                    upset_series,
                    intersection_plot_elements=10,
                    facecolor="navy",
                    shading_color="lightgrey",
                    with_lines=True,
                    show_counts=True,
                    show_percentages=True,
                    totals_plot_elements=2,
                )

                upset.plot()
                upset_file = os.path.join(mode_outdir, f'{condition_column}_{dynamic_prefix}_UpSet_all_SE_SS{ss}_L{length_to_use}.pdf')
                plt.savefig(upset_file, dpi=300, bbox_inches="tight")
                print(f"UpSet plot for all conditions saved at {upset_file}")
                plt.close()

            elif mode == "binary":
                combinations = list(itertools.combinations(split_dfs.keys(), 2))
                for combo in combinations:
                    condition1, condition2 = combo
                    combo_label = f"{condition1}_vs_{condition2}"
                    print(f"Processing binary combination: {condition1} vs {condition2}")
                    common_splice_sites_1, skipped_1, processed_1 = get_label_SS(
                        mode_outdir, directory, split_dfs[condition1], length=length_to_use, primeSS=ss, label=condition1
                    )
                    common_splice_sites_2, skipped_2, processed_2 = get_label_SS(
                        mode_outdir, directory, split_dfs[condition2], length=length_to_use, primeSS=ss, label=condition2
                    )

                    print(f"{condition1} splice sites: {len(common_splice_sites_1) if common_splice_sites_1 else 0}")
                    print(f"{condition2} splice sites: {len(common_splice_sites_2) if common_splice_sites_2 else 0}")

                    if not common_splice_sites_1 or not common_splice_sites_2:
                        print(f"Skipping Venn diagram for {condition1} vs {condition2} due to empty sets.")
                        continue

                    if common_splice_sites_1 == common_splice_sites_2:
                        print(f"Warning: {condition1} and {condition2} have identical splice sites.")
                    elif common_splice_sites_1.isdisjoint(common_splice_sites_2):
                        print(f"Warning: {condition1} and {condition2} have no overlapping splice sites.")

                    # Venn diagram
                    plt.figure(figsize=(10, 10))
                    venn2(
                        subsets=(common_splice_sites_1, common_splice_sites_2),
                        set_labels=(condition1, condition2)
                    )
                    venn_file = os.path.join(mode_outdir, f'{condition_column}_{dynamic_prefix}_Venn_{condition1}_{condition2}_SE_SS{ss}_L{length_to_use}.pdf')
                    plt.title(f"Venn Diagram: {condition1} vs {condition2} (SS={ss}, L={length_to_use})")
                    plt.savefig(venn_file)
                    print(f"Venn Diagram for {condition1} and {condition2} saved at {venn_file}")
                    plt.close()

                    common_combined = common_splice_sites_1.intersection(common_splice_sites_2)

                    df_sub1 = create_SE_df(
                        mode_outdir, common_combined, directory, split_dfs[condition1], length=length_to_use, primeSS=ss, label=condition1, subtype_column=subtype_column
                    )
                    df_sub2 = create_SE_df(
                        mode_outdir, common_combined, directory, split_dfs[condition2], length=length_to_use, primeSS=ss, label=condition2, subtype_column=subtype_column
                    )
                    # Save common splice sites across all conditions, keeping a clean non-NA dataframe (this is our final DF)
                    if df_sub1 is not None and df_sub2 is not None:
                        concatenated_df = pd.concat([df_sub1, df_sub2], axis=0)
                        concat_output_file = os.path.join(mode_outdir, f'{condition_column}_{dynamic_prefix}_concat_{condition1}_{condition2}_SE_SS{ss}_L{length_to_use}.tsv')
                        concatenated_df.to_csv(concat_output_file, sep='\t', index=True)
                        print(f"Saved concatenated file for mode 'binary' at {concat_output_file}\nShape: {concatenated_df.shape}")        

                    summary_info.append({
                        "Mode": mode,
                        "Subtype/Combination": combo_label,
                        "Files_Processed": processed_1 + processed_2,
                        "Files_Skipped": skipped_1 + skipped_2,
                        "Final_Samples": (processed_1 - skipped_1) + (processed_2 - skipped_2),
                        "Splice_Sites": len(common_combined) if common_combined else 0,
                        "File_Length": length_to_use
                    })

    summary_df = pd.DataFrame(summary_info)
    summary_file = os.path.join(outdir, f'{dynamic_prefix}_summary.csv')
    summary_df.to_csv(summary_file, index=False, sep='\t')
    print(f"Summary saved to {summary_file}")
    return summary_df






def ensure_directories_exist(base_outdir, modes, ss_tup):
    """
    Ensure all required directories exist before processing.
    """
    for mode in modes:
        for ss, _ in ss_tup:
            mode_outdir = os.path.join(base_outdir, f'SS{ss}/common_{mode}')
            os.makedirs(mode_outdir, exist_ok=True)
            print(f"Ensured directory exists: {mode_outdir}")






if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process splicing efficiency data for multiple conditions.")
    parser.add_argument("--manifest", required=True, help="Path to the manifest file (e.g., concat_manifest_with_details.tsv)")
    parser.add_argument("--directory", required=True, help="Directory containing .ss3/.ss5 files")
    parser.add_argument("--outdir", required=True, help="Output directory for results")
    parser.add_argument("--condition_column", required=True, help="Column in the manifest file for conditions (e.g., pam50_subtype, condition)")
    parser.add_argument("--subtype_column", default=None, help="Optional column in the manifest file to include in the output (e.g., pam50_subtype)")
    parser.add_argument("--condition_identifier", default=None, help="Specific condition to process in single mode (e.g., LumA)")
    parser.add_argument("--modes", nargs="+", default=["all", "binary", "single"], help="Modes to process (e.g., all, binary, single)")
    parser.add_argument("--ss_tup", nargs="+", default=["3.altSS,0"], help="Splice site and length tuples (e.g., 3.altSS,53). If the number is 0, the median length will be used.")
    parser.add_argument("--cutoff", type=int, default=10, help="Cutoff value for filtering")
    parser.add_argument("--drop_threshold", type=int, default=10, help="Allowed drop threshold")
    parser.add_argument("--individual_drop", type=int, default=2, help="Allowed individual drop")
    parser.add_argument("--drop_cut", type=int, default=0, help="Minimum number of splice sites allowed after intersection in single mode")
    parser.add_argument("--single_selected_file", default=None, help="File with list of splice sites for single_selected mode (one per line)")
    args = parser.parse_args()

    manifest = args.manifest
    directory = args.directory
    outdir = args.outdir
    condition_column = args.condition_column
    subtype_column = args.subtype_column
    condition_identifier = args.condition_identifier
    modes = args.modes
    ss_tup = [tuple(item.split(",")) for item in args.ss_tup]
    cutoff = args.cutoff
    drop_threshold = args.drop_threshold
    individual_drop = args.individual_drop
    drop_cut = args.drop_cut

    # Load and filter the manifest file
    if condition_identifier:
        # Filter the manifest to include only rows matching the condition and condition identifier
        concat_df = pd.read_csv(manifest, sep="\t")
        concat_df = concat_df[concat_df[condition_column] == condition_identifier]
        print(f"Filtered manifest to include only rows matching condition '{condition_identifier}'. Remaining rows: {len(concat_df)}")
    else:
        concat_df = pd.read_csv(manifest, sep="\t")

    # Extract valid file names from the filtered manifest
    valid_file_names = set(concat_df['file_name'].dropna())
    print(f"Valid file names extracted from manifest: {len(valid_file_names)}")

    # Create split_dfs based on unique values in the condition column
    split_dfs = {value: concat_df[concat_df[condition_column] == value] for value in concat_df[condition_column].unique()}
    print(f"Split manifest into {len(split_dfs)} subsets based on condition column.")

    # Ensure only files matching the filtered manifest are processed
    def filter_directory_files(directory, valid_file_names, primeSS):
        """
        Filter files in the directory to match valid file names and primeSS.
        """
        prime = '.ss' + primeSS
        filtered_files = [
            f for f in os.listdir(directory)
            if f.endswith(prime) and f.split('.')[0] in valid_file_names
        ]
        return filtered_files

    # Process modes
    process_modes(
        directory, outdir, split_dfs, modes, ss_tup, cutoff,
        subtype_column=subtype_column,
        condition_identifier=condition_identifier,
        drop_threshold=drop_threshold,
        individual_drop=individual_drop,
        drop_cut=drop_cut,
        single_selected_file=args.single_selected_file
    )