import h5py
import numpy as np
import argparse
import glob

def concatenate_datasets(output_group, datasets, dataset_name):
    """
    Concatenate datasets (either numerical or strings) and write to the output group.
    """
    if isinstance(datasets[0][0], str):  # Check if data is string type
        combined_array = sum(datasets, [])  # Flatten list of lists
        dt = h5py.special_dtype(vlen=str)  # Define variable-length string dtype
        dset = output_group.create_dataset(dataset_name, (len(combined_array),), dtype=dt)
        dset[:] = combined_array  # Assign data to dataset
    else:  # Handle numerical types
        combined_array = np.concatenate(datasets, axis=0)
        output_group.create_dataset(dataset_name, data=combined_array)

def process_group(output_file, group_name, input_files):
    """
    Process a group, recursively handling subgroups and concatenating datasets.
    """
    combined_data = {}
    for i, file_name in enumerate(input_files):
        with h5py.File(file_name, 'r') as f:
            group = f[group_name] if group_name in f else f['/']
            for dset_name, dataset in group.items():
                full_path = f"{group_name}/{dset_name}" if group_name else dset_name
                if isinstance(dataset, h5py.Dataset):
                    if full_path not in combined_data:
                        combined_data[full_path] = []
                    if dataset.dtype.kind in {'O', 'S', 'U'}:  # Object, string types
                        combined_data[full_path].append(list(dataset[:]))
                    else:
                        combined_data[full_path].append(dataset[:])
                elif isinstance(dataset, h5py.Group) and i == 0:
                    # i==0 ensures that we only process subgroups once
                    process_group(output_file, full_path, input_files)

    # Write combined data to the output file
    with h5py.File(output_file, 'a') as f_out:
        for dset_path, data_list in combined_data.items():
            parent_path, dset_name = dset_path.rsplit('/', 1) if '/' in dset_path else ('/', dset_path)
            output_group = f_out.require_group(parent_path)
            concatenate_datasets(output_group, data_list, dset_name)

def concatenate_h5(files, output_file):
    """
    Concatenate HDF5 files, handling nested groups and datasets.
    """
    expanded_files = []
    for pattern in files:
        # Use glob to expand wildcard patterns
        expanded_files.extend(glob.glob(pattern))
    # Remove duplicates that may appear from overlapping patterns
    # unique_files = list(sorted(set(expanded_files)))
    unique_files = expanded_files
    print(f"Concatenating {len(unique_files)} files to {output_file}")
    print("\n".join(unique_files))
    process_group(output_file, '', unique_files)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Concatenate multiple HDF5 files into a single output file, handling nested groups and datasets. All matrix will be concat along axis 0')
    parser.add_argument('files', nargs='+', help='Input HDF5 filenames supporting shell-style wildcards')
    parser.add_argument('-o', '--output', required=True, help='Output HDF5 filename')

    args = parser.parse_args()

    # Ensure the output file is initially empty
    with h5py.File(args.output, 'w') as f:
        pass

    concatenate_h5(args.files, args.output)
