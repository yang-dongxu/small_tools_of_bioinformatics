#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
HDF5 to Zarr and Zarr to HDF5 converter
This script provides functions to convert data between HDF5 and Zarr formats.
"""

import os
import argparse
import numpy as np
import h5py
import zarr
from numcodecs import Blosc

def convert_hdf5_to_zarr(hdf5_file, zarr_store, compression_level=5):
    """
    Convert an HDF5 file to Zarr format
    
    Parameters:
    -----------
    hdf5_file : str
        Path to the HDF5 file
    zarr_store : str
        Path for the output Zarr store
    compression_level : int, optional
        Compression level (0-9) for Zarr datasets
    """
    compressor = Blosc(cname='zstd', clevel=compression_level, shuffle=Blosc.SHUFFLE)
    
    def _copy_attributes(obj, target):
        """Copy attributes from obj to target"""
        for key, value in obj.attrs.items():
            target.attrs[key] = value
    
    def _process_group(group, parent_group):
        """Process HDF5 group recursively"""
        for key, item in group.items():
            if isinstance(item, h5py.Group):
                # Create zarr group
                new_group = parent_group.create_group(key)
                _copy_attributes(item, new_group)
                _process_group(item, new_group)
            elif isinstance(item, h5py.Dataset):
                # Copy dataset
                shape = item.shape
                dtype = item.dtype
                chunks = item.chunks if item.chunks else None
                
                # Create zarr array
                zarr_array = parent_group.create_dataset(
                    key,
                    shape=shape,
                    dtype=dtype,
                    chunks=chunks,
                    compressor=compressor
                )
                
                # Copy data in chunks to handle large datasets
                if shape:  # Only chunk if dataset is not scalar
                    chunk_size = 10_000_000  # 10M elements at a time
                    for slices in _chunk_iterator(shape, chunk_size):
                        zarr_array[slices] = item[slices]
                else:
                    # Handle scalar datasets
                    zarr_array[()] = item[()]
                
                _copy_attributes(item, zarr_array)
    
    def _chunk_iterator(shape, chunk_size):
        """Iterator that yields slice tuples for processing large datasets in chunks"""
        if not shape:
            yield ()
            return
            
        total_size = np.prod(shape)
        if total_size <= chunk_size:
            yield tuple(slice(0, s) for s in shape)
            return
        
        # For simplicity, we'll chunk along the first dimension
        chunk_size_dim0 = max(1, int(chunk_size / np.prod(shape[1:])))
        for i in range(0, shape[0], chunk_size_dim0):
            end = min(i + chunk_size_dim0, shape[0])
            slices = [slice(i, end)]
            slices.extend(slice(0, s) for s in shape[1:])
            yield tuple(slices)
    
    # Create zarr hierarchy
    root = zarr.open(zarr_store, mode='w')
    
    with h5py.File(hdf5_file, 'r') as f:
        _copy_attributes(f, root)
        _process_group(f, root)
    
    print(f"Conversion from HDF5 to Zarr completed: {hdf5_file} → {zarr_store}")


def convert_zarr_to_hdf5(zarr_store, hdf5_file, compression_level=5):
    """
    Convert a Zarr store to HDF5 format
    
    Parameters:
    -----------
    zarr_store : str
        Path to the Zarr store
    hdf5_file : str
        Path for the output HDF5 file
    compression_level : int, optional
        Compression level (0-9) for HDF5 datasets
    """
    def _copy_attributes(obj, target):
        """Copy attributes from obj to target"""
        for key, value in obj.attrs.items():
            target.attrs[key] = value
    
    def _process_group(group, parent_group):
        """Process Zarr group recursively"""
        for key in group.keys():
            item = group[key]
            
            if isinstance(item, zarr.Group):
                # Create hdf5 group
                new_group = parent_group.create_group(key)
                _copy_attributes(item, new_group)
                _process_group(item, new_group)
            elif isinstance(item, zarr.Array):
                # Copy dataset
                shape = item.shape
                dtype = item.dtype
                chunks = item.chunks if item.chunks else None
                
                # Create hdf5 dataset
                hdf5_dataset = parent_group.create_dataset(
                    key,
                    shape=shape,
                    dtype=dtype,
                    chunks=chunks,
                    compression='gzip',
                    compression_opts=compression_level
                )
                
                # Copy data in chunks to handle large datasets
                if shape:  # Only chunk if dataset is not scalar
                    chunk_size = 10_000_000  # 10M elements at a time
                    for slices in _chunk_iterator(shape, chunk_size):
                        hdf5_dataset[slices] = item[slices]
                else:
                    # Handle scalar datasets
                    hdf5_dataset[()] = item[()]
                
                _copy_attributes(item, hdf5_dataset)
    
    def _chunk_iterator(shape, chunk_size):
        """Iterator that yields slice tuples for processing large datasets in chunks"""
        if not shape:
            yield ()
            return
            
        total_size = np.prod(shape)
        if total_size <= chunk_size:
            yield tuple(slice(0, s) for s in shape)
            return
        
        # For simplicity, we'll chunk along the first dimension
        chunk_size_dim0 = max(1, int(chunk_size / np.prod(shape[1:])))
        for i in range(0, shape[0], chunk_size_dim0):
            end = min(i + chunk_size_dim0, shape[0])
            slices = [slice(i, end)]
            slices.extend(slice(0, s) for s in shape[1:])
            yield tuple(slices)
    
    # Open zarr store
    root = zarr.open(zarr_store, mode='r')
    
    # Create HDF5 file
    with h5py.File(hdf5_file, 'w') as f:
        _copy_attributes(root, f)
        _process_group(root, f)
    
    print(f"Conversion from Zarr to HDF5 completed: {zarr_store} → {hdf5_file}")


def main():
    parser = argparse.ArgumentParser(description='Convert between HDF5 and Zarr formats')
    parser.add_argument('input', help='Input file or directory path')
    parser.add_argument('output', help='Output file or directory path')
    parser.add_argument('--compression', type=int, default=5, 
                        help='Compression level (0-9, default=5)')
    parser.add_argument('--format', choices=['auto', 'h5-to-zarr', 'zarr-to-h5'], 
                        default='auto', help='Conversion direction')

    args = parser.parse_args()
    
    # Auto-detect format if not specified
    if args.format == 'auto':
        # Check if input is HDF5 or Zarr
        if os.path.isfile(args.input) and args.input.lower().endswith(('.h5', '.hdf5', '.hdf')):
            args.format = 'h5-to-zarr'
        elif os.path.isdir(args.input) and any(os.path.exists(os.path.join(args.input, m)) 
                                              for m in ['.zgroup', '.zarray', '.zattrs']):
            args.format = 'zarr-to-h5'
        else:
            try:
                # Try to open as HDF5
                with h5py.File(args.input, 'r') as _:
                    args.format = 'h5-to-zarr'
            except:
                try:
                    # Try to open as Zarr
                    zarr.open(args.input, mode='r')
                    args.format = 'zarr-to-h5'
                except:
                    raise ValueError(f"Cannot auto-detect format for {args.input}. Please specify --format explicitly.")
    
    # Perform conversion
    if args.format == 'h5-to-zarr':
        convert_hdf5_to_zarr(args.input, args.output, args.compression)
    else:  # zarr-to-h5
        convert_zarr_to_hdf5(args.input, args.output, args.compression)


if __name__ == "__main__":
    main()
