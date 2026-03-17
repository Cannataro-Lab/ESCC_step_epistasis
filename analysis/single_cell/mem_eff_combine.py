import pandas as pd
import numpy as np
import gzip
from collections import defaultdict

#get all the gene names to help set up the data structure
def get_all_genes(cd45neg_path, cd45pos_path):
    """
    First pass: collect all unique genes from both files.
    """
    genes = set()
    
    # Read genes from CD45neg file
    opener = gzip.open if cd45neg_path.endswith('.gz') else open
    with opener(cd45neg_path, 'rt') as f:
        next(f)  # Skip header
        for line in f:
            gene = line.split()[0]
            genes.add(gene)
    
    # Read genes from CD45pos file
    opener = gzip.open if cd45pos_path.endswith('.gz') else open
    with opener(cd45pos_path, 'rt') as f:
        next(f)  # Skip header
        for line in f:
            gene = line.split()[0]
            genes.add(gene)
    
    # Sort genes alphabetically for consistent ordering
    sorted_genes = sorted(genes)
    print(f"  Found {len(sorted_genes)} unique genes")
    return sorted_genes

def get_cell_names(file_path):
    """Extract cell names from header line."""
    opener = gzip.open if file_path.endswith('.gz') else open
    with opener(file_path, 'rt') as f:
        header = f.readline().strip()
        cell_names = header.split()[1:]  # Skip first column (gene column name)
    return cell_names

def combine_umi_files_efficient(cd45neg_path, cd45pos_path, output_path, chunk_size=1000):
    """    
    Parameters:
    -----------
    cd45neg_path : str
        Path to CD45 negative UMI file
    cd45pos_path : str
        Path to CD45 positive UMI file
    output_path : str
        Path for combined output file
    chunk_size : int
        Number of genes to process at once (default: 1000)
    """
    
    # Step 1: Get all unique genes (memory efficient - just gene names)
    all_genes = get_all_genes(cd45neg_path, cd45pos_path)
    
    # Step 2: Get cell names from both files
    cd45neg_cells = get_cell_names(cd45neg_path)
    cd45pos_cells = get_cell_names(cd45pos_path)
    
    all_cells = cd45neg_cells + cd45pos_cells
    
    # Check for duplicate cell barcodes
    if len(all_cells) != len(set(all_cells)):
        print("WARNING: Duplicate cell barcodes detected!")
    
    # Step 3: Load both files into dictionaries for random access
    cd45neg_data = load_file_to_dict(cd45neg_path, cd45neg_cells)
    cd45pos_data = load_file_to_dict(cd45pos_path, cd45pos_cells)
    
    # Step 4: Write combined file gene by gene
    opener = gzip.open if output_path.endswith('.gz') else open
    
    with opener(output_path, 'wt') as out:
        # Write header
        out.write('\t'.join([''] + all_cells) + '\n')
        
        # Process genes in chunks for progress reporting
        total_umi = 0
        for i, gene in enumerate(all_genes):
            if (i + 1) % chunk_size == 0:
                print(f"  Processed {i + 1}/{len(all_genes)} genes...")
            
            # Get counts for this gene from both files
            neg_counts = cd45neg_data.get(gene, [0] * len(cd45neg_cells))
            pos_counts = cd45pos_data.get(gene, [0] * len(cd45pos_cells))
            
            # Combine counts
            combined_counts = neg_counts + pos_counts
            total_umi += sum(combined_counts)
            
            # Write row
            out.write(gene + '\t' + '\t'.join(map(str, combined_counts)) + '\n')
 

def load_file_to_dict(file_path, cell_names):
    """
    Load UMI file into a dictionary: gene -> list of counts.
    """
    data = {}
    opener = gzip.open if file_path.endswith('.gz') else open
    
    with opener(file_path, 'rt') as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split()
            gene = parts[0]
            counts = [int(x) for x in parts[1:]]
            
            # Verify we have the right number of values
            if len(counts) != len(cell_names):
                raise ValueError(f"Gene {gene}: expected {len(cell_names)} values, got {len(counts)}")
            
            data[gene] = counts
    
    return data



def read_specific_genes(file_path, genes_to_read, cell_names):
    """Read only specific genes from file. Helper function"""
    genes_set = set(genes_to_read)
    data = {}
    
    opener = gzip.open if file_path.endswith('.gz') else open
    with opener(file_path, 'rt') as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split()
            gene = parts[0]
            if gene in genes_set:
                counts = [int(x) for x in parts[1:]]
                data[gene] = counts
    
    if data:
        df = pd.DataFrame(data, index=cell_names).T
    else:
        df = pd.DataFrame(index=[], columns=cell_names)
    
    return df


if __name__ == "__main__":
    combine_umi_files_efficient(
        cd45neg_path='GSE160269_CD45neg_UMIs.txt.gz',
        cd45pos_path='GSE160269_CD45pos_UMIs.txt.gz',
        output_path='all_combined_UMIs.txt.gz',
        chunk_size=1000
    )