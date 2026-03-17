#!/usr/bin/env python3
"""
Pseudotime Re-Analysis Script

Performs trajectory analysis with different root anchoring strategies:
1. Anchored on normal epithelial cells (overall trajectory)
2. Anchored on normal cells for each cell type separately

This script assumes you've already run the main analysis pipeline and have:
- Preprocessed AnnData object
- Cell type classifications
- Tissue type annotations (cancer vs normal)

Usage:
------
python pseudotime_reanalysis.py
"""

import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.sparse import issparse
import os
import warnings
warnings.filterwarnings('ignore')

# Set plotting parameters
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white', figsize=(8, 6))
plt.rcParams['figure.figsize'] = (10, 8)
sns.set_style("whitegrid")


# ============================================================================
# CONFIGURATION
# ============================================================================

# Input data paths
UMI_MATRIX_PATH = 'all_combined_UMIs.txt.gz'
METADATA_PATH = 'GSE160269_series_matrix.txt'
OUTPUT_DIR = './pseudotime_reanalysis_output'

# Genes to visualize along pseudotime
GENES_OF_INTEREST = [
    'NOTCH1', 'FAT1', 'TP53',
    'RB1', 'FBXW7', 'NFE2L2',
    'PIK3CA', 'NOTCH2'
]

# Trajectory parameters
SUBSAMPLE_FRACTION = 0.5  # Use 50% of cells for trajectory (memory efficiency)
RANDOM_SEED = 42

# Checkpointing
USE_CHECKPOINT = True  # Set to False to force reprocessing from scratch

# ============================================================================


class PseudotimeReanalyzer:
    """Re-analyze trajectories with different root anchoring strategies"""
    
    def __init__(self, umi_path, metadata_path, output_dir):
        """Initialize the reanalyzer"""
        self.umi_path = umi_path
        self.metadata_path = metadata_path
        self.output_dir = output_dir
        self.adata = None
        
        import os
        os.makedirs(output_dir, exist_ok=True)
        
    
    def load_and_preprocess(self, use_checkpoint=True):
        """Load data and run minimal preprocessing needed for trajectory
        
        Parameters:
        -----------
        use_checkpoint : bool
            If True, try to load from saved checkpoint instead of reprocessing
        """
        
        # Check if checkpoint exists
        checkpoint_path = f'{self.output_dir}/processed_data.h5ad'
        
        if use_checkpoint and os.path.exists(checkpoint_path):
            print(f"\n  ✓ Found checkpoint: {checkpoint_path}")
            print("  Loading from checkpoint (skip preprocessing)...")
            self.adata = sc.read_h5ad(checkpoint_path)
            
            print(f"\n  Loaded: {self.adata.n_obs} cells, {self.adata.n_vars} HVGs")
            print(f"  Cell types: {self.adata.obs['cell_type'].value_counts().to_dict()}")
            print(f"  Tissue types: {self.adata.obs['tissue_type'].value_counts().to_dict()}")
            print("\n  To reprocess from scratch, delete the checkpoint or set use_checkpoint=False")
            
            return self
        
        print("  No checkpoint found, processing from scratch...")
        
        # Load UMI matrix
        if self.umi_path.endswith('.gz'):
            import gzip
            with gzip.open(self.umi_path, 'rt') as f:
                umi_df = pd.read_csv(f, sep='\t', index_col=0)
        else:
            umi_df = pd.read_csv(self.umi_path, sep='\t', index_col=0)
        
        
        # Parse metadata
        sample_metadata = self._parse_metadata()
        
        # Create AnnData with sparse matrix
        from scipy.sparse import csr_matrix
        sparse_matrix = csr_matrix(umi_df.T.values)
        self.adata = sc.AnnData(X=sparse_matrix)
        self.adata.var_names = umi_df.index
        self.adata.obs_names = umi_df.columns
        
        # Map metadata to cells
        cell_metadata = self._map_metadata_to_cells(sample_metadata)
        self.adata.obs = cell_metadata
        
        del umi_df
        
        # Basic preprocessing
        sc.pp.filter_cells(self.adata, min_genes=200)
        sc.pp.filter_genes(self.adata, min_cells=3)
        
        # Calculate QC metrics
        self.adata.var['mt'] = self.adata.var_names.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(
            self.adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True
        )
        
        # Filter by mito
        self.adata = self.adata[self.adata.obs.pct_counts_mt < 20, :].copy()
        
        # Normalize and log
        sc.pp.normalize_total(self.adata, target_sum=1e4)
        sc.pp.log1p(self.adata)
        
        # Store raw
        self.adata.raw = self.adata
        
        # HVG and scale
        sc.pp.highly_variable_genes(self.adata, n_top_genes=2000)
        self.adata = self.adata[:, self.adata.var.highly_variable].copy()
        sc.pp.scale(self.adata, max_value=10)
        
        # PCA and neighbors
        sc.tl.pca(self.adata, svd_solver='arpack')
        sc.pp.neighbors(self.adata, n_neighbors=10, n_pcs=40)
        
        # UMAP
        sc.tl.umap(self.adata)
        
        # Clustering
        sc.tl.leiden(self.adata, resolution=0.5)
        
        # Cell type classification
        self._classify_cell_types()
        
        
        # Save processed data for checkpointing and export to R
        checkpoint_path = f'{self.output_dir}/processed_data.h5ad'
        print(f"\n  Saving processed data to: {checkpoint_path}")
        self.adata.write(checkpoint_path)
        
        return self
    
    def _parse_metadata(self):
        """Parse metadata from GEO series matrix"""
        
        metadata_dict = {}
        with open(self.metadata_path, 'r') as f:
            for line in f:
                if line.startswith('!Sample_title'):
                    samples = [s.strip('"') for s in line.strip().split('\t')[1:]]
                    metadata_dict['sample_id'] = samples
                elif line.startswith('!Sample_geo_accession'):
                    geo_ids = [g.strip('"') for g in line.strip().split('\t')[1:]]
                    metadata_dict['geo_accession'] = geo_ids
                elif 'tissue: ' in line and line.startswith('!Sample_characteristics'):
                    tissues = [t.strip('"').replace('tissue: ', '') for t in line.strip().split('\t')[1:]]
                    metadata_dict['tissue_annotation'] = tissues
                elif 'cell population: ' in line:
                    cd45 = [c.strip('"').replace('cell population: ', '') for c in line.strip().split('\t')[1:]]
                    metadata_dict['cd45_status'] = cd45
                elif 'patient: ' in line:
                    patients = [p.strip('"').replace('patient: ', '') for p in line.strip().split('\t')[1:]]
                    metadata_dict['patient'] = patients
        
        meta_df = pd.DataFrame(metadata_dict)
        
        # Determine tissue type from sample ID (T=Tumor, N=Normal)
        def determine_tissue_type(sample_id):
            if 'T-' in sample_id or sample_id.endswith('T'):
                return 'cancer'
            elif 'N-' in sample_id or sample_id.endswith('N'):
                return 'normal'
            else:
                return 'unknown'
        
        meta_df['tissue_type'] = meta_df['sample_id'].apply(determine_tissue_type)
        
        return meta_df
    
    def _map_metadata_to_cells(self, sample_metadata):
        """Map sample metadata to cells"""
        cell_meta = pd.DataFrame(index=self.adata.obs_names)
        cell_meta['sample_id'] = 'unknown'
        cell_meta['patient'] = 'unknown'
        cell_meta['tissue_type'] = 'unknown'
        cell_meta['tissue_annotation'] = 'unknown'
        cell_meta['cd45_status'] = 'unknown'
        
        import re
        cell_barcodes = self.adata.obs_names.tolist()
        
        for barcode in cell_barcodes:
            for _, sample_row in sample_metadata.iterrows():
                sample_id = sample_row['sample_id']
                
                # Try exact match
                if sample_id in barcode:
                    cell_meta.loc[barcode, 'sample_id'] = sample_id
                    cell_meta.loc[barcode, 'patient'] = sample_row['patient']
                    cell_meta.loc[barcode, 'tissue_type'] = sample_row['tissue_type']
                    cell_meta.loc[barcode, 'tissue_annotation'] = sample_row['tissue_annotation']
                    cell_meta.loc[barcode, 'cd45_status'] = sample_row['cd45_status']
                    break
                
                # Try P###T or P###N pattern
                pattern_match = re.match(r'(P\d+[TN])', sample_id)
                if pattern_match:
                    short_pattern = pattern_match.group(1)
                    if short_pattern in barcode:
                        cell_meta.loc[barcode, 'sample_id'] = sample_id
                        cell_meta.loc[barcode, 'patient'] = sample_row['patient']
                        cell_meta.loc[barcode, 'tissue_type'] = sample_row['tissue_type']
                        cell_meta.loc[barcode, 'tissue_annotation'] = sample_row['tissue_annotation']
                        cell_meta.loc[barcode, 'cd45_status'] = sample_row['cd45_status']
                        break
        
        return cell_meta
    
    def _classify_cell_types(self):
        """Simple cell type classification using marker genes"""
        marker_genes = {
            'T_cell': ['CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B'],
            'B_cell': ['CD19', 'CD79A', 'MS4A1'],
            'Myeloid': ['CD14', 'CD68', 'LYZ', 'CD163'],
            'NK_cell': ['NCAM1', 'NKG7', 'GNLY', 'KLRD1'],
            'Epithelial': ['EPCAM', 'KRT18', 'KRT19', 'CDH1'],
            'Fibroblast': ['COL1A1', 'COL1A2', 'DCN', 'LUM'],
            'Endothelial': ['PECAM1', 'VWF', 'CDH5'],
            'Malignant': ['MKI67', 'TOP2A', 'CCNA2']
        }
        
        for cell_type, markers in marker_genes.items():
            available_markers = [m for m in markers if m in self.adata.raw.var_names]
            if len(available_markers) > 0:
                sc.tl.score_genes(
                    self.adata, 
                    available_markers, 
                    score_name=f'{cell_type}_score',
                    use_raw=True
                )
        
        score_cols = [col for col in self.adata.obs.columns if col.endswith('_score')]
        if score_cols:
            self.adata.obs['cell_type'] = self.adata.obs[score_cols].idxmax(axis=1)
            self.adata.obs['cell_type'] = self.adata.obs['cell_type'].str.replace('_score', '')
    
    def trajectory_normal_epithelial_root(self, genes_of_interest):
        """
        Trajectory anchored on normal epithelial cells
        """

        
        # Subsample for memory efficiency
        np.random.seed(RANDOM_SEED)
        n_cells = self.adata.n_obs
        n_sample = int(n_cells * SUBSAMPLE_FRACTION)
        sample_indices = np.random.choice(n_cells, size=n_sample, replace=False)
        adata_subset = self.adata[sample_indices, :].copy()
                
        # Find normal epithelial cells
        normal_epi = adata_subset.obs[
            (adata_subset.obs['tissue_type'] == 'normal') & 
            (adata_subset.obs['cell_type'] == 'Epithelial')
        ]
        
        if len(normal_epi) == 0:
            print("No normal epithelial cells found!")
            print("Falling back to first epithelial cell...")
            epi_cells = adata_subset.obs[adata_subset.obs['cell_type'] == 'Epithelial']
            if len(epi_cells) > 0:
                root_idx = epi_cells.index[0]
                adata_subset.uns['iroot'] = np.where(adata_subset.obs_names == root_idx)[0][0]
            else:
                print("Bad labels y'all. Using first cell.")
                adata_subset.uns['iroot'] = 0
        else:
            root_idx = normal_epi.index[0]
            adata_subset.uns['iroot'] = np.where(adata_subset.obs_names == root_idx)[0][0]
        
        # Run trajectory
        sc.tl.paga(adata_subset, groups='leiden')
        
        sc.tl.diffmap(adata_subset)
        sc.tl.dpt(adata_subset)
        
        # Plot 1: Pseudotime on UMAP
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        
        sc.pl.umap(
            adata_subset,
            color='dpt_pseudotime',
            ax=axes[0],
            show=False,
            title='Pseudotime (Normal Epithelial Root)',
            cmap='viridis'
        )
        
        sc.pl.umap(
            adata_subset,
            color='cell_type',
            ax=axes[1],
            show=False,
            title='Cell Types'
        )
        
        sc.pl.umap(
            adata_subset,
            color='tissue_type',
            ax=axes[2],
            show=False,
            title='Tissue Type'
        )
        
        plt.tight_layout()
        plt.savefig(f'{self.output_dir}/01_pseudotime_normal_epithelial_root.png', 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
        # Plot 2: Genes along pseudotime
        available_genes = [g for g in genes_of_interest if g in adata_subset.raw.var_names]
        
        if available_genes:
            self._plot_genes_along_pseudotime(
                adata_subset, 
                available_genes,
                title_suffix='Normal Epithelial Root'
            )
            plt.savefig(f'{self.output_dir}/02_genes_pseudotime_normal_epithelial_root.png', 
                       dpi=300, bbox_inches='tight')
            plt.close()
                
        return self
    
    def trajectory_by_cell_type(self, genes_of_interest):
        """
        Trajectory for each cell type, anchored on normal cells
	Could have done the above with this, but alas, didn't think to
        """
        
        # Get unique cell types
        cell_types = self.adata.obs['cell_type'].unique()
        
        for cell_type in cell_types:
            
            # Subset to this cell type
            adata_celltype = self.adata[self.adata.obs['cell_type'] == cell_type].copy()
            
            if adata_celltype.n_obs < 50:
                print(f"    Cell paucity skipping this type...")
                continue
                        
            # Subsample for memory
            np.random.seed(RANDOM_SEED)
            n_sample = min(int(adata_celltype.n_obs * SUBSAMPLE_FRACTION), adata_celltype.n_obs)
            if n_sample < adata_celltype.n_obs:
                sample_indices = np.random.choice(adata_celltype.n_obs, size=n_sample, replace=False)
                adata_subset = adata_celltype[sample_indices, :].copy()
            else:
                adata_subset = adata_celltype.copy()
                        
            # Find normal cells of this type
            normal_cells = adata_subset.obs[adata_subset.obs['tissue_type'] == 'normal']
            
            if len(normal_cells) == 0:
                print(f"    No normal {cell_type} cells found, skipping...")
                continue
                        
            # Set root to first normal cell
            root_idx = normal_cells.index[0]
            adata_subset.uns['iroot'] = np.where(adata_subset.obs_names == root_idx)[0][0]
            
            # Recompute neighbors for this subset
            sc.pp.neighbors(adata_subset, n_neighbors=10, n_pcs=40)
            sc.tl.umap(adata_subset)
            sc.tl.leiden(adata_subset, resolution=0.5)
            #loads of repeat code I could remove from above function
            # Run trajectory
            sc.tl.paga(adata_subset, groups='leiden')
            
            sc.tl.diffmap(adata_subset)
            sc.tl.dpt(adata_subset)
            
            # Plot 1: Pseudotime on UMAP
            fig, axes = plt.subplots(1, 2, figsize=(12, 5))
            
            sc.pl.umap(
                adata_subset,
                color='dpt_pseudotime',
                ax=axes[0],
                show=False,
                title=f'Pseudotime: {cell_type}\n(Normal Root)',
                cmap='viridis'
            )
            
            sc.pl.umap(
                adata_subset,
                color='tissue_type',
                ax=axes[1],
                show=False,
                title='Tissue Type'
            )
            
            plt.tight_layout()
            safe_name = cell_type.replace('/', '_').replace(' ', '_')
            plt.savefig(
                f'{self.output_dir}/03_pseudotime_{safe_name}_normal_root.png', 
                dpi=300, bbox_inches='tight'
            )
            plt.close()
            
            # Plot 2: Genes along pseudotime
            available_genes = [g for g in genes_of_interest if g in adata_subset.raw.var_names]
            
            if available_genes:
                self._plot_genes_along_pseudotime(
                    adata_subset, 
                    available_genes,
                    title_suffix=f'{cell_type} (Normal Root)'
                )
                plt.savefig(
                    f'{self.output_dir}/04_genes_pseudotime_{safe_name}_normal_root.png', 
                    dpi=300, bbox_inches='tight'
                )
                plt.close()
                    
        return self
    
    def _plot_genes_along_pseudotime(self, adata_subset, genes, title_suffix=''):
        """Plot gene expression along pseudotime trajectory"""
        n_genes = len(genes)
        fig, axes = plt.subplots(n_genes, 1, figsize=(10, 4*n_genes))
        
        if n_genes == 1:
            axes = [axes]
        
        for idx, gene in enumerate(genes):
            # Extract gene expression
            if issparse(adata_subset.raw.X):
                gene_expr = adata_subset.raw[:, gene].X.toarray().flatten()
            else:
                gene_expr = adata_subset.raw[:, gene].X.flatten()
            
            pseudotime = adata_subset.obs['dpt_pseudotime'].values
            
            # Sort by pseudotime
            sort_idx = np.argsort(pseudotime)
            
            # Color by tissue type
            tissue_colors = adata_subset.obs['tissue_type'].astype(str).map({
                'cancer': 'red', 
                'normal': 'blue',
                'unknown': 'gray'
            })
            # Fill any unmapped values with gray
            colors = tissue_colors.fillna('gray').values
            
            # Scatter plot
            axes[idx].scatter(
                pseudotime,
                gene_expr,
                alpha=0.3,
                s=10,
                c=colors
            )
            
            # Add smoothed line
            from scipy.ndimage import gaussian_filter1d
            smooth_expr = gaussian_filter1d(gene_expr[sort_idx], sigma=20)
            axes[idx].plot(pseudotime[sort_idx], smooth_expr, 'k-', linewidth=2, alpha=0.7)
            
            axes[idx].set_xlabel('Pseudotime')
            axes[idx].set_ylabel('Expression')
            axes[idx].set_title(f'{gene} Expression Along Trajectory\n{title_suffix}')
            axes[idx].grid(True, alpha=0.3)
            
            # Add legend
            from matplotlib.patches import Patch
            legend_elements = [
                Patch(facecolor='red', alpha=0.6, label='Cancer'),
                Patch(facecolor='blue', alpha=0.6, label='Normal')
            ]
            axes[idx].legend(handles=legend_elements, loc='upper right')
        
        plt.tight_layout()

def main():
    """Main analysis function"""
    
    print("\n" + "=" * 70)
    print("STARTING PSEUDOTIME RE-ANALYSIS")
    print("=" * 70)
    
    # Initialize analyzer
    analyzer = PseudotimeReanalyzer(
        umi_path=UMI_MATRIX_PATH,
        metadata_path=METADATA_PATH,
        output_dir=OUTPUT_DIR
    )
    
    # Run analyses
    (analyzer
     .load_and_preprocess(use_checkpoint=USE_CHECKPOINT)
     .trajectory_normal_epithelial_root(GENES_OF_INTEREST)
     .trajectory_by_cell_type(GENES_OF_INTEREST)
    )
    
    print("\nPseudotime re-analyses complete!")


if __name__ == '__main__':
    main()
