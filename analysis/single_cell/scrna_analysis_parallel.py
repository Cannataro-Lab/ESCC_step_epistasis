#!/usr/bin/env python3
"""
Parallelized Memory-Optimized scRNA-seq Analysis Pipeline

Key features:
1. Maintains sparse matrices throughout
2. Gene-level parallelization for correlations and DE testing
3. Cell-level parallelization for QC and scoring
4. Configurable thread count in one location
5. Subsampled trajectory analysis (50% of cells)
"""

import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.sparse import issparse, csr_matrix
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import partial
import warnings
import os
warnings.filterwarnings('ignore')

# Set plotting parameters
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white', figsize=(8, 6))
plt.rcParams['figure.figsize'] = (10, 8)
sns.set_style("whitegrid")


# Set this to match your CPU core count
N_THREADS = 12  # <-- MODIFY THIS VALUE



# ============================================================================


class scRNAAnalyzerParallel:
    """Parallelized memory-optimized single-cell RNA-seq analysis pipeline"""
    
    def __init__(self, umi_path, metadata_path, output_dir='./scrna_analysis_output'):
        """
        Initialize analyzer with data paths
        
        Parameters:
        -----------
        umi_path : str
            Path to combined UMI counts matrix (genes x cells)
        metadata_path : str
            Path to GEO series matrix file with sample annotations
        output_dir : str
            Directory for output files and figures
        """
        self.umi_path = umi_path
        self.metadata_path = metadata_path
        self.output_dir = output_dir
        self.adata = None
        self.genes_of_interest = []
        self.n_threads = N_THREADS
        
        import os
        os.makedirs(output_dir, exist_ok=True)
        
    
    def load_data(self):
        """Load UMI matrix and metadata, create AnnData object with sparse matrix"""
        
        # Load UMI matrix
        if self.umi_path.endswith('.gz'):
            import gzip
            with gzip.open(self.umi_path, 'rt') as f:
                umi_df = pd.read_csv(f, sep='\t', index_col=0)
        else:
            umi_df = pd.read_csv(self.umi_path, sep='\t', index_col=0)
        
        print(f"  Matrix shape: {umi_df.shape[0]} genes x {umi_df.shape[1]} cells")
        
        # Parse sample-level metadata
        sample_metadata = self._parse_metadata()
        
        # Create AnnData object with SPARSE matrix
        sparse_matrix = csr_matrix(umi_df.T.values)
        self.adata = sc.AnnData(X=sparse_matrix)
        self.adata.var_names = umi_df.index
        self.adata.obs_names = umi_df.columns
        self.adata.var_names_make_unique()
        
        # Map sample metadata to individual cells based on cell barcodes
        cell_metadata = self._map_metadata_to_cells(sample_metadata)
        self.adata.obs = cell_metadata
        
        # Basic stats
        total_counts = self.adata.X.sum()
        print(f"  Total UMI counts: {total_counts:,.0f}")
        print(f"  Mean UMIs per cell: {total_counts / self.adata.n_obs:.1f}")
        print(f"  Cancer samples: {(self.adata.obs['tissue_type'] == 'cancer').sum()}")
        print(f"  Normal samples: {(self.adata.obs['tissue_type'] == 'normal').sum()}")
        
        del umi_df
        
        return self
    
    def _parse_metadata(self):
        """Parse GEO series matrix metadata file"""
        print("  Parsing sample-level metadata...")
        
        metadata_dict = {}
        with open(self.metadata_path, 'r') as f:
            for line in f:
                if line.startswith('!Sample_title'):
                    samples = line.strip().split('\t')[1:]
                    samples = [s.strip('"') for s in samples]
                    metadata_dict['sample_id'] = samples
                    
                elif line.startswith('!Sample_geo_accession'):
                    geo_ids = line.strip().split('\t')[1:]
                    geo_ids = [g.strip('"') for g in geo_ids]
                    metadata_dict['geo_accession'] = geo_ids
                    
                elif 'tissue: ' in line and line.startswith('!Sample_characteristics'):
                    tissues = line.strip().split('\t')[1:]
                    tissues = [t.strip('"').replace('tissue: ', '') for t in tissues]
                    metadata_dict['tissue_annotation'] = tissues
                    
                elif 'cell population: ' in line:
                    cd45 = line.strip().split('\t')[1:]
                    cd45 = [c.strip('"').replace('cell population: ', '') for c in cd45]
                    metadata_dict['cd45_status'] = cd45
                    
                elif 'patient: ' in line:
                    patients = line.strip().split('\t')[1:]
                    patients = [p.strip('"').replace('patient: ', '') for p in patients]
                    metadata_dict['patient'] = patients
        
        # Create dataframe
        meta_df = pd.DataFrame(metadata_dict)
        
        # Add derived columns based on sample_id naming convention
        # T in name = Tumor, N in name = Normal
        # Examples: P1T-I (tumor), P128N-E (normal)
        def determine_tissue_type(sample_id):
            """Determine tissue type from sample ID naming convention"""
            # Check for T or N in the sample ID
            # Pattern: P[number]T-[I/E] = Tumor
            # Pattern: P[number]N-[I/E] = Normal
            if 'T-' in sample_id or sample_id.endswith('T'):
                return 'cancer'
            elif 'N-' in sample_id or sample_id.endswith('N'):
                return 'normal'
            # Fallback to annotation if available
            else:
                return 'unknown'
        
        meta_df['tissue_type'] = meta_df['sample_id'].apply(determine_tissue_type)
        
        # Count and report
        n_cancer = (meta_df['tissue_type'] == 'cancer').sum()
        n_normal = (meta_df['tissue_type'] == 'normal').sum()
        n_unknown = (meta_df['tissue_type'] == 'unknown').sum()
        
        print(f"  Sample breakdown:")
        print(f"    Tumor samples: {n_cancer}")
        print(f"    Normal samples: {n_normal}")
        if n_unknown > 0:
            print(f"    Unknown samples: {n_unknown}")
        
        return meta_df
    
    def _map_metadata_to_cells(self, sample_metadata):
        """
        Map sample-level metadata to individual cell barcodes
        
        The combined UMI file concatenates CD45neg then CD45pos cells.
        Cell barcodes may contain sample/patient identifiers.
        
        Key pattern: T = Tumor, N = Normal
        Examples: P1T-I (tumor), P128N-E (normal)
        """
        n_cells = self.adata.n_obs
        
        # Create cell-level metadata dataframe
        cell_meta = pd.DataFrame(index=self.adata.obs_names)
        
        # Initialize columns with defaults
        cell_meta['sample_id'] = 'unknown'
        cell_meta['patient'] = 'unknown'
        cell_meta['tissue_type'] = 'unknown'
        cell_meta['tissue_annotation'] = 'unknown'
        cell_meta['cd45_status'] = 'unknown'
        cell_meta['geo_accession'] = 'unknown'
        
        cell_barcodes = self.adata.obs_names.tolist()
        
        # Strategy: Try to match cell barcodes to sample IDs
        matched_count = 0
                
        import re
        
        for idx, barcode in enumerate(cell_barcodes):
            matched = False
            
            # Try to match against each sample
            # IMPORTANT: Sort samples to try normal samples FIRST
            # This prevents tumor samples from matching too broadly
            samples_sorted = sample_metadata.sort_values('tissue_type', ascending=False)  # 'normal' before 'cancer'
            
            for _, sample_row in samples_sorted.iterrows():
                sample_id = sample_row['sample_id']
                patient_id = sample_row['patient']
                tissue_type = sample_row['tissue_type']
                
                # Strategy 1: Try exact sample ID match first (most specific)
                if sample_id in barcode:
                    cell_meta.loc[barcode, 'sample_id'] = sample_id
                    cell_meta.loc[barcode, 'patient'] = patient_id
                    cell_meta.loc[barcode, 'tissue_type'] = tissue_type
                    cell_meta.loc[barcode, 'tissue_annotation'] = sample_row['tissue_annotation']
                    cell_meta.loc[barcode, 'cd45_status'] = sample_row['cd45_status']
                    cell_meta.loc[barcode, 'geo_accession'] = sample_row['geo_accession']
                    matched_count += 1
                    matched = True
                    break
                
                # Strategy 2: Try pattern matching P###T or P###N
                # Extract pattern like "P1T", "P128N" from sample_id
                pattern_match = re.match(r'(P\d+[TN])', sample_id)
                if pattern_match:
                    short_pattern = pattern_match.group(1)
                    if short_pattern in barcode:
                        cell_meta.loc[barcode, 'sample_id'] = sample_id
                        cell_meta.loc[barcode, 'patient'] = patient_id
                        cell_meta.loc[barcode, 'tissue_type'] = tissue_type
                        cell_meta.loc[barcode, 'tissue_annotation'] = sample_row['tissue_annotation']
                        cell_meta.loc[barcode, 'cd45_status'] = sample_row['cd45_status']
                        cell_meta.loc[barcode, 'geo_accession'] = sample_row['geo_accession']
                        matched_count += 1
                        matched = True
                        break
            
            # If still not matched, try less specific patient matching
            # But ONLY if no other match was found
            if not matched:
                for _, sample_row in samples_sorted.iterrows():
                    patient_id = sample_row['patient']
                    
                    # Only match patient ID if it won't create ambiguity
                    # Check if this patient has both T and N samples
                    patient_samples = sample_metadata[sample_metadata['patient'] == patient_id]
                    
                    # If patient has only one tissue type, safe to match by patient alone
                    if len(patient_samples['tissue_type'].unique()) == 1:
                        if patient_id in barcode:
                            cell_meta.loc[barcode, 'sample_id'] = sample_row['sample_id']
                            cell_meta.loc[barcode, 'patient'] = patient_id
                            cell_meta.loc[barcode, 'tissue_type'] = sample_row['tissue_type']
                            cell_meta.loc[barcode, 'tissue_annotation'] = sample_row['tissue_annotation']
                            cell_meta.loc[barcode, 'cd45_status'] = sample_row['cd45_status']
                            cell_meta.loc[barcode, 'geo_accession'] = sample_row['geo_accession']
                            matched_count += 1
                            matched = True
                            break
            
            # Progress indicator
            if (idx + 1) % 50000 == 0:
                print(f"    Processed {idx + 1}/{n_cells} cells...")
        
        match_rate = matched_count / n_cells * 100
        print(f"\n  Matched {matched_count}/{n_cells} cells ({match_rate:.1f}%)")
        
        # Report tissue type distribution
        tissue_counts = cell_meta['tissue_type'].value_counts()
        print(f"\n  Cell distribution by tissue type:")
        for tissue_type, count in tissue_counts.items():
            print(f"    {tissue_type}: {count:,} cells ({count/n_cells*100:.1f}%)")
        
        # Show sample breakdown for matched cells
        if matched_count > 0:
            print(f"\n  Sample breakdown (cells per sample):")
            sample_counts = cell_meta[cell_meta['sample_id'] != 'unknown']['sample_id'].value_counts()
            
            # Group by tissue type
            normal_samples = []
            cancer_samples = []
            
            for sample_id, count in sample_counts.items():
                tissue = cell_meta[cell_meta['sample_id'] == sample_id]['tissue_type'].iloc[0]
                if tissue == 'normal':
                    normal_samples.append((sample_id, count))
                else:
                    cancer_samples.append((sample_id, count))
            
            if normal_samples:
                print(f"\n    Normal samples ({len(normal_samples)}):")
                for sample_id, count in sorted(normal_samples, key=lambda x: x[1], reverse=True)[:5]:
                    print(f"      {sample_id}: {count:,} cells")
            
            if cancer_samples:
                print(f"\n    Cancer samples ({len(cancer_samples)}) - showing top 5:")
                for sample_id, count in sorted(cancer_samples, key=lambda x: x[1], reverse=True)[:5]:
                    print(f"      {sample_id}: {count:,} cells")
        
        return cell_meta
    
    def preprocess(self, min_genes=200, min_cells=3, max_pct_mito=20):
        """
        Preprocess data with CELL-LEVEL PARALLELIZATION for QC
        
        Parameters:
        -----------
        min_genes : int
            Minimum number of genes per cell
        min_cells : int
            Minimum number of cells expressing a gene
        max_pct_mito : float
            Maximum percentage of mitochondrial genes
        """        
        # Calculate QC metrics in parallel by cell chunks
        self.adata.var['mt'] = self.adata.var_names.str.startswith('MT-')
        
        # Split cells into chunks for parallel processing
        n_cells = self.adata.n_obs
        chunk_size = max(1, n_cells // self.n_threads)
        cell_chunks = [
            (i, min(i + chunk_size, n_cells)) 
            for i in range(0, n_cells, chunk_size)
        ]
        
        def compute_qc_chunk(start_idx, end_idx, adata):
            """Compute QC metrics for a chunk of cells"""
            chunk = adata[start_idx:end_idx, :]
            
            # Calculate metrics for this chunk
            n_genes = (chunk.X > 0).sum(axis=1).A1
            total_counts = chunk.X.sum(axis=1).A1
            
            # Mitochondrial content
            mt_genes = chunk.var['mt'].values
            if mt_genes.sum() > 0:
                mt_counts = chunk.X[:, mt_genes].sum(axis=1).A1
                pct_mt = (mt_counts / total_counts) * 100
            else:
                pct_mt = np.zeros(len(n_genes))
            
            return {
                'n_genes_by_counts': n_genes,
                'total_counts': total_counts,
                'pct_counts_mt': pct_mt
            }
        
        # Process chunks in parallel
        with ThreadPoolExecutor(max_workers=self.n_threads) as executor:
            futures = [
                executor.submit(compute_qc_chunk, start, end, self.adata)
                for start, end in cell_chunks
            ]
            
            results = [f.result() for f in futures]
        
        # Combine results
        self.adata.obs['n_genes_by_counts'] = np.concatenate([r['n_genes_by_counts'] for r in results])
        self.adata.obs['total_counts'] = np.concatenate([r['total_counts'] for r in results])
        self.adata.obs['pct_counts_mt'] = np.concatenate([r['pct_counts_mt'] for r in results])
        
        # Filter cells and genes
        sc.pp.filter_cells(self.adata, min_genes=min_genes)
        sc.pp.filter_genes(self.adata, min_cells=min_cells)
        
        self.adata = self.adata[self.adata.obs.pct_counts_mt < max_pct_mito, :].copy()
        
        # Normalize (keeps sparse)
        sc.pp.normalize_total(self.adata, target_sum=1e4)
        sc.pp.log1p(self.adata)
        
        # Store raw counts
        self.adata.raw = self.adata
        
        # Identify highly variable genes
        sc.pp.highly_variable_genes(self.adata, n_top_genes=2000)
        
        # Scale data
        self.adata = self.adata[:, self.adata.var.highly_variable].copy()
        sc.pp.scale(self.adata, max_value=10)
        
        return self
    
    def classify_cell_types(self, use_markers=True):
        """
        Classify cell types for marker scoring
        
        Parameters:
        -----------
        use_markers : bool
            Whether to use canonical marker genes for annotation
        """        
        # PCA
        sc.tl.pca(self.adata, svd_solver='arpack')
        
        # Neighborhood graph
        sc.pp.neighbors(self.adata, n_neighbors=10, n_pcs=40)
        
        # UMAP
        sc.tl.umap(self.adata)
        
        # Leiden clustering
        sc.tl.leiden(self.adata, resolution=0.5)
        
        if use_markers:
            marker_genes = {#these are cannonical and consistent with orig paper
                'T_cell': ['CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B'],
                'B_cell': ['CD19', 'CD79A', 'MS4A1'],
                'Myeloid': ['CD14', 'CD68', 'LYZ', 'CD163'],
                'NK_cell': ['NCAM1', 'NKG7', 'GNLY', 'KLRD1'],
                'Epithelial': ['EPCAM', 'KRT18', 'KRT19', 'CDH1'],
                'Fibroblast': ['COL1A1', 'COL1A2', 'DCN', 'LUM'],
                'Endothelial': ['PECAM1', 'VWF', 'CDH5'],
                'Malignant': ['MKI67', 'TOP2A', 'CCNA2']
            }
                        
            def score_cell_type_parallel(cell_type, markers, adata):
                """Score a cell type using marker genes"""
                available_markers = [m for m in markers if m in adata.raw.var_names]
                if len(available_markers) > 0:
                    # Extract marker expression
                    marker_expr = adata.raw[:, available_markers].X
                    if issparse(marker_expr):
                        marker_expr = marker_expr.toarray()
                    
                    # Calculate mean expression score
                    score = marker_expr.mean(axis=1)
                    return cell_type, score, len(available_markers), len(markers)
                return cell_type, None, 0, len(markers)
            
            # Score all cell types in parallel, stolen shamelessly
            with ThreadPoolExecutor(max_workers=self.n_threads) as executor:
                futures = [
                    executor.submit(score_cell_type_parallel, ct, markers, self.adata)
                    for ct, markers in marker_genes.items()
                ]
                
                results = [f.result() for f in futures]
            
            # Add scores to adata
            for cell_type, score, n_avail, n_total in results:
                if score is not None:
                    self.adata.obs[f'{cell_type}_score'] = score
            
            # Assign cell types based on highest score
            score_cols = [col for col in self.adata.obs.columns if col.endswith('_score')]
            if score_cols:
                self.adata.obs['cell_type'] = self.adata.obs[score_cols].idxmax(axis=1)
                self.adata.obs['cell_type'] = self.adata.obs['cell_type'].str.replace('_score', '')
                print(f"\n  Cell type distribution:")
                print(self.adata.obs['cell_type'].value_counts())
        
        self.adata.obs['cluster'] = self.adata.obs['leiden']
        
        # Plot
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        
        sc.pl.umap(self.adata, color='cluster', ax=axes[0], show=False, title='Clusters')
        if 'cell_type' in self.adata.obs.columns:
            sc.pl.umap(self.adata, color='cell_type', ax=axes[1], show=False, title='Cell Types')
        sc.pl.umap(self.adata, color='tissue_type', ax=axes[2], show=False, title='Tissue Type')
        
        plt.tight_layout()
        plt.savefig(f'{self.output_dir}/01_cell_type_classification.png', dpi=300, bbox_inches='tight')
        plt.close()
                
        return self
    
    def differential_expression_cancer_vs_normal(self, genes_of_interest=None):
        """
         differential expression between cancer and normal
        Uses gene-level parallelization for statistical testing
        
        Parameters:
        -----------
        genes_of_interest : list
            List of genes to specifically analyze
        """
        
        if genes_of_interest:
            self.genes_of_interest = genes_of_interest
        
        # Get cancer vs normal cell indices
        cancer_mask = self.adata.obs['tissue_type'] == 'cancer'
        normal_mask = self.adata.obs['tissue_type'] == 'normal'
        
        cancer_indices = np.where(cancer_mask)[0]
        normal_indices = np.where(normal_mask)[0]
        
        
        def test_gene_de(gene_idx, gene_name, adata_raw, cancer_idx, normal_idx):
            """Test differential expression for a single gene; helper function"""
            # Extract gene expression
            if issparse(adata_raw.X):
                gene_expr = adata_raw.X[:, gene_idx].toarray().flatten()
            else:
                gene_expr = adata_raw.X[:, gene_idx].flatten()
            
            cancer_expr = gene_expr[cancer_idx]
            normal_expr = gene_expr[normal_idx]
            
            # Statistical test
            try:
                statistic, pvalue = stats.mannwhitneyu(
                    cancer_expr, normal_expr, alternative='two-sided'
                )
                
                # Log fold change
                cancer_mean = np.log1p(cancer_expr.mean())
                normal_mean = np.log1p(normal_expr.mean())
                logfc = cancer_mean - normal_mean
                
                return {
                    'gene': gene_name,
                    'logfoldchange': logfc,
                    'pvalue': pvalue,
                    'statistic': statistic
                }
            except Exception as e:
                return {
                    'gene': gene_name,
                    'logfoldchange': 0,
                    'pvalue': 1.0,
                    'statistic': 0
                }
        
        # Test genes in parallel
        gene_indices = list(range(self.adata.raw.n_vars))
        gene_names = list(self.adata.raw.var_names)
        
        with ThreadPoolExecutor(max_workers=self.n_threads) as executor:
            futures = [
                executor.submit(
                    test_gene_de, idx, name, self.adata.raw, 
                    cancer_indices, normal_indices
                )
                for idx, name in zip(gene_indices, gene_names)
            ]
            
            results = []
            for i, future in enumerate(as_completed(futures), 1):
                results.append(future.result())
        
        # Create results dataframe
        de_df = pd.DataFrame(results)
        
        # Multiple testing correction; shamelessly stolen
        from statsmodels.stats.multitest import multipletests
        _, pvals_adj, _, _ = multipletests(de_df['pvalue'], method='fdr_bh')
        de_df['pvalue_adj'] = pvals_adj
        
        # Sort by p-value; don't end me townsend
        de_df = de_df.sort_values('pvalue')
        
        # Save results
        de_df.to_csv(f'{self.output_dir}/DE_cancer_vs_normal_parallel.csv', index=False)
        
        # Volcano plot
        self._plot_volcano(
            de_df, 
            title='Cancer vs Normal',
            genes_to_label=genes_of_interest
        )
        plt.savefig(f'{self.output_dir}/02_volcano_cancer_vs_normal.png', 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
        # Store for later
        self._de_cancer_normal_genes = genes_of_interest
        self._de_results = de_df
        
        return self
    
    def differential_expression_by_cell_type(self, genes_of_interest=None):
        """
        Differential expression between cell types        
        Parameters:
        -----------
        genes_of_interest : list
            List of genes to specifically analyze
        """
        
        if not 'cell_type' in self.adata.obs.columns:
            print("  Warning: Cell types not classified. Skipping.")
            return self
        
        # Scanpy's rank_genes_groups is already optimized and uses parallelization
        sc.tl.rank_genes_groups(
            self.adata,
            groupby='cell_type',
            method='wilcoxon',
            use_raw=True
        )
        
        # Extract and save results
        result = self.adata.uns['rank_genes_groups']
        cell_types = result['names'].dtype.names
        
        for ct in cell_types:
            de_df = pd.DataFrame({
                'gene': result['names'][ct],
                'logfoldchange': result['logfoldchanges'][ct],
                'pvalue': result['pvals'][ct],
                'pvalue_adj': result['pvals_adj'][ct],
            })
            de_df.to_csv(f'{self.output_dir}/DE_celltype_{ct}.csv', index=False)
        
        # Heatmap
        sc.pl.rank_genes_groups_heatmap(
            self.adata, 
            n_genes=10, 
            groupby='cell_type',
            show=False,
            cmap='RdBu_r'
        )
        plt.savefig(f'{self.output_dir}/04_celltype_markers_heatmap.png', 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
        self._de_celltype_genes = genes_of_interest
        
        return self
    
    def visualize_genes_low_dimensional(self, genes_of_interest):
        """
        Visualize gene expression in low-dimensional space (UMAP)
        
        Parameters:
        -----------
        genes_of_interest : list
            List of genes to visualize
        """
        
        available_genes = [g for g in genes_of_interest if g in self.adata.raw.var_names]
        missing_genes = [g for g in genes_of_interest if g not in self.adata.raw.var_names]
        
        if missing_genes:
            print(f"  Warning: {len(missing_genes)} genes not found: {missing_genes[:5]}")
        
        if not available_genes:
            print("  No genes available for visualization")
            return self        
        n_genes = len(available_genes)
        n_cols = min(4, n_genes)
        n_rows = (n_genes + n_cols - 1) // n_cols
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 5*n_rows))
        if n_genes == 1:
            axes = [axes]
        else:
            axes = axes.flatten()
        
        for idx, gene in enumerate(available_genes):
            sc.pl.umap(
                self.adata,
                color=gene,
                use_raw=True,
                ax=axes[idx],
                show=False,
                title=gene,
                cmap='viridis'
            )
        
        for idx in range(len(available_genes), len(axes)):
            axes[idx].axis('off')
        
        plt.tight_layout()
        plt.savefig(f'{self.output_dir}/06_genes_umap.png', dpi=300, bbox_inches='tight')
        plt.close()
                
        return self
    
    def analyze_gene_coexpression(self, genes_of_interest):
        """
        gene co-expression analysis
        Uses gene-pair parallelization
        
        Parameters:
        -----------
        genes_of_interest : list
            List of genes to analyze
        """
        
        available_genes = [g for g in genes_of_interest if g in self.adata.raw.var_names]
        
        if len(available_genes) < 2:
            print("  Need at least 2 genes for co-expression analysis")
            return self
        
        
        n_genes = len(available_genes)
        total_pairs = (n_genes * (n_genes - 1)) // 2 + n_genes  # Including diagonal
        
        def compute_correlation(i, j, gene_i, gene_j, adata_raw):
            """Compute correlation for a gene pair"""
            if issparse(adata_raw.X):
                expr_i = adata_raw[:, gene_i].X.toarray().flatten()
                expr_j = adata_raw[:, gene_j].X.toarray().flatten()
            else:
                expr_i = adata_raw[:, gene_i].X.flatten()
                expr_j = adata_raw[:, gene_j].X.flatten()
            
            corr = np.corrcoef(expr_i, expr_j)[0, 1]
            return i, j, corr
        
        # Generate all pairs
        pairs = [(i, j) for i in range(n_genes) for j in range(i, n_genes)]
                
        # Compute correlations in parallel
        corr_matrix = np.zeros((n_genes, n_genes))
        
        with ThreadPoolExecutor(max_workers=self.n_threads) as executor:
            futures = [
                executor.submit(
                    compute_correlation, i, j, 
                    available_genes[i], available_genes[j], 
                    self.adata.raw
                )
                for i, j in pairs
            ]
            
            completed = 0
            for future in as_completed(futures):
                i, j, corr = future.result()
                corr_matrix[i, j] = corr
                corr_matrix[j, i] = corr
                
                completed += 1
                if completed % 10 == 0 or completed == len(pairs):
                    print(f"    Progress: {completed}/{len(pairs)} pairs")
        
        # Convert to DataFrame
        corr_df = pd.DataFrame(
            corr_matrix,
            columns=available_genes,
            index=available_genes
        )
        
        # Plot correlation heatmap
        plt.figure(figsize=(10, 8))
        sns.heatmap(
            corr_df,
            annot=True,
            cmap='RdBu_r',
            center=0,
            vmin=-1,
            vmax=1,
            square=True,
            linewidths=1
        )
        plt.title('Gene Co-expression Correlation Matrix')
        plt.tight_layout()
        plt.savefig(f'{self.output_dir}/07_gene_correlation.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Save correlation matrix
        corr_df.to_csv(f'{self.output_dir}/gene_correlation_matrix.csv')
        
        # Scatter plots for highly correlated pairs
        high_corr_pairs = []
        for i in range(len(available_genes)):
            for j in range(i+1, len(available_genes)):
                corr_val = corr_matrix[i, j]
                if abs(corr_val) > 0.5:
                    high_corr_pairs.append((available_genes[i], available_genes[j], corr_val))
        
        if high_corr_pairs:            
            n_pairs = min(6, len(high_corr_pairs))
            fig, axes = plt.subplots(2, 3, figsize=(15, 10))
            axes = axes.flatten()
            
            for idx, (gene1, gene2, corr) in enumerate(high_corr_pairs[:n_pairs]):
                if issparse(self.adata.raw.X):
                    expr1 = self.adata.raw[:, gene1].X.toarray().flatten()
                    expr2 = self.adata.raw[:, gene2].X.toarray().flatten()
                else:
                    expr1 = self.adata.raw[:, gene1].X.flatten()
                    expr2 = self.adata.raw[:, gene2].X.flatten()
                
                axes[idx].scatter(expr1, expr2, alpha=0.5, s=10)
                axes[idx].set_xlabel(gene1)
                axes[idx].set_ylabel(gene2)
                axes[idx].set_title(f'r = {corr:.3f}')
                
                z = np.polyfit(expr1, expr2, 1)
                p = np.poly1d(z)
                x_sorted = np.sort(expr1)
                axes[idx].plot(x_sorted, p(x_sorted), "r--", alpha=0.8)
                
                del expr1, expr2
            
            for idx in range(n_pairs, 6):
                axes[idx].axis('off')
            
            plt.tight_layout()
            plt.savefig(f'{self.output_dir}/08_gene_coexpression_scatter.png', 
                       dpi=300, bbox_inches='tight')
            plt.close()
                
        return self
    
    def trajectory_analysis(self, root_cell_type=None):
        """
        Trajectory inference using random 50% of cells
        
        Parameters:
        -----------
        root_cell_type : str
            Cell type to use as trajectory root (optional)
        """
        
        try:
            n_cells = self.adata.n_obs
            n_sample = n_cells // 2
            
            np.random.seed(42)
            sample_indices = np.random.choice(n_cells, size=n_sample, replace=False)
            adata_subset = self.adata[sample_indices, :].copy()
            
            sc.tl.paga(adata_subset, groups='leiden')
            
            sc.pl.paga(
                adata_subset,
                threshold=0.03,
                show=False,
                title='PAGA Trajectory (50% cells)'
            )
            plt.savefig(f'{self.output_dir}/09_paga_trajectory.png', 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
            
            if root_cell_type and 'cell_type' in adata_subset.obs.columns:
                if root_cell_type in adata_subset.obs['cell_type'].values:
                    root_cells = adata_subset.obs['cell_type'] == root_cell_type
                    root_idx = np.where(root_cells)[0][0]
                    adata_subset.uns['iroot'] = root_idx
                    print(f"  Using {root_cell_type} as root")
            else:
                adata_subset.uns['iroot'] = 0
            
            sc.tl.diffmap(adata_subset)
            sc.tl.dpt(adata_subset)
            
            fig, axes = plt.subplots(1, 3, figsize=(18, 5))
            
            sc.pl.umap(
                adata_subset,
                color='dpt_pseudotime',
                ax=axes[0],
                show=False,
                title='Pseudotime',
                cmap='viridis'
            )
            
            if 'cell_type' in adata_subset.obs.columns:
                sc.pl.umap(adata_subset, color='cell_type', ax=axes[1], show=False)
            
            sc.pl.umap(adata_subset, color='tissue_type', ax=axes[2], show=False)
            
            plt.tight_layout()
            plt.savefig(f'{self.output_dir}/10_pseudotime.png', 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
            if self.genes_of_interest:
                available_genes = [g for g in self.genes_of_interest 
                                 if g in adata_subset.raw.var_names]
                
                if available_genes:
                    self._plot_genes_along_pseudotime_sparse(adata_subset, available_genes)
                    plt.savefig(f'{self.output_dir}/11_genes_pseudotime.png', 
                               dpi=300, bbox_inches='tight')
                    plt.close()
            
            del adata_subset
            
        except Exception as e:
            print(f"  Warning: Trajectory failed: {str(e)}")
        
        return self
    
    def perform_dense_operations(self):
        """
        Perform operations requiring dense arrays (e.g., slow and mem intense)
        Could be parallelized by gene, but memory cost usually not worth it
        """
        
        if hasattr(self, '_de_cancer_normal_genes') and self._de_cancer_normal_genes:
            available_genes = [g for g in self._de_cancer_normal_genes 
                             if g in self.adata.raw.var_names]
            if available_genes:
                self._plot_gene_expression_by_condition(
                    available_genes, 
                    groupby='tissue_type',
                    title='Expression by Tissue Type'
                )
                plt.savefig(f'{self.output_dir}/03_genes_cancer_vs_normal.png', 
                           dpi=300, bbox_inches='tight')
                plt.close()
        
        if hasattr(self, '_de_celltype_genes') and self._de_celltype_genes:
            available_genes = [g for g in self._de_celltype_genes 
                             if g in self.adata.raw.var_names]
            if available_genes and 'cell_type' in self.adata.obs.columns:
                self._plot_gene_expression_by_condition(
                    available_genes,
                    groupby='cell_type',
                    title='Expression by Cell Type'
                )
                plt.savefig(f'{self.output_dir}/05_genes_by_celltype.png', 
                           dpi=300, bbox_inches='tight')
                plt.close()        
        return self
    
    def plot_per_sample_cnv_summary(self, cnv_high_quantile=0.70):
        """
        Create per-sample summary of expression heterogeneity among epithelial cells
        and plot median heterogeneity by group.

        Groups:
          - Normal_Epithelial
          - Cancer_Epithelial
          - Cancer_Malignant_CNVhigh

        Note:
        This uses per-cell expression heterogeneity (SD across log-normalized genes
        in adata.raw) as a lightweight proxy for the R cnv_score-like quantity.
        """
        if 'cell_type' not in self.adata.obs.columns:
            print("  Warning: cell_type not found. Run classify_cell_types() first.")
            return self

        # Keep epithelial cells only
        epi_mask = self.adata.obs['cell_type'] == 'Epithelial'
        if epi_mask.sum() == 0:
            print("  Warning: no epithelial cells found.")
            return self

        adata_epi = self.adata[epi_mask].copy()

        # Compute per-cell expression heterogeneity from raw log-normalized data
        X = adata_epi.raw.X
        if issparse(X):
            cell_mean = X.mean(axis=1).A1
            cell_mean_sq = X.power(2).mean(axis=1).A1
            cell_sd = np.sqrt(np.maximum(cell_mean_sq - cell_mean**2, 0))
        else:
            cell_sd = X.std(axis=1)

        adata_epi.obs['cnv_score'] = cell_sd

        # Base groups from tissue type
        adata_epi.obs['expr_group'] = 'Other'
        adata_epi.obs.loc[
            adata_epi.obs['tissue_type'] == 'normal', 'expr_group'
        ] = 'Normal_Epithelial'
        adata_epi.obs.loc[
            adata_epi.obs['tissue_type'] == 'cancer', 'expr_group'
        ] = 'Cancer_Epithelial'

        # Define CNV-high malignant-like cells within each tumor sample
        tumor_mask = adata_epi.obs['tissue_type'] == 'cancer'
        tumor_obs = adata_epi.obs.loc[tumor_mask, ['sample_id', 'cnv_score']].copy()

        sample_thresh = (
            tumor_obs.groupby('sample_id')['cnv_score']
            .quantile(cnv_high_quantile)
            .rename('cnv_thresh')
        )

        adata_epi.obs = adata_epi.obs.join(sample_thresh, on='sample_id')

        malignant_mask = (
            (adata_epi.obs['tissue_type'] == 'cancer') &
            (adata_epi.obs['cnv_score'] >= adata_epi.obs['cnv_thresh'])
        )
        adata_epi.obs.loc[malignant_mask, 'expr_group'] = 'Cancer_Malignant_CNVhigh'

        # Per-sample summary
        ps = (
            adata_epi.obs.loc[
                adata_epi.obs['expr_group'].isin([
                    'Normal_Epithelial',
                    'Cancer_Epithelial',
                    'Cancer_Malignant_CNVhigh'
                ])
            ]
            .groupby(['sample_id', 'patient', 'tissue_type', 'expr_group'], dropna=False)
            .agg(
                n_cells=('cnv_score', 'size'),
                cnv_mean=('cnv_score', 'mean'),
                cnv_median=('cnv_score', 'median')
            )
            .reset_index()
        )

        # Save summary table
        summary_dir = os.path.join(self.output_dir, 'outputs_sc_valid')
        plot_dir = os.path.join(summary_dir, 'plots')
        os.makedirs(summary_dir, exist_ok=True)
        os.makedirs(plot_dir, exist_ok=True)

        ps.to_csv(os.path.join(summary_dir, 'per_sample_summary.csv'), index=False)

        # Plot
        order = [
            'Normal_Epithelial',
            'Cancer_Epithelial',
            'Cancer_Malignant_CNVhigh'
        ]
        labels = [
            'Normal epithelial',
            'ESCC epithelial',
            'ESCC malignant'
        ]

        ps = ps[ps['expr_group'].isin(order)].copy()
        ps['expr_group'] = pd.Categorical(ps['expr_group'], categories=order, ordered=True)

        fig, ax = plt.subplots(figsize=(7, 5))

        grouped = [
            ps.loc[ps['expr_group'] == g, 'cnv_median'].dropna().values
            for g in order
        ]

        ax.boxplot(grouped, positions=np.arange(1, len(order) + 1),
                   widths=0.55, showfliers=False)

        rng = np.random.default_rng(42)
        for i, g in enumerate(order, start=1):
            y = ps.loc[ps['expr_group'] == g, 'cnv_median'].dropna().values
            x = rng.normal(i, 0.06, size=len(y))
            ax.scatter(x, y, alpha=0.6, s=18)

        ax.set_xticks(np.arange(1, len(order) + 1))
        ax.set_xticklabels(labels)
        ax.set_xlabel("Cell classification")
        ax.set_ylabel("Median expression heterogeneity")

        # theme_classic-ish
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, 'CNVscore_by_group_per_sample.png'),
                    dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(plot_dir, 'CNVscore_by_group_per_sample.pdf'),
                    bbox_inches='tight')
        plt.close()

        print("\n  Saved:")
        print(f"    {os.path.join(summary_dir, 'per_sample_summary.csv')}")
        print(f"    {os.path.join(plot_dir, 'CNVscore_by_group_per_sample.png')}")

        return self
    def _plot_volcano(self, de_results, title='Volcano Plot', genes_to_label=None):
        """Create volcano plot"""
        plt.figure(figsize=(10, 8))
        
        de_results['-log10_pval'] = -np.log10(de_results['pvalue'] + 1e-300)
        
        colors = []
        for _, row in de_results.iterrows():
            if row['pvalue_adj'] < 0.05 and abs(row['logfoldchange']) > 1:
                colors.append('red')
            elif row['pvalue_adj'] < 0.05:
                colors.append('orange')
            else:
                colors.append('gray')
        
        plt.scatter(
            de_results['logfoldchange'],
            de_results['-log10_pval'],
            c=colors,
            alpha=0.6,
            s=20
        )
        
        plt.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
        plt.axvline(x=1, color='black', linestyle='--', alpha=0.5)
        plt.axvline(x=-1, color='black', linestyle='--', alpha=0.5)
        
        if genes_to_label:
            for gene in genes_to_label:
                if gene in de_results['gene'].values:
                    gene_data = de_results[de_results['gene'] == gene].iloc[0]
                    plt.annotate(
                        gene,
                        xy=(gene_data['logfoldchange'], gene_data['-log10_pval']),
                        xytext=(5, 5),
                        textcoords='offset points',
                        fontsize=10,
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7)
                    )
        
        plt.xlabel('Log2 Fold Change')
        plt.ylabel('-Log10 P-value')
        plt.title(title)
        plt.grid(True, alpha=0.3)
    
    def _plot_gene_expression_by_condition(self, genes, groupby='tissue_type', title=''):
        """Create violin plots"""
        available_genes = [g for g in genes if g in self.adata.raw.var_names]
        
        if not available_genes:
            return
        
        n_genes = len(available_genes)
        n_cols = min(3, n_genes)
        n_rows = (n_genes + n_cols - 1) // n_cols
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 5*n_rows))
        if n_genes == 1:
            axes = [axes]
        else:
            axes = axes.flatten()
        
        for idx, gene in enumerate(available_genes):
            if issparse(self.adata.raw.X):
                gene_expr = self.adata.raw[:, gene].X.toarray().flatten()
            else:
                gene_expr = self.adata.raw[:, gene].X.flatten()
            
            plot_df = pd.DataFrame({
                'expression': gene_expr,
                'group': self.adata.obs[groupby].values
            })
            
            sns.violinplot(data=plot_df, x='group', y='expression', ax=axes[idx])
            axes[idx].set_title(gene)
            axes[idx].set_xlabel('')
            axes[idx].set_ylabel('Log(Expression + 1)')
            axes[idx].tick_params(axis='x', rotation=45)
            
            del gene_expr, plot_df
        
        for idx in range(n_genes, len(axes)):
            axes[idx].axis('off')
        
        plt.suptitle(title, fontsize=14, y=1.02)
        plt.tight_layout()
    
    def _plot_genes_along_pseudotime_sparse(self, adata_subset, genes):
        """Plot gene expression along pseudotime"""
        available_genes = [g for g in genes if g in adata_subset.raw.var_names]
        
        if not available_genes:
            return
        
        n_genes = len(available_genes)
        fig, axes = plt.subplots(n_genes, 1, figsize=(10, 4*n_genes))
        
        if n_genes == 1:
            axes = [axes]
        
        for idx, gene in enumerate(available_genes):
            if issparse(adata_subset.raw.X):
                gene_expr = adata_subset.raw[:, gene].X.toarray().flatten()
            else:
                gene_expr = adata_subset.raw[:, gene].X.flatten()
            
            pseudotime = adata_subset.obs['dpt_pseudotime'].values
            sort_idx = np.argsort(pseudotime)
            
            axes[idx].scatter(
                pseudotime,
                gene_expr,
                alpha=0.3,
                s=10,
                c=adata_subset.obs['tissue_type'].map({'cancer': 'red', 'normal': 'blue'})
            )
            
            from scipy.ndimage import gaussian_filter1d
            smooth_expr = gaussian_filter1d(gene_expr[sort_idx], sigma=20)
            axes[idx].plot(pseudotime[sort_idx], smooth_expr, 'k-', linewidth=2, alpha=0.7)
            
            axes[idx].set_xlabel('Pseudotime')
            axes[idx].set_ylabel('Expression')
            axes[idx].set_title(f'{gene} Along Trajectory')
            axes[idx].grid(True, alpha=0.3)
            
            del gene_expr
        
        plt.tight_layout()
    


def main():
    """
    Main analysis function with parallelization
    
    CONFIGURATION:
    - Edit N_THREADS at the top of this file to match your CPU cores
    """
    
    # ========== DATA CONFIGURATION ==========
    UMI_MATRIX_PATH = 'all_combined_UMIs.txt.gz'
    METADATA_PATH = 'GSE160269_series_matrix.txt'
    OUTPUT_DIR = 'outputs'
    

    GENES_OF_INTEREST = [ #marker genes, not of main interest to manuscript but like might be good to check
        'TP53', 'CDKN2A', 'PIK3CA', 'NOTCH1',
        'CD3D', 'CD8A', 'CD4',
        'CD19', 'MS4A1',
        'CD68', 'CD14',
        'EPCAM', 'KRT18',
        'MKI67', 'PCNA'
    ]
    GENES_OF_INTEREST = [ #manuscript relevant (e.g., epistatic gene set)
        'NOTCH1', 'FAT1', 'TP53',
        'RB1', 'FBXW7', 'NFE2L2',
        'PIK3CA', 'NOTCH2'
    ]
    
    # ========== RUN PARALLELIZED ANALYSIS ==========
    
    analyzer = scRNAAnalyzerParallel(
        umi_path=UMI_MATRIX_PATH,
        metadata_path=METADATA_PATH,
        output_dir=OUTPUT_DIR
    )
    
    (analyzer
     .load_data()
     .preprocess(min_genes=200, min_cells=3, max_pct_mito=20)
     .classify_cell_types(use_markers=True)
     .differential_expression_cancer_vs_normal(GENES_OF_INTEREST)
     .differential_expression_by_cell_type(GENES_OF_INTEREST)
     .visualize_genes_low_dimensional(GENES_OF_INTEREST)
     .analyze_gene_coexpression(GENES_OF_INTEREST)
     .trajectory_analysis(root_cell_type='Epithelial')
     .perform_dense_operations()
     .plot_per_sample_cnv_summary(cnv_high_quantile=0.70)
    )
    
    

if __name__ == '__main__':
    main()
