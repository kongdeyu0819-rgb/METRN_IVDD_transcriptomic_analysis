#!/usr/bin/env python3
"""
Step 1: FDR Calculation — Differential Expression Analysis for METRN across GEO datasets
========================================================================================
Approach: Download GEO series matrix files directly (more reliable than GEOparse),
then perform DEA with FDR correction.

Datasets:
  - GSE124272: Peripheral blood, AIS vs healthy controls (discovery)
  - GSE150408: Peripheral blood, AIS vs controls (validation)  
  - GSE23130:  Annulus fibrosus/disc tissue, AIS vs controls (tissue)
"""

import os
import sys
import time
import warnings
import gzip
import urllib.request
import urllib.error

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings('ignore')

OUT_DIR = r'D:\workburry数据\AGENT\医学科研专家\METRN_AIS_科创变更补充\03_分析代码'
GEO_CACHE = os.path.join(OUT_DIR, 'geo_cache')
os.makedirs(GEO_CACHE, exist_ok=True)

# GEO FTP base URL for series matrix files
GEO_FTP = 'https://ftp.ncbi.nlm.nih.gov/geo/series'

DATASETS = {
    'GSE124272': {
        'tissue': 'Peripheral blood mononuclear cells',
        'phenotype': 'AIS vs healthy controls',
        'role': 'Discovery',
        'gpl': 'GPL570',   # Affymetrix Human Genome U133 Plus 2.0 Array
    },
    'GSE150408': {
        'tissue': 'Peripheral blood',
        'phenotype': 'AIS vs non-scoliosis controls',
        'role': 'Validation',
        'gpl': 'GPL570',   # Affymetrix Human Genome U133 Plus 2.0 Array (verify)
    },
    'GSE23130': {
        'tissue': 'Annulus fibrosus (disc tissue)',
        'phenotype': 'AIS disc vs control disc',
        'role': 'Tissue-level exploratory',
        'gpl': 'GPL96',    # Affymetrix Human Genome U133A Array (verify)
    },
}


def gse_to_ftp_path(gse_id):
    """Convert GSE ID to FTP directory path."""
    # GSE IDs are grouped: GSE124272 -> GSE124nnn
    prefix = gse_id[:-3] + 'nnn'
    return f"{GEO_FTP}/{prefix}/{gse_id}/matrix/"


def download_series_matrix(gse_id):
    """
    Download GEO series matrix file.
    Returns path to local file, or None if failed.
    """
    cache_file = os.path.join(GEO_CACHE, f"{gse_id}_series_matrix.txt.gz")
    cache_file_txt = os.path.join(GEO_CACHE, f"{gse_id}_series_matrix.txt")
    
    # Check cache first
    if os.path.exists(cache_file):
        print(f"  [Cache] Found {cache_file}")
        return cache_file
    if os.path.exists(cache_file_txt):
        print(f"  [Cache] Found {cache_file_txt}")
        return cache_file_txt
    
    # Try to download
    ftp_path = gse_to_ftp_path(gse_id)
    # Try common file naming patterns
    patterns = [
        f"{gse_id}_series_matrix.txt.gz",
        f"{gse_id}_Homo_sapiens_Series_matrix.txt.gz",
    ]
    
    for pattern in patterns:
        url = ftp_path + pattern
        print(f"  [Download] Trying: {url}")
        try:
            urllib.request.urlretrieve(url, cache_file, 
                                       reporthook=lambda b, bs, sz: None)
            # Check if file is valid
            if os.path.exists(cache_file) and os.path.getsize(cache_file) > 1000:
                print(f"  [Download] Success: {cache_file}")
                return cache_file
        except Exception as e:
            print(f"  [Download] Failed: {e}")
            continue
    
    # Try alternative: use GEO REST API
    print(f"  [Download] Trying GEO REST API...")
    rest_url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse_id}&targ=self&form=text&view=full"
    try:
        req = urllib.request.Request(rest_url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=30) as response:
            data = response.read().decode('utf-8')
            # Check if this is a SOFT file
            if '!dataset_table_begin' in data or 'ID_REF' in data:
                # Parse SOFT format
                pass
    except Exception as e:
        print(f"  [Download] REST API failed: {e}")
    
    print(f"  [Download] All download attempts failed for {gse_id}")
    return None


def parse_series_matrix(filepath):
    """
    Parse GEO series matrix file.
    Returns: (expr_df, sample_info)
      - expr_df: DataFrame with probes as rows, samples as columns
      - sample_info: dict with sample metadata
    """
    print(f"  [Parse] Reading {filepath}")
    
    # Handle gzipped files
    open_func = gzip.open if filepath.endswith('.gz') else open
    mode = 'rt' if filepath.endswith('.gz') else 'r'
    
    sample_info = {}
    expr_data = {}
    probe_ids = []
    is_data_section = False
    header = None
    
    with open_func(filepath, mode, encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.strip()
            
            # Parse sample metadata (lines starting with !)
            if line.startswith('!Sample_'):
                parts = line.split('\t', 1)
                if len(parts) == 2:
                    key = parts[0]
                    vals = parts[1].split('\t')
                    if 'title' in key.lower() or 'characteristics' in key.lower() or 'source' in key.lower():
                        sample_info[key] = vals
            
            # Data section starts after '!dataset_table_end' or specific marker
            if '!dataset_table_begin' in line or 'ID_REF' in line:
                is_data_section = True
                header = line.split('\t')
                continue
            
            if is_data_section and header is not None:
                # Data rows
                parts = line.split('\t')
                if len(parts) >= 2:
                    probe_id = parts[0]
                    values = [float(v) if v else np.nan for v in parts[1:]]
                    expr_data[probe_id] = values
    
    if not expr_data:
        # Try alternative parsing for SOFT format
        print("  [Parse] Trying alternative parsing (SOFT format)...")
        return parse_soft_format(filepath)
    
    # Convert to DataFrame
    columns = header[1:] if header and len(header) > 1 else [f"Sample_{i}" for i in range(len(list(expr_data.values())[0]))]
    expr_df = pd.DataFrame(expr_data).T
    expr_df.columns = columns[:expr_df.shape[1]]
    
    print(f"  [Parse] Parsed {expr_df.shape[0]} probes x {expr_df.shape[1]} samples")
    return expr_df, sample_info


def parse_soft_format(filepath):
    """Parse GEO SOFT format file."""
    open_func = gzip.open if filepath.endswith('.gz') else open
    mode = 'rt' if filepath.endswith('.gz') else 'r'
    
    expr_data = {}
    sample_ids = []
    in_table = False
    
    with open_func(filepath, mode, encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.strip()
            
            if line.startswith('^SERIES'):
                continue
            elif line.startswith('!subset_description'):
                continue
            
            # Sample IDs
            if line.startswith('!Sample_geo_accession'):
                sample_ids = line.split('\t')[1:]
                continue
            
            # Data table
            if line.startswith('ID_REF') or line.startswith('ID'):
                in_table = True
                continue
            
            if in_table and line and not line.startswith('!'):
                parts = line.split('\t')
                if len(parts) > 1:
                    probe_id = parts[0]
                    values = [float(v) if v and v != 'NA' else np.nan for v in parts[1:]]
                    if len(values) == len(sample_ids):
                        expr_data[probe_id] = values
    
    if not expr_data:
        print("  [Parse] ERROR: Could not parse file")
        return None, {}
    
    columns = sample_ids if sample_ids else [f"Sample_{i}" for i in range(len(list(expr_data.values())[0]))]
    expr_df = pd.DataFrame(expr_data).T
    expr_df.columns = columns[:expr_df.shape[1]]
    
    print(f"  [Parse] Parsed (SOFT) {expr_df.shape[0]} probes x {expr_df.shape[1]} samples")
    return expr_df, {}


def label_samples_from_series_matrix(filepath, gse_id):
    """
    Extract sample labels from series matrix file metadata.
    Returns dict: {sample_id: 'AIS' or 'Control' or 'Unknown'}
    """
    open_func = gzip.open if filepath.endswith('.gz') else open
    mode = 'rt' if filepath.endswith('.gz') else 'r'
    
    sample_ids = []
    sample_titles = []
    sample_chars = []
    sample_groups = []
    
    with open_func(filepath, mode, encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.strip()
            if line.startswith('!Sample_geo_accession'):
                sample_ids = line.split('\t')[1:]
            elif line.startswith('!Sample_title'):
                sample_titles = line.split('\t')[1:]
            elif line.startswith('!Sample_characteristics_ch1'):
                sample_chars = line.split('\t')[1:]
            elif line.startswith('!Sample_source_name_ch1'):
                sample_groups = line.split('\t')[1:]
    
    labels = {}
    for i, sid in enumerate(sample_ids):
        title = sample_titles[i] if i < len(sample_titles) else ''
        chars = sample_chars[i] if i < len(sample_chars) else ''
        group = sample_groups[i] if i < len(sample_groups) else ''
        
        all_text = (title + ' ' + chars + ' ' + group).lower()
        
        # Label based on dataset
        if gse_id == 'GSE124272':
            # Known: this dataset has AIS patients and healthy controls
            # Titles typically contain "AIS" or "control"
            if 'ais' in all_text or 'scoliosis' in all_text or 'patient' in all_text:
                label = 'AIS'
            elif 'control' in all_text or 'healthy' in all_text or 'normal' in all_text:
                label = 'Control'
            else:
                label = 'Unknown'
        
        elif gse_id == 'GSE150408':
            if 'ais' in all_text or 'scoliosis' in all_text:
                label = 'AIS'
            elif 'control' in all_text or 'healthy' in all_text or 'non' in all_text:
                label = 'Control'
            else:
                label = 'Unknown'
        
        elif gse_id == 'GSE23130':
            if 'scoliosis' in all_text or 'ais' in all_text:
                label = 'AIS'
            elif 'control' in all_text or 'normal' in all_text:
                label = 'Control'
            else:
                label = 'Unknown'
        
        labels[sid] = {
            'label': label,
            'title': title,
            'characteristics': chars,
            'source': group,
        }
    
    return labels


def annotate_probes_from_gpl(gse_id, expr_df):
    """
    Annotate probes using Bioconductor annotation packages or GEO GPL.
    For now, use a simple approach: map common Affymetrix probe IDs to genes.
    
    This is a simplified version - in practice, use org.Hs.eg.db or similar.
    """
    # For GPL570 (Affy U133 Plus 2.0), we can use a pre-built mapping
    # For now, return empty dict and handle annotation externally
    
    # Try to download GPL annotation from GEO
    gpl_cache = os.path.join(GEO_CACHE, f"{gse_id}_gpl_annotation.txt")
    
    # Common approach: use biomaRt or pre-built mapping
    # Since we can't easily access biomaRt from Python without extra deps,
    # we'll output probe IDs and let user annotate via R/Bioconductor
    
    print("  [Annotate] Probe annotation requires external mapping (see R script)")
    print("  [Annotate] Will proceed with probe IDs; gene symbols added via R")
    
    return {}


def run_dea_simple(expr_df, labels):
    """
    Simplified DEA: t-test + FDR.
    expr_df: rows=probes, cols=samples
    labels: dict {sample_id: {label: 'AIS'/'Control'}}
    """
    ais_samples = [s for s, info in labels.items() if info['label'] == 'AIS']
    ctrl_samples = [s for s, info in labels.items() if info['label'] == 'Control']
    
    ais_samples = [s for s in ais_samples if s in expr_df.columns]
    ctrl_samples = [s for s in ctrl_samples if s in expr_df.columns]
    
    print(f"  DEA: AIS n={len(ais_samples)}, Control n={len(ctrl_samples)}")
    
    if len(ais_samples) < 2 or len(ctrl_samples) < 2:
        print(f"  ERROR: Insufficient samples")
        return pd.DataFrame()
    
    results = []
    for probe in expr_df.index:
        ais_vals = pd.to_numeric(expr_df.loc[probe, ais_samples], errors='coerce').dropna().values
        ctrl_vals = pd.to_numeric(expr_df.loc[probe, ctrl_samples], errors='coerce').dropna().values
        
        if len(ais_vals) < 2 or len(ctrl_vals) < 2:
            continue
        
        mean_ais = np.mean(ais_vals)
        mean_ctrl = np.mean(ctrl_vals)
        log2fc = mean_ais - mean_ctrl
        
        # Welch's t-test
        try:
            t_stat, p_value = stats.ttest_ind(ais_vals, ctrl_vals, equal_var=False)
        except Exception:
            p_value = 1.0
        
        results.append({
            'probe_id': probe,
            'log2FC': log2fc,
            'mean_AIS': mean_ais,
            'mean_Control': mean_ctrl,
            'n_AIS': len(ais_vals),
            'n_Control': len(ctrl_vals),
            'raw_P': p_value,
        })
    
    if not results:
        return pd.DataFrame()
    
    result_df = pd.DataFrame(results)
    
    # FDR correction
    valid_p = np.clip(result_df['raw_P'].values, 1e-300, 1.0)
    reject, fdr, _, _ = multipletests(valid_p, method='fdr_bh', alpha=0.05)
    result_df['FDR'] = fdr
    result_df['significant_FDR005'] = reject
    result_df['significant_rawP005'] = result_df['raw_P'] < 0.05
    
    result_df = result_df.sort_values('raw_P').reset_index(drop=True)
    
    print(f"  Total probes: {len(result_df)}")
    print(f"  Significant (raw P<0.05): {(result_df['raw_P']<0.05).sum()}")
    print(f"  Significant (FDR<0.05): {(result_df['FDR']<0.05).sum()}")
    
    return result_df


def main():
    print("=" * 70)
    print("STEP 1: FDR Calculation — GEO Differential Expression Analysis")
    print("=" * 70)
    
    all_metrn_summary = []
    
    for gse_id in DATASETS:
        config = DATASETS[gse_id]
        print(f"\n{'='*50}")
        print(f"Processing {gse_id}: {config['tissue']}")
        print(f"{'='*50}")
        
        # Download series matrix
        matrix_file = download_series_matrix(gse_id)
        if matrix_file is None:
            print(f"  SKIPPING {gse_id}: Could not download data")
            continue
        
        # Parse
        parsed = parse_series_matrix(matrix_file)
        if parsed is None or parsed[0] is None:
            print(f"  SKIPPING {gse_id}: Could not parse data")
            continue
        expr_df, sample_info = parsed
        
        # Label samples
        labels = label_samples_from_series_matrix(matrix_file, gse_id)
        known = sum(1 for v in labels.values() if v['label'] != 'Unknown')
        print(f"  Labeled samples: {known}/{len(labels)}")
        for lbl in ['AIS', 'Control']:
            n = sum(1 for v in labels.values() if v['label'] == lbl)
            print(f"    {lbl}: {n}")
        
        # Show unlabeled samples for debugging
        unknown_samples = [s for s, v in labels.items() if v['label'] == 'Unknown']
        if unknown_samples:
            print(f"  Unknown samples (first 3):")
            for s in unknown_samples[:3]:
                print(f"    {s}: title='{labels[s]['title']}'")
        
        # Run DEA
        result_df = run_dea_simple(expr_df, labels)
        if result_df.empty:
            continue
        
        # For now, we can't annotate probes in Python easily
        # Save results with probe IDs; annotation will be done in R
        result_df['gene_symbol'] = ''  # placeholder
        
        # Save
        csv_path = os.path.join(OUT_DIR, f'FDR_results_{gse_id}.csv')
        result_df.to_csv(csv_path, index=False)
        print(f"  Saved: {csv_path}")
        
        # Try to find METRN probe by known Affymetrix IDs
        # For GPL570 (U133 Plus 2.0), METRN is often probed by:
        #  (will be annotated in R)
        print(f"  NOTE: Probe-to-gene annotation will be completed in R script")
        print(f"  METRN results will be extracted after annotation")
    
    print(f"\n{'='*50}")
    print("Step 1 complete. Proceed to R script for probe annotation and METRN extraction.")
    print(f"{'='*50}")


if __name__ == '__main__':
    main()
