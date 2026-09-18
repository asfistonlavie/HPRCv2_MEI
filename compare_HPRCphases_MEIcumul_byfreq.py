#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MEI discovery cumulative analysis comparing phased and unphased haplotype data

This script generates publication-quality figures showing cumulative 
Mobile Element Insertion (MEI) discovery across two phases:
- Phase 1: HAP1 haplotypes 
- Phase 2: HAP2 haplotypes 

The analysis categorizes MEIs by population frequency:
- Rare: <1% population frequency
- Polymorphic: 1-95% population frequency  
- Fixed: >95% population frequency

Requirements:
    Python 3.8+ with numpy, pandas, matplotlib

Usage:
    python mei_cumulative_analysis.py -g <genomes_list.txt> -m <mei_data.txt>
    python mei_cumulative_analysis.py --help

Input files:
    -g, --genomes: List of priority genome identifiers (one per line, e.g., HG01123)
    -m, --mei-data: MEI genotyping and frequency data (tab-separated file)

Output:
    - MEI_cumulative_phases_comparison.png  (PNG raster, 300 DPI)
    - MEI_cumulative_phases_comparison.pdf  (PDF vector, for print)
    - MEI_cumulative_phases_comparison.svg  (SVG vector, for web/editing)
    - MEI_cumulative_phases_comparison.html (HTML report with plots and statistics)
    - Console output: Statistical summary

Author: Anna-Sophie Fiston-Lavier
Date: 2026
Version: 1
"""

import sys
import re
import argparse
from collections import defaultdict
from pathlib import Path
from datetime import datetime
from typing import Tuple, Dict, Set, List, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import ScalarFormatter


# ============================================================================
# CONSTANTS AND CONFIGURATION
# ============================================================================

# Publication-ready color scheme
COLORS = {
    # Phase 2 (Total data) - saturated colors
    'phase2_fixed': '#8B0000',      # Dark red
    'phase2_polymorphic': '#E6B800', # Dark yellow/gold
    'phase2_rare': '#2E75B6',        # Deep blue
    
    # Phase 1 (Prioritized) - pastel colors with hatching
    'phase1_fixed': '#FFCCCC',       # Light red
    'phase1_polymorphic': '#FFF3CC', # Light yellow
    'phase1_rare': '#CCEEFF',        # Light blue
}

# Figure dimensions
FIGURE_WIDTH_INCHES = 7.2
FIGURE_HEIGHT_INCHES = 4.8

# Output base name
OUTPUT_BASE = "MEI_cumulative_phases_comparison"

# Line styles
LINE_WIDTH_PHASE2 = 3.0
LINE_WIDTH_PHASE1 = 2.0
MARKER_SIZE = 3
MARKER_INTERVAL = 10

# Hatching pattern for Phase 1
HATCH_PATTERN = '///'
HATCH_ALPHA = 0.9

# Font settings
FONT_FAMILY = 'Arial'
FONT_SIZE_TITLE = 11
FONT_SIZE_AXES = 9
FONT_SIZE_LABEL = 9
FONT_SIZE_LEGEND = 8
FONT_SIZE_ANNOTATION = 8

# Frequency thresholds (as proportions)
RARE_THRESHOLD = 0.01    # <1%
FIXED_THRESHOLD = 0.95   # >95%


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def show_help() -> None:
    """Display detailed help information about the script."""
    help_text = """
================================================================================
MEI CUMULATIVE DISCOVERY ANALYSIS - HELP
================================================================================

DESCRIPTION:
    This script analyzes Mobile Element Insertion (MEI) discovery rates across
    phased haplotype data. It generates cumulative discovery curves comparing
    prioritized haplotypes (Phase 1) versus all available haplotypes (Phase 2).

INPUT FILES:

    1. Genomes list file (-g, --genomes):
       - Text file containing one genome identifier per line
       - Format: HG01123, HG01258, etc.
       - These identifiers are used to select priority haplotypes (HAP1 & HAP2)
       
       Example content:
           HG01123
           HG01258
           HG01372

    2. MEI data file (-m, --mei-data):
       - Tab-separated file containing MEI genotyping and frequency data
       - Expected columns: chr, start, end, type, family, class, strand,
                           length, status, freq, count, total, samples
       - The 'samples' column should contain haplotype information
         (e.g., HG01123_HAP1, HG01123_HAP2)
       
       Expected format example:
           chr1    1000    1200    LINE1   L1HS    +   200    PASS    0.15    3    20    HG01123_HAP1,HG01258_HAP2

OUTPUT FILES:
    All outputs are saved in the current working directory:
    
    - PNG: High-resolution raster image (300 DPI) for quick viewing
    - PDF: Vector format suitable for print publication
    - SVG: Editable vector format for web and figure editing
    - HTML: Interactive report with embedded statistics and download links

STATISTICAL CATEGORIES:
    - Rare: MEI with population frequency < 1%
    - Polymorphic: MEI with population frequency between 1% and 95%
    - Fixed: MEI with population frequency > 95%

PHASES:
    - Phase 1: Cumulative discovery using only HAP1 and HAP2 from prioritized genomes
    - Phase 2: Cumulative discovery using all haplotypes in the dataset

USAGE EXAMPLES:

    Basic usage:
        python mei_cumulative_analysis.py -g genomes.txt -m mei_data.txt

    Using long option names:
        python mei_cumulative_analysis.py --genomes genomes.txt --mei-data mei_data.txt

    Display this help:
        python mei_cumulative_analysis.py --help

DEPENDENCIES:
    - Python 3.8 or higher
    - numpy
    - pandas
    - matplotlib

    Install with: pip install numpy pandas matplotlib

AUTHOR:
    Anna-Sophie Fiston-Lavier
    Website: https://annasfistonlavier.com/

VERSION: 2.1.0
================================================================================
"""
    print(help_text)


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="MEI cumulative discovery analysis for phased haplotype data",
        add_help=False,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        "-g", "--genomes",
        type=str,
        required=False,
        help="Path to file containing list of priority genome identifiers (one per line)"
    )
    
    parser.add_argument(
        "-m", "--mei-data",
        type=str,
        required=False,
        help="Path to MEI genotyping and frequency data file (tab-separated)"
    )
    
    parser.add_argument(
        "-h", "--help",
        action="store_true",
        help="Show this help message and exit"
    )
    
    args = parser.parse_args()
    
    # Show help if requested or if no arguments provided
    if args.help or (args.genomes is None and args.mei_data is None):
        show_help()
        sys.exit(0)
    
    # Validate required arguments
    if args.genomes is None:
        print("Error: Missing required argument: --genomes (-g)")
        print("Use --help for usage information")
        sys.exit(1)
    
    if args.mei_data is None:
        print("Error: Missing required argument: --mei-data (-m)")
        print("Use --help for usage information")
        sys.exit(1)
    
    return args


def normalize_id(identifier: Optional[str]) -> str:
    """Normalize genome identifier by removing whitespace and converting to uppercase."""
    if identifier is None:
        return ""
    identifier = str(identifier).replace("\ufeff", "")
    return re.sub(r"\s+", "", identifier).upper().strip()


def parse_haplotype(raw_haplotype: str) -> Tuple[str, str]:
    """Extract sample name and allele label from haplotype string."""
    if raw_haplotype is None:
        return "", ""
    
    s = str(raw_haplotype).replace("\ufeff", "")
    s = re.sub(r"\s+", "", s).upper().strip()
    
    match_hap = re.match(r"^(HG\d{5})_HAP([12])", s)
    if match_hap:
        return match_hap.group(1), f"HAP{match_hap.group(2)}"
    
    match_num = re.match(r"^(HG\d{5})_([12])", s)
    if match_num:
        return match_num.group(1), f"HAP{match_num.group(2)}"
    
    match_fuzzy = re.match(r"^(HG\d{5})(?:.*HAP([12]))?", s)
    if match_fuzzy:
        sample = match_fuzzy.group(1)
        allele = match_fuzzy.group(2)
        return sample, f"HAP{allele}" if allele else ""
    
    match_sample = re.match(r"^(HG\d{5})", s)
    if match_sample:
        return match_sample.group(1), ""
    
    return s, ""


def categorize_frequency(frequency: float) -> str:
    """Categorize MEI based on population frequency."""
    if frequency < RARE_THRESHOLD:
        return "rare"
    elif frequency < FIXED_THRESHOLD:
        return "common"
    else:
        return "fixed"


def setup_publication_style() -> None:
    """Configure matplotlib settings for publication style."""
    plt.rcParams['font.family'] = FONT_FAMILY
    plt.rcParams['font.size'] = FONT_SIZE_AXES
    plt.rcParams['axes.linewidth'] = 0.8
    plt.rcParams['xtick.major.width'] = 0.8
    plt.rcParams['ytick.major.width'] = 0.8
    plt.rcParams['xtick.direction'] = 'out'
    plt.rcParams['ytick.direction'] = 'out'
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    plt.rcParams['svg.fonttype'] = 'none'


# ============================================================================
# DATA LOADING AND PROCESSING
# ============================================================================

def load_priority_genomes(filepath: str) -> List[str]:
    """Load list of priority genomes from text file."""
    try:
        with open(filepath, "r", encoding="utf-8") as f:
            raw_genomes = [line.strip() for line in f if line.strip()]
        return sorted({normalize_id(x) for x in raw_genomes})
    except FileNotFoundError:
        print(f"Error: File '{filepath}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading '{filepath}': {e}")
        sys.exit(1)


def load_mei_data(filepath: str) -> pd.DataFrame:
    """Load MEI frequency data from tab-separated file."""
    columns = [
        "chr", "start", "end", "type", "family", "class", "strand",
        "length", "status", "freq", "count", "total", "samples"
    ]
    
    try:
        df = pd.read_csv(filepath, sep="\t", header=None, names=columns)
        df["category"] = df["freq"].apply(categorize_frequency)
        return df
    except FileNotFoundError:
        print(f"Error: File '{filepath}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading '{filepath}': {e}")
        sys.exit(1)


def build_haplotype_index(mei_df: pd.DataFrame) -> Tuple[Dict, Dict]:
    """Build haplotype-level index of MEI occurrences."""
    haplotype_data = defaultdict(lambda: {"rare": set(), "common": set(), "fixed": set()})
    haplotype_meta = {}
    
    for idx, row in mei_df.iterrows():
        samples_str = row["samples"]
        if pd.isna(samples_str):
            continue
        
        category = row["category"]
        mei_id = idx
        
        samples = [s for s in str(samples_str).split(",") if s.strip()]
        
        for raw_haplotype in samples:
            sample_name, allele = parse_haplotype(raw_haplotype)
            if not sample_name:
                continue
            
            haplotype_key = f"{sample_name}_{allele}" if allele else f"{sample_name}_UNKNOWN"
            haplotype_meta[haplotype_key] = (sample_name, allele)
            haplotype_data[haplotype_key][category].add(mei_id)
    
    return haplotype_data, haplotype_meta


def create_haplotype_dataframe(haplotype_data: Dict, haplotype_meta: Dict, 
                               priority_genomes: List[str]) -> pd.DataFrame:
    """Create DataFrame from haplotype index with priority flags."""
    records = []
    for haplotype_key, data in haplotype_data.items():
        sample_name, allele = haplotype_meta.get(haplotype_key, parse_haplotype(haplotype_key))
        records.append({
            "haplotype": haplotype_key,
            "sample_name": sample_name,
            "allele": allele,
            "rare": data["rare"],
            "common": data["common"],
            "fixed": data["fixed"],
            "total_count": len(data["rare"]) + len(data["common"]) + len(data["fixed"])
        })
    
    hap_df = pd.DataFrame(records)
    hap_df["sample_name"] = hap_df["sample_name"].apply(normalize_id)
    hap_df["allele"] = hap_df["allele"].apply(normalize_id)
    hap_df["is_priority"] = (
        hap_df["sample_name"].isin(priority_genomes) & 
        hap_df["allele"].isin(["HAP1", "HAP2"])
    )
    
    return hap_df


def compute_cumulative_curves(hap_df: pd.DataFrame) -> pd.DataFrame:
    """Compute cumulative MEI discovery curves for a haplotype set."""
    cumulative_rare = []
    cumulative_common = []
    cumulative_fixed = []
    
    seen_rare = set()
    seen_common = set()
    seen_fixed = set()
    
    for _, row in hap_df.iterrows():
        seen_rare.update(row["rare"])
        seen_common.update(row["common"])
        seen_fixed.update(row["fixed"])
        cumulative_rare.append(len(seen_rare))
        cumulative_common.append(len(seen_common))
        cumulative_fixed.append(len(seen_fixed))
    
    result_df = hap_df.copy()
    result_df["cumul_rare"] = cumulative_rare
    result_df["cumul_common"] = cumulative_common
    result_df["cumul_fixed"] = cumulative_fixed
    result_df["stack_fixed"] = result_df["cumul_fixed"]
    result_df["stack_common"] = result_df["cumul_fixed"] + result_df["cumul_common"]
    result_df["stack_rare"] = result_df["cumul_fixed"] + result_df["cumul_common"] + result_df["cumul_rare"]
    
    return result_df


def sort_haplotypes_for_plotting(hap_df: pd.DataFrame, prioritize: bool = True) -> pd.DataFrame:
    """Sort haplotypes for cumulative plotting."""
    allele_order = {"HAP1": 0, "HAP2": 1, "UNKNOWN": 2, "": 3}
    hap_df = hap_df.copy()
    hap_df["allele_order"] = hap_df["allele"].map(allele_order).fillna(4).astype(int)
    
    if prioritize:
        sort_cols = ["is_priority", "sample_name", "allele_order", "total_count"]
        ascending = [False, True, True, True]
    else:
        sort_cols = ["sample_name", "allele_order", "total_count"]
        ascending = [True, True, True]
    
    return hap_df.sort_values(by=sort_cols, ascending=ascending).reset_index(drop=True)


# ============================================================================
# PLOTTING FUNCTIONS
# ============================================================================

def create_figure() -> Tuple[plt.Figure, plt.Axes]:
    """Create figure and axes for publication-quality plot."""
    fig, ax = plt.subplots(figsize=(FIGURE_WIDTH_INCHES, FIGURE_HEIGHT_INCHES))
    return fig, ax


def plot_phase2_areas(ax: plt.Axes, hap_df_total: pd.DataFrame) -> None:
    """Plot filled areas for Phase 2 (total data)."""
    x = np.arange(len(hap_df_total))
    ax.fill_between(x, 0, hap_df_total["stack_fixed"].values, 
                     color=COLORS['phase2_fixed'], zorder=1, alpha=0.9)
    ax.fill_between(x, hap_df_total["stack_fixed"].values, 
                     hap_df_total["stack_common"].values, 
                     color=COLORS['phase2_polymorphic'], zorder=1, alpha=0.9)
    ax.fill_between(x, hap_df_total["stack_common"].values, 
                     hap_df_total["stack_rare"].values, 
                     color=COLORS['phase2_rare'], zorder=1, alpha=0.9)


def plot_phase2_lines(ax: plt.Axes, hap_df_total: pd.DataFrame) -> None:
    """Plot line boundaries for Phase 2."""
    x = np.arange(len(hap_df_total))
    ax.plot(x, hap_df_total["stack_fixed"].values, 
             color=COLORS['phase2_fixed'], linewidth=LINE_WIDTH_PHASE2, zorder=4)
    ax.plot(x, hap_df_total["stack_common"].values, 
             color=COLORS['phase2_polymorphic'], linewidth=LINE_WIDTH_PHASE2, zorder=4)
    ax.plot(x, hap_df_total["stack_rare"].values, 
             color=COLORS['phase2_rare'], linewidth=LINE_WIDTH_PHASE2, zorder=4)


def plot_phase1_areas(ax: plt.Axes, hap_df_priority: pd.DataFrame) -> None:
    """Plot hatched areas for Phase 1 (prioritized data)."""
    if len(hap_df_priority) == 0:
        return
    x = np.arange(len(hap_df_priority))
    ax.fill_between(x, 0, hap_df_priority["stack_fixed"].values, 
                     color=COLORS['phase1_fixed'], hatch=HATCH_PATTERN, 
                     zorder=2, alpha=HATCH_ALPHA, edgecolor='none')
    ax.fill_between(x, hap_df_priority["stack_fixed"].values, 
                     hap_df_priority["stack_common"].values, 
                     color=COLORS['phase1_polymorphic'], hatch=HATCH_PATTERN, 
                     zorder=2, alpha=HATCH_ALPHA, edgecolor='none')
    ax.fill_between(x, hap_df_priority["stack_common"].values, 
                     hap_df_priority["stack_rare"].values, 
                     color=COLORS['phase1_rare'], hatch=HATCH_PATTERN, 
                     zorder=2, alpha=HATCH_ALPHA, edgecolor='none')


def plot_phase1_lines(ax: plt.Axes, hap_df_priority: pd.DataFrame) -> None:
    """Plot line boundaries for Phase 1 with markers."""
    if len(hap_df_priority) == 0:
        return
    x = np.arange(len(hap_df_priority))
    ax.plot(x, hap_df_priority["stack_fixed"].values, 
             color=COLORS['phase1_fixed'], linewidth=LINE_WIDTH_PHASE1, 
             linestyle='--', marker='o', markersize=MARKER_SIZE, 
             markevery=MARKER_INTERVAL, zorder=5)
    ax.plot(x, hap_df_priority["stack_common"].values, 
             color=COLORS['phase1_polymorphic'], linewidth=LINE_WIDTH_PHASE1, 
             linestyle='--', marker='o', markersize=MARKER_SIZE, 
             markevery=MARKER_INTERVAL, zorder=5)
    ax.plot(x, hap_df_priority["stack_rare"].values, 
             color=COLORS['phase1_rare'], linewidth=LINE_WIDTH_PHASE1, 
             linestyle='--', marker='o', markersize=MARKER_SIZE, 
             markevery=MARKER_INTERVAL, zorder=5)


def add_category_labels(ax: plt.Axes, hap_df_total: pd.DataFrame, 
                        hap_df_priority: pd.DataFrame) -> None:
    """Add category labels and value annotations to the plot."""
    if len(hap_df_total) == 0:
        return
    
    x_max = len(hap_df_total) - 1
    label_x = x_max + x_max * 0.08
    
    fixed_y = hap_df_total["stack_fixed"].values[-1] / 2
    polymorphic_y = hap_df_total["stack_fixed"].values[-1] + (hap_df_total["stack_common"].values[-1] - hap_df_total["stack_fixed"].values[-1]) / 2
    rare_y = hap_df_total["stack_common"].values[-1] + (hap_df_total["stack_rare"].values[-1] - hap_df_total["stack_common"].values[-1]) / 2
    
    ax.text(label_x, rare_y, "Rare", 
             color='black', va='center', ha='left', fontsize=FONT_SIZE_LABEL)
    ax.text(label_x, polymorphic_y, "Polymorphic", 
             color='black', va='center', ha='left', fontsize=FONT_SIZE_LABEL)
    ax.text(label_x, fixed_y, "Fixed", 
             color='black', va='center', ha='left', fontsize=FONT_SIZE_LABEL)
    
    annot_x = x_max + 2
    total_rare = hap_df_total['cumul_rare'].values[-1]
    total_common = hap_df_total['cumul_common'].values[-1]
    total_fixed = hap_df_total['cumul_fixed'].values[-1]
    
    if total_rare > 0:
        ax.text(annot_x, rare_y, f"{int(total_rare)}", 
                 color='black', va='center', ha='left', fontsize=FONT_SIZE_ANNOTATION)
    if total_common > 0:
        ax.text(annot_x, polymorphic_y, f"{int(total_common)}", 
                 color='black', va='center', ha='left', fontsize=FONT_SIZE_ANNOTATION)
    if total_fixed > 0:
        ax.text(annot_x, fixed_y, f"{int(total_fixed)}", 
                 color='black', va='center', ha='left', fontsize=FONT_SIZE_ANNOTATION)
    
    if len(hap_df_priority) > 0:
        p1_x_max = len(hap_df_priority) - 1
        p1_annot_x = p1_x_max + 2
        p1_rare = hap_df_priority["cumul_rare"].values[-1]
        p1_common = hap_df_priority["cumul_common"].values[-1]
        p1_fixed = hap_df_priority["cumul_fixed"].values[-1]
        
        p1_fixed_y = hap_df_priority["stack_fixed"].values[-1] / 2
        p1_polymorphic_y = hap_df_priority["stack_fixed"].values[-1] + (hap_df_priority["stack_common"].values[-1] - hap_df_priority["stack_fixed"].values[-1]) / 2
        p1_rare_y = hap_df_priority["stack_common"].values[-1] + (hap_df_priority["stack_rare"].values[-1] - hap_df_priority["stack_common"].values[-1]) / 2
        
        if p1_rare > 0:
            ax.text(p1_annot_x, p1_rare_y, f"{int(p1_rare)}", 
                     color='black', va='center', ha='left', fontsize=FONT_SIZE_ANNOTATION)
        if p1_common > 0:
            ax.text(p1_annot_x, p1_polymorphic_y, f"{int(p1_common)}", 
                     color='black', va='center', ha='left', fontsize=FONT_SIZE_ANNOTATION)
        if p1_fixed > 0:
            ax.text(p1_annot_x, p1_fixed_y, f"{int(p1_fixed)}", 
                     color='black', va='center', ha='left', fontsize=FONT_SIZE_ANNOTATION)


def add_legend(ax: plt.Axes) -> None:
    """Add publication-style legend to the plot."""
    legend_elements = [
        Patch(facecolor=COLORS['phase2_fixed'], label="Fixed (Phase 2)"),
        Patch(facecolor=COLORS['phase2_polymorphic'], label="Polymorphic (Phase 2)"),
        Patch(facecolor=COLORS['phase2_rare'], label="Rare (Phase 2)"),
        Patch(facecolor=COLORS['phase1_fixed'], hatch=HATCH_PATTERN, 
              edgecolor='white', linewidth=0.5, label="Fixed (Phase 1)"),
        Patch(facecolor=COLORS['phase1_polymorphic'], hatch=HATCH_PATTERN, 
              edgecolor='white', linewidth=0.5, label="Polymorphic (Phase 1)"),
        Patch(facecolor=COLORS['phase1_rare'], hatch=HATCH_PATTERN, 
              edgecolor='white', linewidth=0.5, label="Rare (Phase 1)"),
    ]
    ax.legend(handles=legend_elements, loc="upper left", ncol=1, 
              fontsize=FONT_SIZE_LEGEND, frameon=False)


def configure_axes(ax: plt.Axes, hap_df_total: pd.DataFrame, 
                   hap_df_priority: pd.DataFrame) -> None:
    """Configure axes labels, limits, and appearance."""
    ax.set_xlabel("Nth additional haplotype", fontsize=FONT_SIZE_AXES, color='black')
    ax.set_ylabel("Cumulative MEI discoveries", fontsize=FONT_SIZE_AXES, color='black')
    
    max_y = max(
        hap_df_total["stack_rare"].max() if len(hap_df_total) > 0 else 0,
        hap_df_priority["stack_rare"].max() if len(hap_df_priority) > 0 else 0
    )
    
    x_max = len(hap_df_total) if len(hap_df_total) > 0 else 1
    ax.set_xlim(-5, x_max * 1.35)
    ax.set_ylim(0, max_y * 1.1)
    
    ax.grid(axis="y", alpha=0.3, linestyle='--', linewidth=0.5)
    ax.tick_params(axis='both', colors='black')
    
    if max_y > 10000:
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style='sci', axis='y', scilimits=(4, 4))


# ============================================================================
# OUTPUT FUNCTIONS
# ============================================================================

def save_figure(fig: plt.Figure, output_base: str = OUTPUT_BASE) -> Tuple[str, str]:
    """Save figure in multiple formats: PNG, PDF, and SVG."""
    png_path = f"{output_base}.png"
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    print(f"  {png_path}")
    
    pdf_path = f"{output_base}.pdf"
    fig.savefig(pdf_path, bbox_inches="tight", format='pdf')
    print(f"  {pdf_path}")
    
    svg_path = f"{output_base}.svg"
    fig.savefig(svg_path, bbox_inches="tight", format='svg')
    print(f"  {svg_path}")
    
    return png_path, svg_path


def get_statistics_dict(hap_df_total: pd.DataFrame, hap_df_priority: pd.DataFrame,
                        priority_genomes: List[str]) -> dict:
    """Collect statistics into a dictionary for HTML report."""
    stats = {
        'total_haplotypes': len(hap_df_total),
        'priority_haplotypes': len(hap_df_priority),
        'priority_genomes_requested': len(priority_genomes),
        'phase1': {'rare': 0, 'polymorphic': 0, 'fixed': 0, 'total': 0},
        'phase2': {'rare': 0, 'polymorphic': 0, 'fixed': 0, 'total': 0},
        'missing_genomes': [],
        'timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'rare_threshold': RARE_THRESHOLD * 100,
        'fixed_threshold': FIXED_THRESHOLD * 100
    }
    
    if len(hap_df_priority) > 0:
        stats['phase1']['rare'] = int(hap_df_priority['cumul_rare'].max())
        stats['phase1']['polymorphic'] = int(hap_df_priority['cumul_common'].max())
        stats['phase1']['fixed'] = int(hap_df_priority['cumul_fixed'].max())
        stats['phase1']['total'] = sum(stats['phase1'].values())
    
    if len(hap_df_total) > 0:
        stats['phase2']['rare'] = int(hap_df_total['cumul_rare'].max())
        stats['phase2']['polymorphic'] = int(hap_df_total['cumul_common'].max())
        stats['phase2']['fixed'] = int(hap_df_total['cumul_fixed'].max())
        stats['phase2']['total'] = sum(stats['phase2'].values())
    
    present_samples = set(hap_df_total["sample_name"].unique())
    missing = sorted(set(priority_genomes) - present_samples)
    stats['missing_genomes'] = missing[:20]
    stats['missing_count'] = len(missing)
    
    return stats


def generate_html_report(stats: dict, png_path: str, svg_path: str, 
                         output_base: str = OUTPUT_BASE) -> str:
    """Generate an HTML report with embedded plots and statistics."""
    
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>MEI Cumulative Discovery Analysis Report</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        body {{
            font-family: 'Arial', 'Helvetica', sans-serif;
            background-color: #f5f5f5;
            color: #333;
            line-height: 1.6;
        }}
        
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
        }}
        
        .header {{
            background: linear-gradient(135deg, #1a1a2e 0%, #16213e 100%);
            color: white;
            padding: 40px;
            border-radius: 12px;
            margin-bottom: 30px;
            text-align: center;
        }}
        
        .header h1 {{
            font-size: 28px;
            margin-bottom: 10px;
            font-weight: 600;
        }}
        
        .header p {{
            font-size: 14px;
            opacity: 0.8;
        }}
        
        .card {{
            background: white;
            border-radius: 12px;
            padding: 25px;
            margin-bottom: 25px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        
        .card h2 {{
            font-size: 20px;
            margin-bottom: 20px;
            padding-bottom: 10px;
            border-bottom: 2px solid #e0e0e0;
            color: #1a1a2e;
        }}
        
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin-bottom: 20px;
        }}
        
        .stat-box {{
            background: #f8f9fa;
            padding: 20px;
            border-radius: 8px;
            text-align: center;
            border-left: 4px solid;
        }}
        
        .stat-box.phase1 {{ border-left-color: #E6B800; }}
        .stat-box.phase2 {{ border-left-color: #2E75B6; }}
        
        .stat-number {{
            font-size: 32px;
            font-weight: bold;
            margin: 10px 0;
        }}
        
        .stat-label {{
            font-size: 14px;
            color: #666;
            text-transform: uppercase;
            letter-spacing: 1px;
        }}
        
        .category-table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 15px;
        }}
        
        .category-table th,
        .category-table td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #e0e0e0;
        }}
        
        .category-table th {{
            background-color: #f8f9fa;
            font-weight: 600;
            color: #1a1a2e;
        }}
        
        .category-table tr:hover {{
            background-color: #f8f9fa;
        }}
        
        .badge {{
            display: inline-block;
            padding: 4px 12px;
            border-radius: 20px;
            font-size: 12px;
            font-weight: 500;
        }}
        
        .badge.rare {{ background: #CCEEFF; color: #0066cc; }}
        .badge.polymorphic {{ background: #FFF3CC; color: #8B6914; }}
        .badge.fixed {{ background: #FFCCCC; color: #8B0000; }}
        
        .plot-container {{
            text-align: center;
            margin: 20px 0;
        }}
        
        .plot-container img {{
            max-width: 100%;
            height: auto;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}
        
        .download-links {{
            display: flex;
            gap: 15px;
            justify-content: center;
            margin-top: 20px;
            flex-wrap: wrap;
        }}
        
        .btn {{
            display: inline-block;
            padding: 10px 20px;
            background: #1a1a2e;
            color: white;
            text-decoration: none;
            border-radius: 6px;
            font-size: 14px;
            transition: background 0.3s;
        }}
        
        .btn:hover {{
            background: #16213e;
        }}
        
        .warning {{
            background: #fff3cd;
            border-left: 4px solid #ffc107;
            padding: 15px;
            margin-top: 20px;
            border-radius: 6px;
        }}
        
        .footer {{
            text-align: center;
            padding: 20px;
            font-size: 12px;
            color: #999;
        }}
        
        .footer a {{
            color: #1a1a2e;
            text-decoration: none;
        }}
        
        .footer a:hover {{
            text-decoration: underline;
        }}
        
        .color-indicator {{
            display: inline-block;
            width: 20px;
            height: 20px;
            border-radius: 4px;
            margin-right: 8px;
            vertical-align: middle;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>MEI Cumulative Discovery Analysis</h1>
            <p>Comparative analysis of phased (Phase 1) and unphased (Phase 2) haplotype data</p>
            <p><small>Generated on {stats['timestamp']}</small></p>
        </div>
        
        <div class="card">
            <h2>Summary Statistics</h2>
            <div class="stats-grid">
                <div class="stat-box phase1">
                    <div class="stat-label">Phase 1 (Prioritized)</div>
                    <div class="stat-number">{stats['priority_haplotypes']}</div>
                    <div>haplotypes from {stats['priority_genomes_requested']} genomes</div>
                </div>
                <div class="stat-box phase2">
                    <div class="stat-label">Phase 2 (All data)</div>
                    <div class="stat-number">{stats['total_haplotypes']}</div>
                    <div>total haplotypes analyzed</div>
                </div>
            </div>
            
            <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 20px;">
                <div>
                    <h3 style="font-size: 16px; margin-bottom: 15px;">Phase 1 - Prioritized Haplotypes</h3>
                    <table class="category-table">
                        <thead>
                            <tr><th>Category</th><th>Count</th><th>Frequency threshold</th></tr>
                        </thead>
                        <tbody>
                            <tr><td><span class="badge rare">Rare</span></td><td><strong>{stats['phase1']['rare']}</strong></td><td>&lt;{stats['rare_threshold']:.0f}%</td></tr>
                            <tr><td><span class="badge polymorphic">Polymorphic</span></td><td><strong>{stats['phase1']['polymorphic']}</strong></td><td>{stats['rare_threshold']:.0f}-{stats['fixed_threshold']:.0f}%</td></tr>
                            <tr><td><span class="badge fixed">Fixed</span></td><td><strong>{stats['phase1']['fixed']}</strong></td><td>&gt;{stats['fixed_threshold']:.0f}%</td></tr>
                            <tr style="font-weight: bold; background-color: #f0f0f0;"><td>TOTAL</td><td>{stats['phase1']['total']}</td><td></td></tr>
                        </tbody>
                    </table>
                </div>
                <div>
                    <h3 style="font-size: 16px; margin-bottom: 15px;">Phase 2 - All Haplotypes</h3>
                    <table class="category-table">
                        <thead>
                            <tr><th>Category</th><th>Count</th><th>Frequency threshold</th></tr>
                        </thead>
                        <tbody>
                            <tr><td><span class="badge rare">Rare</span></td><td><strong>{stats['phase2']['rare']}</strong></td><td>&lt;{stats['rare_threshold']:.0f}%</td></tr>
                            <tr><td><span class="badge polymorphic">Polymorphic</span></td><td><strong>{stats['phase2']['polymorphic']}</strong></td><td>{stats['rare_threshold']:.0f}-{stats['fixed_threshold']:.0f}%</td></tr>
                            <tr><td><span class="badge fixed">Fixed</span></td><td><strong>{stats['phase2']['fixed']}</strong></td><td>&gt;{stats['fixed_threshold']:.0f}%</td></tr>
                            <tr style="font-weight: bold; background-color: #f0f0f0;"><td>TOTAL</td><td>{stats['phase2']['total']}</td><td></td></tr>
                        </tbody>
                    </table>
                </div>
            </div>
        </div>
        
        <div class="card">
            <h2>Cumulative MEI Discovery Plot</h2>
            <div class="plot-container">
                <img src="{png_path}" alt="MEI Cumulative Discovery Plot">
            </div>
            <div class="download-links">
                <a href="{png_path}" download class="btn">Download PNG</a>
                <a href="{svg_path}" download class="btn">Download SVG</a>
                <a href="{output_base}.pdf" download class="btn">Download PDF</a>
            </div>
        </div>
        
        <div class="card">
            <h2>Color Legend</h2>
            <div style="display: flex; gap: 30px; flex-wrap: wrap; justify-content: center;">
                <div><span class="color-indicator" style="background: #8B0000;"></span> Fixed (Phase 2)</div>
                <div><span class="color-indicator" style="background: #E6B800;"></span> Polymorphic (Phase 2)</div>
                <div><span class="color-indicator" style="background: #2E75B6;"></span> Rare (Phase 2)</div>
                <div><span class="color-indicator" style="background: #FFCCCC; border: 1px solid #ccc;"></span> Fixed (Phase 1, hatched)</div>
                <div><span class="color-indicator" style="background: #FFF3CC; border: 1px solid #ccc;"></span> Polymorphic (Phase 1, hatched)</div>
                <div><span class="color-indicator" style="background: #CCEEFF; border: 1px solid #ccc;"></span> Rare (Phase 1, hatched)</div>
            </div>
        </div>
"""
    
    if stats['missing_count'] > 0:
        html_content += f"""
        <div class="card">
            <div class="warning">
                <strong>Warning:</strong> {stats['missing_count']} requested genome(s) were not found in the data:
                <ul style="margin-top: 10px; margin-left: 20px;">
                    {''.join(f'<li>{g}</li>' for g in stats['missing_genomes'][:10])}
                    {f'<li>... and {stats["missing_count"] - 10} more</li>' if stats['missing_count'] > 10 else ''}
                </ul>
            </div>
        </div>
"""
    
    html_content += f"""
        <div class="card">
            <h2>Methodology</h2>
            <ul style="margin-left: 20px;">
                <li><strong>Phase 1:</strong> Cumulative MEI discovery using only prioritized haplotypes (HAP1 & HAP2 from selected genomes)</li>
                <li><strong>Phase 2:</strong> Cumulative MEI discovery using all available haplotypes in the dataset</li>
                <li><strong>Rare:</strong> MEI with population frequency &lt; {stats['rare_threshold']:.0f}%</li>
                <li><strong>Polymorphic:</strong> MEI with population frequency between {stats['rare_threshold']:.0f}% and {stats['fixed_threshold']:.0f}%</li>
                <li><strong>Fixed:</strong> MEI with population frequency &gt; {stats['fixed_threshold']:.0f}%</li>
                <li>Phase 1 areas are shown with hatching (<code>{HATCH_PATTERN}</code>) to distinguish them from Phase 2</li>
                <li>Markers on Phase 1 lines are plotted every {MARKER_INTERVAL} haplotype</li>
            </ul>
        </div>
        
        <div class="footer">
            <p>Anna-Sophie Fiston-Lavier team| <a href="https://annasfistonlavier.com/">https://annasfistonlavier.com/</a></p>
        </div>
    </div>
</body>
</html>
"""
    
    html_path = f"{output_base}.html"
    with open(html_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    print(f"  {html_path}")
    return html_path


# ============================================================================
# STATISTICS AND REPORTING
# ============================================================================

def print_statistics(hap_df_total: pd.DataFrame, hap_df_priority: pd.DataFrame,
                     priority_genomes: List[str]) -> None:
    """Print statistical summary to console."""
    print("\n" + "="*60)
    print("STATISTICAL SUMMARY")
    print("="*60)
    print(f"Total haplotypes analyzed: {len(hap_df_total)}")
    print(f"Priority haplotypes (Phase 1): {len(hap_df_priority)}")
    print(f"Priority genomes requested: {len(priority_genomes)}")
    
    print("\nPhase 1 (Prioritized haplotypes only):")
    if len(hap_df_priority) > 0:
        rare_max = hap_df_priority['cumul_rare'].max()
        common_max = hap_df_priority['cumul_common'].max()
        fixed_max = hap_df_priority['cumul_fixed'].max()
        print(f"  Rare (<{RARE_THRESHOLD*100:.0f}%): {int(rare_max)}")
        print(f"  Polymorphic ({RARE_THRESHOLD*100:.0f}-{FIXED_THRESHOLD*100:.0f}%): {int(common_max)}")
        print(f"  Fixed (>{FIXED_THRESHOLD*100:.0f}%): {int(fixed_max)}")
        print(f"  TOTAL: {int(rare_max + common_max + fixed_max)}")
    else:
        print("  No priority haplotypes found")
    
    print("\nPhase 2 (All haplotypes):")
    if len(hap_df_total) > 0:
        rare_max = hap_df_total['cumul_rare'].max()
        common_max = hap_df_total['cumul_common'].max()
        fixed_max = hap_df_total['cumul_fixed'].max()
        print(f"  Rare (<{RARE_THRESHOLD*100:.0f}%): {int(rare_max)}")
        print(f"  Polymorphic ({RARE_THRESHOLD*100:.0f}-{FIXED_THRESHOLD*100:.0f}%): {int(common_max)}")
        print(f"  Fixed (>{FIXED_THRESHOLD*100:.0f}%): {int(fixed_max)}")
        print(f"  TOTAL: {int(rare_max + common_max + fixed_max)}")
    
    present_samples = set(hap_df_total["sample_name"].unique())
    missing = sorted(set(priority_genomes) - present_samples)
    if missing:
        print(f"\nWarning: {len(missing)} requested genomes not found in data:")
        for m in missing[:10]:
            print(f"    - {m}")
        if len(missing) > 10:
            print(f"    ... and {len(missing) - 10} more")


# ============================================================================
# MAIN PIPELINE
# ============================================================================

def main():
    """Main execution pipeline for MEI cumulative discovery analysis."""
    args = parse_arguments()
    
    genomes_file = args.genomes
    mei_file = args.mei_data
    
    setup_publication_style()
    
    print("\n" + "="*60)
    print("MEI CUMULATIVE DISCOVERY ANALYSIS")
    print("="*60)
    
    # Load data
    print("\n[1/5] Loading data...")
    print(f"  Priority genomes: {genomes_file}")
    priority_genomes = load_priority_genomes(genomes_file)
    print(f"  -> {len(priority_genomes)} genomes loaded")
    
    print(f"  MEI data: {mei_file}")
    mei_df = load_mei_data(mei_file)
    print(f"  -> {len(mei_df)} MEI records loaded")
    
    # Process data
    print("\n[2/5] Processing haplotype data...")
    haplotype_data, haplotype_meta = build_haplotype_index(mei_df)
    print(f"  -> {len(haplotype_data)} unique haplotypes identified")
    
    hap_df = create_haplotype_dataframe(haplotype_data, haplotype_meta, priority_genomes)
    print(f"  -> {len(hap_df)} haplotypes in final DataFrame")
    
    # Compute cumulative curves
    print("\n[3/5] Computing cumulative curves...")
    print("  Phase 2 (all haplotypes)...")
    hap_df_total = sort_haplotypes_for_plotting(hap_df, prioritize=True)
    hap_df_total = compute_cumulative_curves(hap_df_total)
    print(f"    -> {len(hap_df_total)} haplotypes sorted")
    
    print("  Phase 1 (priority haplotypes only)...")
    hap_df_priority = hap_df[hap_df["is_priority"]].copy()
    hap_df_priority = sort_haplotypes_for_plotting(hap_df_priority, prioritize=False)
    hap_df_priority = compute_cumulative_curves(hap_df_priority)
    print(f"    -> {len(hap_df_priority)} priority haplotypes")
    
    # Generate figure
    print("\n[4/5] Generating figure...")
    fig, ax = create_figure()
    
    plot_phase2_areas(ax, hap_df_total)
    plot_phase1_areas(ax, hap_df_priority)
    plot_phase2_lines(ax, hap_df_total)
    plot_phase1_lines(ax, hap_df_priority)
    add_category_labels(ax, hap_df_total, hap_df_priority)
    add_legend(ax)
    configure_axes(ax, hap_df_total, hap_df_priority)
    
    # Save outputs
    print("\n[5/5] Saving outputs...")
    print("  Figures:")
    png_path, svg_path = save_figure(fig)
    
    # Print statistics to console
    print_statistics(hap_df_total, hap_df_priority, priority_genomes)
    
    # Generate HTML report
    print("\n  HTML report:")
    stats = get_statistics_dict(hap_df_total, hap_df_priority, priority_genomes)
    html_path = generate_html_report(stats, png_path, svg_path)
    
    plt.close(fig)
    
    print("\n" + "="*60)
    print("ANALYSIS COMPLETE")
    print("="*60)
    print(f"\nOutput files generated:")
    print(f"  {OUTPUT_BASE}.png / .pdf / .svg")
    print(f"  {html_path}")
    print("\n")


if __name__ == "__main__":
    main()
