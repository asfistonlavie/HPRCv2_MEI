# HPRCv2_MEI - Mobile element detection and analysis in human pangenomic data (HPRC)

# MEI Cumulative Discovery Analysis

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://python.org)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey)]()

A Python pipeline for analyzing cumulative Mobile Element Insertion (MEI) discovery rates across phased haplotype data, generating publication-ready figures. This tool analyzes MEI genotyping data to generate cumulative discovery curves comparing two experimental phases:

- **Phase 1**: HPRC1 haplotypes (HAP1 and HAP2 from selected genomes)
- **Phase 2**: All available haplotypes in the dataset

The analysis categorizes MEIs by population frequency into three classes:

- **Rare**: <1% population frequency
- **Polymorphic**: 1-95% population frequency
- **Fixed**: >95% population frequency


## Table of Contents
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Input files](#input-files)
- [Output files](#output-files)
- [Output interpretation](#output-interpretation)
- [License & Authors](#license-authors)


## Requirements
- Python 3.8 or higher
- Required packages:
  - numpy ≥ 1.21.0
  - pandas ≥ 1.3.0
  - matplotlib ≥ 3.4.0


## Installation
### 1. Clone the repository

```bash
git clone https://github.com/annasophie/mei-cumulative-analysis.git
cd mei-cumulative-analysis
```

2. Install dependencies

Using pip:
```bash
pip install numpy pandas matplotlib
```
Or using conda:
```bash

conda install numpy pandas matplotlib
```
3. Verify installation
```bash

python mei_cumulative_analysis.py --help
```
## Usage
Basic usage
```bash

python mei_cumulative_analysis.py -g genomes.txt -m mei_data.txt
```
With long option names
```bash

python mei_cumulative_analysis.py --genomes genomes.txt --mei-data mei_data.txt
```
Display help
```bash

python mei_cumulative_analysis.py --help
```
## input-files
1. Genomes list file (-g, --genomes)

Text file containing one genome identifier per line:
```text

HG01123
HG01258
HG01372
HG01489
HG01595

These identifiers are used to select priority haplotypes (HAP1 and HAP2) for Phase 1 analysis.

Important notes:

    Identifiers are case-insensitive (HG01123 = hg01123)

    Spaces and special characters are automatically removed

    Duplicate entries are automatically deduplicated
```
2. MEI data file (-m, --mei-data)

```text
Tab-separated file containing MEI genotyping and frequency data with the following columns:
Column	Description	Type	Example
chr	Chromosome	string	chr1
start	Start position	integer	1000
end	End position	integer	1200
type	MEI type	string	LINE1
family	Repeat family	string	L1HS
class	Repeat class	string	LINE
strand	Strand orientation	string	+
length	Insert length	integer	200
status	Call status	string	PASS
freq	Population frequency	float	0.15
count	Allele count	integer	3
total	Total alleles	integer	20
samples	Haplotype samples	string	HG01123_HAP1,HG01258_HAP2

Critical requirement: The samples column must contain haplotype information in the format {SAMPLE_ID}_HAP1 or {SAMPLE_ID}_HAP2 (e.g., HG01123_HAP1).
```
Example line:
```text

chr1	1000	1200	LINE1	L1HS	LINE	+	200	PASS	0.15	3	20	HG01123_HAP1,HG01258_HAP2,HG01372_HAP1
```
## output-files

All outputs are saved in the current working directory with the base name MEI_cumulative_phases_comparison:
File	Format	Description	Use case
MEI_cumulative_phases_comparison.png	PNG	High-resolution raster (300 DPI)	Quick viewing, presentations
MEI_cumulative_phases_comparison.pdf	PDF	Vector format	Print publication
MEI_cumulative_phases_comparison.svg	SVG	Editable vector	Web, figure editing
MEI_cumulative_phases_comparison.html	HTML	Interactive report	Data exploration, sharing
Example workflow
Step 1: Prepare input files
```bash

# Create genomes list
cat > genomes.txt << EOF
HG01123
HG01258
HG01372
HG01489
HG01595
EOF

# Verify your MEI data format
head -n 5 your_mei_data.txt

# Check column count (should be 13 columns)
awk -F'\t' '{print NF; exit}' your_mei_data.txt
```
Step 2: Run analysis
```bash

python mei_cumulative_analysis.py -g genomes.txt -m your_mei_data.txt
```
Step 3: Monitor progress
```text

======================================================================
MEI CUMULATIVE DISCOVERY ANALYSIS
======================================================================

[1/5] Loading data...
  Priority genomes: genomes.txt
  -> 5 genomes loaded
  MEI data: your_mei_data.txt
  -> 15000 MEI records loaded

[2/5] Processing haplotype data...
  -> 250 unique haplotypes identified
  -> 250 haplotypes in final DataFrame

[3/5] Computing cumulative curves...
  Phase 2 (all haplotypes)...
    -> 250 haplotypes sorted
  Phase 1 (priority haplotypes only)...
    -> 10 priority haplotypes

[4/5] Generating figure...

[5/5] Saving outputs...
  Figures:
    MEI_cumulative_phases_comparison.png
    MEI_cumulative_phases_comparison.pdf
    MEI_cumulative_phases_comparison.svg
  HTML report:
    MEI_cumulative_phases_comparison.html
```
Step 4: View results
```bash

# On macOS
open MEI_cumulative_phases_comparison.html

# On Linux
xdg-open MEI_cumulative_phases_comparison.html

# On Windows
start MEI_cumulative_phases_comparison.html
```
## output-interpretation
Cumulative curves figure
```text
The generated figure displays:

X-axis: "Nth additional haplotype" - The order in which haplotypes are added
Y-axis: "Cumulative MEI discoveries" - Running total of unique MEIs found

Two visualization layers:

    Phase 2 (solid colors): Cumulative discovery using all haplotypes

        Dark red: Fixed variants (>95% frequency)

        Gold: Polymorphic variants (1-95% frequency)

        Deep blue: Rare variants (<1% frequency)

    Phase 1 (hatched patterns): Cumulative discovery using only prioritized haplotypes

        Light red with hatching: Fixed variants

        Light yellow with hatching: Polymorphic variants

        Light blue with hatching: Rare variants

Labels on the right:

    Category names (Fixed, Polymorphic, Rare) centered in their respective regions

    Numerical values showing total counts per category

    Separate annotations for Phase 1 and Phase 2 when applicable
```
## statistics-output

The console provides:
```text

======================================================================
STATISTICAL SUMMARY
======================================================================
Total haplotypes analyzed: 250
Priority haplotypes (Phase 1): 10
Priority genomes requested: 5

Phase 1 (Prioritized haplotypes only):
  Rare (<1%): 1250
  Polymorphic (1-95%): 3200
  Fixed (>95%): 450
  TOTAL: 4900

Phase 2 (All haplotypes):
  Rare (<1%): 3100
  Polymorphic (1-95%): 5200
  Fixed (>95%): 890
  TOTAL: 9190

HTML report sections

The interactive HTML report includes:

    Header section: Title, description, and generation timestamp

    Summary statistics cards: Key metrics for both phases

    Category tables: Detailed counts with color-coded badges

    Embedded figure: Interactive plot with download buttons

    Color legend: Visual guide for interpretation

    Methodology: Detailed analysis description

    Warnings: Missing genomes or data issues

    Footer: Author information and website link

Customization
Adjusting frequency thresholds
```
Edit these constants in the script:
```python

# In mei_cumulative_analysis.py
RARE_THRESHOLD = 0.01    # <1% (default)
FIXED_THRESHOLD = 0.95   # >95% (default)

# Example for stricter rare definition
RARE_THRESHOLD = 0.005   # <0.5%
FIXED_THRESHOLD = 0.99   # >99%
```
Changing figure dimensions
```python

# Nature single column: 89mm (3.5 inches)
FIGURE_WIDTH_INCHES = 3.5
FIGURE_HEIGHT_INCHES = 2.8

# Nature double column: 183mm (7.2 inches)
FIGURE_WIDTH_INCHES = 7.2
FIGURE_HEIGHT_INCHES = 4.8

# Full page width
FIGURE_WIDTH_INCHES = 10.0
FIGURE_HEIGHT_INCHES = 6.7
```
Modifying color scheme
```python

COLORS = {
    # Phase 2 (saturated)
    'phase2_fixed': '#8B0000',      # Dark red
    'phase2_polymorphic': '#E6B800', # Gold
    'phase2_rare': '#2E75B6',        # Deep blue
    
    # Phase 1 (pastel with hatching)
    'phase1_fixed': '#FFCCCC',       # Light red
    'phase1_polymorphic': '#FFF3CC', # Light yellow
    'phase1_rare': '#CCEEFF',        # Light blue
}

# Alternative colorblind-friendly scheme
COLORS = {
    'phase2_fixed': '#D55E00',       # Vermilion
    'phase2_polymorphic': '#009E73', # Blue-green
    'phase2_rare': '#0072B2',        # Blue
    'phase1_fixed': '#F0E442',       # Yellow (hatched)
    'phase1_polymorphic': '#56B4E9', # Sky blue (hatched)
    'phase1_rare': '#CC79A7',        # Pink (hatched)
}
```
Changing hatching pattern
```python

# Available hatch patterns
HATCH_PATTERN = '///'   # Diagonal lines (default)
# HATCH_PATTERN = 'xxx'  # Cross-hatch
# HATCH_PATTERN = '---'  # Horizontal lines
# HATCH_PATTERN = '|||'  # Vertical lines
# HATCH_PATTERN = '+++'  # Grid pattern
# HATCH_PATTERN = 'ooo'  # Dots
```
Modifying font settings
```python

# Change font family
FONT_FAMILY = 'Helvetica'     # or 'Arial', 'Times New Roman', 'DejaVu Sans'

# Adjust font sizes
FONT_SIZE_TITLE = 11
FONT_SIZE_AXES = 9
FONT_SIZE_LABEL = 9
FONT_SIZE_LEGEND = 8
FONT_SIZE_ANNOTATION = 8
```
Changing output filename
```python

# At the top of the script
OUTPUT_BASE = "my_custom_analysis_name"

# Or modify in the function call
save_figure(fig, output_base="my_custom_name")

Troubleshooting
Common errors and solutions
Error: "Missing required argument: --genomes"

Cause: One or both required input files not specified.
```
Solution: Always provide both input files:
```bash

python mei_cumulative_analysis.py -g genomes.txt -m mei_data.txt

Error: "File not found"

Cause: Specified file path does not exist.
```
Solution: Verify file paths:
```bash

# Check if files exist
ls -la genomes.txt
ls -la path/to/mei_data.txt

# Use absolute paths if needed
python mei_cumulative_analysis.py -g /full/path/genomes.txt -m /full/path/mei_data.txt

Warning: "requested genomes not found in data"

Cause: Genomes in your list don't appear in the MEI data.
```
Solution: Check available genomes:
```bash

# Extract unique genome IDs from MEI data
cut -f13 mei_data.txt | tr ',' '\n' | cut -d'_' -f1 | sort -u

# Compare with your genomes list
comm -23 <(sort genomes.txt) <(cut -f13 mei_data.txt | tr ',' '\n' | cut -d'_' -f1 | sort -u)

Empty plot or missing Phase 1

Cause: No priority haplotypes found with HAP1/HAP2 suffixes.
```
Solution: Verify haplotype format:
```bash

# Check samples column format
cut -f13 mei_data.txt | head -5

# Should show format like: HG01123_HAP1,HG01258_HAP2

Matplotlib font warnings

Cause: Arial or Helvetica fonts not available on your system.
```
Solution: Use default matplotlib fonts:
```python

# In mei_cumulative_analysis.py
FONT_FAMILY = 'DejaVu Sans'  # Default matplotlib font

Memory error with large datasets

Cause: Dataset too large for available RAM.
```
Solution: Process in chunks or increase memory:
```python

# Add chunking to load_mei_data function
def load_mei_data(filepath, chunksize=10000):
    chunks = []
    for chunk in pd.read_csv(filepath, sep="\t", header=None, 
                             names=columns, chunksize=chunksize):
        chunk["category"] = chunk["freq"].apply(categorize_frequency)
        chunks.append(chunk)
    return pd.concat(chunks, ignore_index=True)

Performance guidelines
Dataset size	Expected memory	Expected time	Recommendation
<10,000 MEIs	<500 MB	<1 min	Standard execution
10,000-50,000 MEIs	500 MB - 2 GB	1-5 min	Monitor memory usage
50,000-100,000 MEIs	2-4 GB	5-15 min	Use high-memory instance
>100,000 MEIs	>4 GB	>15 min	Consider downsampling
Validation checks
```
Run these checks before analysis:
```bash

# 1. Check genomes file format
file genomes.txt
wc -l genomes.txt

# 2. Check MEI data columns
head -1 mei_data.txt | awk -F'\t' '{print NF}'
# Should output 13

# 3. Check for missing values
cut -f13 mei_data.txt | grep -c "^$"

# 4. Validate haplotype format
cut -f13 mei_data.txt | grep -E "HG[0-9]{5}_(HAP[12]|[12])" | wc -l

```
## license-authors

This project is licensed under the MIT License - see the LICENSE file for details.
Authors : Anna-Sophie Fiston-Lavier, Capucine Mayoud, Shadi Yacoub
Acknowledgments : Human Pangenome Reference Consortium (HPRC) for MEI data standards and resources
