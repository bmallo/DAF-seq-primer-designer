# DAF-seq Primer Designer

A Python pipeline for designing PCR primers optimized for DAF-seq treated DNA. This tool addresses the unique challenge that cytosine deamination creates: converted cytosines become uracils (and subsequently thymines after PCR), reducing primer binding efficiency when primers contain guanines.

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
- [Configuration](#configuration)
- [Output](#output)
- [Algorithm Details](#algorithm-details)
- [Troubleshooting](#troubleshooting)
- [Advanced Usage](#advanced-usage)

---

## Features

- **DAF-seq primer design**: Minimizes guanines in primers and penalizes G positions near the 3' end
- **Flexible region specification**: Define target regions by coordinates or DNA sequence
- **Comprehensive primer evaluation**: 
  - Tm calculation (Primer3-compatible)
  - Secondary structure checking (hairpins, self-dimers, primer-dimers)
  - Specificity checking (no off-target binding within provided sequence region)
- **Diversity enforcement**: Prevents overlapping primers in final results
- **Customizable scoring**: Adjustable penalties for all primer properties
- **Silent or verbose mode**: Run quietly for batch processing or with detailed output

---

## Installation

### Requirements

- Python 3.7+
- Required packages:
  ```bash
  pip install biopython pyyaml
  ```

### Setup

1. Clone or download all Python files to a directory:
   ```
   primer_designer/
   ├── config.py
   ├── sequence_parser.py
   ├── primer_generator.py
   ├── primer_evaluator.py
   ├── primer_pairer.py
   ├── primer_specificity.py
   └── primer_designer.py
   ```

2. Ensure all files are in the same directory

3. Test installation:
   ```bash
   python primer_designer.py --create-example-config
   ```

---

## Quick Start

### 1. Generate a configuration file

```bash
python primer_designer.py --create-example-config
```

This creates `primer_design_config.yaml` with default settings.

### 2. Edit the configuration

Open `primer_design_config.yaml` and modify:

```yaml
sequence_file: "path/to/your/template.fasta"
target_region:
  core_sequence: "ATCGATCGATCG"  # Or use core_start/core_end
  min_amplicon: 2000
  max_amplicon: 7000
output_file: "my_primers.txt"
output_dir: "results"
```

### 3. Run primer design

```bash
# Silent mode (default)
python primer_designer.py primer_design_config.yaml

# Verbose mode
python primer_designer.py primer_design_config.yaml -v
```

### 4. View results

Results are saved to `results/my_primers.txt` with tab-delimited primer information.

---

## Usage

### Command Line Options

```bash
python primer_designer.py [config_file] [options]

Arguments:
  config_file              YAML configuration file

Options:
  -o, --output FILE       Specify output file (overrides config)
  -v, --verbose           Print detailed progress
  -q, --quiet             Force silent mode (overrides config)
  --create-example-config Generate example configuration file
  -h, --help              Show help message
```

### Examples

```bash
# Basic usage
python primer_designer.py config.yaml

# Specify output file
python primer_designer.py config.yaml -o custom_results.txt

# Verbose output
python primer_designer.py config.yaml -v

# Force silent mode
python primer_designer.py config.yaml -q
```

---

## Configuration

### Complete Configuration Reference

```yaml
# Input sequence
sequence_file: "template.fasta"

# Target region definition
target_region:
  # Option 1: Define by coordinates (0-based, will be converted to 1-based in output)
  core_start: 1000
  core_end: 3000
  
  # Option 2: Define by sequence (automatically located)
  # core_sequence: "ATCGATCGATCG"
  
  # Amplicon size constraints
  min_amplicon: 2000
  max_amplicon: 7000
  
  # Optional: Restrict primer search regions
  # Can use coordinates OR sequences for each
  # forward_region: [500, 1100]
  # reverse_region: [2900, 3500]
  # forward_region_sequence: "GATTACAGAT"
  # reverse_region_sequence: "CTAGCTAGCT"

# Primer constraints
primer_constraints:
  # Length
  min_length: 18
  max_length: 25
  
  # Guanine restrictions (key for DAF-seq primers!)
  max_guanines: 2
  no_gc_3prime_bases: 3  # No G or C in last 3 bases
  
  # Tm range
  tm_min: 58.0
  tm_max: 62.0
  
  # GC content
  min_gc_content: 20.0
  max_gc_content: 60.0
  
  # Other constraints
  max_poly_base: 4  # Max consecutive identical bases
  
  # Specificity checking
  check_specificity: true
  allow_multiple_exact_matches: false
  specificity_min_match_length: 15
  specificity_allow_mismatches: 2
  
  # Tm calculation parameters (salt concentrations)
  Na: 50.0    # mM
  K: 0.0      # mM
  Mg: 0.0     # mM
  dNTPs: 0.0  # mM
  primer_conc: 250.0  # nM

# Pairing constraints
pairing_constraints:
  max_tm_difference: 2.0  # Max Tm difference between primers
  max_primer_reuse: 3     # Max times a primer can appear in results
  min_primer_spacing: 50  # Min distance between primers (for diversity)

# Scoring weights (higher = more penalty)
scoring_weights:
  guanine_penalty: 10.0
  guanine_3prime_penalty: 50.0  # Extra penalty for G near 3' end
  tm_deviation_penalty: 1.0
  gc_deviation_penalty: 0.5
  length_deviation_penalty: 0.5
  hairpin_penalty: 5.0
  self_dimer_penalty: 5.0
  poly_base_penalty: 3.0

# Output settings
num_primer_sets: 10
output_file: "my_primers.txt"
output_dir: "results"
verbose: false  # Set to true for detailed output
```

### Key Parameters Explained

#### Guanine Restrictions
- **`max_guanines`**: Maximum G's allowed per primer (default: 2)
  - Fewer G's = better for DAF-seq-treated DNA
  - 0 G's is ideal but may be too restrictive
  
- **`no_gc_3prime_bases`**: Number of bases at 3' end that cannot contain G or C
  - Critical for PCR efficiency
  - Default: 3 bases

#### Specificity Checking
- **`check_specificity`**: Enable/disable off-target checking
  - `true`: Reject primers with multiple binding sites (recommended)
  - `false`: Skip specificity check (faster, less stringent)

- **`min_primer_spacing`**: Minimum distance between primers in final set
  - Ensures diversity in primer positions
  - Increase for more distinct primer locations
  - Decrease if you need more primer options

#### Tm Calculation
- Salt concentrations affect Tm calculations
- Default values (Na=50mM) are typical for PCR
- Adjust based on your reaction conditions

---

## Output

### Output File Format

Results are saved as a tab-delimited text file with the following columns:

```
Rank                    Primer set ranking (1 = best)
Pair_Score              Composite score (lower = better)
Amplicon_Size           PCR product size (bp)
Tm_Difference           ΔTm between forward and reverse primers (°C)
Forward_Sequence        5'→3' primer sequence
Forward_Position        1-based genomic position
Forward_Length          Primer length (bp)
Forward_Tm              Melting temperature (°C)
Forward_GC%             GC content (%)
Forward_Num_Gs          Number of guanines
Forward_G_Positions     G positions from 3' end
Forward_Score           Individual primer score
Reverse_Sequence        5'→3' primer sequence
Reverse_Position        1-based genomic position
Reverse_Length          Primer length (bp)
Reverse_Tm              Melting temperature (°C)
Reverse_GC%             GC content (%)
Reverse_Num_Gs          Number of guanines
Reverse_G_Positions     G positions from 3' end
Reverse_Score           Individual primer score
Primer_Dimer            Risk of primer-dimer formation (Yes/No)
Covers_Core             Amplicon covers entire core region (Yes/No)
```

### Reading Results

- **All positions are 1-based** (compatible with genome browsers and analysis tools)
- Primers are ranked by score (lower = better)
- Focus on primer sets with:
  - Low Pair_Score
  - Few or no guanines
  - Tm_Difference < 2°C
  - No primer dimers

### Example Output

```
Rank	Pair_Score	Amplicon_Size	Tm_Difference	Forward_Sequence	Forward_Position	...
1	15.23	2847	0.45	ATATATATATATATATATAT	1245	...
2	18.67	3102	1.23	TATATATATATATAATATAT	1389	...
3	21.04	2654	0.89	AATATATATATATATAATAT	1156	...
```

---

## Algorithm Details

### Primer Design Pipeline

1. **Sequence Parsing**
   - Load FASTA file
   - Locate core region (by coordinates or sequence)
   - Determine primer search regions

2. **Primer Generation**
   - Generate all k-mers (18-25 bp) using sliding window
   - Pre-filter by G content and 3' restrictions
   - Calculate basic properties (GC%, poly-bases)

3. **Primer Evaluation**
   - Calculate Tm (Primer3 nearest-neighbor method)
   - Check secondary structures (hairpins, self-dimers)
   - Score primers based on multiple criteria

4. **Specificity Checking**
   - Search for off-target binding sites
   - Check 3' end specificity (most critical)
   - Reject non-specific primers

5. **Primer Pairing**
   - Create all valid forward/reverse combinations
   - Filter by amplicon size and Tm difference
   - Check primer-dimer formation
   - Score pairs

6. **Diversity Selection**
   - Select top N pairs
   - Enforce minimum spacing between primers
   - Limit primer reuse

### Scoring System

Primers are scored with penalties for deviations from optimal:

- **Guanine penalties**: Position-weighted (3' end = highest penalty)
- **Tm deviation**: From optimal Tm (midpoint of range)
- **GC deviation**: From optimal GC% (midpoint of range)
- **Length deviation**: From optimal length (midpoint of range)
- **Secondary structures**: Hairpins and self-dimers
- **Primer-dimer risk**: For pairs

**Lower scores are better** (0 = perfect primer, though unattainable)

---

## Troubleshooting

### No primers found

**Symptoms:** "No primers found with current constraints"

**Solutions:**
1. Increase `max_guanines` (try 3 instead of 2)
2. Widen Tm range (`tm_min: 55`, `tm_max: 65`)
3. Expand primer length range
4. Check that your target regions make sense

### No specific primers found

**Symptoms:** "No specific primers found!"

**Solutions:**
1. Disable specificity checking: `check_specificity: false`
2. Allow multiple matches: `allow_multiple_exact_matches: true`
3. Increase `specificity_allow_mismatches: 3`
4. Check for repetitive sequences in your template

### Fewer than N primer sets returned

**Symptoms:** "Only X diverse primer pairs found"

**Solutions:**
1. Reduce `min_primer_spacing` (try 25 instead of 50)
2. Increase `max_primer_reuse` (try 5 instead of 3)
3. Relax other constraints to generate more candidate primers

### All primers have guanines at 3' end

**Symptoms:** No primers pass filters

**Solutions:**
- Your target region may have high GC content
- Try increasing `no_gc_3prime_bases` search range
- Consider a different target region

### Tm calculation errors

**Symptoms:** Error during Tm calculation

**Solutions:**
- Check that BioPython is installed: `pip install biopython`
- Verify primer sequences are valid DNA (ATCG only)
- Check salt concentration parameters are reasonable

---

## Advanced Usage

### Batch Processing

Process multiple regions with a shell script:

```bash
#!/bin/bash
for region in region1 region2 region3; do
    # Modify config for each region
    sed "s/REGION_NAME/$region/g" template_config.yaml > ${region}_config.yaml
    
    # Run primer design
    python primer_designer.py ${region}_config.yaml -q
done
```

### Custom Scoring

Adjust scoring weights for your specific needs:

```yaml
scoring_weights:
  guanine_penalty: 20.0        # Strongly avoid Gs
  guanine_3prime_penalty: 100.0  # Extremely avoid 3' Gs
  tm_deviation_penalty: 0.5    # Less strict on Tm
```

### Region Definition Strategies

**Strategy 1: Precise coordinates**
```yaml
target_region:
  core_start: 1000
  core_end: 3000
  forward_region: [500, 1100]
  reverse_region: [2900, 3500]
```

**Strategy 2: Sequence-based (flexible)**
```yaml
target_region:
  core_sequence: "ATCGATCGATCG"
  # Automatically finds and flanks with primers
```

**Strategy 3: Mixed approach**
```yaml
target_region:
  core_start: 1000
  core_end: 3000
  forward_region_sequence: "GATTACA"  # Let it find the sequence
  reverse_region: [2900, 3500]         # Use exact coordinates
```

### Working with Large Genomes

For large templates (>100 kb):

1. Extract your region of interest first:
   ```bash
   samtools faidx genome.fa chr1:1000000-1100000 > region.fasta
   ```

2. Update coordinates in config to match extracted region

3. Run primer designer on smaller region

### Interpreting Scores

**Excellent primers** (score < 20):
- 0-1 guanines
- Optimal Tm
- No secondary structures
- Perfect specificity

**Good primers** (score 20-40):
- 1-2 guanines
- Tm within range
- Minor structural issues
- Specific

**Acceptable primers** (score 40-60):
- 2 guanines
- Tm at range edges
- Some structural concerns
- Still specific

**Problematic primers** (score > 60):
- Consider relaxing constraints or choosing different region

---

## File Descriptions

- **`config.py`**: Configuration parsing and data structures
- **`sequence_parser.py`**: FASTA reading and region extraction
- **`primer_generator.py`**: K-mer generation and basic filtering
- **`primer_evaluator.py`**: Tm calculation, structure checking, scoring
- **`primer_pairer.py`**: Primer pairing and compatibility checking
- **`primer_specificity.py`**: Off-target detection and specificity filtering
- **`primer_designer.py`**: Main pipeline script

---

## Citation

If you use this tool in your research, please cite:

```
[Your citation information here]
```

---

## Support

For issues, questions, or feature requests:
- Check the [Troubleshooting](#troubleshooting) section
- Review example configurations
- Ensure all dependencies are installed

---

## License

[Your license information here]

---

## Changelog

### Version 1.0
- Initial release
- Core primer design functionality
- DAF-seq-aware scoring
- Specificity checking
- Diversity enforcement

---

## Acknowledgments

- BioPython for sequence analysis tools
- Primer3 team for Tm calculation methods
- [Any other acknowledgments]
