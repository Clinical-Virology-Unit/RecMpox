# RecMpox

RecMpox flags **recombination in mpox consensus genomes** by classifying each genome at **diagnostic SNPs** (and optionally indels) between two references (e.g. Clade Ia vs Ib, or IIa vs IIb). It aligns refs with [Squirrel](https://github.com/aineniamh/squirrel), finds diagnostic sites, then classifies each consensus at those positions. Outputs include a **TSV**, an **interactive HTML report** (sortable table, bar chart, diagnostic sites and recombination tracts per sample), and **all_sequences.fasta**. Potential recombinants are flagged when both refs contribute ≥5% of diagnostic sites; the report shows recombination tracts and breakpoints per genome.

## Installation

### 1. Create conda environment with required tools

```bash
conda create -n recmpox python=3.9 -y
conda activate recmpox
conda install -c bioconda -c conda-forge minimap2 samtools -y
conda install -c bioconda -c conda-forge squirrel -y
conda env config vars set PYTHONNOUSERSITE=1
conda deactivate && conda activate recmpox

```

### 2. Install RecMpox:

```bash
git clone https://github.com/DaanJansen94/RecMpox.git
cd RecMpox
pip install .
```

### 3. Re-installation (when updates are available):

```bash
conda activate RecMpox  # Make sure you're in the right environment
cd RecMpox
git pull  # Get the latest updates from GitHub
pip uninstall RecMpox
pip install .
```

## Usage

### Basic usage

```bash
# Use built-in references: Ia vs Ib or IIa vs IIb (sequences downloaded automatically)
recmpox -i fasta/ -o output -ref Ia,Ib -t 4
recmpox -i fasta/ -o output -ref IIa,IIb -t 4

# Input can be: FASTA file, directory of .fa/.fasta/.fna, or NCBI accession(s)
recmpox -i consensus.fa -o output -ref Ia,Ib
recmpox -i OZ375330.1,PX739443.1 -o output -ref IIa,IIb
recmpox -i accessions.txt -o output -ref Ia,Ib   # one accession per line or comma-separated
```

**Note**: Either `-ref` (e.g. `Ia,Ib` or `IIa,IIb`) or both `-ref1` and `-ref2` are required. With `-ref`, default references are used (Ia=OZ254474.1, Ib=PP601219.1, IIa=OZ287284.1, IIb=NC_063383.1).

### Command-line options

#### Required
- `-i, --input`: Input: FASTA file, directory of `.fa`/`.fasta`/`.fna`, `.txt` file of accessions (one per line or comma-separated), or NCBI accession(s). Accessions are downloaded and used as queries.

#### Reference (use one of)
- `-ref`: Reference pair: `Ia,Ib` or `IIa,IIb`. Uses built-in defaults.
- `-ref1`, `-ref2`: Custom references (path or NCBI accession). Use with `-ref1_g`/`-ref2_g` for labels (e.g. `-ref1_g Ia -ref2_g Ib`).

#### Optional
- `-o, --output`: Output directory (default: `output/`)
- `-ref1_g`, `-ref2_g`: Genotype labels for TSV/HTML (default from `-ref` or accession)
- `-include-indels`: Include diagnostic indels (default: SNPs only)
- `-min-indel-size`: Min indel length (bp) when using `-include-indels` (default: 100)
- `-t, --threads`: Number of threads
- `-q, --quiet`: Log to file only

### Examples

```bash
# Built-in Ia vs Ib
recmpox -i fasta/ -o output -ref Ia,Ib -t 4

# Built-in IIa vs IIb
recmpox -i fasta/ -o output -ref IIa,IIb -t 4

# Custom references
recmpox -i fasta/ -o output -ref1 NC_003310.1 -ref2 PP601219.1 -ref1_g Ia -ref2_g Ib -t 4

# Mixed clades (e.g. Ia vs IIb)
recmpox -i fasta/ -o output -ref1 ACC1 -ref2 ACC2 -ref1_g Ia -ref2_g IIb

# Include diagnostic indels
recmpox -i fasta/ -o output -ref Ia,Ib -include-indels
```

## Output files

- **recmpox_results.tsv**: Per-genome counts (n_ref1, n_ref2, n_other), percentages (pct_ref1, pct_ref2, pct_other), and recombinant call (no recombinant / potential recombinant).
- **recmpox_results.html**: Interactive report (summary, sortable table, stacked bar chart, diagnostic SNP positions, diagnostic sites per sample, recombination tracts and breakpoints per sample). Split into multiple files + index when >100 genomes.
- **all_sequences.fasta**: Ref1, ref2, and all query sequences (aligned).
- **potential_recombinants_diagnostic_sites.tsv**: Diagnostic site classification per potential recombinant (when any exist).
- **.recmpox.log**: Log file (in output directory).

Intermediate files (e.g. diagnostic_snps.txt, Squirrel outputs) are written under `work/`.

## Interpretation

- **No recombinant**: One ref dominates (minor ref &lt; 5% of diagnostic sites).
- **Potential recombinant**: Both refs contribute ≥5% (minor ref % ≥ 5%). The HTML report shows recombination tracts (beginning/end of each tract) and breakpoints between tracts; a single tract means the genome is entirely one clade (no recombination).
- **High pct_other**: Many Ns, gaps, or non-ref bases at diagnostic sites (poor coverage or alignment).

## Citation

If you use RecMpox in your research, please cite:

```
Jansen, D., Laumen, J., Siebenmann, E., & Vercauteren, K. (2025). LassaSeq: A Command-Line Tool for Downloading, Processing and Analyzing Lassa Virus Sequences for Phylogenetic Analysis (v0.1.2). Zenodo. https://doi.org/10.5281/zenodo.14936276
```

## License
This project is licensed under the GNU General Public License v3.0 (GPL-3.0) - see the LICENSE file for details.

## Contributing
Contributions are welcome! Please feel free to submit a Pull Request.

## Support
If you encounter any problems or have questions, please open an issue on GitHub.
