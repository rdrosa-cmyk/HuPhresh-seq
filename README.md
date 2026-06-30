# Human Peptidome PhIP-seq Library (JP)

## Overview
This pipeline generates a genome-wide human PhIP-seq peptide library and converts it
into synthesis-ready DNA oligos. It tiles two protein sources plus a control:

- **GENCODE v50** protein-coding translations (the canonical human peptidome)
- **Ribo-seq dark-proteome ncORFs** (non-canonical / dark-proteome ORFs)
- **GFP (P42212)** as a negative control

Design follows the i-PanSeq library approach
(https://www.protocols.io/view/i-panseq-library-design-rm7vz9od4gx1/v1).

## Folders
| Path | Contents |
|------|----------|
| `peptidome_proteomes/` | Raw input proteomes (GENCODE FASTA + GTF, Ribo-seq ORFs) |
| `peptidome_library_construction/GENCODE/` | GENCODE intermediates, tiles, CD-HIT outputs |
| `peptidome_library_construction/Riboseq/` | Ribo-seq ncORF intermediates, tiles, CD-HIT outputs |
| `peptidome_library_construction/NegativeControls/` | GFP control intermediates and tiles |
| `peptidome_library_construction/` | Final merged peptides and encoded oligos |
| `JP_peptidome.ipynb` | The full pipeline (run cells top to bottom) |

## Input Data
| Source | File | Notes |
|--------|------|-------|
| GENCODE peptidome | `1_peptidome_gencode.v50.pc_translations.fa.gz` | Filtered to `transcript_type == protein_coding` using the GTF |
| GENCODE annotation | `1_gencode.v50.annotation.gtf.gz` | Used to select protein_coding transcript IDs |
| Ribo-seq ncORFs | `2_darkproteome_Ribo-seq_ORFs.primary.faa` | Dark-proteome / non-canonical ORFs |
| Negative control | GFP `P42212` (UniProt) | Non-human, no expected human reactivity |

## Key Design Choice: long vs short proteins
Each source is split into two streams **before** tiling and clustering:

- **Long (≥64 aa):** tiled normally; deduplicated with similarity-based CD-HIT.
- **Short (5–63 aa):** N-terminally padded to exactly 64 aa with an **XTEN linker**,
  then deduplicated only at 100% full-length identity.

**Rationale:** short proteins share the same XTEN padding sequence. If padded
peptides were clustered with the standard 95% local settings, unrelated short
proteins would collapse together merely because they share the XTEN region.
Keeping them separate and deduplicating only at exact full-length identity
prevents this artifact. Proteins shorter than 5 aa are removed.

XTEN linker (from the Bodansky AEGIS paper, via John):
`ATPESSGSETPGTSESATPESSGSEAPGTSTEPSEGSAPGSEPATSGSETPGTSESATP`

## Pipeline Steps

### 1. Process GENCODE translations (Cell 0)
Parse the GTF for `protein_coding` transcript IDs and keep only matching records
from the translation FASTA.

### 2. GENCODE cleanup + split (Cell 1)
Strip `*`, remove <5 aa, split into long-unpadded (≥64 aa) and short-XTEN-padded
(5–63 aa → 64 aa). Cell 2 plots the sub-64 aa length distribution.

### 3. GENCODE tiling + CD-HIT (Cell 3)
Using `pepsyn`:
- `x2ggsg` — replace X runs with a length-preserving GGSG linker
- `tile -l 64 -p 42` — 64-aa tiles, 42-aa overlap (22-aa step)
- `ctermpep -l 64` — C-terminal peptides generated separately to preserve protein ends
- `disambiguateaa` — resolve ambiguous residues (B, J, Z, X)

CD-HIT settings:

| Stream | Identity | Coverage |
|--------|----------|----------|
| Long tiles | 95% (`-c 0.95 -G 0`) | `-A 53` (local, ≥53 aa) |
| Long C-terminal | 95% (`-c 0.95 -G 0`) | full length (`-aL 1.0 -aS 1.0`) |
| Short XTEN-padded | 100% (`-c 1.0 -G 1`) | full length (`-aL 1.0 -aS 1.0`) |

The 53-aa alignment threshold prevents consecutive tiles (42-aa overlap) from
collapsing into each other. C-terminal peptides use full-length coverage instead
(following Mohan et al. 2018 / i-PanSeq) so the protein's terminal residues — and
any C-terminal epitopes — are never lost to a partial-overlap collapse.

### 4. Record CD-HIT collapses (Cell 4)
Parse `.clstr` files; mark representative headers that absorbed members with `#`
(plus `collapsed_n` and `cluster` metadata) and write a full cluster-membership TSV.

### 5–9. Ribo-seq ncORFs (Cells 5–9)
Same cleanup → split → tiling → CD-HIT → collapse-record workflow as GENCODE,
applied to the dark-proteome ORFs.

### 10–11. GFP control (Cells 10–11)
Copy and tile GFP (238 aa, naturally long); deduplicate only at 100% full-length
identity so control representation is preserved.

### 12. Merge all peptide sets (Cell 12)
Combine all eight sub-libraries into one amino-acid FASTA. Each record gets a
compact, **whitespace-free unique ID**:

```
<source>|<protein>|<specific_id>|<tile>
```
e.g. `GENCODE_long_w64s42|OR4F5|ENSP00000493376.2|0-64`

This is required because the downstream `pepsyn` encoding truncates FASTA headers
at the first whitespace; putting a space-free ID first keeps a real per-peptide
identifier on every oligo. The cell validates header shapes and **raises on any
duplicate ID** (read-to-peptide mapping must be 1:1). Where available it prefers
the `marked_representatives` FASTAs so CD-HIT cluster metadata is preserved.

### 13. Oligo generation (Cell 13)
Convert peptides to DNA with `pepsyn`:

| Step | Tool | Description |
|------|------|-------------|
| 1 | `stripstop` | Remove stop codons |
| 2 | `pad -l 64` | Pad inserts to 64 aa (not relevant for this set but kept for future compatibility) |
| 3 | `disambiguateaa` | Resolve ambiguous residues |
| 4 | `revtrans --codon-freq-threshold 0.01` | Reverse translate using E. coli codons (>1% frequency) |
| 5 | `prefix` | Add 5' adapter (19 bp) |
| 6 | `suffix` | Add 3' FLAG tag (24 bp) + adapter (19 bp) |
| 7 | `recodesite` | Remove EcoRI / HindIII / XhoI sites from the coding region |

### 14. GC / homopolymer repair — DNAChisel (Cell 14)
`pepsyn` already produces E. coli codon-optimized oligos, but a small minority fall
outside the synthesis-friendly and PCR window. This cell reads the **raw** pepsyn
oligos, re-encodes **only flagged** oligos, and writes a separate **post-repair**
FASTA (`encoded_oligos.fasta`), leaving the raw file (`raw_encoded_oligos.fasta`)
intact for auditing/diffing. It writes the full library and reports how many
flagged oligos (if any) could not be brought within the parameters.

An oligo is flagged if: global **GC < 35% or > 65%**, OR any **homopolymer run ≥ 10**,
OR it contains an EcoRI/HindIII/XhoI site.

Repair (`DnaOptimizationProblem`):
- `EnforceTranslation` over the 192-bp coding region (peptide unchanged)
- `AvoidChanges` locks the 5' adapter and 3' FLAG+adapter
- `EnforceGCContent(0.35, 0.65)` — global GC as a **hard band constraint** (never minimized)
- `AvoidPattern` — no homopolymer run ≥ 10; no EcoRI/HindIII/XhoI sites
- objective `CodonOptimize(e_coli, match_codon_usage)` — keeps codons E. coli-natural

GC is enforced as a band rather than minimized so that codons are nudged only
enough to enter the window, preserving E. coli expression.

### 15–16. GC analysis (Cells 15–16)
Overall oligo GC distribution (Cell 15) and GC split by source — GENCODE vs
Ribo-seq vs GFP (Cell 16).

### 17. Quality control (Cell 17)
On the final encoded oligos:
- **ID uniqueness** — every oligo maps back to exactly one peptide
- **Stop codons** — none in the coding region
- **Restriction sites** — no EcoRI/HindIII/XhoI anywhere in the oligo
- **Length** — every oligo is 254 bp
- **Flanking sequences untouched** — 5' adapter and 3' FLAG+adapter match exactly
- **GC band** — every oligo is within 35–65% GC
- **Homopolymer** — no single-base run ≥ 10 bp
- **Record count** — oligo count matches the merged peptide library (nothing lost or added)

## Oligo Structure
```
5'─[Adapter 19 bp]─[Peptide 192 bp]─[FLAG 24 bp]─[Adapter 19 bp]─3'
   CCTATACTTCCAAGGCGCA   (variable)   GACTACAAAGATGACGACGATAAA  GGTGACTCTCTGTCTTGGC
```
Total oligo length: **254 bp** (64-aa insert × 3 + 62 bp fixed flanks).

## Output Files
| File | Description |
|------|-------------|
| `peptidome_library_construction/human_peptidome_JP_library.merged_peptides.fasta` | Merged peptides (all sources) with unique IDs |
| `peptidome_library_construction/human_peptidome_JP_library.raw_encoded_oligos.fasta` | Raw pepsyn DNA oligos (pre-repair) |
| `peptidome_library_construction/human_peptidome_JP_library.encoded_oligos.fasta` | Final synthesis-ready DNA oligos (post DNAChisel repair) |
| `GENCODE/`, `Riboseq/`, `NegativeControls/` | Per-source tiles, CD-HIT outputs, and `cdhit_collapse_records/` (marked FASTAs + cluster TSVs) |

## Tools Required
- **pepsyn** — peptide tiling and reverse translation
- **cd-hit** — sequence clustering
- **DNAChisel** — GC / homopolymer repair (Cell 14)
- **Python 3** with `matplotlib`, `numpy`

Create the environment from the included spec:
```bash
conda env create -f environment.yml
conda activate peptidome
```
The bash cells call `pepsyn` and `cd-hit` directly, so run the notebook from a
kernel where those are on `PATH`.

## Repository Contents
Only the pipeline and docs are version-controlled. The data files (multi-GB,
fully regenerable) are excluded via `.gitignore`:
- `JP_peptidome.ipynb` — the full pipeline
- `README.md` — this file
- `environment.yml` — conda environment (pepsyn, cd-hit, dnachisel, etc.)
- `.gitignore` — excludes `peptidome_proteomes/` and `peptidome_library_construction/`

To reproduce: create the environment, download the raw inputs (GENCODE v50
translations + GTF, Ribo-seq ORFs, UniProt GFP P42212) into `peptidome_proteomes/`,
then run the notebook.

## Usage
Run cells 0–17 sequentially in `JP_peptidome.ipynb`. If the XTEN linker, adapters,
or tiling parameters change, rerun the whole pipeline from Cell 0.
