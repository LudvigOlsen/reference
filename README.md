# *reference* Kmer counter

Ultra‑fast **k‑mer counter** for any reference genome. Count globally or in specified windows with blacklist filtering.

Written in Rust for speed and parallelism. Outputs standard NumPy/SciPy formats for easy use in python.

Counting the occurence of reference k‑mers inside thousands of genomic windows is a common step in many bioinformatical pipelines. E.g. useful as features in machine learning models or for normalization and bias correction.

**! WORK-IN-PROGRESS !**

---

## Key features

| capability                   | details                                                                                                     |
| ---------------------------- | ----------------------------------------------------------------------------------------------------------- |
| **Any k (1‑27)**             | pass one or more values with `-k 3,5,11`                                                                    |
| **Multiple windowing modes** | fixed length (`--by-size 10_000`), BED intervals (`--by-bed sites.bed`), or single genome‑wide (`--global`) |
| **Blacklist masking**        | exclude repeats/artefacts with one or several BEDs                                                          |
| **Canonical kmers**          | merge reverse complements with `-c`                                                                         |
| **Dense *or* sparse output** | dense `.npy` for small k or SciPy‑compatible COO `.npz` for large k (`--save-sparse`)                       |
| **Multithreaded**            | set `-t <N>` to fill N cores                                                                                |
| **Runs on 2‑bit file**       | use 2bit reference file (e.g., `hg38.2bit`)                                                                 |

---

## Installation

```bash
$ cargo install --git https://github.com/ludvigolsen/reference
# or clone + build
$ git clone https://github.com/ludvigolsen/reference
$ cd reference && cargo build --release
```

---

## Quick‑start example

The below example ran in ~25s (or <14s with 12 threads):

```bash
reference \
  --ref-2bit hg38.2bit \                  # reference genome in 2‑bit
  --output-dir results \                  # where to write files
  --kmer-sizes 3,5 \                      # count 3-mers and 5‑mers
  --by-size 1000000 \                     # tiling 1Mb windows
  --n-threads 4 \                         # use 4 CPU cores (max. one per chromosome)
  --blacklist encode_blacklist.bed        # mask ENCODE blacklist

// Optional:
  --save-sparse                           # write SciPy COO for large k‑mers
  --canonical                             # collapse to lexicographically lowest canons
```

**Output explained**

(Shown for k=3)

| file                    | description                                       |
| ----------------------- | ------------------------------------------------- |
| `results/k3_motifs.txt` | one motif per column in the count array           |
| `results/k3_counts.npy` | dense matrix `[windows × 4^3]` of `uint64` counts |
| `results/bins.bed`      | coordinates of every window + % blacklist overlap |

for sparse arrays (instead of `*_counts.npy`):

| file                           | description                                                 |
| ------------------------------ | ----------------------------------------------------------- |
| `results/k3_counts_sparse.npz` | SciPy COO archive (`data`, `row`, `col`, `shape`, `format`) |

---

## Working with the output in Python

Dense outputs are plain `.npy` arrays:

```python
import numpy as np
counts_3 = np.load("results/k3_counts.npy")                  # shape (N_windows, N_kmers)
motifs_3 = np.loadtxt("results/k3_motifs.txt", dtype="U")    # len == N_kmers
```

Sparse outputs can be loaded with SciPy:

```python
import scipy.sparse
import numpy as np

coo_3 = scipy.sparse.load_npz("results/k3_counts_sparse.npz")  # (N_windows, N_kmers)

# Convert once for fast row slicing
csr_3 = coo_3.tocsr()

# Iterate over windows
for i in range(csr_3.shape[0]):
    sparse_row = csr.getrow(i)
    dense_row = sparse_row.toarray().ravel()
    # downstream analysis here ...
```

---


## Command‑line reference

```bash
reference --help
```

| option                      | purpose                                                 |
| --------------------------- | ------------------------------------------------------- |
| `-r`, `--ref-2bit <path>`   | two‑bit reference genome                                |
| `-k`, `--kmer-sizes <list>` | k values (1–27)                                         |
| `-c`, `--canonical`         | merge reverse complements (to lexicographically lowest) |
| **Window selection**        |                                                         |
| `--by-size <bp>`            | fixed‑length windows                                    |
| `--by-bed <BED>`            | custom intervals                                        |
| `--global`                  | one big window per chromosome set                       |
| **Filtering**               |                                                         |
| `-b`, `--blacklist <BED>`   | mask repeats/artefacts                                  |
| `--blacklist-min-size <bp>` | drop tiny blacklist entries                             |
| **Chromosome selection**    |                                                         |
| `--chromosomes <list>`      | chromosomes to process (default: chr1-22)               |
| `--chromosomes-file <path>` | file with chromosomes to process                        |
| **Output**                  |                                                         |
| `--save-sparse`             | write SciPy‑loadable COO                                |
| **Performance**             |                                                         |
| `-t`, `--n-threads <N>`     | CPU threads                                             |


