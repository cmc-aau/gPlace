# gPlace

**gPlace** places genomes on the GTDB reference trees **exactly as GTDB-Tk's pplacer does**, bit for bit, while using a fraction of the memory and time. It is a re-implementation of `pplacer -c <refpkg>` (pplacer v1.1.alpha19, the build GTDB-Tk ships): the WAG + FastTree CAT reference packages used by GTDB-Tk, and every other model pplacer supports (LG, JTT, GTR for DNA; CAT or gamma rates; model or empirical frequencies).

* **Identical output:** same edges, likelihoods, like-weight ratios and branch lengths. Every number in the jplace file is textually identical to pplacer's (`%.12g`), as is the reference tree string.
* **Bounded memory:** under 44 GiB for the largest GTDB-Tk R226 class-level tree, which needs ~104 GB with pplacer.
* **Fast:** all 127,343 GlobDB MAGs on the 8 R226 bacterial class-level trees placed in 28 minutes on one server (peak 43.8 GiB).

## Features

* **Drop-in placement:** `gplace place` takes the same inputs as pplacer (aligned queries + refpkg) and writes the same jplace.
* **Pooled batches:** `gplace place_batch` places many alignments at once; alignments sharing a reference package are placed in a single pass.
* **Memory limits and huge trees:** `--max_memory 64G` keeps a run within a memory limit. gPlace places without tiles when that fits, and otherwise tile by tile, with the largest tiles that fit, and identical results. This places genomes on the full GTDB bacterial tree (273,290 nodes), which would need ~490 GiB untiled. `--tile_edges N` sets the tile size directly.
* **Every pplacer model:** WAG, LG, JTT (protein) and GTR (DNA), with FastTree CAT or gamma rate variation and model or empirical frequencies, read from the reference package as pplacer does. `--fix_gamma_underflow` optionally avoids pplacer's gamma underflow on large trees (likelihood `-inf`).
* **Self-verifying installation:** `gplace check_install` checks the compiled engine bit for bit against values taken from GTDB-Tk's pplacer binary; `gplace test` places bundled genomes and compares them with pplacer's own output.

## Documentation

Documentation for gPlace can be found [here](https://cmc-aau.github.io/gPlace/).

## Installation

gPlace needs Python >= 3.12 and an x86-64 CPU with AVX2 and FMA (any recent Intel/AMD server).

```bash
pip install gtdb-gplace
```
On x86-64 Linux this installs a prebuilt wheel: no compiler needed, one wheel for every Linux
distribution with glibc >= 2.17 (CentOS 7 and newer) and every Python >= 3.12. The wheel is
built in the manylinux2014 container and only published after it passes `gplace check_install`
and `gplace test` there and on several distributions.

From source (other platforms, or development), which needs a C compiler (gcc) and `make`:
```bash
git clone https://github.com/cmc-aau/gPlace.git
cd gPlace
pip install .          # also compiles the placement engine (builds GSL 1.16 once, ~1-3 minutes)

# verify
gplace check_install
gplace test
```

When installing from source, the engine is compiled during installation. Optional environment variables for the build:

| Variable | Use |
|---|---|
| `GPLACE_GSL_TARBALL` | path to `gsl-1.16.tar.gz` for offline installs (otherwise downloaded from ftp.gnu.org and SHA-256 checked) |
| `GPLACE_GSL_PREFIX` | reuse an existing static, `-fPIC` GSL 1.16 build |
| `GPLACE_BUILD_DIR` | scratch directory for the GSL build (default `~/.cache/gplace/build`, reused between installs) |
| `GPLACE_BUILD_JOBS` | parallel `make` jobs (default 8) |

`gplace build_engine` rebuilds the engine without reinstalling (advanced).

## Usage

You can view the general help menu at any time:
```bash
gplace -h
```

### Place one alignment
```bash
gplace place --msa_file user_msa_file.fasta \
             --refpkg gtdbtk.package.2.refpkg \
             --out_file pplacer.class_level.bac120.json \
             --cpus 40
```
This is equivalent to GTDB-Tk's `pplacer -m WAG -j 40 -c gtdbtk.package.2.refpkg -o pplacer.class_level.bac120.json user_msa_file.fasta`.

### Place many alignments
```bash
gplace place_batch --batchfile batch.tsv --out_dir gplace_out --cpus 40
```
The file provided to `--batchfile` is a three-column, tab-separated file: query alignment, reference package, output jplace.
```
/runs/part00/tree_2/user_msa_file.fasta    /gtdbtk/split/class_level/pplacer/gtdbtk.package.2.refpkg    /runs/part00/tree_2/pplacer.class_level.bac120.json
/runs/part01/tree_2/user_msa_file.fasta    /gtdbtk/split/class_level/pplacer/gtdbtk.package.2.refpkg    /runs/part01/tree_2/pplacer.class_level.bac120.json
...
```
All rows that use the same reference package are placed together, so the reference likelihoods are computed once. `scripts/make_gtdbtk_batchfile.py` writes this file for a directory of GTDB-Tk class-level alignments.

### Other models and DNA
The model comes from the reference package's `phylo_model` JSON, as in pplacer (`-m` is ignored for reference packages): `subs_model` WAG, LG, JTT or GTR, `ras_model` Price-CAT or gamma, `empirical_frequencies`. Reference packages made with `taxit create` from FastTree, RAxML or PhyML runs work unchanged; DNA queries use pplacer's nucleotide alphabet (IUPAC codes included).

pplacer's gamma model fails on large trees: once the slow rates of a fast-evolving site underflow to zero, pplacer stops rescaling that site and the likelihood becomes `-inf` (on the GTDB-Tk class-level tree 8 with gamma, for 342 of 343 genomes). gPlace reproduces this by default and warns; `--fix_gamma_underflow` fixes it, at the cost of no longer matching pplacer (its results elsewhere change only in the last digits). See the documentation's models page.

### Compare jplace files
```bash
gplace compare --jplace_a gplace.jplace --jplace_b pplacer.jplace
```

### Choosing `--cpus` and planning a job
Placement is limited by memory bandwidth rather than core count: on a 64-vCPU AMD EPYC server, 32-40 processes were fastest and 64 were about 1.5x slower. Peak memory barely depends on `--cpus`: the large arrays are shared between processes.

```bash
# fastest --cpus on this machine, plus memory/time for the full job and a SLURM request
gplace benchmark --refpkg gtdbtk.package.2.refpkg --msa_file user_msa_file.fasta

# memory estimate only
gplace benchmark --refpkg gtdbtk.package.2.refpkg --queries 50000
```
`--cpus auto` uses the value saved by `gplace benchmark` on this machine; without one, it measures a few worker counts on the first part of the run and keeps the fastest. Under SLURM it stays within the CPUs allocated to the job.

### Memory limits and the full GTDB tree
```bash
# stay within 64 GB: untiled if it fits, otherwise tile by tile (identical results)
gplace place --msa_file gtdbtk.bac120.user_msa.fasta --refpkg gtdb_r226_bac120.refpkg \
             --out_file gplace.jplace --cpus 40 --max_memory 64G

# or choose the tile size: tiles of at most 20,000 edges
gplace place ... --tile_edges 20000

# which tile size --max_memory would pick, without placing anything
gplace benchmark --refpkg gtdb_r226_bac120.refpkg --queries 127343 --max_memory 64G

# fastest --cpus and time estimate for the tiled job (benchmarks tile by tile, within 64 GB)
gplace benchmark --refpkg gtdb_r226_bac120.refpkg --msa_file user_msa_file.fasta --queries 127343 --max_memory 64G
```
With `--max_memory`, gPlace estimates the peak memory from the reference tree and the number of queries, keeps 5% of the limit as headroom, and uses the largest tiles that fit (larger tiles are faster). With `place_batch` the choice is made per reference package. Tiled placement rebuilds each tile in every pass, so it takes 1.7-2.5x longer than untiled placement on a tree that fits in memory.

## Testing

Check that gPlace is installed correctly and reproduces pplacer with the `check_install` and `test` subcommands:

```bash
gplace check_install
gplace test --out_dir /tmp/gplace_test --cpus 4
gplace test --models --cpus 4      # also LG, JTT, GTR, empirical frequencies and gamma
```

Unit tests:
```bash
python -m unittest discover -s tests
```

## How it works

pplacer caches three likelihood vectors (distal, proximal, supernode) of *sites × 20* doubles (protein CAT; × rates for gamma, × 4 for DNA) for every node of the reference tree, which is ~104 GB for the largest GTDB-Tk R226 class-level tree. gPlace keeps the arithmetic identical and changes only what is stored:

* distal/proximal vectors are kept for internal nodes only, before pplacer's final power-of-two rescaling, which is re-applied exactly when a vector is used;
* leaf vectors are rebuilt on demand in pplacer's operation order, and supernodes are computed during the prefilter;
* SIMD lanes only process *different* numbers side by side, with no reordered sums and no fused multiply-add;
* the engine links the same numerical libraries as pplacer's statically linked binary: GSL 1.16, glibc 2.19's `hypot`, and correctly-rounded `exp`/`log`/`pow`/`sin` (glibc 2.19's are correctly rounded, and CORE-MATH gives the same results, also inside GSL's gamma functions);
* two engines are built: 20 states (protein) and 4 states (DNA).

Two behaviours of pplacer are reproduced deliberately:

* After each evaluation the tree keeps the *last* evaluated branch lengths, not the optimiser's best point.
* A branch-length assertion failing inside a GSL callback becomes a NaN, which triggers GSL's error path and pplacer's retry from default branch lengths.

See the documentation for the full list of rules that had to be matched.

## Validation

On the 127,343 GlobDB genomes classified with GTDB-Tk R226, across the 8 bacterial class-level reference packages, gPlace reproduces GTDB-Tk's jplace files for **127,343 / 127,343 genomes (100%)**: best edge, full placement list and every printed value, with all 112 tree strings byte-identical. One `place_batch` run took 28 min 25 s with 40 processes under a 64 GB limit (peak 43.8 GiB).

On the full GTDB R226 bacterial tree (273,290 nodes), `--max_memory 64G` placed all 127,343 genomes tile by tile in 2 h 19 min (6 tiles, peak 59.2 GiB; untiled would need ~490 GiB). Every genome that is also in GTDB-Tk's full-tree classification tree (3,459) lands on the same edge with the same branch lengths. Tiled runs on the class-level trees are identical to GTDB-Tk's jplace files in every printed value. See the documentation's validation page.

The other models (LG, JTT, empirical frequencies, GTR, gamma) were compared with the stock pplacer binary on its own test data, the bundled test package, simulated 16S-like DNA (up to 5,000 references) and GTDB-Tk's class-level tree 8: every placement record identical, untiled and tiled, including pplacer's `-inf` values on tree 8 with gamma.

## Citing gPlace

A manuscript is in preparation. Please also cite pplacer, GTDB-Tk and GTDB (see the documentation's references page).

## License

gPlace is released under the GNU General Public License (Version 3). It includes code from CORE-MATH (MIT license) and from the GNU C Library 2.19 (`e_hypot.c`, fdlibm notice preserved); see `gplace/engine/csrc/`.
