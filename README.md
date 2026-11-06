<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/src/_static/gPlace-stacked-dark.svg">
    <img src="docs/src/_static/gPlace-stacked.svg" alt="gPlace" width="220">
  </picture>
</p>

<p align="center">
  <a href="https://github.com/cmc-aau/gPlace/actions/workflows/tests.yml"><img alt="Unit tests" src="https://img.shields.io/endpoint?url=https://gist.githubusercontent.com/pchaumeil/08bdf557b640135e66112abccdc87bcc/raw/gplace-tests.json"></a>
  <a href="https://github.com/cmc-aau/gPlace/actions/workflows/tests.yml"><img alt="pplacer parity" src="https://img.shields.io/endpoint?url=https://gist.githubusercontent.com/pchaumeil/08bdf557b640135e66112abccdc87bcc/raw/gplace-parity.json"></a>
</p>

> [!WARNING]
> **Work in progress.** gPlace is under active development and its results are not yet fully tested. It reproduces pplacer's placements on the data sets described under [Validation](#validation), but it has not been tested widely beyond them, and its newest features (`--full_model` for IQ-TREE and RAxML-NG fits, macOS builds) have had the least use. Options and outputs may still change. Before relying on gPlace for published results, check it on your own data, for example by comparing its jplace files with pplacer's on a subset of your genomes (`gplace compare`). Please report problems on the [issue tracker](https://github.com/cmc-aau/gPlace/issues).

**gPlace** is a replacement for [pplacer](https://github.com/matsen/pplacer) that gives **exactly pplacer's results**, bit for bit, using a fraction of its memory and time. It places aligned sequences on a fixed reference tree as pplacer v1.1.alpha19 does, from the same reference packages and alignments, with every model pplacer supports (WAG, LG, JTT for protein and GTR for DNA; CAT or gamma rates; model or empirical frequencies), and writes the same jplace file.

* **Identical output:** same edges, likelihoods, like-weight ratios and branch lengths. Every number in the jplace file is textually identical to pplacer's (`%.12g`), as is the reference tree string.
* **Bounded memory:** pplacer holds three likelihood vectors for every node of the reference tree; gPlace holds two for each internal node, and with `--max_memory` it places tile by tile within a limit. The full GTDB R226 bacterial tree (136,646 sequences) fits in 64 GB, where pplacer used 655 GiB.
* **Fast:** on the full GTDB R226 bacterial tree, gPlace placed 10,000 genomes in 14 minutes on the server where pplacer took 17 hours, with identical results.

## Features

* **pplacer's command line:** replace `pplacer` with `gplace place` and keep the rest: a reference package (`-c`) or reference files (`-t`/`-r`/`-s`/`-d`), query alignments in FASTA or Stockholm (gzipped or not) at the end, `-o`, `--out-dir`, `-j`. The jplace files are pplacer's, with pplacer's names. See [Replace a pplacer command](#replace-a-pplacer-command).
* **pplacer's behaviour:** reference packages with a taxonomy get pplacer's `classification` field (edge painting, or `--mrca_class`); identical query sequences are placed once and written as one record; queries that don't overlap the reference alignment stop the run as in pplacer, or are dropped with `--discard_nonoverlapped`.
* **pplacer's options:** its search settings (`--max_strikes`, `--keep_at_most`, `--no_pre_mask`, ...), posterior probabilities (`-p`), model refinement (`--always_refine`, or reference sequences in the query file) and reference files with RAxML 7, PhyML 3.0 or JSON models, all with pplacer's results, under pplacer's names (`--keep-at-most`, ...) or gPlace's (`--keep_at_most`, ...).
* **Every pplacer model:** WAG, LG, JTT (protein) and GTR (DNA), with FastTree CAT or gamma rate variation and model or empirical frequencies, read from the reference package as pplacer does. `--fix_gamma_underflow` optionally avoids pplacer's gamma underflow on large trees (likelihood `-inf`).
* **Trees from IQ-TREE and RAxML-NG:** reference packages built from their fits (`taxit create --tree-stats`), or the `.iqtree` file / RAxML-NG log given directly with `--stats`, place as pplacer places them; gPlace says which parts of the fitted model pplacer does not use. `--full_model` places with the fit's complete model instead (IQ-TREE: every reversible matrix, fitted frequencies, mean or median `+G`, `+I`, `+R`, protein mixtures such as C20 or LG4X; RAxML-NG: every DNA model and protein matrix, `+FC`/`+FU`, mean or median `+G`, `+I`, `+R`, LG4M, LG4X), checked against the fitting program's own likelihoods.
* **Memory limits and huge trees:** `--max_memory 64G` keeps a run within a memory limit. gPlace places without tiles when that fits, and otherwise tile by tile, with the largest tiles that fit, and identical results. This places sequences on trees whose likelihood vectors do not fit in memory, such as the full GTDB bacterial tree (273,290 nodes, which would need ~490 GiB untiled for 127,343 genomes). `--tile_edges N` sets the tile size directly.
* **Many alignments at once:** `gplace place_batch` places many alignments in one run; alignments sharing a reference package are placed in a single pass.
* **No pplacer, guppy or taxtastic needed:** `gplace tog` writes the same tree as `guppy tog` (byte for byte, Newick or phyloXML), and `gplace taxit create/check/info` builds and checks reference packages as taxtastic's `taxit` does.
* **Self-verifying installation:** `gplace check_install` checks the compiled engine bit for bit against values taken from the pplacer v1.1.alpha19 release binary; `gplace test` places bundled sequences and compares them with pplacer's own output.

## Documentation

Documentation for gPlace can be found [here](https://cmc-aau.github.io/gPlace/).

## Installation

gPlace needs Python >= 3.12 on 64-bit Linux (x86-64 or ARM64) or macOS (Apple Silicon or Intel), the platforms pplacer runs on. It uses SIMD kernels on x86-64 CPUs with AVX2 and FMA and on ARM64 (NEON, including Apple Silicon), and a portable engine with the same results on other x86-64 CPUs.

```bash
pip install gtdb-gplace
```
On x86-64 and ARM64 Linux this installs a prebuilt wheel: no compiler needed, one wheel per architecture for every Linux
distribution with glibc >= 2.17 (CentOS 7 and newer) and every Python >= 3.12. The wheel is
built in the manylinux2014 container and only published after it passes `gplace check_install`
and `gplace test` there and on several distributions. On macOS, pip compiles gPlace from source:
install Apple's command line tools first (`xcode-select --install`, which provides clang and `make`).

From source (macOS, or development), which needs a C compiler (gcc 4.8 or newer, e.g. CentOS 7's, or clang, including Apple's) and `make` for the engine and the third-party C it links (GSL 1.16, CORE-MATH, glibc 2.19's `hypot` and `log1p`):
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
| `GPLACE_GSL_TARBALL` | path to `gsl-1.16.tar.gz` for offline installs (otherwise downloaded from a GNU mirror and SHA-256 checked) |
| `GPLACE_GSL_PREFIX` | reuse an existing static, `-fPIC` GSL 1.16 build |
| `GPLACE_BUILD_DIR` | scratch directory for the GSL build (default `~/.cache/gplace/build`, reused between installs) |
| `GPLACE_BUILD_JOBS` | parallel `make` jobs (default 8) |

`gplace build_engine` rebuilds the engine without reinstalling (advanced).

## Usage

You can view the general help menu at any time:
```bash
gplace -h
```

### Replace a pplacer command
Replace `pplacer` with `gplace place`; the rest of the command line stays the same:
```bash
pplacer      -j 32 -c my.refpkg -o queries.jplace queries.fasta
gplace place -j 32 -c my.refpkg -o queries.jplace queries.fasta    # the same jplace file
```
gPlace takes pplacer's options under pplacer's names (`-c`, `-o`, `-j`, `--out-dir`, `-t`/`-r`/`-s`/`-d`, `-m`, `-p`, `--keep-at-most`, ...) or its own (`--refpkg`, `--out_file`, `--cpus`, `--out_dir`, `--keep_at_most`, ...; the query alignment can also be given with `--msa_file`):
```bash
gplace place --refpkg my.refpkg --msa_file queries.fasta --out_file queries.jplace --cpus 32
```
As in pplacer, query files without `-o` are written to `<name>.jplace` in `--out-dir` (`queries.fasta.gz` gives `queries.jplace`), and several query files give one jplace file each. pplacer's `--mmap-file`, `--verbosity` and `--timing` are accepted and change nothing. The options gPlace does not implement (`--map-identity`, `--map-mrca`, `--groups`, `--fig-*`, `--fantasy`, `--write-masked`, `--write-pre-masked`, `--pretend`, `--check-like`) stop the run with a message.

### Place many alignments
```bash
gplace place_batch --batchfile batch.tsv --out_dir gplace_out --cpus 40
```
The file provided to `--batchfile` is a three-column, tab-separated file: query alignment, reference package, output jplace.
```
runs/sample1.fasta    refpkgs/tree_a.refpkg    out/sample1.tree_a.jplace
runs/sample2.fasta    refpkgs/tree_a.refpkg    out/sample2.tree_a.jplace
runs/sample3.fasta    refpkgs/tree_b.refpkg    out/sample3.tree_b.jplace
...
```
All rows that use the same reference package are placed together, so the reference likelihoods are computed once.

### Other models and DNA
The model comes from the reference package's `phylo_model` JSON, as in pplacer (`-m` is ignored for reference packages): `subs_model` WAG, LG, JTT or GTR, `ras_model` Price-CAT or gamma, `empirical_frequencies`. Reference packages made with `taxit create` from FastTree, RAxML or PhyML runs work unchanged; DNA queries use pplacer's nucleotide alphabet (IUPAC codes included).

pplacer's gamma model fails on large trees: once the slow rates of a fast-evolving site underflow to zero, pplacer stops rescaling that site and the likelihood becomes `-inf` (on the GTDB R226 class-level tree 8 with gamma, for 342 of 343 genomes). gPlace reproduces this by default and warns; `--fix_gamma_underflow` fixes it, at the cost of no longer matching pplacer (its results elsewhere change only in the last digits). See the documentation's models page.

### Trees from IQ-TREE or RAxML-NG
```bash
# as pplacer would place them (identical results); gPlace warns about the parts of the fit pplacer ignores
gplace place --tree run.treefile --ref_aln ref.fasta --stats run.iqtree --msa_file queries.fasta --out_file q.jplace
# with the complete fitted model (results differ from pplacer's): IQ-TREE reads run.ckp.gz next to run.iqtree,
# RAxML-NG reads run.raxml.bestModel next to run.raxml.log
gplace place --tree run.treefile --ref_aln ref.fasta --stats run.iqtree --full_model --msa_file queries.fasta --out_file q.jplace
gplace place --tree run.raxml.bestTree --ref_aln ref.fasta --stats run.raxml.log --full_model --msa_file queries.fasta --out_file q.jplace
```
`--full_model` covers every reversible DNA model and protein matrix of IQ-TREE 3, fitted frequencies, `+G` (mean or median), `+I`, `+I+G`, `+R`, `+I+R` and protein mixtures (C10-C60, EX2/EX3/EHO, UL2/UL3, LG4M, LG4X, `MIX{...}`), and every DNA model and protein matrix of RAxML-NG 2 with `+FC`/`+FE`/`+FU`/`+FO`, `+G` (mean or median), `+I`/`+IC`, `+R`, LG4M and LG4X. **Not supported yet:** `+I` inside a mixture, partitioned models, GTR20 and user-defined matrices, DNA mixtures, and models outside reversible DNA/protein substitution (codon, binary, multi-state, morphology, non-reversible, heterotachy `+H`, PoMo, `+ASC`, site frequency models). See the FAQ.

### Compare jplace files
```bash
gplace compare --jplace_a gplace.jplace --jplace_b pplacer.jplace
```

### Placement tree and reference packages (no pplacer or taxtastic needed)
```bash
gplace tog -o placements.tree gplace.jplace                  # = guppy tog
gplace taxit create -P my.refpkg -l bac120 -f msa.faa -t tree.nwk -s fasttree.log    # = taxit create
gplace taxit check my.refpkg
```

### Choosing `--cpus` and planning a job
Placement is limited by memory bandwidth rather than core count: on a 64-vCPU AMD EPYC server, 32-40 processes were fastest and 64 were about 1.5x slower. Peak memory barely depends on `--cpus`: the large arrays are shared between processes.

```bash
# fastest --cpus on this machine, plus memory/time for the full job and a SLURM request
gplace benchmark --refpkg my.refpkg --msa_file queries.fasta

# memory estimate only
gplace benchmark --refpkg my.refpkg --queries 50000
```
`--cpus auto` uses the value saved by `gplace benchmark` on this machine; without one, it measures a few worker counts on the first part of the run and keeps the fastest. Under SLURM it stays within the CPUs allocated to the job.

### Memory limits and very large trees
```bash
# the full GTDB bacterial tree within 64 GB: untiled if it fits, otherwise tile by tile (identical results)
gplace place --msa_file queries.fasta --refpkg gtdb_r226_bac120.refpkg \
             --out_file gplace.jplace --cpus 40 --max_memory 64G

# or choose the tile size: tiles of at most 20,000 edges
gplace place ... --tile_edges 20000

# which tile size --max_memory would pick, without placing anything
gplace benchmark --refpkg gtdb_r226_bac120.refpkg --queries 127343 --max_memory 64G

# fastest --cpus and time estimate for the tiled job (benchmarks tile by tile, within 64 GB)
gplace benchmark --refpkg gtdb_r226_bac120.refpkg --msa_file queries.fasta --queries 127343 --max_memory 64G
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

pplacer caches three likelihood vectors (distal, proximal, supernode) of *sites × 20* doubles (protein CAT; × rates for gamma, × 4 for DNA) for every node of the reference tree, which is ~104 GB for the largest GTDB R226 class-level tree. gPlace keeps the arithmetic identical and changes only what is stored:

* distal/proximal vectors are kept for internal nodes only, before pplacer's final power-of-two rescaling, which is re-applied exactly when a vector is used;
* leaf vectors are rebuilt on demand in pplacer's operation order, and supernodes are computed during the prefilter;
* SIMD lanes only process *different* numbers side by side, with no reordered sums and no fused multiply-add;
* the engine links the same numerical libraries as pplacer's statically linked binary: GSL 1.16, glibc 2.19's `hypot` and `log1p`, and correctly-rounded `exp`/`log`/`pow`/`sin` (glibc 2.19's are correctly rounded, and CORE-MATH gives the same results, also inside GSL's gamma functions);
* two engines are built: 20 states (protein) and 4 states (DNA).

Two behaviours of pplacer are reproduced deliberately:

* After each evaluation the tree keeps the *last* evaluated branch lengths, not the optimiser's best point.
* A branch-length assertion failing inside a GSL callback becomes a NaN, which triggers GSL's error path and pplacer's retry from default branch lengths.

See the documentation for the full list of rules that had to be matched.

## Validation

gPlace was compared with pplacer v1.1.alpha19 (`pplacer -m WAG -c <refpkg>`) on 127,343 GlobDB r226 genomes, each placed on one of the 8 GTDB R226 bacterial class-level reference trees: **127,343 / 127,343 genomes identical** in best edge, full placement list and every printed value, with all 112 tree strings byte-identical. One `place_batch` run took 28 min 25 s with 40 processes under a 64 GB limit (peak 43.8 GiB).

On the full GTDB R226 bacterial tree (273,290 nodes), 10,000 random GlobDB genomes placed with `--max_memory 64G` (5 tiles) were identical to pplacer's output in every printed value: on the same 1 TB server, gPlace took 13 min 47 s and 55.7 GiB, pplacer 17 h 13 min and 654.9 GiB. `--max_memory 64G` also placed all 127,343 genomes on this tree in 2 h 19 min (6 tiles, peak 59.2 GiB; untiled would need ~490 GiB). Tiled runs on the class-level trees are identical to pplacer's jplace files in every printed value. See the documentation's validation page, and its comparison with other placement tools (BSCAMPP, Splicer, APPLES-2, EPA-ng).

The other models (LG, JTT, empirical frequencies, GTR, gamma) were compared with the stock pplacer binary on its own test data, the bundled test package, simulated 16S-like DNA (up to 5,000 references) and the GTDB R226 class-level tree 8: every placement record identical, untiled and tiled, including pplacer's `-inf` values on tree 8 with gamma.

## Citing gPlace

gPlace reproduces pplacer's placements. If you use gPlace, please cite pplacer, and the source of your reference tree (for example GTDB); see the documentation's references page.

## License

gPlace is released under the GNU General Public License (Version 3). It includes code from CORE-MATH (MIT license) and from the GNU C Library 2.19 (`e_hypot.c`, `s_log1p.c`, fdlibm notices preserved); see `gplace/engine/csrc/`.
