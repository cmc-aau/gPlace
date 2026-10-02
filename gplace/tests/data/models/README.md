# Test sets for pplacer's other models

Every `*.pplacer.jplace` was written by the stock pplacer v1.1.alpha19 (commit 807f6f3, the build
GTDB-Tk uses): `pplacer -j 4 -c <refpkg> -o <out> <queries>`.

* `mini_<model>.phylo_model.json`: the bundled `mini.refpkg` (its tree and alignment) with another
  model: LG, JTT, empirical frequencies (`_emp`), gamma (`_g<alpha>x<n_cats>`). Queries:
  `../queries.fasta`.
* `mini_WAG_g0.05x4.fix.jplace`: the same with `--fix_gamma_underflow` (gPlace only; not pplacer
  output). At alpha 0.05 the blocked rescaling changes the last digits.
* `like_gtr.refpkg`, `like_wag.refpkg`: pplacer's own likelihood test data
  (pplacer `tests/data/like/{gtr,wag}`, GPL-3.0, https://github.com/matsen/pplacer): HIV DNA with
  GTR + gamma and actin protein with WAG + gamma, PhyML trees and models. The queries
  (`*_queries.fasta`) are 3 mutated copies of each reference (2-10% substitutions, a few ambiguity
  codes, one fragment), named Q0_/Q1_/Q2_<reference>.
* `dna_cat.refpkg`: the HIV alignment with a FastTree 2.2 `-nt -gtr` tree and taxit refpkg
  (GTR + CAT, empirical frequencies); queries `like_gtr_queries.fasta`.
