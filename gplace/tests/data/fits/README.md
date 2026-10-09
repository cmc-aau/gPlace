# IQ-TREE and RAxML-NG fits (tests/test_fitted.py)

Fits on the fixed trees of the bundled test packages (branch lengths and model parameters optimised):
IQ-TREE 3.1.4 `iqtree3 -s <aln> -te <tree> -m <model>` on `mini.refpkg` (protein; `iq_prot_*`) and on
`models/like_gtr.refpkg` (HIV DNA; `iq_dna_*`), RAxML-NG 2.0.3 `raxml-ng --evaluate --msa <aln> --tree <tree>
--model <model>` (`rx_*`; the protein tree without its internal labels, which RAxML-NG cannot read).
`*.pplacer.jplace`: stock pplacer v1.1.alpha19 (`pplacer -c <package>`) on the package that taxtastic's
`taxit create --aln-fasta <aln> --tree-file <tree> --tree-stats <fit>` built from the fit, with the queries of
`queries.fasta` (protein) or `models/like_gtr_queries.fasta` (DNA).

`full_*.{ckp.gz,iqtree,treefile}`: IQ-TREE 3.1.4 fits for `place --full_model` (Q.pfam+G4, HKY+F+I+R2 on the DNA set,
EX2, LG4X, LG+G4 with `--gamma-median`), same commands. `full_model_expected.json`: gPlace's best placement of the
first queries and IQ-TREE's log-likelihood of the tree with each query attached there (model fixed, `-blfix`, on the
columns gPlace uses); `iqtree_gamma`: IQ-TREE's gamma rates (`RateGamma::computeRates`, compiled from its source).
`full_partitioned_dna.*`: a partitioned IQ-TREE fit (HIV DNA, two halves, `-p` with GTR+G4) that `--full_model`
does not support yet (test of the error).

`rx_full_*.raxml.{bestModel,bestTree,log}`: RAxML-NG 2.0.3 fits for `place --full_model` (HKY+FC+I+G4 and
GTR+I+R3 on the DNA set; LG4X, LG+G4a (median gamma rates) and WAG+FC+IC+G4 on the protein set), `raxml-ng
--evaluate --msa <aln> --tree <tree> --model <model> --seed 1`. `rx_full_expected.json`: gPlace's best placement
of each query and RAxML-NG's log-likelihood of the tree with the query attached there (`raxml-ng --evaluate
--opt-model off --opt-branches off`, the bestModel values, on the columns gPlace uses); `gamma`: RAxML-NG's
discrete gamma rates (coraxlib's `corax_compute_gamma_cats`, compiled from its source) for a few shapes.
