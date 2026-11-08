# Test data for posterior probabilities (pplacer's -p)

Every `*.pplacer.jplace` was written by the stock pplacer v1.1.alpha19 (commit 807f6f3, the statically
linked release binary); `tests/test_posterior.py` places the same inputs with the same settings, untiled and
tiled, and compares.

* `mini_q4.fasta`: the first 4 queries of `../queries.fasta`, on `../mini.refpkg` (WAG + CAT):
  `mini_p` (`-p`), `mini_p_unif` (`-p --unif-prior`), `mini_p_inform`
  (`-p --inform-prior --prior-lower 0.1`), `mini_p_err1e-6` (`-p --pp-rel-err 0.000001`: the
  integrations fail to converge, inner ones inside the outer one, until pplacer gives up),
  `mini_keep3_p_strikes6` (`--keep-at-most 3 -p --max-strikes 6`: option order).
* `wag_g_q6.fasta`: the first 6 queries of `../models/like_wag_queries.fasta`, on
  `../models/like_wag.refpkg` (WAG + gamma): `wag_g_p` (`-p`).
* `mini_short.refpkg`: `../mini.refpkg` with edge 93 at 4e-6, edge 11 at 0 and leaf edge 45 at
  2e-6; `short_q6.fasta`: 6 queries of `../queries.fasta` placed on or next to them;
  `mini_short_inform` (`-p --inform-prior`): edges shorter than 1e-5 (only the pendant length is
  integrated) and a prior that is 0 at the ML pendant length (upper limit found by halving).
