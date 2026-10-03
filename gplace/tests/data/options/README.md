# Test data for pplacer's options and inputs

Every `*.pplacer.jplace` was written by the stock pplacer v1.1.alpha19 (commit 807f6f3, the build
GTDB-Tk uses) with `-j 4` and the options in its `invocation` field (local paths removed);
`tests/test_options.py` places the same inputs with the same settings and compares.

* `mini.*`, `gtr_g.*`, `wag_g.*`, `dna_cat.*`: search settings on `../mini.refpkg`,
  `../models/like_gtr.refpkg`, `../models/like_wag.refpkg`, `../models/dna_cat.refpkg`.
* `hiv.*`, `actin.*`, `actin_model.json`: pplacer's own likelihood test data
  (pplacer `tests/data/like/{gtr,wag}`, GPL-3.0, https://github.com/matsen/pplacer): PhyML 3.0 trees,
  alignments, stats files and JSON models.
* `nogamma_stats.txt`: `hiv.phy_phyml_stats.txt` with "Discrete gamma model: No".
* `RAxML_info.gtr_fe`, `RAxML_info.wag723`: RAxML 7.2.6 info files in the `-f e` and the 7.2.3
  formats pplacer parses, written for these tests with the PhyML parameters.
* `RAxML8_info.gtr` (RAxML 8.2.12) and `phyml33_stats.txt` (PhyML 3.3): real files that pplacer
  v1.1.alpha19 cannot parse.
