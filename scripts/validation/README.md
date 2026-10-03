# Validation against GTDB-Tk's pplacer binary

GTDB-Tk's pplacer v1.1.alpha19 (`pplacer v1.1.alpha19-0-g807f6f3`) is a statically linked
x86-64 binary with debug symbols. These gdb Python scripts read values out of the running binary,
and were used to establish gPlace's exactness rules. Run them with the binary from the GTDB-Tk
conda environment. pplacer forks its workers even with `-j 1`, so the scripts set
`follow-fork-mode child`.

| Script | What it records | How it was used |
|---|---|---|
| `gdb_trace_evaluations.py` | every pendant/distal branch length set by `Three_tax` and every `mat_log_like3_c` result (hex) | eval-by-eval comparison with gPlace's trace (`xp_trace_open`); 1,113/1,113 evaluations identical for a test genome |
| `gdb_dump_eigen.py` | the matrix passed to `gsl_eigen_symmv` and its eigenvalues/eigenvectors | showed the input matrix was identical and the solver differed, which led to the GSL 1.16 + glibc 2.19 `hypot` requirement; the dump is `gplace/engine/data/stock_wag_eigen.txt` |
| `gdb_dump_vectors.py` | `dediagonalize` arguments and the first three-taxon vectors | localised differences to the model set-up |
| `gdb_call_libm.py` | `hypot`/`exp`/`log` evaluated inside the binary on random inputs | `gplace/engine/data/stock_libm_vectors.tsv` |
| `gdb_call_log1p.py` | `log1p` evaluated inside the binary (bits passed exactly): every argument GSL passes for 21-300 gamma categories, plus random and special inputs | the `q` rows of `stock_libm_vectors.tsv`; showed glibc 2.19's `s_log1p.c` gives the binary's values |

Examples:

```bash
# trace one genome (outputs TRACE_OUT, hex floats)
TRACE_OUT=stock_trace.txt gdb -q -batch -x gdb_trace_evaluations.py --args \
    pplacer -m WAG -j 1 -c gtdbtk.package.8.refpkg -o q1.json q1.fasta

# eigendecomposition dump
DUMP_DIR=dump gdb -q -batch -x gdb_dump_eigen.py --args pplacer -m WAG -j 1 -c <refpkg> -o x.json q.fasta

# libm values (inputs.txt / stock_vals.txt written in the current directory)
gdb -q -batch -x gdb_call_libm.py --args pplacer --version
```

Symbol names (e.g. `camlThree_tax__set_pend_bl_10608`) are specific to this binary build.
