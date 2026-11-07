<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/src/_static/gPlace-stacked-dark.svg">
    <img src="https://raw.githubusercontent.com/cmc-aau/gPlace/main/docs/src/_static/gPlace-stacked.svg" alt="gPlace" width="220">
  </picture>
</p>

<p align="center">
  <a href="https://github.com/cmc-aau/gPlace/actions/workflows/tests.yml"><img alt="Unit tests" src="https://img.shields.io/endpoint?url=https://gist.githubusercontent.com/pchaumeil/08bdf557b640135e66112abccdc87bcc/raw/gplace-tests.json"></a>
  <a href="https://github.com/cmc-aau/gPlace/actions/workflows/tests.yml"><img alt="pplacer parity" src="https://img.shields.io/endpoint?url=https://gist.githubusercontent.com/pchaumeil/08bdf557b640135e66112abccdc87bcc/raw/gplace-parity.json"></a>
</p>

> [!WARNING]
> **Work in progress.** gPlace is still being tested. Before relying on its results, compare them with pplacer's on some of your own data (`gplace compare`), and please report problems on the [issue tracker](https://github.com/cmc-aau/gPlace/issues).

**gPlace** is a replacement for [pplacer](https://github.com/matsen/pplacer): it places aligned sequences on a fixed reference tree and gives exactly pplacer's results, bit for bit, using a fraction of its memory and time. It takes pplacer's command line: replace `pplacer` with `gplace place`.

## Documentation

**https://cmc-aau.github.io/gPlace/**: [installing gPlace](https://cmc-aau.github.io/gPlace/installing/index.html), the [quick start](https://cmc-aau.github.io/gPlace/quickstart.html), the commands, memory limits and large trees, models, validation, and how gPlace works.

## Citing

If you use gPlace, please cite pplacer, and the source of your reference tree (for example GTDB); see the documentation's [References](https://cmc-aau.github.io/gPlace/references.html) page.

## License

gPlace is released under the GNU General Public License (Version 3). It includes code from CORE-MATH (MIT license) and from the GNU C Library 2.19 (`e_hypot.c`, `s_log1p.c`, fdlibm notices preserved); see `gplace/engine/csrc/`.
