# Logo and social preview images

The logo files are in `../src/_static/`:

* `gPlace-mark.svg` (large sizes), `gPlace-mark-medium.svg`, `gPlace-icon-32.svg` and `gPlace-icon-16.svg`
  (drawn for 32 and 16 px), `gPlace-lockup.svg` (mark and wordmark): the logo as designed.
* `gPlace-stacked.svg` / `gPlace-stacked-dark.svg`: the mark with the wordmark under it (README; docs sidebar).
* `gPlace-lockup-dark.svg`, `gPlace-stacked-dark.svg`: for dark backgrounds, the wordmark's "g" in off-white.
  The mark itself reads on dark backgrounds as it is.
* `gPlace-mark-white.svg` and `gPlace-mark-white-{1024,512,256,180}.png`: the mark centred on a white square
  (avatars, slides; 180 px is the Apple touch icon of the docs).
* `gPlace-favicon-16.png`, `gPlace-favicon-32.png`: the docs' browser-tab icons, from the 16 and 32 px icons.

Here: `gplace-social-preview.png` (logo and one line) and `gplace-social-preview-minimal.png` (logo only),
GitHub social preview cards, 1280 x 640; upload one in the repository's Settings > General > Social preview.

To redraw them after a change to the logo or the wording:

    python docs/social/render_social.py          # the two cards, from social_preview*.html
    python docs/social/render_icons.py           # the white-square PNGs and the favicons

These files are not part of the published documentation (Sphinx builds `docs/src` only).
