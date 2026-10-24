"""PNG icons of the gPlace logo, from the SVGs in docs/src/_static:

* gPlace-mark-white-{1024,512,256,180}.png: the mark centred on a white square (gPlace-mark-white.svg);
  180 px is the size Apple uses for home-screen icons (apple-touch-icon);
* gPlace-favicon-16.png / -32.png: the browser-tab icons, each from the icon drawn for that size
  (gPlace-icon-16.svg, gPlace-icon-32.svg), on a transparent background.

    python docs/social/render_icons.py        (needs playwright with Chromium)
"""
import os

from playwright.sync_api import sync_playwright

STATIC = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src', '_static')   # where the icons live
JOBS = [('gPlace-mark-white.svg', f'gPlace-mark-white-{n}.png', n) for n in (1024, 512, 256, 180)] + \
       [('gPlace-icon-16.svg', 'gPlace-favicon-16.png', 16), ('gPlace-icon-32.svg', 'gPlace-favicon-32.png', 32)]

with sync_playwright() as p:
    browser = p.chromium.launch()
    for source, png, size in JOBS:
        with open(os.path.join(STATIC, source)) as fh:
            svg = fh.read()
        page = browser.new_page(viewport={'width': size, 'height': size}, device_scale_factor=1)
        sized = svg.replace('<svg ', f'<svg style="display:block;width:{size}px;height:{size}px" ', 1)
        page.set_content(f'<html><body style="margin:0;background:transparent">{sized}</body></html>')
        page.wait_for_timeout(200)
        page.screenshot(path=os.path.join(STATIC, png), clip={'x': 0, 'y': 0, 'width': size, 'height': size},
                        omit_background=True)
        page.close()
        print(f'wrote {png}')
    browser.close()
