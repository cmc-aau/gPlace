"""Render the GitHub social preview cards (1280 x 640, GitHub's recommended size):
social_preview.html -> gplace-social-preview.png, social_preview_minimal.html -> gplace-social-preview-minimal.png.

    python docs/social/render_social.py        (needs playwright with Chromium, and sphinx_rtd_theme for the fonts)

The cards use the logo's lockup and the documentation's typeface (Lato), taken from the installed sphinx_rtd_theme.
"""
import os

import sphinx_rtd_theme
from playwright.sync_api import sync_playwright

HERE = os.path.dirname(os.path.abspath(__file__))
FONTS = 'file://' + os.path.join(os.path.dirname(sphinx_rtd_theme.__file__), 'static', 'css', 'fonts')

with sync_playwright() as p:
    browser = p.chromium.launch()
    page = browser.new_page(viewport={'width': 1280, 'height': 640}, device_scale_factor=1)
    for source, png in (('social_preview.html', 'gplace-social-preview.png'),
                        ('social_preview_minimal.html', 'gplace-social-preview-minimal.png')):
        with open(os.path.join(HERE, source)) as fh:
            html = fh.read().replace('FONTS/', FONTS + '/')
        tmp = os.path.join(HERE, '.render_' + source)          # next to the card, so its relative paths resolve
        with open(tmp, 'w') as fh:
            fh.write(html)
        try:
            page.goto('file://' + tmp)
            page.wait_for_timeout(800)
            page.screenshot(path=os.path.join(HERE, png))
        finally:
            os.remove(tmp)
        print('wrote', png)
    browser.close()
