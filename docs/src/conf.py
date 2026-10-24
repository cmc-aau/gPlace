# Configuration file for the Sphinx documentation builder.
import os
import sys
from datetime import datetime

sys.path.insert(0, os.path.abspath('../..'))
from gplace import __author__, __version__, __title__, __maintainer__, __url__

project = __title__
copyright = f'{datetime.now().year}, {__maintainer__}'
author = __author__
pygments_style = 'sphinx'
github_url = __url__
release = __version__
version = __version__

extensions = ['sphinxarg.ext', 'sphinx.ext.napoleon', 'sphinx.ext.autodoc', 'myst_parser']
templates_path = ['_templates']
exclude_patterns = []
source_suffix = {'.rst': 'restructuredtext', '.md': 'markdown'}

html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'logo_only': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': True,
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False,
}
html_static_path = ['_static']
html_context = {
    'display_github': True,
    'github_user': 'cmc-aau',
    'github_repo': 'gPlace',
    'github_version': 'main/docs/src/',
}
html_css_files = ['theme_overrides.css']
html_logo = '_static/gPlace-stacked-dark.svg'     # the mark with the name under it, for the dark sidebar ('g' in off-white)
# favicons (16 and 32 px, each drawn for its size) and the Apple touch icon: _templates/layout.html
html_baseurl = 'https://cmc-aau.github.io/gPlace/'
