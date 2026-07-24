# Configuration file for the Sphinx documentation builder.
#
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
from datetime import date

sys.path.insert(0, os.path.abspath('../..'))

# -- Project information -----------------------------------------------------

project = 'mofstructure'
author = 'Dinga Wonanke'
copyright = f'2024-{date.today().year}, {author}'

# read the version from the installed package so the docs cannot drift from it
try:
    from importlib.metadata import version as _package_version
    release = _package_version('mofstructure')
except Exception:
    release = '0.0.0'
version = release

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx',
    'sphinxcontrib.mermaid',
    'sphinx.ext.autosectionlabel',
    'sphinx_copybutton',
    'sphinx.ext.imgconverter',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
master_doc = 'index'

# autosectionlabel collides across documents without this
autosectionlabel_prefix_document = True

autodoc_member_order = 'bysource'
autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': True,
}

napoleon_google_docstring = True
napoleon_numpy_docstring = True

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
}

# -- HTML output -------------------------------------------------------------

html_theme = 'alabaster'
html_logo = '_static/logo.png'
html_static_path = ['_static']
html_css_files = ['style.css']
# adds the mobile navigation bar; the stylesheet degrades gracefully without it
html_js_files = ['nav.js']

html_allow_unicode = True
html_use_index = True
html_use_search = True
html_show_sourcelink = True
html_show_copyright = True
html_show_powered_by = False

# xcode is chosen by measurement, not taste. Against the #edf2f4 code
# background, one of its 61 coloured tokens falls below the 4.5:1 contrast
# threshold, and that one sits at 4.49. Pygments' default style puts 20 tokens
# below it, including every string literal at 3.13:1.
pygments_style = 'xcode'

# Palette: #2b2d42 #8d99ae #edf2f4 #f4f3ee #000000, with #fb8500 for
# navigation, hover and active state. See _static/style.css for the measured
# contrast constraints on where orange and grey may be used as text.
#
# sidebar_width is load bearing: alabaster gives div.bodywrapper a left margin
# of exactly this value, so the stylesheet keeps the sidebar's padding inside
# it with box-sizing: border-box. Without that the sidebar overlaps the text.
html_theme_options = {
    'description': 'Topology, porosity and building-unit analysis of MOFs',
    'github_user': 'bafgreat',
    'github_repo': 'mofstructure',
    'github_button': True,
    'github_type': 'star',
    'fixed_sidebar': True,
    'page_width': '1200px',
    'sidebar_width': '270px',

    'extra_nav_links': {
        'Source on GitHub': 'https://github.com/bafgreat/mofstructure',
        'Package on PyPI': 'https://pypi.org/project/mofstructure/',
        "Author's website": 'https://www.dingawonanke.com',
    },

    'body_bg': '#f4f3ee',
    'body_text': '#2b2d42',
    'link': '#2b2d42',
    'link_hover': '#fb8500',
    'sidebar_header': '#fb8500',
    'sidebar_text': '#8d99ae',
    'sidebar_link': '#edf2f4',
    'sidebar_link_underscore': '#fb8500',
    'sidebar_hr': 'rgba(251, 133, 0, 0.35)',
    'sidebar_search_button': '#fb8500',
    'anchor': '#fb8500',
    'anchor_hover_fg': '#fb8500',
    'gray_1': '#2b2d42',
    'gray_2': '#edf2f4',
    'gray_3': '#8d99ae',
    'pre_bg': '#edf2f4',
    'note_bg': '#edf2f4',
    'note_border': '#2b2d42',
    'warn_bg': '#edf2f4',
    'warn_border': '#fb8500',
}

todo_include_todos = True
