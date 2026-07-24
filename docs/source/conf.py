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

html_allow_unicode = True
html_use_index = True
html_use_search = True
html_show_sourcelink = True
html_show_copyright = True
html_show_powered_by = False

# Palette: #000000 #14213d #fb8500 (navigation) #fca311 #e5e5e5 #ffffff
# The colours alabaster renders itself are set here; anything it does not
# expose as an option is handled in _static/style.css.
#
# sidebar_width is load bearing: alabaster gives div.bodywrapper a left margin
# of exactly this value, so the stylesheet keeps the sidebar's padding inside
# it with box-sizing: border-box. Widening the padding without that makes the
# sidebar overlap the body text.
html_theme_options = {
    'description': 'Topology, porosity and building-unit analysis of MOFs',
    'github_user': 'bafgreat',
    'github_repo': 'mofstructure',
    'github_button': True,
    'github_type': 'star',
    'fixed_sidebar': True,
    'page_width': '1200px',
    'sidebar_width': '270px',

    'body_text': '#1f2a44',
    'link': '#14213d',
    'link_hover': '#fb8500',
    'sidebar_header': '#fb8500',
    'sidebar_text': '#c8cedb',
    'sidebar_link': '#eaeef5',
    'sidebar_link_underscore': '#fb8500',
    'sidebar_hr': 'rgba(251, 133, 0, 0.3)',
    'sidebar_search_button': '#fb8500',
    'anchor': '#fb8500',
    'anchor_hover_fg': '#fb8500',
    'gray_1': '#14213d',
    'gray_2': '#f4f5f7',
    'gray_3': '#5a6478',
    'pre_bg': '#f4f5f7',
    'note_bg': '#f7f8fa',
    'note_border': '#14213d',
    'warn_bg': '#fff8f0',
    'warn_border': '#fb8500',
}

todo_include_todos = True
