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

# Palette: #000000 #14213d #fca311 #e5e5e5 #ffffff
# The colours alabaster itself renders are set here; everything alabaster does
# not expose as an option is handled in _static/style.css.
html_theme_options = {
    'description': 'Topology, porosity and building-unit analysis of MOFs',
    'github_user': 'bafgreat',
    'github_repo': 'mofstructure',
    'github_button': True,
    'github_type': 'star',
    'fixed_sidebar': True,
    'page_width': '1180px',
    'sidebar_width': '260px',

    'body_text': '#14213d',
    'link': '#14213d',
    'link_hover': '#fca311',
    'sidebar_header': '#fca311',
    'sidebar_text': '#e5e5e5',
    'sidebar_link': '#ffffff',
    'sidebar_hr': '#fca311',
    'gray_1': '#14213d',
    'gray_2': '#e5e5e5',
    'gray_3': '#000000',
    'pre_bg': '#e5e5e5',
    'note_bg': '#e5e5e5',
    'note_border': '#14213d',
    'warn_bg': '#fdf0d5',
    'warn_border': '#fca311',
}

todo_include_todos = True
