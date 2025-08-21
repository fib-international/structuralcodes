"""Configuration for the docs of structuralcodes."""

import datetime
import importlib
import inspect
import os
import sys

sys.path.insert(1, '..')

project = 'StructuralCodes documentation'
copyright = f'{datetime.datetime.today().year}, fib International'
author = 'fib International'

# General configuration
extensions = [
    'myst_parser',
    'sphinx_copybutton',
    'sphinx_design',
    'sphinx.ext.napoleon',
    'sphinx.ext.autodoc',
    'sphinx.ext.linkcode',
    'sphinx.ext.intersphinx',
]

myst_enable_extensions = [
    'colon_fence',
    'deflist',
    'dollarmath',
]

templates_path = ['_templates']
exclude_patterns = [
    '_build',
    'Thumbs.db',
    '.DS_Store',
    'venv',
    'README.md',
    'node_modules',
    'src',
    'static',
    'build',
    'docs',
]

# Options for HTML output
html_title = project
html_theme = 'furo'
html_static_path = ['_static']
html_theme_options = {
    'footer_icons': [
        {
            'name': 'GitHub',
            'url': 'https://github.com/fib-international/structuralcodes',
            'html': """
                <svg stroke="currentColor" fill="currentColor" stroke-width="0" viewBox="0 0 16 16">
                    <path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0 0 16 8c0-4.42-3.58-8-8-8z"></path>
                </svg>
            """,
            'class': '',
        },
    ],
    'light_logo': 'structuralcodeslogo-300.png',
    'dark_logo': 'structuralcodeslogo-300.png',
    'source_repository': 'https://github.com/fib-international/structuralcodes',
    'source_branch': 'main',
    'source_directory': 'docs/',
    'sidebar_hide_name': False,
}

# Options for intersphinx
intersphinx_mapping = {
    'shapely': ('https://shapely.readthedocs.io/en/stable/', None),
}

# Autodoc settings
autodoc_type_aliases = {
    'ArrayLike': 'ArrayLike',
    'npt.ArrayLike': 'ArrayLike',
}


# Function to resolve links to source
def linkcode_resolve(domain, info):
    """Resolve link to source code using sphinx.ext.linkcode."""
    if domain != 'py':
        return None
    if not info['module']:
        return None
    mod = importlib.import_module(info['module'])
    if '.' in info['fullname']:
        objname, attrname = info['fullname'].split('.')
        obj = getattr(mod, objname)
        try:
            obj = getattr(obj, attrname)
        except AttributeError:
            return None
    else:
        obj = getattr(mod, info['fullname'])

    try:
        file = inspect.getsourcefile(obj)
        lines = inspect.getsourcelines(obj)
    except TypeError:
        return None

    filename = os.path.relpath(file, os.path.abspath('..'))
    start, end = lines[1], lines[1] + len(lines[0]) - 1

    return (
        f'https://github.com/fib-international/structuralcodes/tree/main/{filename}'
        f'#L{start}-L{end}'
    )
