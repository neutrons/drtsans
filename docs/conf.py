# -*- coding: utf-8 -*-
#
# drt-sans documentation build configuration file, created by
# sphinx-quickstart on Tue Sep  3 10:08:07 2019.
#
# This file is execfile()d with the current directory set to its
# containing dir.
#
# Note that not all possible configuration values are present in this
# autogenerated file.
#
# All configuration values have a default; values that are commented out
# serve to show the default.

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))

# The full version, including alpha/beta/rc tags.
import os
import sys
from os.path import basename
from io import StringIO
from docutils.parsers.rst import Directive
from docutils import nodes, statemachine

sys.path.insert(0, os.path.abspath(".."))

from drtsans import __version__ as release

# -- General configuration ------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.doctest",
    "sphinx.ext.todo",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.githubpages",
    "sphinx.ext.graphviz",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = ".rst"

# The master toctree document.
master_doc = "index"

# General information about the project.
project = "drt-sans"
copyright = "2022, SCSE@ORNL"
author = "SCSE@ORNL"

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# 'release' is imported at the top of the file
#
# The short X.Y version.
version = ".".join(release.split(".")[:2])

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = "en"

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This patterns also effect to html_static_path and html_extra_path
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "legacy"]

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True


# -- Options for HTML output ----------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "alabaster"

autosummary_generate = True

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
# html_theme_options = {'show_related': True,
#                       'page_width': '100%'}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']

# Custom sidebar templates, must be a dictionary that maps document names
# to template names.
#
# This is required for the alabaster theme
# refs: http://alabaster.readthedocs.io/en/latest/installation.html#sidebars
# html_sidebars = {
#     '**': [
#         'about.html',
#         'navigation.html',
#         # 'localtoc.html',
#         'relations.html',  # needs 'show_related': True theme option to display
#         'searchbox.html',
#     ]
# }


# -- Options for HTMLHelp output ------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = "drt-sansdoc"


# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',
    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',
    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',
    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, "drt-sans.tex", "drt-sans Documentation", "Jose Borreguero", "manual"),
]


# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [(master_doc, "drt-sans", "drt-sans Documentation", [author], 1)]


# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (
        master_doc,
        "drt-sans",
        "drt-sans Documentation",
        author,
        "drt-sans",
        "One line description of project.",
        "Miscellaneous",
    ),
]

# Napoleon settings
napoleon_google_docstring = False
napoleon_numpy_docstring = True


# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/1.22/", None),
    "mantid": ("https://docs.mantidproject.org/nightly/", None),
    "matplotlib": ("https://matplotlib.org/", None),
}

# -- Custom directives ----------------------------------------------------


class ExecDirective(Directive):
    r"""
    Execute the specified python code and insert the output into the document

    Credit goes to:
    https://stackoverflow.com/questions/27875455/displaying-dictionary-data-in-sphinx-documentation/29789910#29789910
    """
    has_content = True

    def run(self):
        old_stdout, sys.stdout = sys.stdout, StringIO()

        tab_width = self.options.get("tab-width", self.state.document.settings.tab_width)
        source = self.state_machine.input_lines.source(self.lineno - self.state_machine.input_offset - 1)

        try:
            exec("\n".join(self.content))
            text = sys.stdout.getvalue()
            lines = statemachine.string2lines(text, tab_width, convert_whitespace=True)
            self.state_machine.insert_input(lines, source)
            return []
        except Exception:
            error_message = f"Unable to execute python code at {basename(source)}:{self.lineno}:"
            return [
                nodes.error(None, nodes.paragraph(text=error_message), nodes.paragraph(text=str(sys.exc_info()[1])))
            ]
        finally:
            sys.stdout = old_stdout


def setup(app):
    app.add_directive("exec", ExecDirective)
