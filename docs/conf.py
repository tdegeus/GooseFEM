import os
import subprocess
import sys

project = "GooseFEM"
copyright = "2017-2021, Tom de Geus"
author = "Tom de Geus"

subprocess.call("cd ..; python setup.py build --build-type Release -vv", shell=True)
mybuild = os.listdir("../_skbuild")[0]
sys.path.insert(0, os.path.abspath(f"../_skbuild/{mybuild}/cmake-install/python"))

doxydir = "_doxygen"
os.makedirs(doxydir, exist_ok=True)
subprocess.call(f"cmake .. -B{doxydir} -DBUILD_DOCS=1", shell=True)
subprocess.call(f"cd {doxydir}; make html", shell=True)
subprocess.call(f"python -m breathe.apidoc -m -f -p GooseFEM -o api {doxydir}/xml", shell=True)

extensions = [
    "breathe",
    "sphinx.ext.mathjax",
    "sphinx.ext.todo",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx_tabs.tabs",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
html_theme = "furo"

breathe_projects = {
    "GooseFEM": f"{doxydir:s}/xml/",
}
