"""
DiTing: A Pipeline to Infer and Compare Biogeochemical Pathways From Metagenomic and Metatranscriptomic Data
"""

__author__ = "Chunxu Xue; Heyu Lin"
__contact__ = "xuechunxu@outlook.com; heyu.lin@qut.edu.au"
from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("diting")
except PackageNotFoundError:
    __version__ = "unknown"
