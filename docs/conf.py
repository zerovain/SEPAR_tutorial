import os  
import sys  
sys.path.insert(0, os.path.abspath('..'))  

project = 'SEPAR'  
copyright = '2025'  
author = 'ZL'  

extensions = [  
    'myst_parser',  
    'sphinx.ext.mathjax',  
    'sphinx_rtd_theme',  
]  

# 支持 markdown  
source_suffix = {  
    '.rst': 'restructuredtext',  
    '.md': 'markdown',  
}  

html_theme = 'sphinx_rtd_theme'  

master_doc = 'index'  
