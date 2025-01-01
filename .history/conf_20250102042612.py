import os  
import sys  
sys.path.insert(0, os.path.abspath('..'))  

project = 'SEPAR'  
copyright = '2025'  
author = 'Zhang Lei'  

extensions = [  
    'myst_parser',  # 用于支持 Markdown  
    'sphinx.ext.mathjax',  
]  

# 支持 markdown  
source_suffix = {  
    '.rst': 'restructuredtext',  
    '.md': 'markdown',  
}  

html_theme = 'sphinx_rtd_theme'  