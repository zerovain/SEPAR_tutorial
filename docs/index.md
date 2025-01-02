# SEPAR Tutorial

SEPAR (Spatial gene Expression PAttern Recognition) is a computational method designed for analyzing spatial transcriptomics data. This tutorial will guide you through the installation and basic usage of SEPAR.

## Installation

### Prerequisites

Before installing SEPAR, ensure you have Python 3.7 or later installed. SEPAR requires the following package dependencies:

```bash
pandas>=2.0.3
numpy>=1.23.5
scanpy>=1.9.6
anndata>=0.8.0
matplotlib>=3.6.1
scipy>=1.10.0
scikit-learn>=1.2.0
tqdm>=4.64.1
```

### Setting Up the Environment

1. We recommend using conda to create a new environment:

```bash
conda create -n separ python=3.8
conda activate separ
```

2. Install the required packages:

```bash
conda install pandas numpy scipy matplotlib scikit-learn tqdm
conda install -c conda-forge scanpy anndata
```

### Installing SEPAR

Clone the SEPAR repository from GitHub:

```bash
git clone https://github.com/zerovain/SEPAR.git
cd SEPAR
```

## Quick Start

Here's a minimal example to get you started with SEPAR:

```python
import scanpy as sc
from SEPAR_model import SEPAR

# Load your data (example with anndata format)
adata = sc.read_h5ad('your_data.h5ad')

# Initialize SEPAR
separ = SEPAR(adata, n_cluster=8)

# Preprocess the data
separ.preprocess(min_counts=0, min_cells=0, n_top_genes=3000)

# Compute spatial graph
separ.compute_graph()

# Compute weights
separ.compute_weight()

# Run SEPAR algorithm
separ.separ_algorithm(r=30, alpha=1.0, beta=0.1, gamma=0.1)

# Perform clustering
separ.clustering()
```

## Detailed Usage

### Data Preprocessing

The `preprocess` function handles basic data preprocessing:

```python
separ.preprocess(
    min_counts=0,     # Minimum counts required for a cell to pass filtering
    min_cells=0,      # Minimum cells required for a gene to pass filtering
    n_top_genes=3000, # Number of highly variable genes to select
    normalize=True    # Whether to normalize the data
)
```

### Spatial Graph Construction

Compute the spatial neighborhood graph:

```python
separ.compute_graph(
    rad_cutoff1=None,     # Radius cutoff for neighborhood graph
    radius_rate=1.2       # Rate to adjust the radius if radius cutoff is not given
)
```

### Pattern Recognition

Run the SEPAR algorithm to identify spatial patterns:

```python
separ.separ_algorithm(
    r=30,           # Number of patterns to identify
    alpha=1.0,      # Weight for graph regularization
    beta=0.1,       # Weight for sparsity penalty
    gamma=0.1,      # Weight for pattern orthogonality
)
```

### Clustering and Visualization

Perform spatial domain clustering and visualize results:

```python
# Clustering
separ.clustering(n_cluster=8)

```

## Example Datasets

For demonstration purposes, you can use publicly available spatial transcriptomics datasets, such as:
- 10x Visium datasets
- Stereo-seq datasets
- osmFISH datasets
- MISAR-seq datasets

## Output Interpretation

SEPAR generates several key outputs:
1. Spatial patterns (accessible via `separ.Wpn`)
2. Gene loadings (accessible via `separ.Hpn`)
3. Clustering results (accessible via `separ.labelres`)
4. Pattern-specific genes (can be identified using `identify_pattern_specific_genes()`)

## Citation

If you use SEPAR in your research, please cite:

```bibtex
[Citation information will be added upon publication]
```

## Support

For questions and issues, please visit our [GitHub repository](https://github.com/zerovain/SEPAR) or contact us through the issues page.

## License

SEPAR is released under the MIT License. See the LICENSE file in the repository for more details.
