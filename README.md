# importABCatlas

## 🧠 About the Tool

**`importABCatlas`** provides a streamlined interface for accessing and analyzing data from the **Allen Brain Cell (ABC) Atlas** using the **Seurat** ecosystem. Built with accessibility and ease-of-use in mind, this package is intended to empower researchers—particularly those more comfortable in R environments—to explore high-resolution single-cell data from the ABC Atlas without needing Python or high-performance computing (HPC) infrastructure.

The **Allen Brain Cell Atlas** is an open-access, multimodal platform developed by the **Allen Institute** as part of the Brain Knowledge Platform. Its mission is to provide a unified, high-resolution view of the molecular and spatial diversity of brain cell types across mammalian species. This groundbreaking initiative continues to expand, offering researchers unprecedented capabilities to:

- 🔬 Identify more cell types in the brain  
- 🧭 Investigate spatial organization of these cell types  
- 🧬 Explore gene expression and co-expression patterns  
- 🧠 Refine the molecular definitions of brain regions

Despite the wealth of data available in the ABC Atlas, programmatic access is currently supported primarily via **Python APIs**. However, nearly half of the neuroscience community—including many wet lab scientists—prefer **R** for downstream single-cell transcriptomic analyses. This creates a gap in accessibility and ease-of-use.

**`importABCatlas`** aims to bridge that gap by:

- Providing seamless access to ABC Atlas data **entirely in R**
- Leveraging the **Seurat** framework for downstream analysis
- Supporting flexible subsetting and filtering of metadata
- Enabling usage on **non-HPC machines**, making it viable for most labs worldwide
- Offering compatibility with other popular R-based frameworks like **SingleCellExperiment** and **monocle3** by allowing users to extract data from Seurat objects

Whether you're investigating a small cell population in a low-resource environment or exploring atlas-wide transcriptomic patterns, `importABCatlas` is designed to provide a simple, robust, and R-native workflow.

---

## 🚧 Current Status

This package is currently **under active development**. Several core functions are already implemented, including support for:

- ✅ Setting up the environment and loading ABC Atlas metadata  
- ✅ Filtering metadata based on user-specified criteria  
- ✅ Downloading and formatting gene expression data into **Seurat** objects  

---

## 🚀 Usage Instructions

You can extract gene expression data from the **Allen Brain Cell Atlas** adult human brain dataset (WHB) in just a few steps:

---

### 🔧 Step 1: Clone the Repository

```bash
git clone https://github.com/me37uday/importABCatlas.git
```

### 🧬 Step 2: Load the Package in R

```
# Set the base download path in the github repo where the raw data will be stored
setwd("/path/to/importABCatlas/")
download_base <- "abc_download_root"

# Load required library (if not already installed: install.packages("devtools"))
library(devtools)

# Load the package from your local clone
devtools::load_all()
```
Replace "/path/to/importABCatlas/" with the actual path to your cloned repository.

### 📥 Step 3: Load Metadata

```
mem.maxVsize(Inf)
WHB_metadata <- load_data_WHB()
```
This returns a list with the following:

- WHB_metadata$cell_metadata: metadata for all WHB cells
- WHB_metadata$gene_data: reference gene annotation
- WHB_metadata$unique_values: helper list showing column names and unique values (used for filtering)

### 🧪 Step 4: Define Filtering Criteria
Create a named vector to specify the cell populations of interest. Example:
```
filters <- c(
  feature_matrix_label = "WHB-10Xv3-Nonneurons",
  cluster = "VendV_17"
)
```
You can explore valid column names and their possible values using:
```
WHB_metadata$unique_values
```
### 🧱 Step 5: Extract Seurat Object
```
obj <- fetch_data(
  metadata  = WHB_metadata$cell_metadata,
  filters   = filters,
  gene_data = WHB_metadata$gene_data,
  genes     = WHB_metadata$gene_data$gene_symbol
)
```
This returns a Seurat object containing the filtered cell population.

### 📤 Export for Other Frameworks (Optional)
To use the data in other single-cell analysis frameworks like SingleCellExperiment, monocle3, etc., you can extract the raw counts and metadata:
```
counts <- obj@assays$RNA@counts
meta   <- obj[[]]
```

### 🧱 Current Dataset Support
Currently, the package supports the WHB (Whole Human Brain) dataset. Additional datasets are in development:

- 🧠 Whole Mouse Brain (WMB)
- ⏳ Aging Mouse Brain (AgingMouse)
- 🔬 Human-Mammalian Brain Atlas (HMBA)
- 🧬 Human Postmortem-Derived Brain Sequencing (PMDBS)

All required raw files will be downloaded to the location defined by download_base.

⚙️ System requirements and package dependencies are currently being tested. Detailed installation and configuration instructions will be available in the final documentation.

Thank you for your patience and support as the package evolves.
