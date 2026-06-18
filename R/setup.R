#' Setup the Python environment
#' 
#' This function sets up the Python environment required to run the data loading functions.
#' @export
setup_environment <- function() {
  library(reticulate)
  
  # Attempt to find Python 3.10 in the system path
  python_path <- Sys.which("python3.12")
  
  # Check if Python 3.10 was found
  if (python_path == "") {
    stop("Python 3.10 not found in the system PATH. Please ensure Python 3.10 is installed and available in the PATH.")
  }
  #use_python(python_path, required = TRUE)
  virtualenv_dir <- "r-reticulate-env"
  
  if (!virtualenv_exists(virtualenv_dir)) {
    virtualenv_create(virtualenv_dir, python = python_path)
    virtualenv_install(virtualenv_dir, packages = c("pandas", "pathlib", "numpy", "anndata"))
    
    # Install abc_atlas_access from GitHub
    virtualenv_install(virtualenv_dir, packages = "git+https://github.com/alleninstitute/abc_atlas_access", ignore_installed = TRUE)
  }
  
  use_virtualenv(virtualenv_dir, required = TRUE)

  print("Finished setting up Python environment")
}
