dada2nephele R package
==============================

  - [Usage](#usage)
  - [Requirements](#requirements)
  - [Installation](#installation)

#### Usage

  - **Function reference:** [dada2nephele R package manual](doc/Reference_Manual_dada2nephele.md)
  - **R:** To use the dada2nephele R package on its own outside of Nephele, see the example script, [standalone.R](scripts/standalone.R) and the ``dada2compute`` [function help](doc/Reference_Manual_dada2nephele.md#dada2compute).
  - **Python:** With rpy2, use ``trycomputewrapper``. See [trycomputewrapper
    help](doc/Reference_Manual_dada2nephele.md#dada2compute).
  
    - Will need to pass the location of the reference database see refdb and refdb\_species
      arguments.
    
  - The global constants are set in an R option ``dparams`` when the package loads. This includes default values for parameters not set by the user as well as those that are user options. See ``.onLoad``
    [help](doc/Reference_Manual_dada2nephele.md#onload).
  - [Nephele user docs](doc/user_doc.md)

#### Requirements

- **R version**: We use R version **4.0.x** and Bioconductor version **3.12**. 
- **Databases**
  - [dada2 formatted SILVA
    dbs](https://benjjneb.github.io/dada2/training.html)
  - HOMD databases: in s3 nephele-dev-nephele-refdbs/dada2_homd.tar.gz
      - Original files downloaded from [HOMD db
        page](http://www.homd.org/index.php?name=seqDownload&file&type=R)
      - Used HOMD 16S rRNA RefSeq Version 15.22 (Starts from position 9)
        & eHOMD 16S rRNA RefSeq Version 15.22 Taxonomy File for MOTHUR to
        format the DADA2 database.  Script: [format_homd_dada2.py](../format_homd_dada2.py) (:py:mod:`format_homd_dada2.py <nephele2.pipelines.DADA2.format_homd_dada2>`)
  - [DECIPHER modified SILVA r132 SSU
    db](http://www2.decipher.codes/Downloads.html)
- R dependencies installed from [dependencies.R](scripts/dependencies.R), installation below.

#### Installation

  - To install all dependencies, run the commands in
    [dependencies.R](scripts/dependencies.R). You may need to change the
    path to the dada2nephele directory.

  - To build the package along with the user and package documentation,
    run the commands in [build.R](scripts/build.R). This will produce
    the standard man files for within-R help, as well as the [package
    manual](doc/Reference_Manual_dada2nephele.md) and the [user
    docs](doc/user_doc.rst).
    
    ``` r
    source("scripts/build.R")
    ```

##### Installing dada2nephele

  - You can install this package, dada2nephele, with remotes within R
    from a locally cloned repository (may need to change directory):
    
    ``` r
    remotes::install_local("/path/to/nephele2/pipelines/DADA2/dada2nephele", 
        dependencies = TRUE, force = TRUE)
    ```

  - or you can install from the command line (**this does not install
    dependencies, which should already be installed**):
    
    ``` bash
    R CMD INSTALL /path/to/nephele2/pipelines/DADA2/dada2nephele
    ```

  - or to use remotes to install dada2nephele from the NIAID github,
    you will need to generate a [GitHub personal access
    token](https://help.github.com/articles/creating-a-personal-access-token-for-the-command-line/).
    In R:
    
    ``` r
    # change token to token string
    Sys.setenv(GITHUB_PAT = "token")
    # change ref to whichever branch
    remotes::install_github("bcbb/nephele2/pipelines/DADA2/dada2nephele", 
        host = "https://github.niaid.nih.gov/api/v3", ref = "dada2nephele", 
        dependencies = TRUE)
    ```
