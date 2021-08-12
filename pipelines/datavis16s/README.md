datavis16s R package README
================

-   [Usage](#usage)
-   [Dependencies](#dependencies)
-   [Installation](#installation)

R package for graphs for Nephele 16S pipelines.

## Usage

- **Function reference:** [datavis16s R package manual](doc/Reference_Manual_datavis16s.md) (:ref:`pipelines.datavis16s.readme:R package manual`)

- **Nephele User docs:** [GitHub](doc/user_doc.md), [html for Nephele2 website](https://github.niaid.nih.gov/bcbb/nephele2/blob/master/pipelines/datavis16s/doc/datavis16s_pipeline.html), (:ref:`sphinx <pipelines.datavis16s.readme:User docs>`)

- There are 4 main functions for making graphs: `adivboxplot`, `morphheatmap`, `pcoaplot`, and `rarefactioncurve`, as well as `allgraphs` which makes all 4. See the [manual](doc/Reference_Manual_datavis16s.md) for the arguments for these functions.

**Python with rpy2**

- Will need to import rpy2 library:

- ```python
  from rpy2.robjects.packages import importr  ## to import R package
  from rpy2.robjects.vectors import IntVector  ## to handle return values
  import rpy2.rinterface  ## to start R instance
  ```

- The generic wrapper function to be called from rpy2 is [trygraphwrapper](doc/Reference_Manual_datavis16s.md#trygraphwrapper). See the [function help](doc/Reference_Manual_datavis16s.md#trygraphwrapper) for examples.
  - It returns 0 for success, 1 for warnings, and raises `rpy2.rinterface.RRuntimeError` on major error (like problems with input file), which you can catch.  
  - Must pass the full paths for the output directory, mapping file, biom file[1](#fn1).

- `trygraphwrapper`, by default, prints the R sessionInfo to the logfile before it runs the function. If you do not want it to do this (e.g. if you are calling the function multiple times in the same script), you can pass `info = FALSE`.

- To call the function, import the R library and call `trygraphwrapper`:

  - ```python
    datavis16s = importr('datavis16s')
    exit_code = datavis16s.trygraphwrapper(datafile="/path/to/outputs/out.biom", 
                                           outdir="/path/to/outputs/", 
                                           mapfile = "/path/to/inputs/mapfile.txt", 
                                           FUN = functionname, otherarguments_for_functionname)
    ```

**R**

-   You can use `readindata` to create an ampvis2 object, and pass that instead of the mapping file and biom file.

**Sampling depth**
-   The sampling depth argument, `sampdepth`, is optional for all functions.
-   If specified, it is used to remove samples with read counts below `sampdepth`.
-   For functions with `rarefy` argument, setting to TRUE will rarefy the OTU table to `sampdepth` reads. `adivboxplot` will use the smallest sample size to rarefy if `sampdepth` is not specified.
-   Currently, we use 10000 as the default.

1Can optionally pass tab-delimited text file to each of the functions instead of the biom file. See \[`readindata`\]\[doc/Reference\_Manual\_datavis16s.md\#readindata\] for more details.

### Outputs

-   `trygraphwrapper` creates a subdirectory “graphs” in output directory, and passes that as the output directory to the individual functions. The individual graph functions just use the specified output directory passed to them. - Creates html files in outdir as well as subdirectory “lib” which contains the external js,css,etc files for the graphs. See [user doc](https://github.niaid.nih.gov/bcbb/nephele2/blob/master/pipelines/datavis16s/doc/user_doc.md) for more info.

## Dependencies

-   For the Nephele AMI, debian package libxml2-dev is needed.

-   Requires [rpy2](https://rpy2.bitbucket.io) to use with python.

-   To install all R dependencies, run the commands in [dependencies.R](scripts/datavis16s_dependencies.R). You may need to change the path to the datavis16s directory.

-   To build the package along with the user and package documentation, run the commands in [build.R](scripts/build.R). This will produce the standard man files for within-R help, as well as the [package manual](doc/Reference_Manual_datavis16s.md) and the [user docs](doc/user_doc.md)

    ``` r
    source("scripts/build.R")
    ```
- See [Tools and References](doc/user_doc.md#tools-and-references) in the user doc for dependency citations.

## Installation

-   You can install this package, datavis16s, with devtools within R from a locally cloned repository (may need to change directory):

    ``` r
    devtools::install_local("/path/to/nephele2/pipelines/datavis16s", 
        dependencies = TRUE, force = TRUE)
    ```

-   or you can install from the command line (**this does not install dependencies, which should already be installed**):

    ``` bash
    R CMD INSTALL --no-help /path/to/nephele2/pipelines/datavis16s
    ```

-   or to use devtools to install datavis16s from the NIAID github, you will need to generate a [GitHub personal access token](https://help.github.com/articles/creating-a-personal-access-token-for-the-command-line/). In R:

    ``` r
    # change token to token string
    Sys.setenv(GITHUB_PAT = "token")
    # change ref to whichever branch
    devtools::install_github("bcbb/nephele2/pipelines/datavis16s", 
        host = "https://github.niaid.nih.gov/api/v3", ref = "datavis16s", 
        dependencies = TRUE)
    ```
