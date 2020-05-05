dada2nephele R package
==============================

  - [About](#about)
  - [Requirements](#requirements)
  - [Installation](#installation)

#### About

  - **Function reference:** [dada2nephele R package manual](doc/Reference_Manual_dada2nephele.md)
  - **R:** To use the dada2nephele R package on its own outside of Nephele, see the example script, [standalone.R](scripts/standalone.R) and the ``dada2compute`` [function help](doc/Reference_Manual_dada2nephele.md#dada2compute).
  - **Python:** With rpy2, use ``trycomputewrapper``. See [trycomputewrapper
    help](doc/Reference_Manual_dada2nephele.md#dada2compute).

    - Will need to pass the location of the reference database see refdb and refdb\_species
      arguments.
  - The global constants are set in an R option ``dparams`` when the package loads. This includes default values for parameters not set by the user as well as those that are user options. See ``.onLoad``
    [help](doc/Reference_Manual_dada2nephele.md#onload).
  - [Nephele user docs](doc/user_doc.md)
  - Third party tools used for this package are cited in the [references of the user docs](doc/user_doc.md#tools-and-references).  Of note is the [biomformat R library v1.4](https://github.com/joey711/biomformat-oldfork) from which the ``make_biom`` and ``write_biom`` functions in [biomformat.R](R/biomformat.R) are modified and licensed under GPL-2.

#### Requirements

- **R version**: We use R version **3.5.x** and Bioconductor version **3.8**.
- **Databases**
  - [dada2 formatted SILVA
    dbs](https://benjjneb.github.io/dada2/training.html)
  - HOMD databases: in s3 nephele-dev-nephele-refdbs/dada2_homd.tar.gz
      - Original files downloaded from [HOMD db
        page](http://www.homd.org/index.php?name=seqDownload&file&type=R)
      - Used HOMD 16S rRNA RefSeq Version 15.1 (Starts from position 9)
        & eHOMD 16S rRNA RefSeq Version 15.1 Taxonomy File for MOTHUR to
        format the DADA2 database.  Script: [format_homd_dada2.py](../format_homd_dada2.py) (:py:mod:`format_homd_dada2.py <nephele2.pipelines.DADA2.format_homd_dada2>`)
  - [DECIPHER modified SILVA r132 SSU
    db](http://www2.decipher.codes/Downloads.html)
- R dependencies installed from [dependencies.R](scripts/dependencies.R), installation below.

#### Installation

To install all dependencies, run the commands in [dependencies.R](scripts/dependencies.R). You may need to change the path to the dada2nephele directory

##### Installing dada2nephele

  - You can install this package, dada2nephele, with devtools within R
    from a locally cloned repository (may need to change directory):

    ``` r
    devtools::install_local("/path/to/nephele2/pipelines/DADA2/dada2nephele",
        dependencies = TRUE, force = TRUE)
    ```

  - or you can install from the command line (**this does not install
    dependencies, which should already be installed**):

    ``` bash
    R CMD INSTALL /path/to/nephele2/pipelines/DADA2/dada2nephele
    ```

  - or to use devtools to install dada2nephele from the NIAID github,
    you will need to generate a [GitHub personal access
    token](https://help.github.com/articles/creating-a-personal-access-token-for-the-command-line/).
