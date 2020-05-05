
##### None of these rules are absolute  
There's usually a case to break them in some scenario. Comment when you do though!


Let's go through the prototype /demo.py:  

```python  
#!/usr/bin/env python3

"""Simplistic demo Pipe."""
```

The reason for this shebang line is we don't always know where python3 might be, use env. Gratuitous comment. 

```python  
import argparse
import traceback
import sh
from nephele2.pipelines import pipeline_error
from nephele2.pipelines.pipebase import PipeBase
```

Argparse is used to interface with the pipe.  
Traceback is used to print up exception details (into logs usually).  
sh is a handy module for interacting with other software.  
pipelines.pipeline_error If something goes wrong in a pipeline, we raise a pipeline exception (else raise whatever exception we have).  
PipeBase is a base class for python pipes that can contain anything that is generally useful for pipelines. For example a lot of pipelines will need a reference database, so there's a get_reference_DBs function, and so on.

```python  
class TestPipe(PipeBase):
    """Start pipe class definition."""
```
Inherit PipeBase. Another gratuitous comment. 

```python  
    class Conf:
        """
        Putting this inside the pipe class defn.
        the names in these pipes can be a bit grim. Try to keep them in one place in here.
        Caps denote fnames are constants.
        """
        PLOT_DIR = 'loads_of_plots/'
        FILE_A = 'a_plot.png'
        FILE_B = PLOT_DIR + 'x_plot.png'
        FILE_C = PLOT_DIR + 'y_plot.png'
        # END OF Conf #
```

I'm embedding some constants which are relevent this this pipe, inside the TestPipe class. 

Calling into super with args here, and __init__ in pipebase inits a logfile, and writes these arguments into the top of the logfile.


```python
    @staticmethod
    def bad_thing(do_bad_thing):
        """simulates an error"""
        if do_bad_thing:
            raise pipeline_error.UnknownPipeError(msg='A bad thing happened in do_bad_thing() '\
                                                  'resulting from arg:' + str(do_bad_thing))
```





1. Pre-pipeline - python
   * make sure the fastq/fasta files are there and unzipped, etc. 
   * make sure the mapping file exists
   * make sure arguments make sense
   * make sure we can access the current working directory
     	  * either pass the full path
	  * or change to the working directory before executing script
   * pass the output directory
   

2. During the pipeline - R
   * Logging
     * pass the logfile name in the python - write R function to accept logfile name
     * log to stdout/stderr - try this
   * Catching errors - use try-catch for errors
   * put results and plots in output directory
   * return value 0


4. Post-pipeline - python
   * Compressing output
   
