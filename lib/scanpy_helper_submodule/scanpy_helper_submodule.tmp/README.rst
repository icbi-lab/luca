Scanpy helper functions
=======================

A bunch of custom helper functions to help with cell type annotation 
and dataset integration.

Installation
------------

Development install using 

```bash
git clone git@github.com:grst/scanpy_helpers.git
cd scanpy_helpers
flit install -s
```

This is still evolving, no stable version identifiers are available. 
Consider adding as a git submodule to your project. 

Modules
-------

`de` - Differential expression
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Wrapper functions around differential expression testing methods and convenience
methods for integration with `scanpy`/`anndata`. 

`annotation` - Cell-type annotation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Helper functions for clustering, annotating and subclustering.

`integration` - Altas-level data integration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Helper functions for standardizing data annotations across datasets and 
integrating them. 

