This folder contains annotations of single cell data that have
been generated (semi-)automatically using the code in this repo
like doublet or cell-type annotations.

Since results cannot be necessary reproduced on every system
we store the results here and treat them like precious manually
curated input data. This ensures that while the preprocessing
might be not entirely reproducible, at least the results
can be regenerated from these intermediate data.

* "preliminary" annotations are coarse grained and based on "classical"
  preprocessing. The final annotations are based on scVI preprocessing and
  have more fine-grained annotations.
* The X_scVI files contain the latent representation derived using scVI.
  Since scVI is not necessarily reproducible on different systems,
  this is the basis for making the cell-type annotations reproducible.
