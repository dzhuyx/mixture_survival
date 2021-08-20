# Synthetic Taiwan LHID-2010 dataset
This folder contains a synthetic dataset created for runnning data analyses on. Code and relevant parameters used for synthetic data creation are also provided.

`synthetic_Taiwan.rda` is the created synthetic dataset.
`synthetic_Taiwan.R` is the code used to create `synthetic_Taiwan.rda`.
To run `synthetic_Taiwan.R`, two rda files are needed—`Taiwan_X_spec.rda`, which contains distribution parameters of the covariate matrices, and `Taiwan_par_spec.rda`, which contains the model coefficients for generating event time.

Codes used for data analyses can be run on `synthetic_Taiwan.rda` to create synthetic results.
