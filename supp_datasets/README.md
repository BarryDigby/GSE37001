This directory contains results of differentially expressed genes, isofroms, introns and exons given by Dominissini et. al.

The only change made to the datasets was changing the value `1.79769e+308` to `Inf` in the isoforms file. These values were tedious to work with in R (certain functions failed working on the column, and I could not view the dataframe) hence their replacement with `Inf`. For posterity, the original has been uploaded in addition to the `Inf` isoforms file. 
