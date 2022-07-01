# MetaExact2x2

R package for exact inference for meta-analysis of 2x2 tables.

## Installation

To install this package from GitHub:

```R
library("devtools")
devtools::install_github("UW-GAC/MetaExact2x2")
```

## Citation

## Example

```R
library(MetaExact2x2)

# load the example data
data(example, package = 'MetaExact2x2')
dim(example)
head(example)
```

Try running on all 48 studies:
```R
meta_exact_2x2(example)
# Error in meta_exact_2x2(example) : 
#   The polynomial roots of study 41 could not be solved. Either remove that study or import roots using the lambdas input.
```

Run on the first 40 studies:
```R
meta_exact_2x2(example[1:40,])
#                   est      ci.lb     ci.ub         p
# cMLE_Blaker 0.3836173 -0.1523135 0.9459545 0.1714682
# cMLE_Fisher 0.3836173 -0.1398335 0.9070682 0.1508932
```

Provide polynomial roots calculated with Mathematica for studies 41 and 42 and re-run:
```R
lambdas <- vector("list", length = nrow(example))
lambdas[[41]] <- c(-1.26469, -1.22776, -1.19768, -1.17118, -1.14699, -1.12446,
                    -1.10319, -1.08291, -1.06342, -1.04459, -1.02629, -1.00842,
                    -0.990897, -0.973644, -0.956585, -0.939643, -0.922737, -0.905774,
                    -0.88864, -0.871185, -0.853192, -0.83431, -0.813869, -0.790103)
lambdas[[42]] <- c(-3.35687, -3.25556, -3.1744, -3.10398, -3.04065, -2.98251, -2.92841,
                    -2.87759, -2.8295, -2.78374, -2.74002, -2.69807, -2.65772, -2.61878,
                    -2.58114, -2.54467, -2.50927, -2.47487, -2.44138, -2.40873, -2.37688,
                    -2.34576, -2.31533, -2.28555, -2.25638, -2.22778, -2.19972, -2.17217,
                    -2.14511, -2.11849, -2.09231, -2.06654, -2.04115, -2.01614, -1.99146,
                    -1.96712, -1.94309, -1.91935, -1.89589, -1.87269, -1.84974, -1.82702,
                    -1.80451, -1.7822, -1.76008, -1.73812, -1.71632, -1.69465, -1.6731,
                    -1.65165, -1.63029, -1.60898, -1.58772, -1.56646, -1.5452, -1.52389,
                    -1.5025, -1.48098, -1.4593, -1.43737, -1.41513, -1.39248, -1.36926,
                    -1.3453, -1.32027, -1.29369, -1.26456, -1.23027)

meta_exact_2x2(example, lambdas = lambdas)
#                   est      ci.lb     ci.ub         p
# cMLE_Blaker 0.3548394 0.01025874 0.6939258 0.0458027
# cMLE_Fisher 0.3548394 0.02897376 0.6807051 0.0328239
```

