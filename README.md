
# FPLDpdQ

**Probability Density Quantile Estimation for the Five-Parameter Lambda Distribution**

`FPLDpdQ` is an R package for fitting the five-parameter lambda distribution (FPLD) using probability density quantile (pdQ) methods.

The package estimates the three shape parameters by minimising the Hellinger distance between empirical and theoretical pdQ functions. The location and inverse-scale parameters are subsequently estimated using several alternative approaches.

## Installation

The package can be installed from GitHub using:

```r
# install.packages("remotes")
remotes::install_github("YOUR_USERNAME/FPLDpdQ")
```

*The GitHub repository address will be updated after publication.*

## Estimation methods

Five estimation variants are available:

| Method | Description |
|---|---|
| `dens` | Density-based estimation |
| `qm` | Quartile matching |
| `mom` | Moment matching |
| `qr` | Quantile regression |
| `qr` with `weights = TRUE` | Weighted quantile regression |

All five variants use the same pdQ-based shape estimation procedure.

## Usage

```r
library(FPLDpdQ)

set.seed(123)
x <- rnorm(100)

# Density-based estimation (default)
fit_dens <- FPLDpdQ(x, method = "dens")

# Quartile matching
fit_qm <- FPLDpdQ(x, method = "qm")

# Moment matching
fit_mom <- FPLDpdQ(x, method = "mom")

# Quantile regression
fit_qr <- FPLDpdQ(x, method = "qr")

# Weighted quantile regression
fit_wqr <- FPLDpdQ(x, method = "qr", weights = TRUE)

fit_dens
```

## Function arguments

```r
FPLDpdQ(x, J = 50, method = "dens", tau = 0.5, weights = FALSE)
```

- `x`: Numeric sample.
- `J`: Number of evaluation points used in empirical pdQ estimation.
- `method`: Location and scale estimation method.
- `tau`: Quantile regression parameter.
- `weights`: Whether to use weighted quantile regression.

## Output

The function returns a list containing five estimated parameters:

- `lambda1`: Location parameter.
- `lambda2`: Inverse-scale parameter.
- `lambda3`, `lambda4`, `lambda5`: Shape parameters.

## License

MIT License.


