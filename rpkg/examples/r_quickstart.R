# Minimal R example for TANGO.
# Run from rpkg with devtools::load_all(".").

z <- c(2.1, 1.7, -0.4, 0.8, 1.2)
corr <- diag(length(z))

res <- tango_test(z, corr = corr)
print(res)
