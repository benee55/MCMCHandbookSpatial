
library(cmdstanr)   # Use this package to employ Stan.
library(RSpectra)   # Use this package to eigendecompose the Moran operator.
library(parallel)   # Use this package to parallelize.

options(mc.cores = parallel::detectCores())   # How many CPU cores are available?

# The following function multiplies matrices A and B in embarrassingly parallel fashion.

matmult = function(cl, A, B)
{
    if (ncol(A) != nrow(B))
        stop("The supplied matrices are not conformable.")
    ind = splitIndices(nrow(A), length(cl))
    Alist = lapply(ind, function(ii) A[ii, , drop = FALSE])
    ans = clusterApply(cl, Alist, get("%*%"), B)
    do.call(rbind, ans)
}

# Read in the thrush data and the adjacency data.

thrush = read.csv("thrush_data.csv")
y = thrush$thrush
cover = as.numeric(scale(thrush$cover))           # Normalize the two covariates.
elevation = as.numeric(scale(thrush$elevation))
n = length(y)
A = matrix(0, n, n)
adj = readLines("thrush_adjacency.txt", n = -1)
for (i in 1:n)
{
    temp = as.numeric(strsplit(adj[[i]], " ")[[1]])
    A[i, temp] = 1
}

I = diag(1, n)
J = matrix(1, n, n)
P = I - J / n
q = 400                                       # Use 400 basis vectors.

start = proc.time()

cl = makeCluster(spec = 14, type = "PSOCK")   # Use 14 CPU cores.
M = matmult(cl, matmult(cl, P, A), P)         # Build the Moran matrix.
stopCluster(cl)
eig = eigs_sym(M, q)          # Get the first q eigenvectors of M.
B = eig$vectors               # Store the q basis functions as matrix B.
L = diag(rowSums(A), n) - A   # This is the graph Laplacian.
V = solve(t(B) %*% L %*% B)   # V / tau is the prior covariance matrix for gamma.

(stop = proc.time() - start)  # running time for pre-processing

# Compile the Stan model and draw posterior samples.

start = proc.time()

mod = cmdstan_model("hermit_thrush.stan")
fit = mod$sample(data = list(n = n, q = q, y = y, cover = cover, elevation = elevation,
                 B = B, V = V, zero = rep(0, q)),
                 chains = 1,                         # Generate one sample path
                 iter_sampling = 10000,              # of length 10,000.
                 refresh = 1000)

(stop = proc.time() - start)  # running time for sampling


