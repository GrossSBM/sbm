library(nnet)

npc <- 150
Q <- 3
n <- npc * Q
p <- 2
# 
B <- matrix(c(-1, 1, 0, 1, -1, 0), byrow = TRUE, ncol = Q)

set.seed(123);X <- matrix(rnorm(p*n), ncol = p)

softmax_rows <- function(x){
    x_shift <- x - apply(x, 1, max)
    ex <- exp(x_shift)
    ex / rowSums(ex)
}

pZ <- softmax_rows(X %*%B)
apply(pZ, 1, which.max)

Z_factor <- sapply(seq_len(nrow(pZ)), function(idx) {
    sample(x = seq_len(ncol(pZ)), size = 1, prob = pZ[idx, ])
})

Z <- t(sapply(seq_along(Z_factor), function(idx) {
    vec <- rep(0, ncol(pZ))
    vec[Z_factor[idx]] <- 1
    vec
}))

# Exp with Sophie's discovery
taus <- pZ  # + rnorm(npc * ncol(pZ), sd = 3) taus <- softmax_rows(taus)

unscaled_taus <- taus/taus[,ncol(taus)]

(t(X)%*%X)^(-1)%*%t(X) %*% log(unscaled_taus)


indata <- data.frame(Z = as.factor(Z_factor), X)
indata$Z <- relevel(indata$Z, ref = paste(ncol(B)))

fit_multinom <- multinom(Z ~ 0+X, data=indata)
summary(fit_multinom)

epsilon <- 0.001

P <- matrix(c(0.5+epsilon, 0, 0,
              0, 0.5, 0,
              0, 0, 0.5-epsilon), byrow = TRUE, nrow = Q)
M <- 1 * (matrix(runif(n * n), n, n) < Z %*% P %*% t(Z)) ## adjacency matrix

devtools::load_all()
# res_sbm <- estimateSimpleSBM(netMat = M)

best_param <- function(fit) {
    fit[["model_parameters"]][[which.max(fit[["ICL"]])]]
}

bm_fit <- BM_bernoulli(
    membership_type = "SBM", 
    adj = M, 
    verbosity=6,
    autosave='',
    plotting=character(0),
    exploration_factor=1.5,
    exploration_direction=numeric(0),
    explore_min=4,
    explore_max=Inf,
    nodes_covariates=list(node = matrix(X, ncol = p)))
bm_fit$estimate()
bm_fit$memberships[[3]][["B"]]

devtools::load_all()
res_sbm_cov <- estimateSimpleSBM(netMat = M, nodes_covariates = matrix(X, ncol = p))
res_sbm_cov$nodesCovarParam
res_sbm <- estimateSimpleSBM(netMat = M)
res_sbm$nodesCovarParam

reorder_beta <- function(fit) {
max_idx <- which.max(fit$ICL)
fit$memberships[[max_idx]]$B

orderLabels <- order(colMeans(fit$memberships[[max_idx]]$alpha) %*%
fit$model_parameters[[max_idx]]$pi, decreasing = TRUE)
beta_hat <- matrix(fit$memberships[[max_idx]]$B[,orderLabels],ncol = max_idx)

beta_tilde <- beta_hat - beta_hat[,ncol(beta_hat)]
beta_tilde
}

reorder_beta(fit_optim)