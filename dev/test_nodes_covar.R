library(nnet)

npc <- 15
Q <- 3
n <- npc * Q
p <- 2
q <- 3
# 
B <- matrix(c(-1, 1, 0, 1, -1, 0), byrow = TRUE, ncol = Q)

G <- matrix(c(5,-5,0,
              -3, 7, 0,
              1, -1, 0), byrow = TRUE, ncol = Q)

set.seed(123);X <- matrix(rnorm(p*n), ncol = p);Y_col <- matrix(rnorm(q*n), ncol = q)

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

pZ2 <- softmax_rows(Y_col %*%G)
apply(pZ2, 1, which.max)

Z2_factor <- sapply(seq_len(nrow(pZ2)), function(idx) {
    sample(x = seq_len(ncol(pZ2)), size = 1, prob = pZ2[idx, ])
})

Z2 <- t(sapply(seq_along(Z2_factor), function(idx) {
    vec <- rep(0, ncol(pZ2))
    vec[Z2_factor[idx]] <- 1
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

indata2 <- data.frame(Z2 = as.factor(Z2_factor), Y_col)
indata2$Z2 <- relevel(indata2$Z2, ref = paste(ncol(G)))

fit_multinom2 <- multinom(Z2 ~ 0+Y_col, data=indata2)
summary(fit_multinom2)


epsilon <- 0.001

P <- matrix(c(0.5+epsilon, 0, 0,
              0, 0.5, 0,
              0, 0, 0.5-epsilon), byrow = TRUE, nrow = Q)
M <- 1 * (matrix(runif(n * n), n, n) < Z %*% P %*% t(Z)) ## adjacency matrix

M_lbm <- 1 * (matrix(runif(n * n), n, n) < Z %*% P %*% t(Z2))
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

bm_fit_lbm <- BM_bernoulli(
    membership_type = "LBM", 
    adj = M_lbm, 
    verbosity=6,
    autosave='',
    plotting=character(0),
    exploration_factor=1.5,
    exploration_direction=numeric(0),
    explore_min=4,
    explore_max=Inf,
    nodes_covariates=list(row = matrix(X, ncol = p), col = matrix(Y_col, ncol = q)))
bm_fit_lbm$estimate()
bm_fit_lbm$memberships[[which.max(bm_fit_lbm$ICL)]]
devtools::load_all()
res_sbm_cov <- estimateSimpleSBM(netMat = M, nodes_covariates = matrix(X, ncol = p))
res_sbm_cov$nodesCovarParam
res_sbm <- estimateSimpleSBM(netMat = M)
res_sbm$nodesCovarParam

index_mat <- seq(1,nrow(M_lbm)*ncol(M_lbm))


set.seed(321);to_remove <- sample(x = index_mat, size = max(1, floor(0.1 * length(index_mat))))

M_lbm_mis <- M_lbm
M_lbm_mis[to_remove] <- NA 

res_lbm <- estimateBipartiteSBM(netMat = M_lbm, dimLabels = c(row = "tutu", col = "lala"))

res_lbm_cov <- estimateBipartiteSBM(netMat = M_lbm, nodes_covariates = list(row = matrix(X, ncol = p), col = matrix(Y_col, ncol = q)), dimLabels = c(row = "tutu", col = "lala"))

res_lbm_mis <- estimateBipartiteSBM(netMat = M_lbm_mis)

res_lbm_cov_mis <- estimateBipartiteSBM(netMat = M_lbm_mis, nodes_covariates = list(row = matrix(X, ncol = p), col = matrix(Y_col, ncol = q)))