devtools::load_all()

set.seed(123)
npc <- 50 # nodes per class
Q <- 3 # classes
n <- npc * Q # nodes


col_covariates <- cbind(1,
                        c(rep(6, npc), rep(0, npc), rep(-1,npc)),
                        c(rep(2,npc), rep(1, npc/2), rep(-1,1.5*npc)),
                        c(rep(1,npc), rep(2, npc/2), rep(-1,1.5*npc)))
row_covariates <- cbind(
    1,
    c(rep(3, npc), rep(-3, npc), rep(0, npc)),
    c(rep(-3, npc), rep(3, npc), rep(0, npc))
)
colnames(row_covariates) <- c("X1", "X2", "X3")
colnames(col_covariates) <- c("X1", "X2", "X3", "X4")


B <- matrix(
    c(
        -2.5, -2.5, 0,
        2.5, -2.5, 0,
        -2.5, 2.5, 0
    ),
    nrow = 3,
    byrow = TRUE
)

G <- matrix(
    c(
        -2.5, -2.5, 0,
        2.5, -2.5, 0,
        -2.5, 2.5, 0,
        -2, 2, 0
    ),
    nrow = 4,
    byrow = TRUE
)

softmax_rows <- function(x){
    x_shift <- x - apply(x, 1, max)
    ex <- exp(x_shift)
    ex / rowSums(ex)
}
Z1 <- round(softmax_rows(row_covariates %*% B))
Z2 <- round(softmax_rows(col_covariates %*% G))


Z1_factor <- as.factor(apply(Z1, 1, which.max))
indata <- data.frame(Z1 = Z1_factor, row_covariates)
indata$Z1 <- relevel(indata$Z1, ref = "3")

new_data <- data.frame(X1 = c(1,1,1), X2 = c(3, -3, 0), X3 = c(-3,3,-3))

require(nnet)

multinom_fit_Z1 <- multinom(Z1 ~ -1 + X1 + X2 + X3, data = indata)
summary(multinom_fit_Z1)
predict(multinom_fit_Z1, newdata = new_data)

Z2_factor <- as.factor(apply(Z2, 1, which.max))
indata2 <- data.frame(Z2 = Z2_factor, col_covariates)
# indata2[["Z2"]] <- relevel(indata[["Z2"]], ref = "3")
multinom_fit_Z2 <- multinom(Z2 ~ -1 + X1 + X2 + X3 + X4, data = indata2)
summary(multinom_fit_Z2)

# Z1<-diag(Q[1])%x%matrix(1,npc[1],1)
# Z2 <- diag(Q[2]) %x% matrix(1, npc[2], 1)


P <- matrix(runif(Q * Q), Q, Q)
M <- 1 * (matrix(runif(n * n), n, n) < Z1 %*% P %*% t(Z1)) ## adjacency matrix

start_Z <- matrix(1/3, nrow = n, ncol = Q)



# fit_only <- dispatcher(membership_name = "LBM", model_name = "bernoulli", membership_init = list(Z1=Z1, Z2=Z2), network = list(adjacency = M), real_EM = TRUE)

# fit_covariates$membership$alpha1
# fit_only$membership$alpha1
library(blockmodels)
fit_no_cov <- estimateSimpleSBM(netMat = M)
fit_cov <- estimateSimpleSBM(netMat = M, nodes_covariates = row_covariates)


sbm_no_covar <- BM_bernoulli(membership_type = "SBM", adj = M, plotting = character(0), ncores = 1L, verbosity = 6)
sbm_no_covar$estimate()

sbm <- BM_bernoulli(membership_type = "SBM", adj = M, plotting = character(0), ncores = 1L, verbosity = 6)
sbm$estimate()

sbm_cov_gradient <- BM_bernoulli(membership_type = "SBM", adj = M, plotting = character(0), ncores = 1L, nodes_covariates = list(node = row_covariates), verbosity = 6)
sbm_cov_gradient$estimate()


sbm_cov_bfgs <- BM_bernoulli(membership_type = "SBM", adj = M, plotting = character(0), ncores = 1L, nodes_covariates = list(nodes = row_covariates), verbosity = 6)
sbm_cov_bfgs$estimate()

(mb_gradient <- microbenchmark("Gradient" = {
sbm_cov <- BM_bernoulli(membership_type = "SBM", adj = M, plotting = character(0), ncores = 1L, nodes_covariates = list(nodes = row_covariates), verbosity = 6)
sbm_cov$estimate()}, times = 3L
))

(mb_bfgs <- microbenchmark("BFGS" = {
sbm_cov <- BM_bernoulli(membership_type = "SBM", adj = M, plotting = character(0), ncores = 1L, nodes_covariates = list(nodes = row_covariates), verbosity = 6)
sbm_cov$estimate()}, times = 3L
))

model <- BM_bernoulli(membership_type = "LBM", adj = M, plotting = character(0), ncores = 1L, nodes_covariates = list(row = row_covariates, col = col_covariates), verbosity = 6)

model$estimate()
max(model$ICL, na.rm = TRUE)

## SBM gaussien multivarié avec covariables de noeuds

Mu1 <- 4 * matrix(runif(Q * Q), Q, Q)
Mu2 <- 4 * matrix(runif(Q * Q), Q, Q)
Noise1 <- matrix(rnorm(n * n, sd = 1), n, n)
Noise2 <- matrix(rnorm(n * n, sd = 1), n, n)
M1_gaussian <- Z1 %*% Mu1 %*% t(Z1) + Noise1
M2_gaussian <- Z1 %*% Mu2 %*% t(Z1) + 0.5 * Noise1 + Noise2

sbm_gaussian_multi_cov <- BM_gaussian_multivariate(
    membership_type = "SBM",
    adj = list(M1_gaussian, M2_gaussian),
    plotting = character(0),
    ncores = 1L,
    nodes_covariates = list(nodes = row_covariates),
    verbosity = 6
)

sbm_gaussian_multi_cov$estimate()
max(sbm_gaussian_multi_cov$ICL, na.rm = TRUE)

## SBM_sym gaussien multivarié avec covariables de noeuds

Mu1_sym <- 4 * matrix(runif(Q * Q), Q, Q)
Mu2_sym <- 4 * matrix(runif(Q * Q), Q, Q)
Mu1_sym[lower.tri(Mu1_sym)] <- t(Mu1_sym)[lower.tri(Mu1_sym)]
Mu2_sym[lower.tri(Mu2_sym)] <- t(Mu2_sym)[lower.tri(Mu2_sym)]

Noise1_sym <- matrix(rnorm(n * n, sd = 1), n, n)
Noise2_sym <- matrix(rnorm(n * n, sd = 1), n, n)
Noise1_sym[lower.tri(Noise1_sym)] <- t(Noise1_sym)[lower.tri(Noise1_sym)]
Noise2_sym[lower.tri(Noise2_sym)] <- t(Noise2_sym)[lower.tri(Noise2_sym)]

M1_gaussian_sym <- Z1 %*% Mu1_sym %*% t(Z1) + Noise1_sym
M2_gaussian_sym <- Z1 %*% Mu2_sym %*% t(Z1) + 0.5 * Noise1_sym + Noise2_sym
M1_gaussian_sym[lower.tri(M1_gaussian_sym)] <- t(M1_gaussian_sym)[lower.tri(M1_gaussian_sym)]
M2_gaussian_sym[lower.tri(M2_gaussian_sym)] <- t(M2_gaussian_sym)[lower.tri(M2_gaussian_sym)]

sbm_sym_gaussian_multi_cov <- BM_gaussian_multivariate(
    membership_type = "SBM_sym",
    adj = list(M1_gaussian_sym, M2_gaussian_sym),
    plotting = character(0),
    ncores = 1L,
    nodes_covariates = list(nodes = row_covariates),
    verbosity = 6
)

sbm_sym_gaussian_multi_cov$estimate()
max(sbm_sym_gaussian_multi_cov$ICL, na.rm = TRUE)

## LBM gaussien multivarié avec covariables de noeuds (avec Z2)

Mu1_lbm <- 4 * matrix(runif(Q * Q), Q, Q)
Mu2_lbm <- 4 * matrix(runif(Q * Q), Q, Q)
Noise1_lbm <- matrix(rnorm(n * n, sd = 1), n, n)
Noise2_lbm <- matrix(rnorm(n * n, sd = 1), n, n)

M1_gaussian_lbm <- Z1 %*% Mu1_lbm %*% t(Z2) + Noise1_lbm
M2_gaussian_lbm <- Z1 %*% Mu2_lbm %*% t(Z2) + 0.5 * Noise1_lbm + Noise2_lbm

devtools::load_all()
lbm_gauss_multi_mb <- microbenchmark(
"Cov + BFGS" = {
lbm_gaussian_multi_cov <- BM_gaussian_multivariate(
    membership_type = "LBM",
    adj = list(M1_gaussian_lbm, M2_gaussian_lbm),
    plotting = character(0),
    ncores = 1L,
    nodes_covariates = list(row = row_covariates, col = col_covariates),
    verbosity = 6
)

lbm_gaussian_multi_cov$estimate()
max(lbm_gaussian_multi_cov$ICL, na.rm = TRUE)
},
"No cov" = {
lbm_gaussian_multi <- BM_gaussian_multivariate(
    membership_type = "LBM",
    adj = list(M1_gaussian_lbm, M2_gaussian_lbm),
    plotting = character(0),
    ncores = 1L,
    verbosity = 6
)

lbm_gaussian_multi$estimate()
max(lbm_gaussian_multi$ICL, na.rm = TRUE)
}, times = 10)
