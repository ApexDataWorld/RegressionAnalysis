# Load necessary library
library(faraway)

# Fit the interaction model with prostate data
out <- lm(lpsa ~ lcavol * lweight, data = prostate)
summary(out)

# Hypothesis test: H0: β1 + 3.6β3 = 0
bethat <- out$coefficients
Vhat <- vcov(out)
cvec <- c(0, 1, 0, 3.6)
est <- t(cvec) %*% bethat
se <- sqrt(t(cvec) %*% Vhat %*% cvec)
tstat <- est / se
pv <- pt(abs(tstat), lower.tail = FALSE, df = out$df.residual)
data.frame(estimate = est, se = se, tstat = tstat, pvalue = pv)

# Reduced model
out_red <- lm(lpsa ~ lweight, data = prostate)
anova(out_red)

# Hypothesis test: H0: β1 + 6β3 = 0
cvec <- c(0, 1, 0, 6)
est <- t(cvec) %*% bethat
se <- sqrt(t(cvec) %*% Vhat %*% cvec)
tstat <- est / se
pv <- pt(abs(tstat), lower.tail = FALSE, df = out$df.residual)
data.frame(estimate = est, se = se, tstat = tstat, pvalue = pv)

# Full model for F-test
out_full <- lm(lpsa ~ ., data = prostate)
summary(out_full)

# Hypothesis test: H0: Kβ = 0
K <- rbind(
  c(0, 1, -1, rep(0, 6)),
  c(0, rep(0, 6), 1, -1)
)
bethat <- out_full$coefficients
Vhat <- vcov(out_full)
est <- K %*% bethat
M <- K %*% Vhat %*% t(K)
Fstat <- t(est) %*% solve(M) %*% est / nrow(K)
pv <- pf(Fstat, lower.tail = FALSE, df1 = nrow(K), df2 = out_full$df.residual)
data.frame(Fstat = Fstat, pvalue = pv)

# Reduced model for F-test
prostate$X1 <- prostate$lcavol + prostate$lweight
prostate$X2 <- prostate$gleason + prostate$pgg45
out_red <- lm(lpsa ~ . - lcavol - lweight - gleason - pgg45, data = prostate)
summary(out_red)

# ANOVA to compare full and reduced models
anova(out_red, out_full)
