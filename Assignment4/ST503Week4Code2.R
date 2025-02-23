# Load necessary libraries
library(faraway)
library(multcomp)

# Fit the model with prostate data
out <- lm(lpsa ~ ., data = prostate)
summary(out)

# 95% Confidence Interval for β_lcavol
alpha <- 0.05
bhat <- out$coefficients[2]
se <- sqrt(vcov(out)[2,2])
tcrit <- qt(alpha/2, lower.tail = FALSE, df = out$df.residual)
CI <- bhat + c(-1,1)*tcrit*se
data.frame(estimate = bhat, se = se, lower = CI[1], upper = CI[2])

# Built-in confidence interval
confint(out, "lcavol", level = 0.95)

# Confidence interval for β_lcavol - β_lweight
cvec <- c(0, 1, -1, rep(0, 6))
bethat <- out$coefficients
Vhat <- vcov(out)
est <- c(t(cvec) %*% bethat)
se <- c(sqrt(t(cvec) %*% Vhat %*% cvec))
tcrit <- qt(alpha/2, lower.tail = FALSE, df = out$df.residual)
CI <- est + c(-1,1)*tcrit*se
data.frame(estimate = est, se = se, critical = tcrit, lower = CI[1], upper = CI[2])

# Using glht() for the same interval
K <- matrix(cvec, nrow = 1)
rownames(K) <- "lcavol - lweight"
out.hyp <- glht(out, linfct = K)
confint(out.hyp)

# Prediction interval for lcavol = 2
out_pred <- lm(lpsa ~ lcavol, data = prostate)
xnew <- c(1, 2)
bethat <- coef(out_pred)
Vhat <- vcov(out_pred)
sigsq <- (sigma(out_pred))^2
pred <- c(t(xnew) %*% bethat)
pse <- sqrt(t(xnew) %*% Vhat %*% xnew + sigsq)
tcrit <- qt(alpha/2, lower.tail = FALSE, df = out_pred$df.residual)
PI <- pred + c(-1,1)*tcrit*pse
data.frame(prediction = pred, se = pse, critical = tcrit, lower = PI[1], upper = PI[2])

# Using predict() function
data.frame(predict(out_pred, newdata = data.frame(lcavol = 2), se.fit = TRUE, interval = "prediction"))
data.frame(predict(out_pred, newdata = data.frame(lcavol = 2), se.fit = TRUE, interval = "confidence"))
