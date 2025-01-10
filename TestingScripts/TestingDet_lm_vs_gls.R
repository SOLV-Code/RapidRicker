# COMPARING 2 DETERMINISTIC RICKER FITS
# V1: RapidRicker fn calcDetModelFit() is currently built around lm()
# V2: Hamazaki Shiny app (deterministic version is built around gls() fro package {nlme} (which has an AR1 version)
# Deterministic Hamazaki app: https://hamachan.shinyapps.io/Spawner_Recruit/

# Reason for this test: just making sure...

# Note: Excel LINEST() gives identical results for Stock A, but different from
#       both R versions for Stock B (closer to gls version)




library(tidyverse)
library(RapidRicker)

# built in data set
gls_test_data


fits.out <- data.frame(Stock = character(), Fit = character(), ln.alpha = numeric(), ln.alpha.c = numeric(),
                       alpha = numeric(), alpha.c= numeric(), beta = numeric(), sigma= numeric())


for(stk.do in c("Stock A", "Stock B")){

sr.stk <-   gls_test_data %>% dplyr::filter(Stock == stk.do)


# RapidRicker fit
rr.fit <- calcDetModelFit(sr_obj = sr.stk,sr.scale = 1, resids = TRUE)
fits.out <- bind_rows(fits.out, bind_cols(Stock = stk.do, Fit = "RapidRicker", rr.fit$pars %>% select(-n_obs)))

# Double check the lm fit
lm.fit <-lm(sr.stk$logRpS ~ sr.stk$Spn)

lm.sigma <- sigma(lm.fit)
lm.lna <- lm.fit$coefficients[1]
lm.lna.c <- lm.lna + (lm.sigma^2  / 2)
lm.a <- exp(lm.lna)
lm.a.c <- exp(lm.lna.c)
lm.b <- - lm.fit$coefficients[2]

fits.out <- bind_rows(fits.out, bind_cols(Stock = stk.do, Fit = "lm()",
                                          ln.alpha = lm.lna, ln.alpha.c = lm.lna.c,
                                          alpha = lm.a, alpha.c= lm.a.c, beta = lm.b, sigma= lm.sigma ))



# Do the gls fit (from Hamazaki code)
require(nlme)

gls.fit <-gls(logRpS ~ Spn,data=sr.stk,method='ML')
names(gls.fit)

gls.sigma <- sigma(gls.fit)
gls.lna <- gls.fit$coefficients[1]
gls.lna.c <- gls.lna + (gls.sigma^2  / 2)
gls.a <- exp(gls.lna)
gls.a.c <- exp(gls.lna.c)
gls.b <- - gls.fit$coefficients[2]

fits.out <- bind_rows(fits.out, bind_cols(Stock = stk.do, Fit = "gls()",
                                          ln.alpha = gls.lna, ln.alpha.c = gls.lna.c,
                                          alpha = gls.a, alpha.c= gls.a.c, beta = gls.b, sigma= gls.sigma ))




# Do the gls fit with AR1 (from Hamazaki code)


gls.ar1.fit <-gls(logRpS ~ Spn,data=sr.stk,correlation=corAR1(form=~1),method='ML')
names(gls.ar1.fit)

gls.ar1.sigma <- sigma(gls.ar1.fit)
gls.ar1.lna <- gls.ar1.fit$coefficients[1]
gls.ar1.lna.c <- gls.ar1.lna + (gls.ar1.sigma^2  / 2)
gls.ar1.a <- exp(gls.ar1.lna)
gls.ar1.a.c <- exp(gls.ar1.lna.c)
gls.ar1.b <- - gls.ar1.fit$coefficients[2]

fits.out <- bind_rows(fits.out, bind_cols(Stock = stk.do, Fit = "gls() with AR1",
                                          ln.alpha = gls.ar1.lna, ln.alpha.c = gls.ar1.lna.c,
                                          alpha = gls.ar1.a, alpha.c= gls.ar1.a.c, beta = gls.ar1.b, sigma= gls.ar1.sigma ))


}

fits.out




# function test


test.lm <- calcDetModelFit(sr_obj = sr.stk,sr.scale = 1, min.obs=15,resids = FALSE, fn.use = "lm", ar1 = FALSE)
test.gls <- calcDetModelFit(sr_obj = sr.stk,sr.scale = 1, min.obs=15,resids = FALSE, fn.use = "gls", ar1 = FALSE)
test.gls.ar1 <- calcDetModelFit(sr_obj = sr.stk,sr.scale = 1, min.obs=15,resids = FALSE, fn.use = "gls", ar1 = TRUE)

test.lm$pars
test.gls$pars
test.gls.ar1$pars
