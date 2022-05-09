# COMPARING 2 DETERMINISTIC RICKER FITS
# V1: RapidRicker fn calcDetModelFit() is currently built around lm()
# V2: Hamazaki Shiny app (deterministic version is built around gls() fro package {nlme} (which has an AR1 version)
# Deterministic Hamazaki app: https://hamachan.shinyapps.io/Spawner_Recruit/

# Reason for this test: identical results for stock A, very different fits for Stock B.

# Note: Excel LINEST() gives identical results for Stock A, but different from
#       both R versions for Stock B (closer to gls version)




library(tidyverse)
library(RapidRicker)

# built in data set
gls.test.data


fits.out <- data.frame(Stock = character(), Fit = character(), ln.alpha = numeric(), ln.alpha.c = numeric(),
                       alpha = numeric(), alpha.c= numeric(), beta = numeric(), sigma= numeric())


for(stk.do in c("Stock A", "Stock B")){

sr.stk <-   gls.test.data %>% dplyr::filter(Stock == stk.do)


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


}

fits.out
