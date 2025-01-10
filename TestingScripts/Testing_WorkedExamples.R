library(RapidRicker)

test.df <- data.frame(ln.alpha = log(seq(1.5,8.5,by=0.5)),beta = 10^(-05),sigma = 0.5 )
test.df


test.df <- calcRickerSmsy(test.df, method = "Scheuerell2016", sr.scale = 1, out.type = "Full")
test.df

test.df <- calcRickerSgen(test.df, method = "Connorsetal2022",sr.scale = 1, out.type = "Full",tracing = FALSE)
test.df










