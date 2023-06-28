library(fitdistrplus)
library(ggplot2)

data <- data.matrix(HP_testverb_Ttest_Cz_vals)
data_vector = as.vector(data)
isNum = is.numeric(data_vector)
descdist(data_vector, discrete = FALSE)
denscomp(ft)