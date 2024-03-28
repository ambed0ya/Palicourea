library(RevGadgets)
library(ggplot2)

# specify the input file
file <- paste0("~/Palicourea/Inflorescence_ASE/palicourea.inflorescence_type_freeK.log")

# read the trace and discard burnin
trace_quant <- readTrace(path = file, burnin = 0.25)

# produce the plot object, showing the posterior distributions of the rates.
p <- plotTrace(trace = trace_quant, vars = paste0("rate[",1:3,"]"))[[1]] +
  # modify legend location using ggplot2
  theme(legend.position = c(0.88,0.85))

ggsave(paste0("~/Palicourea/Inflorescence_ASE/Inflorescence_ASE_rates_freeK.pdf"), p, width = 5, height = 5)


ERM Marginal likelihood
-51.12247
-51.11971

freeK Marginal likelihood
-49.78148
-49.79963