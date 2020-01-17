library(ggplot2)
library(relater)

df <- data.frame()

df <- rbind(df, cbind(read.coal("./example_relate.coal"), type = "Relate trees"))
df   <- rbind(df, cbind(read.coal("./example_true.coal"), type = "True trees"))
df  <- rbind(df, cbind(read.coal("./data/example_truth.coal"), type = "Truth"))

p <- ggplot(df) + geom_step(aes(x = epoch.start, y = 0.5/haploid.coalescence.rate, color = type, group = paste(group1, group2,type))) + scale_x_continuous(trans = "log10", limit = c(1e2, 1e7)) + scale_y_continuous(trans = "log10", limit = c(1e2, 1e7)) + theme_bw()

ggsave(p, file = "plot.pdf", width = 5, height = 5)

