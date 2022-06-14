library(mclust)
library(mixtools)
library(stats)
library(modeest)
library(ggplot2)
library(ggprism)
library(magrittr)

load(file = "CorData.RData")

dat <- CorData[,c(1,3)]
colnames(dat) <- c("CNV_Protein", "CNV_mRNA")
dat <- dat[complete.cases(dat),]
atten <- dat$CNV_mRNA - dat$CNV_Protein


mod2 <- Mclust(atten, G = 2, modelNames = "E")


dat$Class <- as.factor(mod2$classification)
dat$Group <- ifelse(dat$Class==1, "High", "Low")

ggplot(data = dat, aes(x = CNV_mRNA, y = CNV_Protein, color = Group)) +
  geom_point() +
  stat_density2d(data = dat, aes(x = .data$CNV_mRNA, y = .data$CNV_Protein, group=.data$Class),
                 color = "black", lwd = 1) +
  scale_x_continuous(
    limits = c(-1, 1)
  ) +
  scale_y_continuous(
    limits = c(-1, 1)
  ) +
  scale_color_manual(
    values = c("red", "blue")
  )+
  theme_prism() +
  theme(
    panel.background = element_rect(color = "black"),
    legend.position = c(0.1, 0.9)
  )


# two component Gaussian mixture model
# https://www.biostars.org/p/194236/
load("dat4.RData")
sampleExpr <- log2(dat4[[1]]$Expression$plex1_126)
b<-bw.SJ(sampleExpr) # Find bandwith with Shafer-Jones method
d<-density(sampleExpr, bw=b, kernel="gaussian") # Gaussian kernel density estimate
mode<-d$x[which.max(d$y)]
model<-normalmixEM(sampleExpr,mu=c(mode,mode),sigma=NULL)
sigma = min(model$sigma)


#' Plot a Mixture Component
#' 
#' @param x Input data
#' @param mu Mean of component
#' @param sigma Standard deviation of component
#' @param lam Mixture weight of component
plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}

set.seed(1)
mixmdl <- model
## number of iterations= 29
data.frame(x = mixmdl$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), binwidth = 0.1, colour = "black", 
                 fill = "white") +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                colour = "red", lwd = 1) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                colour = "blue", lwd = 1) +
  ylab("Density")



