## devtools::install_github("alicebalard/parasiteLoad", ref = "v2.1")
library(parasiteLoad)
library(optimx)
library(fitdistrplus)
library(parallel)

library(ggplot2)
library(ggplotify)
library(pheatmap)
library(patchwork)

library(reshape2)

# for reproducibility
set.seed(8)

BAIdata <- read.csv("https://raw.githubusercontent.com/alicebalard/parasiteLoad/master/data/WATWMdata.csv", na.strings = c("", " ", NA))
# Keep individuals with hybrid index, sex recorded and pinworms counted
BAIdata <- BAIdata[!is.na(BAIdata$HI) & !is.na(BAIdata$Sex) & !is.na(BAIdata$Aspiculuris.Syphacia),]

# Balard et al. 2020
BALdata <- read.csv("https://raw.githubusercontent.com/alicebalard/Article_IntensityEimeriaHMHZ/master/data/cleanedData.csv", na.strings = c("", " ", NA))
# Keep individuals with hybrid index, sex recorded and pinworms counted
BALdata <- BALdata[!is.na(BALdata$HI) & !is.na(BALdata$Sex) & !is.na(BALdata$Aspiculuris_Syphacia),]


simulateHybrid <- function(size, L1, L2, alpha, theta){
    hybind  <-  sample(c(BAIdata$HI, BALdata$HI), size)
    simdata  <-  data.frame(ID = 1:size, HI =hybind,
                            Sex = factor(rep(c("female", "male"), size/2)),
                            meanHere = parasiteLoad::MeanLoad(
                                                         L1 = L1, L2 = L2,
                                                         alpha = alpha,
                                                         hybridIndex =hybind))
    simdata$load <- sapply(simdata$meanHere, function(x){
        rnegbin(n = 1, mu = x, theta = theta)
    })
    simdata
}



cs <- function (df){
    ## ## THE TESTS ############
    ## Test 1: chi2 wormy VS non wormy, cut at 250, hybrids 12.5% < HI < 87.5%
    ## 250 is the ~ 0.9 quantile, better doing 0.8 for more power
    chisq.test(df$load>quantile(df$load,0.8),
               df$HI<0.125|df$HI>0.875)$p.value
}


kw <- function(df){
    df$genotype <- "Hybrid"
    df$genotype[df$HI < 0.2] <- "Mmd"
    df$genotype[df$HI > 0.6] <- "Mmm"
    kruskal.test(load ~ genotype, data = df)$p.value
}


ml <- function (df) {
    ML <- parasiteLoad::analyse(data = df, response = "load", model = "negbin",
                                group = "Sex", hybridIndex = "HI")
    ML$H0$Gtest$pvalue
}

alphaV <- seq(0, 1, by=0.1)
sizes <- seq(100, 1000, by=100)

power.cs <- mclapply(sizes,  function(s) {
    alphas <- sapply(alphaV, function(x) {
        p.vals <- mclapply(seq(100), function (i) {
            df <- simulateHybrid(s, 100, 100, x, 0.8)
            cs(df)
        }, mc.cores=12)
        length(p.vals[p.vals<0.05])
    })
}, mc.cores=10)


power.kw <- mclapply(sizes,  function(s) {
    alphas <- sapply(alphaV, function(x) {
        p.vals <- mclapply(seq(100), function (i) {
            df <- simulateHybrid(s, 100, 100, x, 0.8)
            kw(df)
        }, mc.cores=12)
        length(p.vals[p.vals<0.05])
    })
}, mc.cores=10)


power.ml <- mclapply(sizes,  function(s) {
    alphas <- sapply(alphaV, function(x) {
        ## remember to change this to 100 for final evaluation
        p.vals <- mclapply(seq(100), function (i) {
            df <- simulateHybrid(s, 100, 100, x, 0.8)
            ml(df)
        }, mc.cores=24)
        length(p.vals[p.vals<0.05])
    })
}, mc.cores=10)

## quick save in my home to avoid re-doing this
saveRDS(list(power.ml, power.kw, power.cs),
        "~/review_power.Rds")

MLC <- do.call(rbind, power.cs)/100
colnames(MLC) <- alphaV
rownames(MLC) <- sizes

MLK <- do.call(rbind, power.kw)/100
colnames(MLK) <- alphaV
rownames(MLK) <- sizes


MLP <- do.call(rbind, power.ml)/100
colnames(MLP) <- alphaV
rownames(MLP) <- sizes

MLC.pm <- as.ggplot(
    pheatmap(MLC, cluster_rows=F, cluster_cols=F,
             main="Power Chisq test",
             legend_breaks = seq(0, 1, by=0.2),
             legend=FALSE
             )
) 


MLK.pm <- as.ggplot(
    pheatmap(MLK, cluster_rows=F, cluster_cols=F,
             main="Power Kruskall-Wallis test",
             legend_breaks = seq(0, 1, by=0.2),
             legend=FALSE)
)

MLP.pm <- as.ggplot(
    pheatmap(MLP, cluster_rows=F, cluster_cols=F,
             main="Power Maximum Likelihood estimate",
             legend_breaks = seq(0, 1, by=0.2),
             legend_labels = seq(0, 1, by=0.2))
)


speci <- cbind("Chisq"=MLC[,1], "Kr-Wa"=MLK[,1], "ML"=MLP[,1])
speci <-  melt(speci)


FP <- ggplot(speci, aes(Var1, value, color=Var2)) +
    geom_point() + 
    scale_y_continuous("False positive proportion")+
    scale_x_continuous("Sample size")+
    geom_line()+    
    theme_minimal() +
    theme(legend.position = "none") 

powerT <- rbind(
    cbind("Chisq"=MLC[c(1), -1], "Kr-Wa"=MLK[c(1), -1], "ML"=MLP[c(1), -1]),
    cbind("Chisq"=MLC[c(6), -1], "Kr-Wa"=MLK[c(6), -1], "ML"=MLP[c(6), -1])
)

powerT <-  reshape2::melt(powerT)
powerT$Ssize <- rep(rep(c(100, 600), each=10), times=3)

Pow <- ggplot(powerT, aes(Var1, value, color=Var2)) +
    geom_point() + 
    scale_y_continuous("Power (prop. true positives)")+
    scale_x_continuous("Effect size")+
    facet_wrap(~Ssize) +
    geom_line()+
    theme_minimal()


pdf("power_Fig.pdf", width=12, height=8)
(MLC.pm + MLK.pm + MLP.pm)/(FP + Pow) +
      plot_annotation(tag_levels = 'a')
dev.off()    
