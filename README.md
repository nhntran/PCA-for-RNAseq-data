Quick And Easy Pipeline for Principal Component Analysis on RNASeq Data using R
================
Tran Nguyen
================

This is the quick and simple pipeline for PCA on RNAseq data using R. Feel free to use the codes for your analysis. If you want to learn more about PCA, check out the excellent "StatQuest with Josh Starmer" video series on Youtube. (<https://www.youtube.com/watch?v=FgakZw6K1QQ>) The materials below were generated based on the "StatQuest with Josh Starmer" tutorials, with some modifications by me.

### Contents:
-   [1. Get and set working directory](#get-and-set-working-directory)
-   [2. Generate the example data](#generate-the-example-data)
-   [3. Important tips on PCA](#important-tips-on-pca)
-   [4. Perform PCA with prcomp()](#perform-pca-with-prcomp)
-   [5. Plot the PCA](#plot-the-pca)
-   [6. Plot results](#plot-results)



## 1. Get and set working directory

``` r
### See the current directory
getwd() # =>[1] "/Users/trannguyen"
setwd("/Users/trannguyen/R_tutorials/PCA_RNASeq/") 
#setwd("../") #=> go back to the parent directory: [1] "/Users/trannguyen"
rm(list=ls()) #Clear the workspace before working to keep things running smoothly

#Loading library
library(ggplot2) 
```

## 2. Generate the example data


-   Reference: StatQuest with Josh Starmer - Youtube

``` r
#The material is from StatQuest: PCA in R - Joh Starmer
#https://www.youtube.com/watch?v=0Jp4gsfOLMs&t=190s

#example of expression data containing 6 samples and 100 genes
# 3 wildtype (wt) and 3 knock-out (ko) samples
edata<-matrix(nrow=100, ncol=6)
colnames(edata)<-c(paste("wt",1:3,sep=""), paste("ko",1:3,sep="")) 
rownames(edata)<-paste("gene",1:100,sep="")
for (i in 1:100){
  wt.values <- rpois(3, lambda=sample(x=10:1000, size=1))
  ko.values <- rpois(3, lambda=sample(x=10:1000, size=1))
  edata[i,] <- c(wt.values, ko.values)
}
```

## 3. Important tips on PCA

1.  Scale all the variables so that all the variables are roughly equivalent to avoid bias towards one of them. Standard practice for scaling: dividing each variable by its standard deviation.
2.  Centering the data.
3.  In theory, there is 1 PC for each variable. But if the number of samples &lt; number of variables,the number of samples puts an upper bound on the number of PCs that have eigenvalue&gt;0.

-   The 1st PC: accounts for the most variation in the data
-   The 2nd PC: accounts for the second most variation, and so on...

## 4. Perform PCA with prcomp()

``` r
pca <-prcomp(t(edata), scale=TRUE)
#prcomp() expects the samples (variables) to be rows, genes to be columns
#=> need to transpose the matrix using t()
#scale = TRUE => center and scale the data, cannot be used if there are zero or constant 

summary(pca)
#Output of prcomp():
#1. PC1, 2, ...: pca$x[,1], pca$x[,2], ...
#2. standard deviation: pca$sdev
#3. the loading scores: Loading score1, 2, ... pca$rotation[,1], pca$rotation[,2], ...
# big absolute loading scores: mainly responsible for pushing samples to either direction
```

## 5. Plot the PCA


``` r
#Setting up the color for the plot
mycolorcollection=  c('dodgerblue3', 'cornflowerblue', 'aquamarine2','limegreen', 'yellow', 'hotpink','darkorange','salmon','chocolate4','darkorchid4','gray32')
palette(mycolorcollection)
par(pch = 19) #set the character to be a filled dot

#plot the principal components - Plot1 - PCA_RNAseq_plot1.png
plot(pca$x[,1], pca$x[,2], col=2, xlab="PC1", ylab="PC2")

#plot all the PCs using Scree Plot - Plot2 - PCA_RNAseq_plot2.png
pca.var <- pca$sdev^2
pcs <- round(pca.var/sum(pca.var)*100,1) #calculate percent variation
barlabel<-pcs
pcs<-t(matrix(pcs)) #convert it to matrix and transpose
colnames(pcs)<-paste("PC",1:6,sep="")
x<-barplot(pcs, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation",
        ylim=c(0,100), col=8)
text(x,barlabel+2.5,format(barlabel,T))

#fancy clustering plot - Plot3 - PCA_RNAseq_plot3.png
#creating the pcadata
pca.data<-data.frame(Sample=rownames(pca$x), X=pca$x[,1],Y=pca$x[,2])

ggplot(data=pca.data, aes(x=X, y=Y, label=Sample)) + geom_text() + 
  xlab(paste("PC1 - ",barlabel[1],"%", sep="")) +
  ylab(paste("PC2 - ", barlabel[2], "%", sep="")) +
  theme_bw() + ggtitle("The PCA Graph")
```

## 6. Plot results


#### Plot 1:

By plotting the first 2 PCs, the data is shown to be clustered into 2 separate groups.

<p align="center">
  <img src="./img/PCA_RNAseq_plot1.png" alt="Size Limit CLI">
</p>


#### Plot 2:

According to the plot, the PC1 accounts for the most variation in the data (90.4%)

<p align="center">
  <img src="./img/PCA_RNAseq_plot2.png" alt="Size Limit CLI">
</p>

#### Plot 3:

By using ggplot2, we can clearly demonstrate which sample belong to each cluster/group.

<p align="center">
  <img src="./img/PCA_RNAseq_plot3.png" alt="Size Limit CLI">
</p>

This document was processed on: 2019-01-15.

``` r
sessionInfo() 
```

    ## R version 3.5.1 (2018-07-02)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS High Sierra 10.13.6
    ## 
    ## Matrix products: default
    ## BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] compiler_3.5.1  magrittr_1.5    tools_3.5.1     htmltools_0.3.6
    ##  [5] yaml_2.2.0      Rcpp_1.0.1      stringi_1.4.3   rmarkdown_1.12 
    ##  [9] knitr_1.22      stringr_1.4.0   xfun_0.6        digest_0.6.18  
    ## [13] evaluate_0.13

``` r
#devtools::session_info()
```
