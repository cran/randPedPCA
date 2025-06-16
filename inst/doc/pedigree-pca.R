## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(randPedPCA)

## -----------------------------------------------------------------------------
# generate a pedigree object from the example dataset provided
ped <- pedigree(pedMeta2$fid,
                pedMeta2$mid,
                pedMeta2$id)
pc01 <- rppca(ped, center=TRUE)
plot(pc01, col = factor(pedMeta2$population))

## -----------------------------------------------------------------------------
li <- sparse2spam(getLInv(ped))
pc02 <- rppca(li, center = TRUE)
plot(pc02, col = factor(pedMeta2$population))

## -----------------------------------------------------------------------------
pc03 <- rppca(pedLInv2, center = TRUE)
plot(pc03, col = factor(pedMeta2$population))

## ----echo=F, message=F, fig.dim = c(6, 3)-------------------------------------
set.seed(123345) # set random seed as Hutch++ uses random numbers
ttv <- hutchpp(pedLInv2,num_queries = 100, center = TRUE) # estimate
opar <- par(mfrow=c(1,2))
plot(rppca(ped), main="Un-centred", col=factor(pedMeta2$population))
plot(rppca(ped, center=TRUE, totVar = ttv), main="Centred", col=factor(pedMeta2$population))
par(opar)

## -----------------------------------------------------------------------------
# True total variance computed via the pedigree's inbreeding values
sum(inbreeding(ped) + 1) # "ped" was defined at the top

# Now estimate the total variance (the trace of A) from the corresponding
#  L inverse matrix using Hutch++
li <- sparse2spam(getLInv(ped)) # generate L inverse and convert to spam format
set.seed(123345) # set random seed as Hutch++ uses random numbers
hutchpp(li) # Hutch++ with default settings

# for higher accuracy increase num_queries (increases running time)
hutchpp(li,num_queries = 100)

## -----------------------------------------------------------------------------
# Get L, the "data matrix" of a pedigree
ll <- getL(ped)
# centre L (because L is upper triangular, we centre the rows)
llc <- apply(ll, 1, function(x) x - mean(x))
# compute additive relationship matrix of the centred data
ac <- llc %*% t(llc)
sum(diag(ac)) # exact value (would be too expensive to compute for a large pedigree)

# Obtain centred estimate from L inverse using Hutch++
li <- sparse2spam(getLInv(ped)) # generate L inverse and convert to spam format
set.seed(123345) # set random seed as Hutch++ uses random numbers
hutchpp(li,num_queries = 100, center = TRUE) # estimate

## -----------------------------------------------------------------------------
summary(pc02)

## -----------------------------------------------------------------------------
pc04 <- rppca(ped, center=FALSE)
summary(pc04)

## -----------------------------------------------------------------------------
pc05 <- rppca(pedLInv2, center=FALSE, totVar=3521.534)
summary(pc05)

## -----------------------------------------------------------------------------
pc07 <- rppca(ped, center=FALSE, totVar=123)
summary(pc07)

## -----------------------------------------------------------------------------
pc06 <- rppca(ped, center=TRUE) 
summary(pc06)

## -----------------------------------------------------------------------------
# No proportions shown by default
pc08 <- rppca(pedLInv2, center=TRUE) 
summary(pc08)

## -----------------------------------------------------------------------------
# Only when estimate is supplied
pc09 <- rppca(pedLInv2, center=TRUE, totVar=2673.6) 
summary(pc09)

## -----------------------------------------------------------------------------
# PC1 and PC2 are plotted by default
plot(pc06, col=factor(pedMeta2$population), main="My pedigree PCA")
# plot PC1 and PC3 instead
plot(pc06, dims=c(1,3), col=factor(pedMeta2$population), main="My pedigree PCA,\ncustom PCs shown")

