## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(randPedPCA)

## -----------------------------------------------------------------------------
# generate a pedigree object from the example dataset provided
ped <- pedigree(pedMeta$fid,
                pedMeta$mid,
                pedMeta$id)
pc01 <- rppca(ped, center=T)
plot(pc01, col = pedMeta$population)

## -----------------------------------------------------------------------------
li <- sparse2spam(getLInv(ped))
pc02 <- rppca(li, center = T)
plot(pc02, col = pedMeta$population)

## -----------------------------------------------------------------------------
pc03 <- rppca(pedLInv, center = T)
plot(pc03, col = pedMeta$population)

## ----echo=F, message=F, fig.dim = c(6, 3)-------------------------------------
set.seed(123345) # set random seed as Hutch++ uses random numbers
ttv <- hutchpp(pedLInv,num_queries = 100, center = T) # estimate
opar <- par(mfrow=c(1,2))
plot(rppca(ped), main="Un-centred", col=pedMeta$population)
plot(rppca(ped, center=T, totVar = ttv), main="Centred", col=pedMeta$population)
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
hutchpp(li,num_queries = 100, center = T) # estimate

## -----------------------------------------------------------------------------
summary(pc02)

## -----------------------------------------------------------------------------
pc04 <- rppca(ped, center=F)
summary(pc04)

## -----------------------------------------------------------------------------
pc05 <- rppca(pedLInv, center=F, totVar=3521.534)
summary(pc05)

## -----------------------------------------------------------------------------
pc07 <- rppca(ped, center=F, totVar=123)
summary(pc07)

## -----------------------------------------------------------------------------
pc06 <- rppca(ped, center=T) 
summary(pc06)

## -----------------------------------------------------------------------------
# No proportions shown by default
pc08 <- rppca(pedLInv, center=T) 
summary(pc08)

## -----------------------------------------------------------------------------
# Only when estimate is supplied
pc09 <- rppca(pedLInv, center=T, totVar=2673.6) 
summary(pc09)

## -----------------------------------------------------------------------------
# PC1 and PC2 are plotted by default
plot(pc06, col=pedMeta$population, main="My pedigree PCA")
# plot PC1 and PC3 instead
plot(pc06, dims=c(1,3), col=pedMeta$population, main="My pedigree PCA,\ncustom PCs shown")

