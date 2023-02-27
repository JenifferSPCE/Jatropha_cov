####################################################################### Repeated Measures Data ########################################################################
########################################################################################################################################################################

##>> Import dataset
dataset <- read.table("C:\\Users\\Jeniffer\\dataset.txt", header = TRUE)

##>> Transform numeric variables in factors
dataset <- transform(dataset, Measurements = factor(Measurements), Plots = factor(Plots), Replications = factor(Replications), Progenies = factor(Progenies), ID = factor(ID))

##>> Loading needed packages
require(asreml)
require(Matrix)

##>> AIC & BIC
calc <- function (asreml){
  summ <- summary(asreml)
  vc <- summ$varcomp
  p <- nrow(vc)
  if ("Fixed" %in% levels(vc$constraint))
    p <- p - table(vc$constraint)["Fixed"]
  if ("Constrained" %in% levels(vc$constraint))
    p <- p - table(vc$constraint)["Constrained"]
  names(p) <- ""
  logREML <- summ$loglik
  AIC <- -2 * logREML + 2 * p
  BIC <- -2 * logREML + p * log(summ$nedf)
  data.frame(p, AIC, BIC, logREML)
}

############################################################################## Modelling the residual effects ################################################################################
##>> Models
m1 <- asreml(fixed = yield ~ Measurements + Measurements:Replications,
             random = ~ ID + Progenies + Measurements:Progenies + Plots,
             residual = ~ Measurements:ID,                        
             data = dataset,
             maxiter = 10000)

m2 <- asreml(fixed = yield ~ Measurements + Measurements:Replications,
             random = ~ ID + Progenies + Measurements:Progenies + Plots,
             residual = ~ dsum(~id(units)|Measurements),                       
             data = dataset,
             maxiter = 10000)

m3 <- asreml(fixed = yield ~ Measurements + Measurements:Replications,
             random = ~ ID + Progenies + Measurements:Progenies + Plots,
             residual = ~ ar1h(Measurements):ID,                       
             data = dataset,
             maxiter = 10000)

rbind(calc(m1), calc(m2), calc(m3))


############################################################################ Modelling the permanent effects ###############################################################################
##>> Models
m3 <- asreml(fixed = yield ~ Measurements + Measurements:Replications,
             random = ~ ID + Progenies + Measurements:Progenies + Plots,
             residual = ~ ar1h(Measurements):ID,                       
             data = dataset,
             maxiter = 10000)

m4 <- asreml(fixed = yield ~ Measurements + Measurements:Replications,
             random = ~ at(Measurements):ID + Progenies + Measurements:Progenies + Plots,
             residual = ~ ar1h(Measurements):ID,                        
             data = dataset,
             maxiter = 10000)

rbind(calc(m1), calc(m2), calc(m3), calc(m4))

########################################################################## Modelling the plots effects ##############################################################################
##>> Models
m5 <-asreml(fixed = yield ~ Measurements + Measurements:Replications,
            random = ~ ID:idh(Measurements) + Progenies + Measurements:Progenies + at(Measurements):Plots,
             residual = ~ ar1h(Measurements):ID,                        
             data = dataset,
              maxiter = 10000)

m6 <- asreml(fixed = yield ~ Measurements + Measurements:Replications,
             random = ~ ID:idh(Measurements) + Progenies + Measurements:Progenies + corh(Measurements):Plots,
             residual = ~ ar1h(Measurements):ID,                        
             data = dataset,
             maxiter = 10000)

m7 <- asreml(fixed = yield ~ Measurements + Measurements:Replications,
             random = ~ ID:idh(Measurements) + Progenies + Measurements:Progenies + ar1h(Measurements):Plots,
             residual = ~ ar1h(Measurements):ID,                        
             data = dataset,
             maxiter = 10000)

rbind(calc(m1), calc(m2), calc(m3), calc(m4), calc(m5), calc(m6), calc(m7))

######################################################################## Modelling the Progenies effects ############################################################################

m8 <- asreml(fixed = yield ~ Measurements + Measurements:Replications,
             random = ~ ID:idh(Measurements) + at(Measurements):Progenies + corh(Measurements):Plots,
             residual = ~ ar1h(Measurements):ID,                     
             data = dataset,
             maxiter = 10000)

m9 <- asreml(fixed = yield ~ Measurements + Measurements:Replications,
             random = ~ ID:idh(Measurements)+ corh(Measurements):Progenies + corh(Measurements):Plots,
             residual = ~ ar1h(Measurements):ID,                        
             data = dataset,
             maxiter = 10000)

m10 <- asreml(fixed = yield ~ Measurements + Measurements:Replications,
              random = ~ ID:idh(Measurements) + fa(Measurements,1):Progenies + corh(Measurements):Plots,
              residual = ~ ar1h(Measurements):ID,                        
              data = dataset,
              maxiter = 10000)

rbind(calc(m1), calc(m2), calc(m3), calc(m4), calc(m5), calc(m6), calc(m7), calc(m8), calc(m9), calc(m10))

summary(m10)$varcomp
x=coef(m10)
effects=x$random

#BLUP
pred = predict(m10, classify = 'Progenies', vcov = T)
predvalor = pred$pvals

write.table(effects,"C:\\Users\\Jeniffer\\Desktop\\effects.txt")



####################################################################################
