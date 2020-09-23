########### dynamic linear model for analyzing fishery policy  #########################
########################################################################################
#### install packages
install.packages('zoo')
install.packages('lm')
install.packages("tseries")
#### build time-lag function
L<-function(x,shift_by){
   stopifnot(is.numeric(shift_by))
   if(shift_by==0){
      out<-x
   } else {
      out<-c(rep(NA,shift_by), head(x,-shift_by))
   }
   out
}
#### build stepAICc function based on https://github.com/biometry/APES/blob/master/Data/Dormann2013/stepAICc.r
extractAICc=function (fit, scale, k = 2, ...) 
{
   res <- logLik(fit)
   
   edf <- attr(res, "df")
   
   n=length(residuals(fit))	
   
   c(edf,-2 * res + k * edf+2*edf*(edf+1)/(n-edf-1))	
   
}
dropterm.AICc=function (object, scope, scale = 0, test = c("none", "Chisq"), k = 2, sorted = FALSE, trace = FALSE, ...) 
{
   
   tl <- attr(terms(object), "term.labels")
   
   if (missing(scope)) 
      
      scope <- drop.scope(object)
   
   else {
      
      if (!is.character(scope)) 
         
         scope <- attr(terms(update.formula(object, scope)), 
                       
                       "term.labels")
      
      if (!all(match(scope, tl, 0L))) 
         
         stop("scope is not a subset of term labels")
      
   }
   
   ns <- length(scope)
   
   ans <- matrix(nrow = ns + 1L, ncol = 2L, dimnames = list(c("<none>", 
                                                              
                                                              scope), c("df", "AIC")))
   
   ans[1, ] <- extractAICc(object, scale, k = k, ...)
   
   env <- environment(formula(object))
   
   n0 <- length(object$residuals)
   
   for (i in seq(ns)) {
      
      tt <- scope[i]
      
      if (trace) {
         
         message("trying -", tt)
         
         utils::flush.console()
         
      }
      
      nfit <- update(object, as.formula(paste("~ . -", tt)), 
                     
                     evaluate = FALSE)
      
      nfit <- eval(nfit, envir = env)
      
      ans[i + 1, ] <- extractAICc(nfit, scale, k = k, ...)
      
      if (length(nfit$residuals) != n0) 
         
         nfit$residuals<-na.omit(nfit$residuals)
      
   }
   
   dfs <- ans[1L, 1L] - ans[, 1L]
   
   dfs[1L] <- NA
   
   aod <- data.frame(Df = dfs, AIC = ans[, 2])
   
   o <- if (sorted) 
      
      order(aod$AIC)
   
   else seq_along(aod$AIC)
   
   test <- match.arg(test)
   
   if (test == "Chisq") {
      
      dev <- ans[, 2L] - k * ans[, 1L]
      
      dev <- dev - dev[1L]
      
      dev[1L] <- NA
      
      nas <- !is.na(dev)
      
      P <- dev
      
      P[nas] <- safe_pchisq(dev[nas], dfs[nas], lower.tail = FALSE)
      
      aod[, c("LRT", "Pr(Chi)")] <- list(dev, P)
      
   }
   
   aod <- aod[o, ]
   
   head <- c("Single term deletions", "\nModel:", deparse(as.vector(formula(object))))
   
   if (scale > 0) 
      
      head <- c(head, paste("\nscale: ", format(scale), "\n"))
   
   class(aod) <- c("anova", "data.frame")
   
   attr(aod, "heading") <- head
   
   aod
   
}
addterm.AICc=function (object, scope, scale = 0, test = c("none", "Chisq"),  k = 2, sorted = FALSE, trace = FALSE, ...) 
{
   
   if (missing(scope) || is.null(scope)) 
      
      stop("no terms in scope")
   
   if (!is.character(scope)) 
      
      scope <- add.scope(object, update.formula(object, scope))
   
   if (!length(scope)) 
      
      stop("no terms in scope for adding to object")
   
   ns <- length(scope)
   
   ans <- matrix(nrow = ns + 1L, ncol = 2L, dimnames = list(c("<none>", 
                                                              
                                                              scope), c("df", "AIC")))
   
   ans[1L, ] <- extractAICc(object, scale, k = k, ...)
   
   n0 <- length(object$residuals)
   
   env <- environment(formula(object))
   
   for (i in seq(ns)) {
      
      tt <- scope[i]
      
      if (trace) {
         
         message("trying +", tt)
         
         utils::flush.console()
         
      }
      
      nfit <- update(object, as.formula(paste("~ . +", tt)), 
                     
                     evaluate = FALSE)
      
      nfit <- eval(nfit, envir = env)
      
      ans[i + 1L, ] <- extractAICc(nfit, scale, k = k, ...)
      
      if (length(nfit$residuals) != n0) 
         
         stop("number of rows in use has changed: remove missing values?")
      
   }
   
   dfs <- ans[, 1L] - ans[1L, 1L]
   
   dfs[1L] <- NA
   
   aod <- data.frame(Df = dfs, AIC = ans[, 2L])
   
   o <- if (sorted) 
      
      order(aod$AIC)
   
   else seq_along(aod$AIC)
   
   test <- match.arg(test)
   
   if (test == "Chisq") {
      
      dev <- ans[, 2L] - k * ans[, 1L]
      
      dev <- dev[1L] - dev
      
      dev[1L] <- NA
      
      nas <- !is.na(dev)
      
      P <- dev
      
      P[nas] <- safe_pchisq(dev[nas], dfs[nas], lower.tail = FALSE)
      
      aod[, c("LRT", "Pr(Chi)")] <- list(dev, P)
      
   }
   
   aod <- aod[o, ]
   
   head <- c("Single term additions", "\nModel:", deparse(as.vector(formula(object))))
   
   if (scale > 0) 
      
      head <- c(head, paste("\nscale: ", format(scale), "\n"))
   
   class(aod) <- c("anova", "data.frame")
   
   attr(aod, "heading") <- head
   
   aod
   
}
stepAICc=function (object, scope, scale = 0, direction = c("both", "backward", "forward"), trace = 1, keep = NULL, steps = 1000, use.start = FALSE, k = 2, ...) 
{
   
   mydeviance <- function(x, ...) {
      
      dev <- deviance(x)
      
      if (!is.null(dev)) 
         
         dev
      
      else extractAICc(x, k = 0)[2L]
      
   }
   
   cut.string <- function(string) {
      
      if (length(string) > 1L) 
         
         string[-1L] <- paste("\n", string[-1L], sep = "")
      
      string
      
   }
   
   re.arrange <- function(keep) {
      
      namr <- names(k1 <- keep[[1L]])
      
      namc <- names(keep)
      
      nc <- length(keep)
      
      nr <- length(k1)
      
      array(unlist(keep, recursive = FALSE), c(nr, nc), list(namr, 
                                                             
                                                             namc))
      
   }
   
   step.results <- function(models, fit, object, usingCp = FALSE) {
      
      change <- sapply(models, "[[", "change")
      
      rd <- sapply(models, "[[", "deviance")
      
      dd <- c(NA, abs(diff(rd)))
      
      rdf <- sapply(models, "[[", "df.resid")
      
      ddf <- c(NA, abs(diff(rdf)))
      
      AIC <- sapply(models, "[[", "AIC")
      
      heading <- c("Stepwise Model Path \nAnalysis of Deviance Table", 
                   
                   "\nInitial Model:", deparse(as.vector(formula(object))), 
                   
                   "\nFinal Model:", deparse(as.vector(formula(fit))), 
                   
                   "\n")
      
      aod <- if (usingCp) 
         
         data.frame(Step = change, Df = ddf, Deviance = dd, 
                    
                    `Resid. Df` = rdf, `Resid. Dev` = rd, Cp = AIC, 
                    
                    check.names = FALSE)
      
      else data.frame(Step = change, Df = ddf, Deviance = dd, 
                      
                      `Resid. Df` = rdf, `Resid. Dev` = rd, AIC = AIC, 
                      
                      check.names = FALSE)
      
      attr(aod, "heading") <- heading
      
      class(aod) <- c("Anova", "data.frame")
      
      fit$anova <- aod
      
      fit
      
   }
   
   Terms <- terms(object)
   
   object$formula <- Terms
   
   if (inherits(object, "lme")) 
      
      object$call$fixed <- Terms
   
   else if (inherits(object, "gls")) 
      
      object$call$model <- Terms
   
   else object$call$formula <- Terms
   
   if (use.start) 
      
      warning("'use.start' cannot be used with R's version of glm")
   
   md <- missing(direction)
   
   direction <- match.arg(direction)
   
   backward <- direction == "both" | direction == "backward"
   
   forward <- direction == "both" | direction == "forward"
   
   if (missing(scope)) {
      
      fdrop <- numeric(0)
      
      fadd <- attr(Terms, "factors")
      
      if (md) 
         
         forward <- FALSE
      
   }
   
   else {
      
      if (is.list(scope)) {
         
         fdrop <- if (!is.null(fdrop <- scope$lower)) 
            
            attr(terms(update.formula(object, fdrop)), "factors")
         
         else numeric(0)
         
         fadd <- if (!is.null(fadd <- scope$upper)) 
            
            attr(terms(update.formula(object, fadd)), "factors")
         
      }
      
      else {
         
         fadd <- if (!is.null(fadd <- scope)) 
            
            attr(terms(update.formula(object, scope)), "factors")
         
         fdrop <- numeric(0)
         
      }
      
   }
   
   models <- vector("list", steps)
   
   if (!is.null(keep)) 
      
      keep.list <- vector("list", steps)
   
   if (is.list(object) && (nmm <- match("nobs", names(object), 
                                        
                                        0)) > 0) 
      
      n <- object[[nmm]]
   
   else n <- length(residuals(object))
   
   fit <- object
   
   bAIC <- extractAICc(fit, scale, k = k, ...)
   
   edf <- bAIC[1L]
   
   bAIC <- bAIC[2L]
   
   if (is.na(bAIC)) 
      
      stop("AIC is not defined for this model, so stepAIC cannot proceed")
   
   nm <- 1
   
   Terms <- terms(fit)
   
   if (trace) {
      
      cat("Start:  AICc=", format(round(bAIC, 2)), "\n", cut.string(deparse(as.vector(formula(fit)))), 
          
          "\n\n", sep = "")
      
      utils::flush.console()
      
   }
   
   models[[nm]] <- list(deviance = mydeviance(fit), df.resid = n - 
                           
                           edf, change = "", AIC = bAIC)
   
   if (!is.null(keep)) 
      
      keep.list[[nm]] <- keep(fit, bAIC)
   
   usingCp <- FALSE
   
   while (steps > 0) {
      
      steps <- steps - 1
      
      AIC <- bAIC
      
      ffac <- attr(Terms, "factors")
      
      if (!is.null(sp <- attr(Terms, "specials")) && !is.null(st <- sp$strata)) 
         
         ffac <- ffac[-st, ]
      
      scope <- factor.scope(ffac, list(add = fadd, drop = fdrop))
      
      aod <- NULL
      
      change <- NULL
      
      if (backward && length(scope$drop)) {
         
         aod <- dropterm.AICc(fit, scope$drop, scale = scale, trace = max(0, 
                                                                          
                                                                          trace - 1), k = k, ...)
         
         rn <- row.names(aod)
         
         row.names(aod) <- c(rn[1L], paste("-", rn[-1L], sep = " "))
         
         if (any(aod$Df == 0, na.rm = TRUE)) {
            
            zdf <- aod$Df == 0 & !is.na(aod$Df)
            
            nc <- match(c("Cp", "AIC"), names(aod))
            
            nc <- nc[!is.na(nc)][1L]
            
            ch <- abs(aod[zdf, nc] - aod[1, nc]) > 0.01
            
            if (any(ch)) {
               
               warning("0 df terms are changing AIC")
               
               zdf <- zdf[!ch]
               
            }
            
            if (length(zdf) > 0L) 
               
               change <- rev(rownames(aod)[zdf])[1L]
            
         }
         
      }
      
      if (is.null(change)) {
         
         if (forward && length(scope$add)) {
            
            aodf <- addterm.AICc(fit, scope$add, scale = scale, 
                                 
                                 trace = max(0, trace - 1), k = k, ...)
            
            rn <- row.names(aodf)
            
            row.names(aodf) <- c(rn[1L], paste("+", rn[-1L], 
                                               
                                               sep = " "))
            
            aod <- if (is.null(aod)) 
               
               aodf
            
            else rbind(aod, aodf[-1, , drop = FALSE])
            
         }
         
         attr(aod, "heading") <- NULL
         
         if (is.null(aod) || ncol(aod) == 0) 
            
            break
         
         nzdf <- if (!is.null(aod$Df)) 
            
            aod$Df != 0 | is.na(aod$Df)
         
         aod <- aod[nzdf, ]
         
         if (is.null(aod) || ncol(aod) == 0) 
            
            break
         
         nc <- match(c("Cp", "AIC"), names(aod))
         
         nc <- nc[!is.na(nc)][1L]
         
         o <- order(aod[, nc])
         
         if (trace) {
            
            print(aod[o, ])
            
            utils::flush.console()
            
         }
         
         if (o[1L] == 1) 
            
            break
         
         change <- rownames(aod)[o[1L]]
         
      }
      
      usingCp <- match("Cp", names(aod), 0) > 0
      
      fit <- update(fit, paste("~ .", change), evaluate = FALSE)
      
      fit <- eval.parent(fit)
      
      if (is.list(fit) && (nmm <- match("nobs", names(fit), 
                                        
                                        0)) > 0) 
         
         nnew <- fit[[nmm]]
      
      else nnew <- length(residuals(fit))
      
      if (nnew != n) 
         
         stop("number of rows in use has changed: remove missing values?")
      
      Terms <- terms(fit)
      
      bAIC <- extractAICc(fit, scale, k = k, ...)
      
      edf <- bAIC[1L]
      
      bAIC <- bAIC[2L]
      
      if (trace) {
         
         cat("\nStep:  AICc=", format(round(bAIC, 2)), "\n", 
             
             cut.string(deparse(as.vector(formula(fit)))), 
             
             "\n\n", sep = "")
         
         utils::flush.console()
         
      }
      
      if (bAIC >= AIC + 1e-07) 
         
         break
      
      nm <- nm + 1
      
      models[[nm]] <- list(deviance = mydeviance(fit), df.resid = n - 
                              
                              edf, change = change, AIC = bAIC)
      
      if (!is.null(keep)) 
         
         keep.list[[nm]] <- keep(fit, bAIC)
      
   }
   
   if (!is.null(keep)) 
      
      fit$keep <- re.arrange(keep.list[seq(nm)])
   
   step.results(models = models[seq(nm)], fit, object, usingCp)
   
}
##############################################################
#### load data and libraries
setwd("C:/Users/xiong/Desktop/MS_policy analysis China")
library(tseries)
#### N_tr
N_tr<-read.csv("N_tr.csv",header=TRUE,sep=',')
str(N_tr)
N_tr$N_tr
## explore normality
shapiro.test(N_tr$N_tr) # p < 0.05, not normal
hist(N_tr$N_tr) # close to normal, use z-score
hist(N_tr$RNT)
## test autocorrelation
N_tr_ts<-ts(N_tr,start = c(1951,1),end = c(2018,1),frequency=1)
summary(N_tr_ts)
acf1<-acf(N_tr_ts[,2],type="p")
acf1
## arima

library(forecast)
rg<-factor(N_tr$Agenda21)
rg
ari0<-auto.arima(N_tr$N_tr, trace=T) # 
summary(ari0)
Box.test(ari0$residuals) ## null: the autocorrelations are all zerios. when p < 0.05, reject it.

## contrast model
plot(N_tr~Year,N_tr)
cont_n_tr<-lm(N_tr~L(N_tr,1)+L(N_tr,2),data=N_tr)
summary(cont_n_tr) # 0.2681
cont_n_tr$rse <- vcovHC(cont_n_tr, type="HC1")
coeftest(cont_n_tr, cont_n_tr$rse) # Welch's t test on coefficients
# the default setting in R: options(contrasts = c(unordered = "contr.treatment", ordered = "contr.poly"))
## full model
options(contrasts = c(unordered = "contr.treatment", ordered = "contr.treatment")) # make sure factor() does not derive levels
full_n_tr<-lm(N_tr~L(N_tr,1)+L(N_tr,2)+factor(Agenda21)+factor(L(FueSub,2))+factor(MarFisPol)+factor(L(EcoRef,1))+factor(L(DouCon,1))+factor(L(VesScrap,1)),data=N_tr,na.action = na.exclude)

summary(full_n_tr) # 0.7901

# examine statistical effects
library(lmtest)
library(sandwich)
full_n_tr$rse <- vcovHC(full_n_tr, type="HC1")
coeftest(full_n_tr, full_n_tr$rse) # Welch's t test on coefficients
coef<-coeftest(full_n_tr, full_n_tr$rse)
write.csv(coef,"Coefficient full_n_tr.csv")
## residual check for independence, homoskedasticity, and normality
plot(full_n_tr$residuals)
acf(full_n_tr$residuals,type="p")
Box.test(full_n_tr$residuals, lag = 6) # autocorrelation not sign.
lmtest::bptest(full_n_tr)
shapiro.test(full_n_tr$residuals)
qqnorm(full_n_tr$residuals)
qqline(full_n_tr$residuals, col = "steelblue", lwd = 2)

plot(opt_n_tr$residuals)
acf(opt_n_tr$residuals,type="p")
Box.test(opt_n_tr$residuals, lag = 6) # autocorrelation not sign.
lmtest::bptest(opt_n_tr)
shapiro.test(opt_n_tr$residuals)
qqnorm(opt_n_tr$residuals)
qqline(opt_n_tr$residuals, col = "steelblue", lwd = 2)

# box-cox transformation
BCNT <- caret::BoxCoxTrans(N_tr$N_tr)
BCNT
## compare two models
anova(cont_n_tr,full_n_tr)
## contrast with coxtest for non-nested models
#library(lmtest)
#coxtest(dlmb00,dlmb2) # F=11.035, p < .001
##### dominance analysis
remove.packages("dominanceanalysis")
install.packages("dominanceanalysis")
library(dominanceanalysis)
dominanceAnalysis(full_n_tr)
## simple DLMs
# Agenda21, DouCon,1, VesScrap,1
N_tr2<-N_tr[,-c(3,5,8)]
str(N_tr2)
r<-c(1,0,1)
for(i in 3:5){
   l<-r[i-2]
   dlm2 <- lm(N_tr~L(N_tr,1)+L(N_tr,2)+L(N_tr,10)+factor(L(N_tr2[,i],l)), data = N_tr2)
   print(names(N_tr2)[i]);print(l)
   dlm2$rse <- vcovHC(dlm2, type="HC1")
   print(coeftest(dlm2, dlm2$rse)) # Welch's t test on coefficients
}
#### H_tr
H_tr<-read.csv("H_tr.csv",header=TRUE,sep=',')
str(H_tr)

##create times-series dataframe
H_tr1_ts<-ts(H_tr,start = c(1951,1),end = c(2018,1),frequency=1)
summary(H_tr1_ts)

## explore normality
shapiro.test(H_tr$H_tr) # p < 0.05, not normal
hist(H_tr$H_tr) # close to normal
## contrast model
acf1<-acf(H_tr1_ts[,2],type="p")
acf1
plot(H_tr~Year,H_tr)
cont_h_tr<-lm(H_tr~L(H_tr,1),data=H_tr[-c(1,2),])
summary(cont_h_tr) # 0.730
fitted(cont_h_tr)
# the default setting in R: options(contrasts = c(unordered = "contr.treatment", ordered = "contr.poly"))
## full model
options(contrasts = c(unordered = "contr.treatment", ordered = "contr.treatment")) # make sure factor() does not derive levels
full_h_tr<-lm(H_tr~L(H_tr,1)+factor(L(Agenda21,1))+factor(L(DouCon,1))+factor(L(EcoRef,3))+factor(L(FisPerReg,1))+factor(FueSub)+factor(MarFisPol)+factor(L(VesScrap,1)),data=H_tr1_ts,na.action = na.exclude)
summary(full_h_tr) # 0.893
fitted(full_h_tr)
# examine statistical effects
library(lmtest)
library(sandwich)
full_h_tr$rse <- vcovHC(full_h_tr, type="HC1")
coeftest(full_h_tr, full_h_tr$rse) # Welch's t test on coefficients
coef<-coeftest(full_h_tr, full_h_tr$rse)
write.csv(coef,"Coefficient full_h_tr.csv")
## compare two models
anova(cont_h_tr,full_h_tr)
##### dominance analysis
library(dominanceanalysis)
dominanceAnalysis(full_h_tr)
## simple DLMs
# Agenda21, DouCon,1, VesScrap,1
coef
str(H_tr)
cpd2<-H_tr[,-5]
str(cpd2)
r<-c(1,1,1,0,0,1)
for(i in 3:8){
   l<-r[i-2]
   dlm2 <- lm(H_tr~L(H_tr,1)+factor(L(cpd2[,i],l)), data = cpd2)
   print(names(cpd2)[i]);print(l)
   dlm2$rse <- vcovHC(dlm2, type="HC1")
   print(coeftest(dlm2, dlm2$rse)) # Welch's t test on coefficients
}
#### H_dt
H_dt<-read.csv("H_dt.csv",header=TRUE,sep=',')
str(H_dt)

##create times-series dataframe
H_dt1_ts<-ts(H_dt,start = c(1951,1),end = c(2018,1),frequency=1)
summary(H_dt1_ts)

## explore normality
shapiro.test(H_dt$H_dt) # p = .229, normal
hist(H_dt$H_dt) # close to normal
## contrast model
acf1<-acf(na.omit(H_dt1_ts[,2]),type="p")
acf1
plot(H_dt~Year,H_dt)
cont_H_dt<-lm(H_dt~factor(FinCris),data=H_dt)
summary(cont_H_dt) # # 0.2612
fitted(cont_H_dt)
cont_H_dt$rse <- vcovHC(cont_H_dt, type="HC1")
coeftest(cont_H_dt, cont_H_dt$rse) # Welch's t test on coefficients
# the default setting in R: options(contrasts = c(unordered = "contr.treatment", ordered = "contr.poly"))
## full model
options(contrasts = c(unordered = "contr.treatment", ordered = "contr.treatment")) # make sure factor() does not derive levels
full_H_dt<-lm(H_dt~factor(EcoRef)+factor(FisPerReg)+factor(FinCris)+factor(L(FueSub,1))+factor(MarFisPol)+factor(L(VesScrap,1))+factor(L(FishAgree,1)),data=H_dt1_ts,na.action = na.exclude)
summary(full_H_dt) # 0.663
fitted(full_H_dt)
# op_H_dt<-stepAICc(full_H_dt, direction="backward"); summary(op_H_dt)
# examine statistical effects
library(lmtest)
library(sandwich)
full_H_dt$rse <- vcovHC(full_H_dt, type="HC1")
coeftest(full_H_dt, full_H_dt$rse) # Welch's t test on coefficients
coef<-coeftest(full_H_dt, full_H_dt$rse)
write.csv(coef,"Coefficient full_H_dt.csv")
## compare two models
anova(cont_H_dt,full_H_dt)
##### dominance analysis
library(dominanceanalysis)
dominanceAnalysis(full_H_dt)
## simple DLMs
# Agenda21, DouCon,1, VesScrap,1
coef
str(H_dt)
cpd2<-H_dt[,-7]
str(cpd2)
r<-c(0,0,1,0,1)
for(i in 3:7){
   l<-r[i-2]
   dlm2 <- lm(H_dt~factor(L(cpd2[,i],l)), data = cpd2)
   print(names(cpd2)[i]);print(l)
   dlm2$rse <- vcovHC(dlm2, type="HC1")
   print(coeftest(dlm2, dlm2$rse)) # Welch's t test on coefficients
}

#### C_bt
C_bt<-read.csv("C_bt.csv",header=TRUE,sep=',')
str(C_bt)

##create times-series dataframe
C_bt1_ts<-ts(C_bt,start = c(1951,1),end = c(2014,1),frequency=1)
summary(C_bt1_ts)

## explore normality
shapiro.test(C_bt$C_bt) # p = .229, normal
hist(C_bt$C_bt) # close to normal
## contrast model
acf1<-acf(na.omit(C_bt1_ts[,2]),type="p")
acf1
plot(C_bt~Year,C_bt)
cont_C_bt<-lm(C_bt~1,data=C_bt[-1,])
summary(cont_C_bt) 
fitted(cont_C_bt)
# the default setting in R: options(contrasts = c(unordered = "contr.treatment", ordered = "contr.poly"))
## full model
options(contrasts = c(unordered = "contr.treatment", ordered = "contr.treatment")) # make sure factor() does not derive levels
full_C_bt<-lm(C_bt~factor(FishAgree)+factor(L(Agenda21,1))+factor(L(EEZ,1))+factor(FueSub)+factor(LawEnf)+factor(L(MarFisPol,1))+factor(OutCon)+factor(SumMor),data=C_bt1_ts,na.action = na.exclude)
summary(full_C_bt) # Adj R-sq: 0.930, p < 0.001
fitted(full_C_bt)
# examine statistical effects
library(lmtest)
library(sandwich)
full_C_bt$rse <- vcovHC(full_C_bt, type="HC1")
coeftest(full_C_bt, full_C_bt$rse) # Welch's t test on coefficients
coef<-coeftest(full_C_bt, full_C_bt$rse)
write.csv(coef,"Coefficient full_C_bt.csv")
## compare two models
anova(cont_C_bt,full_C_bt)
##### dominance analysis
library(dominanceanalysis)
dominanceAnalysis(full_C_bt)
## simple DLMs
# FishAgree, LawEnf, SumMor
coef
str(C_bt)
cpd2<-C_bt
str(cpd2)

for(i in 3:6){
   l<-r[i-2]
   dlm2 <- lm(C_bt~factor(cpd2[,i]), data = cpd2)
   print(names(cpd2)[i]);print(l)
   dlm2$rse <- vcovHC(dlm2, type="HC1")
   print(coeftest(dlm2, dlm2$rse)) # Welch's t test on coefficients
}

#### C_dbt
C_dbt<-read.csv("C_dbt.csv",header=TRUE,sep=',')
str(C_dbt)

##create times-series dataframe
C_dbt1_ts<-ts(C_dbt,start = c(1951,1),end = c(2014,1),frequency=1)
summary(C_dbt1_ts)

## explore normality
shapiro.test(C_dbt$C_dbt) # p = .229, normal
hist(C_dbt$C_dbt) # close to normal
## contrast model
acf1<-acf(na.omit(C_dbt1_ts[,2]),type="p")
acf1
plot(C_dbt~Year,C_dbt)
cont_C_dbt<-lm(C_dbt~1,data=C_dbt[-1,])
summary(cont_C_dbt) 
fitted(cont_C_dbt)
# the default setting in R: options(contrasts = c(unordered = "contr.treatment", ordered = "contr.poly"))
## full model
options(contrasts = c(unordered = "contr.treatment", ordered = "contr.treatment")) # make sure factor() does not derive levels
full_C_dbt<-lm(C_dbt~factor(FishAgree)+factor(FueSub)+factor(MarFisPol)+factor(OutCon),data=C_dbt1_ts,na.action = na.exclude)
summary(full_C_dbt) # Adj R-sq: 0.3829, p < 0.012
fitted(full_C_dbt)
# examine statistical effects
library(lmtest)
library(sandwich)
full_C_dbt$rse <- vcovHC(full_C_dbt, type="HC1")
coeftest(full_C_dbt, full_C_dbt$rse) # Welch's t test on coefficients
coef<-coeftest(full_C_dbt, full_C_dbt$rse)
write.csv(coef,"Coefficient full_C_dbt.csv")
## compare two models
anova(cont_C_dbt,full_C_dbt)
##### dominance analysis
library(dominanceanalysis)
dominanceAnalysis(full_C_dbt)

## simple DLMs
# FishAgree, LawEnf, SumMor
coef
str(C_dbt)
cpd2<-C_dbt[,-c(3,6)]
str(cpd2)
r<-c()
for(i in 3:4){
   dlm2 <- lm(C_dbt~factor(cpd2[,i]), data = cpd2)
   print(names(cpd2)[i])
   dlm2$rse <- vcovHC(dlm2, type="HC1")
   print(coeftest(dlm2, dlm2$rse)) # Welch's t test on coefficients
}

#### HpV_bt
HpV_bt<-read.csv("HpV_bt.csv",header=TRUE,sep=',')
str(HpV_bt)

##create times-series dataframe
HpV_bt1_ts<-ts(HpV_bt,start = c(1951,1),end = c(2018,1),frequency=1)
summary(HpV_bt1_ts)

## explore normality
shapiro.test(HpV_bt$HpV_bt) # p < 0.05, not normal
hist(HpV_bt$HpV_bt) # close to normal
## contrast model
acf1<-acf(HpV_bt1_ts[,2],type="p")
acf1
plot(HpV_bt~Year,HpV_bt)
cont_HpV_bt<-lm(HpV_bt~L(HpV_bt,7),data=HpV_bt)
summary(cont_HpV_bt) # 0.0738
fitted(cont_HpV_bt)
cont_HpV_bt$rse <- vcovHC(cont_HpV_bt, type="HC1")
coeftest(cont_HpV_bt, cont_HpV_bt$rse) # Welch's t test on coefficients
cont_HpV_bt<-lm(HpV_bt~1,data=HpV_bt)
summary(cont_HpV_bt) # 0
# the default setting in R: options(contrasts = c(unordered = "contr.treatment", ordered = "contr.poly"))
## full model
options(contrasts = c(unordered = "contr.treatment", ordered = "contr.treatment")) # make sure factor() does not derive levels
full_HpV_bt<-lm(HpV_bt~factor(L(DouCon,1))+factor(L(NotrawlZone,1))+factor(Agenda21)+factor(L(EcoRef,1))+factor(FisPerReg)+factor(FueSub)+factor(MarFisPol)+factor(L(VesScrap,1)),data=HpV_bt1_ts,na.action = na.exclude)
summary(full_HpV_bt) # 0.452
fitted(full_HpV_bt)
opt_HpV_bt<-stepAICc(full_HpV_bt, direction = "backward")
summary(opt_HpV_bt) # 0.501
opt_HpV_bt$rse <- vcovHC(opt_HpV_bt, type="HC1")
coeftest(opt_HpV_bt, opt_HpV_bt$rse) # Welch's t test on coefficients
# examine statistical effects
library(lmtest)
library(sandwich)
full_HpV_bt$rse <- vcovHC(full_HpV_bt, type="HC1")
coeftest(full_HpV_bt, full_HpV_bt$rse) # Welch's t test on coefficients
coef<-coeftest(full_HpV_bt, full_HpV_bt$rse)
write.csv(coef,"Coefficient full_HpV_bt.csv")
## compare two models
cont_HpV_bt<-lm(HpV_bt~1,data=HpV_bt[-1,])
anova(cont_HpV_bt,full_HpV_bt)
##### dominance analysis
library(dominanceanalysis)
dominanceAnalysis(full_HpV_bt)

## simple DLMs
# Agenda21, DouCon,1, VesScrap,1
coef
str(HpV_bt)
cpd2<-HpV_bt[,-c(8,10)]
str(cpd2)
r<-c(0,1,1,0,0,0)
for(i in 3:8){
   l<-r[i-2]
   dlm2 <- lm(HpV_bt~factor(L(cpd2[,i],l)), data = cpd2)
   print(names(cpd2)[i]);print(l)
   dlm2$rse <- vcovHC(dlm2, type="HC1")
   print(coeftest(dlm2, dlm2$rse)) # Welch's t test on coefficients
}

#### HpV_dbt
HpV_dbt<-read.csv("HpV_dbt.csv",header=TRUE,sep=',')
str(HpV_dbt)

##create times-series dataframe
HpV_dbt1_ts<-ts(HpV_dbt,start = c(1951,1),end = c(2018,1),frequency=1)
summary(HpV_dbt1_ts)

## explore normality
shapiro.test(HpV_dbt$HpV_dbt) # p < 0.05, not normal
hist(HpV_dbt$HpV_dbt) # close to normal
## contrast model
acf1<-acf(na.omit(HpV_dbt1_ts[,2]),type="p")
acf1
plot(HpV_dbt~Year,HpV_dbt)
cont_HpV_dbt<-lm(HpV_dbt~factor(FinCris),data=HpV_dbt)
cont_HpV_dbt$rse <- vcovHC(cont_HpV_dbt, type="HC1")
coeftest(cont_HpV_dbt, cont_HpV_dbt$rse) # Welch's t test on coefficients

summary(cont_HpV_dbt) # 0.220
fitted(cont_HpV_dbt)
# the default setting in R: options(contrasts = c(unordered = "contr.treatment", ordered = "contr.poly"))
## full model
options(contrasts = c(unordered = "contr.treatment", ordered = "contr.treatment")) # make sure factor() does not derive levels
full_HpV_dbt<-lm(HpV_dbt~factor(FinCris)+factor(L(EcoRef,1))+factor(L(EEZ,1))+factor(L(DouCon,1))+factor(L(FishAgree,1))+factor(FisPerReg)+factor(L(FueSub,1))+factor(MarFisPol)+factor(L(VesScrap,1)),data=HpV_dbt1_ts,na.action = na.exclude)
summary(full_HpV_dbt) # 0.392
# examine statistical effects
library(lmtest)
library(sandwich)
full_HpV_dbt$rse <- vcovHC(full_HpV_dbt, type="HC1")
coeftest(full_HpV_dbt, full_HpV_dbt$rse) # Welch's t test on coefficients
coef<-coeftest(full_HpV_dbt, full_HpV_dbt$rse)
write.csv(coef,"Coefficient full_HpV_dbt.csv")
## compare two models
anova(cont_HpV_dbt,full_HpV_dbt)

opt_HpV_dbt<-stepAICc(cont_HpV_dbt, scope= list(upper=full_HpV_dbt,
                                                lower=cont_HpV_dbt),
                      direction = "forward")
summary(opt_HpV_dbt) # 0.501
anova(cont_HpV_dbt,opt_HpV_dbt)

##### dominance analysis
library(dominanceanalysis)
dominanceAnalysis(full_HpV_dbt)

## simple DLMs
# Agenda21, DouCon,1, VesScrap,1
coef
str(HpV_dbt)
cpd2<-HpV_dbt[,-c(4,8,9)]
str(cpd2)
r<-c(1,0,1,1,1,1)
fincris<-as.factor(HpV_dbt$FinCris)
for(i in 3:8){
   l<-r[i-2]
   dlm2 <- lm(HpV_dbt~fincris+factor(L(cpd2[,i],l)), data = cpd2)
   print(names(cpd2)[i]);print(l)
   dlm2$rse <- vcovHC(dlm2, type="HC1")
   print(coeftest(dlm2, dlm2$rse)) # Welch's t test on coefficients
}

#### CpH_bt
CpH_bt<-read.csv("CpH_bt.csv",header=TRUE,sep=',')
str(CpH_bt)

##create times-series dataframe
CpH_bt1_ts<-ts(CpH_bt,start = c(1951,1),end = c(2014,1),frequency=1)
summary(CpH_bt1_ts)

## explore normality
shapiro.test(CpH_bt$CpH_bt) # p = .229, normal
hist(CpH_bt$CpH_bt) # close to normal
## contrast model
acf1<-acf(na.omit(CpH_bt1_ts[,2]),type="p")
acf1
plot(CpH_bt~Year,CpH_bt)
cont_CpH_bt<-lm(CpH_bt~L(CpH_bt,1)+L(CpH_bt,2)+L(CpH_bt,4),data=CpH_bt)
summary(cont_CpH_bt) 
fitted(cont_CpH_bt)
cont_CpH_bt$rse <- vcovHC(cont_CpH_bt, type="HC1")
coeftest(cont_CpH_bt, cont_CpH_bt$rse) # Welch's t test on coefficients
cont_CpH_bt<-lm(CpH_bt~L(CpH_bt,2),data=CpH_bt)
summary(cont_CpH_bt) # 0.245
# the default setting in R: options(contrasts = c(unordered = "contr.treatment", ordered = "contr.poly"))
## full model
options(contrasts = c(unordered = "contr.treatment", ordered = "contr.treatment")) # make sure factor() does not derive levels
full_CpH_bt<-lm(CpH_bt~L(CpH_bt,2)+factor(L(EcoRef,1))+factor(EEZ)+factor(LawEnf)+factor(L(MarFisPol,1))+factor(NotrawlZone)+factor(SumMor),data=CpH_bt,na.action = na.exclude)
summary(full_CpH_bt) # Adj R-sq: 0.410, p < 0.05
fitted(full_CpH_bt)
# examine statistical effects
library(lmtest)
library(sandwich)
full_CpH_bt$rse <- vcovHC(full_CpH_bt, type="HC1")
coeftest(full_CpH_bt, full_CpH_bt$rse) # Welch's t test on coefficients
coef<-coeftest(full_CpH_bt, full_CpH_bt$rse)
write.csv(coef,"Coefficient full_CpH_bt.csv")
# optimal model
opt_CpH_bt<-stepAICc(full_CpH_bt, direction = "backward")
summary(opt_CpH_bt) # Adj R-sq: 0.446
opt_CpH_bt$rse <- vcovHC(opt_CpH_bt, type="HC1")
coeftest(opt_CpH_bt, opt_CpH_bt$rse) # Welch's t test on coefficients
## compare two models
anova(cont_CpH_bt,full_CpH_bt)
anova(cont_CpH_bt, opt_CpH_bt)
##### dominance analysis
library(dominanceanalysis)
dominanceAnalysis(full_CpH_bt)

## simple DLMs
# FishAgree, LawEnf, SumMor
coef
str(CpH_bt)
cpd2<-CpH_bt[,c(2,5)]
str(cpd2)
dlm2 <- lm(CpH_bt~L(CpH_bt,2)+factor(L(cpd2[,2],1)), data = cpd2)
print(names(cpd2)[2])
dlm2$rse <- vcovHC(dlm2, type="HC1")
print(coeftest(dlm2, dlm2$rse)) # Welch's t test on coefficients

#### RN_bt
RN_bt<-read.csv("RN_bt.csv",header=TRUE,sep=',')
str(RN_bt)

##create times-series dataframe
RN_bt1_ts<-ts(RN_bt,start = c(1951,1),end = c(2014,1),frequency=1)
summary(RN_bt1_ts)

## explore normality
shapiro.test(RN_bt$RN_bt) # p = .229, normal
hist(RN_bt$RN_bt) # close to normal
## contrast model
acf1<-acf(na.omit(RN_bt1_ts[,2]),type="p")
acf1
plot(RN_bt~Year,RN_bt)
cont_RN_bt<-lm(RN_bt~L(RN_bt,1)+Year,data=RN_bt)
summary(cont_RN_bt) 
fitted(cont_RN_bt)
cont_RN_bt$rse <- vcovHC(cont_RN_bt, type="HC1")
coeftest(cont_RN_bt, cont_RN_bt$rse) # Welch's t test on coefficients
cont_RN_bt<-lm(RN_bt~1,data=RN_bt)
summary(cont_RN_bt) # 0.245
# the default setting in R: options(contrasts = c(unordered = "contr.treatment", ordered = "contr.poly"))
## full model
options(contrasts = c(unordered = "contr.treatment", ordered = "contr.treatment")) # make sure factor() does not derive levels
full_RN_bt<-lm(RN_bt~factor(L(EcoRef,2))+factor(EEZ)+factor(DouCon)+factor(L(MarFisPol,1))+factor(NotrawlZone)+factor(FishAgree)+factor(FueSub),data=RN_bt,na.action = na.exclude)
summary(full_RN_bt) # Adj R-sq: 0.703, p < 0.001
# examine statistical effects
library(lmtest)
library(sandwich)
full_RN_bt$rse <- vcovHC(full_RN_bt, type="HC1")
coeftest(full_RN_bt, full_RN_bt$rse) # Welch's t test on coefficients
coef<-coeftest(full_RN_bt, full_RN_bt$rse)
write.csv(coef,"Coefficient full_RN_bt.csv")
# optimal model
## compare two models
cont_RN_bt<-lm(RN_bt~1,data=RN_bt[-c(1,2),])
anova(cont_RN_bt,full_RN_bt)
##### dominance analysis
library(dominanceanalysis)
dominanceAnalysis(full_RN_bt)

## simple DLMs
# FishAgree, LawEnf, SumMor
coef
str(RN_bt)
cpd2<-RN_bt[,-c(3,4,6,9)]
str(cpd2)
r<-c(2,0,0)
for(i in 3:5){
   l<-r[i-2]
   dlm2 <- lm(RN_bt~L(RN_bt,2)+factor(cpd2[,i]), data = cpd2)
   print(names(cpd2)[i]);print(l)
   dlm2$rse <- vcovHC(dlm2, type="HC1")
   print(coeftest(dlm2, dlm2$rse)) # Welch's t test on coefficients
}
####
#### RH_bt
RH_bt<-read.csv("RH_bt.csv",header=TRUE,sep=',')
str(RH_bt)

##create times-series dataframe
RH_bt1_ts<-ts(RH_bt,start = c(1951,1),end = c(2014,1),frequency=1)
summary(RH_bt1_ts)

## explore normality
shapiro.test(RH_bt$RH_bt) # p = .229, normal
hist(RH_bt$RH_bt) # close to normal
## contrast model
acf1<-acf(na.omit(RH_bt1_ts[,2]),type="p")
acf1
plot(RH_bt~Year,RH_bt)
cont_RH_bt<-lm(RH_bt~L(RH_bt,1)+L(RH_bt,2)+Year,data=RH_bt)
summary(cont_RH_bt) 
cont_RH_bt$rse <- vcovHC(cont_RH_bt, type="HC1")
coeftest(cont_RH_bt, cont_RH_bt$rse) # Welch's t test on coefficients
cont_RH_bt<-lm(RH_bt~L(RH_bt,1)+L(RH_bt,2),data=RH_bt)
summary(cont_RH_bt) # 0.164
# the default setting in R: options(contrasts = c(unordered = "contr.treatment", ordered = "contr.poly"))
## full model
options(contrasts = c(unordered = "contr.treatment", ordered = "contr.treatment")) # make sure factor() does not derive levels
full_RH_bt<-lm(RH_bt~L(RH_bt,1)+L(RH_bt,2)+factor(L(EcoRef,3))+factor(L(Agenda21,1))+factor(L(DouCon,1))+factor(L(MarFisPol,1))+factor(NotrawlZone)+factor(L(VesScrap,1))+factor(SumMor),data=RH_bt,na.action = na.exclude)
summary(full_RH_bt) # Adj R-sq: 0.525, p < 0.001
# examine statistical effects
library(lmtest)
library(sandwich)
full_RH_bt$rse <- vcovHC(full_RH_bt, type="HC1")
coeftest(full_RH_bt, full_RH_bt$rse) # Welch's t test on coefficients
coef<-coeftest(full_RH_bt, full_RH_bt$rse)
write.csv(coef,"Coefficient full_RH_bt.csv")
# optimal model
## compare two models
cont_RH_bt<-lm(RH_bt~L(RH_bt,1)+L(RH_bt,2),data=RH_bt[-1,])
anova(cont_RH_bt,full_RH_bt)
##### dominance analysis
library(dominanceanalysis)
dominanceAnalysis(full_RH_bt)

## simple DLMs
# FishAgree, LawEnf, SumMor
coef
str(RH_bt)
cpd2<-RH_bt[,-c(4,6)]
str(cpd2)
r<-c(0,3,0,1,1)
for(i in 3:7){
   l<-r[i-2]
   dlm2 <- lm(RH_bt~L(RH_bt,1)+L(RH_bt,2)+factor(cpd2[,i]), data = cpd2)
   print(names(cpd2)[i]);print(l)
   dlm2$rse <- vcovHC(dlm2, type="HC1")
   print(coeftest(dlm2, dlm2$rse)) # Welch's t test on coefficients
}
####
#### RC_bt_eez
RC_bt_eez<-read.csv("RC_bt_eez.csv",header=TRUE,sep=',')
str(RC_bt_eez)

##create times-series dataframe
RC_bt_eez1_ts<-ts(RC_bt_eez,start = c(1951,1),end = c(2014,1),frequency=1)
summary(RC_bt_eez1_ts)

## explore normality
shapiro.test(RC_bt_eez$RC_bt_eez) # p = .229, normal
hist(RC_bt_eez$RC_bt_eez) # close to normal
## contrast model
acf1<-acf(na.omit(RC_bt_eez1_ts[,2]),type="p")
acf1
plot(RC_bt_eez~Year,RC_bt_eez)
cont_RC_bt_eez<-lm(RC_bt_eez~Year,data=RC_bt_eez)
summary(cont_RC_bt_eez) 
cont_RC_bt_eez$rse <- vcovHC(cont_RC_bt_eez, type="HC1")
coeftest(cont_RC_bt_eez, cont_RC_bt_eez$rse) # Welch's t test on coefficients
cont_RC_bt_eez<-lm(RC_bt_eez~1,data=RC_bt_eez)
summary(cont_RC_bt_eez) # 0
# the default setting in R: options(contrasts = c(unordered = "contr.treatment", ordered = "contr.poly"))
## full model
options(contrasts = c(unordered = "contr.treatment", ordered = "contr.treatment")) # make sure factor() does not derive levels
RC_bt_eez$SumMor=as.factor(RC_bt_eez$SumMor)
RC_bt_eez$NotrawlZone=as.factor(RC_bt_eez$NotrawlZone)
RC_bt_eez$MarFisPol=as.factor(RC_bt_eez$MarFisPol)
RC_bt_eez$EcoRef=as.factor(L(RC_bt_eez$EcoRef,1))
RC_bt_eez$DouCon=as.factor(L(RC_bt_eez$DouCon,1))
str(RC_bt_eez)
RC_bt_eez<-RC_bt_eez[-1,-1]

full_RC_bt_eez<-lm(RC_bt_eez~.,data=RC_bt_eez)

#full_RC_bt_eez<-lm(RC_bt_eez~factor(L(EcoRef,1))+factor(L(DouCon,1))+factor(MarFisPol)+factor(NotrawlZone)+factor(SumMor),data=RC_bt_eez,na.action = na.exclude)
summary(full_RC_bt_eez) # Adj R-sq: 0.021
## optimal model
cont_RC_bt_eez<-lm(RC_bt_eez~1,data=RC_bt_eez)
summary(cont_RC_bt_eez)
opt_RC_bt<-stepAICc(full_RC_bt_eez, direction = "backward",trace=T)
opt_RC_bt<-stepAICc(cont_RC_bt_eez, scope= list(upper=full_RC_bt_eez,
                                                lower=cont_RC_bt_eez),
                    direction = "forward", trace=T)
summary(opt_RC_bt) # Adj R-sq: 0.066
# examine statistical effects
library(lmtest)
library(sandwich)
full_RC_bt_eez$rse <- vcovHC(full_RC_bt_eez, type="HC1")
coeftest(full_RC_bt_eez, full_RC_bt_eez$rse) # Welch's t test on coefficients
coef<-coeftest(full_RC_bt_eez, full_RC_bt_eez$rse)
write.csv(coef,"Coefficient full_RC_bt_eez.csv")
# optimal model
## compare two models
cont_RC_bt_eez<-lm(RC_bt_eez~1,data=RC_bt_eez[-1,])
anova(cont_RC_bt_eez,full_RC_bt_eez)

anova(cont_RC_bt_eez,opt_RC_bt)

##### dominance analysis
library(dominanceanalysis)
dominanceAnalysis(full_RC_bt_eez)

## simple DLMs
# FishAgree, LawEnf, SumMor
coef
str(RC_bt_eez)
cpd2<-RC_bt_eez
str(cpd2)
r<-c(1,1,0,0,0)
for(i in 3:7){
   l<-r[i-2]
   dlm2 <- lm(RC_bt_eez~factor(cpd2[,i]), data = cpd2)
   print(names(cpd2)[i]);print(l)
   dlm2$rse <- vcovHC(dlm2, type="HC1")
   print(coeftest(dlm2, dlm2$rse)) # Welch's t test on coefficients
}
####
#### FIB
FIB<-read.csv("FIB.csv",header=TRUE,sep=',')
str(FIB)

##create times-series dataframe
FIB1_ts<-ts(FIB,start = c(1951,1),end = c(2014,1),frequency=1)
summary(FIB1_ts)

## explore normality
shapiro.test(FIB$FIB) # p = .229, normal
hist(FIB$FIB) # close to normal
## contrast model
acf1<-acf(na.omit(FIB1_ts[,2]),type="p")
acf1
plot(FIB~Year,FIB)
cont_FIB<-lm(FIB~L(FIB,1),data=FIB)
summary(cont_FIB) 
cont_FIB$rse <- vcovHC(cont_FIB, type="HC1")
coeftest(cont_FIB, cont_FIB$rse) # Welch's t test on coefficients
summary(cont_FIB) # 0.083
# the default setting in R: options(contrasts = c(unordered = "contr.treatment", ordered = "contr.poly"))
## full model
options(contrasts = c(unordered = "contr.treatment", ordered = "contr.treatment")) # make sure factor() does not derive levels
full_FIB<-lm(FIB~L(FIB,1)+factor(L(EcoRef,2))+factor(L(DouCon,1))+factor(L(VesScrap,1))+factor(L(MarFisPol,1))+factor(NotrawlZone)+factor(L(EEZ,1))+factor(FisPerReg),data=FIB,na.action = na.exclude)
summary(full_FIB) # Adj R-sq: 0.342
# examine statistical effects
library(lmtest)
library(sandwich)
full_FIB$rse <- vcovHC(full_FIB, type="HC1")
coeftest(full_FIB, full_FIB$rse) # Welch's t test on coefficients
coef<-coeftest(full_FIB, full_FIB$rse)
write.csv(coef,"Coefficient full_FIB.csv")
# optimal model
## compare two models
cont_FIB<-lm(FIB~L(FIB,1),data=FIB[-1,])

anova(cont_FIB,full_FIB)

##### dominance analysis
library(dominanceanalysis)
dominanceAnalysis(full_FIB)

## simple DLMs
# FishAgree, LawEnf, SumMor
coef
str(FIB)
cpd2<-FIB[,-c(3,6)]
str(cpd2)
r<-c(2,1,1,1,0)
for(i in 3:7){
   l<-r[i-2]
   dlm2 <- lm(FIB~L(FIB,1)+factor(cpd2[,i]), data = cpd2)
   print(names(cpd2)[i]);print(l)
   dlm2$rse <- vcovHC(dlm2, type="HC1")
   print(coeftest(dlm2, dlm2$rse)) # Welch's t test on coefficients
}

##################################### the end for regression models #############################
#################### plot figures 3 - 6
### Fig. 3 plot for input-policy
setwd("C:/Users/xiong/Desktop/MS_policy analysis China")
## read data
df<-read.csv("China fishery and policy data2.csv",header=TRUE,sep=',')
str(df)
# set up data for each era
df<-df[,c(1:15)]
str(df)

time <- df$Year
time
N_tr <- df$N_tr
N_dv<-df$N_dt
H_tr<-df$H_tr
H_dv<-df$H_dt
RN_tr<-df$RN_bt
RN_tr
RH_tr<-df$RH_bt
tiff("Figure/Fig. 3 input.tif",width = 7.7, height = 7.5, units = 'in', res = 300)
par(mar=c(2,3.5, 1, 3),mfrow=c(2,1))## add extra space to right margin of plot within frame
plot(time,N_tr/1000,pch=16,axes=FALSE, ylim=c(0,80),xlab="", ylab="",
     xlim=c(1950,2020),type="o",cex=.8)
axis(2, col="black",las=1)  ## las=1 makes horizontal labels
lines(time, N_dv/1000, pch=17,  xlab="", ylab="", 
      axes=FALSE, type="o",cex=.8)## Plot the second plot and put axis scale on right
mtext("Total number of Chinese BTs (x 1000)",side=2,font=2,line=2.5)

box()
par(new=TRUE)## Allow a second plot on the same graph
max(na.omit(H_tr)/1000000)
plot(time, H_tr/1000000, pch=1,col="red", xlab="", ylim=c(0,8),xlim=c(1950,2020),ylab="", 
     axes=FALSE, cex=.8, type="o")## Plot the second plot and put axis scale on right
lines(time, H_dv/1000000, pch=2, xlab="", xlim=c(1950,2020),ylab="", 
      axes=FALSE, cex=.8, col="red",type="o")## Plot the second plot and put axis scale on right
mtext("Total hp of Chinese BTs (x 1,000,000 kW)",side=4, font=2,col="red",line=2)

Labels<-c("Economic reform\n1978",
          "Increasing\nreform speed\n1992",
          "Ocean Agenda\n1996",
          "Permit\nprovisions\n2004",
          "Double Control\n2003",
          "Interim measures\nof subsidy funds\n2010",
          "Accelerating\nfishery\nupgrading\n2013")
abline(v=c(1978,1992,1996,2004,2003,2010,2013),lty="dotted",col=c("purple","purple","darkgreen","purple","darkgreen","purple","blue"))
text(x=c(1979,1993,1995,2003,2003,2009,2012),y=c(1,1,2.5,1.5,5,7.3,4),labels = Labels,cex=.8, pos=c(2,2,4,4,4,4,4),col=c("purple","purple","darkgreen","purple","darkgreen","purple","blue"))

axis(4, col="red",col.axis="red",las=1)
axis(1,pretty(range(time)))## Draw the time axis
xtick<-seq(1950, 2020, by=1)
axis(side=1, at=xtick, tck=-0.01,labels = FALSE)
legend("topleft",legend=c("All waters","Distant waters beyond C4S","All waters","Distant waters beyond C4S"),
       cex=1,text.col=c("black","black","red","red"),lwd=1, lty=c(1,1,1,1), pch=c(16,17,1,2),col=c("black","black","red","red"))
mtext("a)",side=3,adj=0,cex=1,font=2)

## b)
str(df)
HpV_bt<-df$HpV_bt
HpV_dbt<-df$HpV_dbt
RN_bt<-df$RN_bt
RH_bt<-df$RH_bt
plot(time,RN_bt,pch=16,axes=FALSE, ylim=c(0,140),xlab="", ylab="",
     xlim=c(1950,2020),type="o",cex=.8)
axis(2, col="black",las=1)  ## las=1 makes horizontal labels
lines(time, RH_bt, pch=17,  xlab="", ylab="", 
      axes=FALSE, type="o",cex=.8)## Plot the second plot and put axis scale on right
mtext("% of BTs in Chinese fishing fleets",font=2,side=2,line=2.5)

box()
par(new=TRUE)## Allow a second plot on the same graph
plot(time,HpV_bt/100,pch=1,axes=FALSE, ylim=c(0,16),xlab="", ylab="",
     xlim=c(1950,2020),type="o",col="red",cex=.8)
lines(time, HpV_dbt/100, pch=2,  xlab="", ylab="", 
      axes=FALSE, cex=.8,col="red", type="o")## Plot the second plot and put axis scale on right
mtext("Mean hp of Chinese BTs (x100 kW)",side=4,font=2,col="red",line=2)
axis(4, col="red",col.axis="red",las=1)
axis(1,pretty(range(time)))## Draw the time axis
xtick<-seq(1950, 2020, by=1)
axis(side=1, at=xtick, tck=-0.01,labels = FALSE)
Labels<-c("No-trawl zone\n1955",
          "Economic reform\n1978",
          "Inshore-fisheries\nprotection\n1981",
          "Agenda 21\n1994",
          "Double Control\n1997",
          "Sino-Japanese\nFishery Agreement\n2000",
          "Sino-South Korean\nFishery Agreement\n2001",
          "Vessel\nbuyback\n2002",
          "Sino-Vietnamese\nFishery Agreement\n2004",
          "Permit\nprovisions\n2013")
abline(v=c(1955,1978,1981,1994,1997,2000,2001,2002,2004,2013),lty="dotted",col=c("darkgreen","purple","darkgreen","darkgreen","blue","darkgreen","darkgreen","purple","blue","purple"))
text(x=c(1954,1979,1982,1995,1996,1999,2000,2001,2003,2012),y=c(10,4,6,15,15,13,10.5,4,0.5,3.4),labels = Labels,cex=.8, pos=c(4,2,2,2,4,4,4,4,4,4),col=c("darkgreen","purple","darkgreen","darkgreen","blue","darkgreen","darkgreen","purple","blue","purple"))

axis(4, col="red",col.axis="red",las=1)
axis(1,pretty(range(time)))## Draw the time axis
xtick<-seq(1950, 2020, by=1)
axis(side=1, at=xtick, tck=-0.01,labels = FALSE)
legend("topleft",legend=c("All waters","Distant waters beyond C4S","All waters","Distant waters beyond C4S"),
       cex=1,text.col=c("black","black","red","red"),lwd=1, lty=c(1,1,1,1), pch=c(16,17,1,2),col=c("black","black","red","red"))
mtext("b)",side=3,adj=0,cex=1,font=2)
dev.off()


# Fig. S2.3 plot for output-policy
str(df)
Cbt<-df$C_bt
Cdv<-df$C_dbt
CpH<-df$CpH_bt #7.7
FIB<-df$FIB # 1.9
max(na.omit(Cdv/1000000))
par(mar=c(2, 3, 1, 3),mfrow=c(2,1))## add extra space to right margin of plot within frame
plot(time, Cbt/1000000,pch=16, xlab = "", ylab="", 
     xlim=c(1950,2020),ylim=c(0,10),type="o",col="black",las=1) ## Plot first set of data and draw its axis
lines(time, Cdv/1000000, pch=17,  xlab="", ylab="", 
      axes=FALSE, type="o", col="black")## Plot the second plot and put axis scale on right

mtext("Catch (Mt)",side=2,col="black",line=2)
par(new=T)
lines(time, CpH, pch=15,  xlab="", ylab="", 
      axes=FALSE, type="o", col="red")## Plot the second plot and put axis scale on right
axis(4, col="red",col.axis="red",las=1)
mtext("Catch per unit effort (t/(kW*year))",side=4,col="red",line=2)  
xtick<-seq(1950, 2020, by=1)
axis(side=1, at=xtick, tck=-0.01,labels = FALSE)
Labels<-c("No-trawl zone\n1980",
          "Developing\ndistant-water\nfisheries\n1985",
          "Summer\nmoratorium\n1995",
          "UNCLOS\nRatification\n1996",
          "Summer\nmoratorium\n1998",
          "Negative\nGrowth\n2000",
          "Protecting Fisheries\n2009",
          "Interim measures\nof subsidy funds\n2010")
abline(v=c(1980,1985,1995,1996,1998,2000,2009,2010),lty="dotted",
       col=c("darkgreen","purple","blue","blue","purple","darkgreen","purple","purple"))

text(x=c(1980,1985,1995,1996,1998,2000,2009,2010),
     y=c(3,6,8,9,6,4,2,6),labels = Labels,pos=c(2,2,2,4,4,4,4,4),cex=.8, col=c("darkgreen","purple","blue","blue","purple","darkgreen","purple","purple"))
legend("topleft",legend=c("All waters","Distant waters beyond C4S","All waters"),
       cex=1,text.col=c("black","black","red"),lwd=1, lty=c(1,1,1), pch=c(16,17,15),col=c("black","black","red"))
mtext("a)",side=3,adj=0,cex=1,font=2)

## b)
par(mar=c(2,3.5,1,1))
plot(time, FIB,pch=16, xlab = "", ylab="", 
     xlim=c(1950,2020),type="o",col="black",las=1) ## Plot first set of data and draw its axis

mtext("Fishing-in-balance index",side=2,col="black",line=2.5) 

Labels<-c("No-trawl zone\n1955",
          "No-trawl zone\n1957",
          "UNCLOS Ratification\n1996",
          "Double Control\n1997",
          "EEZ Law\n1998",
          "Vessel buyback\n2002")
abline(v=c(1955,1957,1996,1997,1998,2002),lty="dotted",
       col="darkgreen")
text(x=c(1955,1957,1997,1998,1998,2002),
     y=c(0.25,1.15,0.8,0.5,0.1,1),labels = Labels,pos=c(4,4,2,2,4,4), col="darkgreen")

xtick<-seq(1950, 2020, by=1)
axis(side=1, at=xtick, tck=-0.01,labels = FALSE)
mtext("b)",side=3,adj=0,cex=1,font=2)


### Fig. 4 get confidence interval
getwd()
tiff("Figure/Fig. 4 input and output.tif",width = 9, height = 6.5, units = 'in', res = 300)
pred_n_tr<-as.data.frame(na.omit(predict(full_n_tr,interval = "confidence")))
pred_n_tr
min<-min(min(pred_n_tr$lwr/1000),N_tr$N_tr/1000)
min
max<-max(max(pred_n_tr$upr/1000),N_tr$N_tr/1000)
max
## plot fit
par(mai=c(0.3,0.7,0.2,0.01),mfrow=c(4,3),oma=c(0,0,0,0),las=0)
plot(N_tr$N_tr/1000,pch = 19, xaxt='n',xlab="",ylab="",ylim=c(min,max),las=1)
axis(1,at=seq(1,68,5),labels=seq(1951,2018,5))## Draw the time axis
lines(x=seq(3,68,1),fitted(cont_n_tr)/1000,col="black",lty="solid")
lines(fitted(full_n_tr)/1000,col="red",lty="solid",lwd=2)
polygon(x=c(seq(3,68,1),rev(seq(3,68,1))),
        y=c(pred_n_tr$lwr/1000, rev(pred_n_tr$upr/1000)),
        col=adjustcolor("red",alpha.f = 0.2),border = NA)
mtext("Difference (x1000)",cex=0.8,side=2,line = 3,las=0)
text(x=2,y=-12,labels = "Adjusted R-squared: 0.268", pos=4, col="black")
text(x=2,y=-15,labels = "Adjusted R-squared: 0.790", pos=4, col="red")
text(x=53,y=-15,labels = "P < 0.001", pos=4, col="red")
mtext("a) Number of bottom trawlers (BTs)",cex=0.8,side = 3,adj = 0, line=0.1)
### Fig. Sc get confidence intervals
pred_h_tr<-as.data.frame(na.omit(predict(full_h_tr,interval = "confidence")))
pred_h_tr
min<-min(min(pred_h_tr$lwr)/10000,H_tr$H_tr/10000)
min
max<-max(max(pred_h_tr$upr)/10000,H_tr$H_tr/10000)
max
## plot fit
plot(H_tr$H_tr/10000,pch = 19, xaxt='n',xlab="",ylab="",ylim=c(min,max),las=1)
axis(1,at=seq(1,68,5),labels=seq(1951,2018,5))## Draw the time axis
lines(x=seq(4,68,1),fitted(cont_h_tr)/10000,col="black",lty="solid")
lines(fitted(full_h_tr)/10000,col="red",lty="solid",lwd=2)
polygon(x=c(seq(4,68,1),rev(seq(4,68,1))),
        y=c(pred_h_tr$lwr/10000, rev(pred_h_tr$upr/10000)),
        col=adjustcolor("red",alpha.f = 0.2),border = NA)
mtext("Difference (x 10,000 kW)",cex=0.8,side=2,line = 3)
mtext("b) Horsepower of BTs",cex=0.8,side = 3,adj = 0, line=0)
text(x=2,y=-25,labels = "Adjusted R-squared: 0.730", pos=4, col="black")
text(x=2,y=-40,labels = "Adjusted R-squared: 0.893", pos=4, col="red")
text(x=53,y=-40,labels = "P < 0.001", pos=4, col="red")
### Fig. Sd get confidence intervals
pred_H_dt<-as.data.frame(na.omit(predict(full_H_dt,interval = "confidence")))
pred_H_dt
min<-min(min(pred_H_dt$lwr)/10000,na.omit(H_dt$H_dt/10000))
min
max<-max(max(pred_H_dt$upr)/10000,na.omit(H_dt$H_dt/10000))
max
## plot fit
plot(H_dt$H_dt/10000,pch = 19, xaxt='n',xlab="",ylab="",ylim=c(min,max),las=1)
axis(1,at=seq(1,68,5),labels=seq(1951,2018,5))## Draw the time axis
lines(x=seq(36,68,1),fitted(cont_H_dt)/10000,col="black",lty="solid")
lines(fitted(full_H_dt)/10000,col="red",lty="solid",lwd=2)
polygon(x=c(seq(36,68,1),rev(seq(36,68,1))),
        y=c(pred_H_dt$lwr/10000, rev(pred_H_dt$upr/10000)),
        col=adjustcolor("red",alpha.f = 0.2),border = NA)
mtext("Difference (x 10,000 kW)",cex=0.8,side=2,line = 3)
mtext("c) Horsepower of distant-water BTs",cex=0.8,side = 3,adj = 0, line=0.1)
text(x=2,y=-15,labels = "Adjusted R-squared: 0.261", pos=4, col="black")
text(x=2,y=-25,labels = "Adjusted R-squared: 0.663", pos=4, col="red")
text(x=54,y=-25,labels = "P < 0.01", pos=4, col="red")
### Fig. Se get confidence intervals
pred_C_bt<-as.data.frame(na.omit(predict(full_C_bt,interval = "confidence")))
pred_C_bt
min<-min(min(pred_C_bt$lwr)/1000000,na.omit(C_bt$C_bt/1000000))
min
max<-max(max(pred_C_bt$upr)/1000000,na.omit(C_bt$C_bt/1000000))
max
## plot fit
plot(C_bt$C_bt/1000000,pch = 19, xaxt='n',xlab="",ylab="",ylim=c(min,max),las=1)
axis(1,at=seq(1,64,5),labels=seq(1951,2014,5))## Draw the time axis
lines(x=seq(2,64,1),fitted(cont_C_bt)/1000000,col="black",lty="solid")
lines(fitted(full_C_bt)/1000000,col="red",lty="solid",lwd=2)
polygon(x=c(seq(2,64,1),rev(seq(2,64,1))),
        y=c(pred_C_bt$lwr/1000000, rev(pred_C_bt$upr/1000000)),
        col=adjustcolor("red",alpha.f = 0.2),border = NA)
mtext("Difference (Mt)",cex=0.8,side=2,line = 3)
mtext("d) Catch by BTs",cex=0.8,side = 3,adj = 0, line=0.1)
text(x=2,y=-0.5,labels = "Adjusted R-squared: 0.930", pos=4, col="red")
text(x=50,y=-0.5,labels = "P < 0.001", pos=4, col="red")
### Fig. Sf get confidence intervals
pred_C_dbt<-as.data.frame(na.omit(predict(full_C_dbt,interval = "confidence")))
pred_C_dbt
min<-min(min(pred_C_dbt$lwr)/1000000,na.omit(C_dbt$C_dbt/1000000))
min
max<-max(max(pred_C_dbt$upr)/1000000,na.omit(C_dbt$C_dbt/1000000))
max
## plot fit
plot(C_dbt$C_dbt/1000000,pch = 19, xaxt='n',xlab="",ylab="",ylim=c(min,max),las=1)
axis(1,at=seq(1,64,5),labels=seq(1951,2014,5))## Draw the time axis
lines(x=seq(36,64,1),fitted(cont_C_dbt)/1000000,col="black",lty="solid")
lines(fitted(full_C_dbt)/1000000,col="red",lty="solid",lwd=2)
polygon(x=c(seq(36,64,1),rev(seq(36,64,1))),
        y=c(pred_C_dbt$lwr/1000000, rev(pred_C_dbt$upr/1000000)),
        col=adjustcolor("red",alpha.f = 0.2),border = NA)
mtext("Difference (Mt)",cex=0.8,side=2,line = 3)
mtext("e) Catch by distant-water BTs",cex=0.8,side = 3,adj = 0, line=0.1)
text(x=2,y=-0.7,labels = "Adjusted R-squared: 0.383", pos=4, col="red")
text(x=50,y=-0.7,labels = "P < 0.05", pos=4, col="red")
### Fig. Sg get confidence intervals
pred_HpV_bt<-as.data.frame(na.omit(predict(full_HpV_bt,interval = "confidence")))
pred_HpV_bt
min<-min(min(pred_HpV_bt$lwr)/100,HpV_bt$HpV_bt/100)
min
max<-max(max(pred_HpV_bt$upr)/100,HpV_bt$HpV_bt/100)
max
## plot fit
plot(HpV_bt$HpV_bt/100,pch = 19, xaxt='n',xlab="",ylab="",ylim=c(min,max),las=1)
axis(1,at=seq(1,68,5),labels=seq(1951,2018,5))## Draw the time axis
lines(x=seq(2,68,1),fitted(cont_HpV_bt)/100,col="black",lty="solid")
lines(fitted(full_HpV_bt)/100,col="red",lty="solid",lwd=2)
polygon(x=c(seq(2,68,1),rev(seq(2,68,1))),
        y=c(pred_HpV_bt$lwr/100, rev(pred_HpV_bt$upr/100)),
        col=adjustcolor("red",alpha.f = 0.2),border = NA)
mtext("Difference (x 100 kW)",cex=0.8,side=2,line = 3)
mtext("f) Mean horsepower of BTs",cex=0.8,side = 3,adj = 0, line=0.1)
text(x=2,y=-1.20,labels = "Adjusted R-squared: 0.452", pos=4, col="red")
text(x=53,y=-1.20,labels = "P < 0.001", pos=4, col="red")
### Fig. Sh get confidence intervals
pred_HpV_dbt<-as.data.frame(na.omit(predict(full_HpV_dbt,interval = "confidence")))
pred_HpV_dbt
min<-min(min(pred_HpV_dbt$lwr)/100,na.omit(HpV_dbt$HpV_dbt/100))
min
max<-max(max(pred_HpV_dbt$upr)/100,na.omit(HpV_dbt$HpV_dbt/100))
max
## plot fit
plot(HpV_dbt$HpV_dbt/100,pch = 19, xaxt='n',xlab="",ylab="",ylim=c(min,max),las=1)
axis(1,at=seq(1,68,5),labels=seq(1951,2018,5))## Draw the time axis
lines(x=seq(36,68,1),fitted(cont_HpV_dbt)/100,col="black",lty="solid")
lines(x=seq(36,68,1),fitted(opt_HpV_dbt)/100,col="blue",lty="solid",lwd=2)
lines(fitted(full_HpV_dbt)/100,col="red",lty="solid",lwd=2)
polygon(x=c(seq(36,68,1),rev(seq(36,68,1))),
        y=c(pred_HpV_dbt$lwr/100, rev(pred_HpV_dbt$upr/100)),
        col=adjustcolor("red",alpha.f = 0.2),border = NA)
mtext("Difference (x 100 kW)",cex=0.8,side=2,line = 3)
mtext("g) Mean horsepower of distant-water BTs",cex=0.8,side = 3,adj = 0, line=0.1)
text(x=2,y=-3,labels = "Adjusted R-squared: 0.220", pos=4, col="black")
text(x=2,y=-5,labels = "Adjusted R-squared: 0.392", pos=4, col="red")
text(x=2,y=-7,labels = "Adjusted R-squared: 0.501", pos=4, col="blue")
text(x=53,y=-5,labels = "P = 0.149", pos=4, col="red")
text(x=53,y=-7,labels = "P < 0.01", pos=4, col="blue")
### Fig. Si get confidence intervals
pred_CpH_bt<-as.data.frame(na.omit(predict(full_CpH_bt,interval = "confidence")))

min<-min(min(pred_CpH_bt$lwr)/1,na.omit(CpH_bt$CpH_bt/1))
min
max<-max(max(pred_CpH_bt$upr)/1,na.omit(CpH_bt$CpH_bt/1))
max
## plot fit
plot(CpH_bt$CpH_bt/1,pch = 19, xaxt='n',xlab="",ylab="",ylim=c(min-0.2,max),las=1)
axis(1,at=seq(1,64,5),labels=seq(1951,2014,5))## Draw the time axis
lines(x=seq(3,64,1),fitted(cont_CpH_bt)/1,col="black",lty="solid")
lines(fitted(full_CpH_bt)/1,col="red",lty="solid",lwd=2)
polygon(x=c(seq(3,64,1),rev(seq(3,64,1))),
        y=c(pred_CpH_bt$lwr, rev(pred_CpH_bt$upr)),
        col=adjustcolor("red",alpha.f = 0.2),border = NA)
mtext("Difference (t/kW)",cex=0.8,side=2,line = 3)
mtext("h) Catch per unit effort of BTs ",cex=0.8,side = 3,adj = 0, line=0.1)
text(x=2,y=-1.5,labels = "Adjusted R-squared: 0.245", pos=4, col="black")
text(x=2,y=-1.8,labels = "Adjusted R-squared: 0.410", pos=4, col="red")
text(x=50,y=-1.8,labels = "P < 0.05", pos=4, col="red")

###
pred_RN_bt<-as.data.frame(na.omit(predict(full_RN_bt,interval = "confidence")))
pred_RN_bt
min<-min(min(pred_RN_bt$lwr),RN_bt$RN_bt)
min
max<-max(max(pred_RN_bt$upr),RN_bt$RN_bt)
max
## plot fit
plot(RN_bt$RN_bt,pch = 19, xaxt='n',xlab="",ylab="",ylim=c(min,max),las=1)
axis(1,at=seq(1,68,5),labels=seq(1951,2018,5))## Draw the time axis
lines(x=seq(1,68,1),fitted(cont_RN_bt),col="black",lty="solid")
lines(fitted(full_RN_bt),col="red",lty="solid",lwd=2)
polygon(x=c(seq(3,68,1),rev(seq(3,68,1))),
        y=c(pred_RN_bt$lwr, rev(pred_RN_bt$upr)),
        col=adjustcolor("red",alpha.f = 0.2),border = NA)
mtext("Difference",cex=0.8,side=2,line = 3,las=0)
text(x=6,y=-30,labels = "Adjusted R-squared: 0.703", pos=4, col="red")
text(x=54,y=-30,labels = "P < 0.001", pos=4, col="red")
mtext("i) Percentage of BTs (by number)",cex=0.8,side = 3,adj = 0, line=0.1)
### Fig. Sxb
pred_RH_bt<-as.data.frame(na.omit(predict(full_RH_bt,interval = "confidence")))
pred_RH_bt
min<-min(min(pred_RH_bt$lwr),RH_bt$RH_bt)
min
max<-max(max(pred_RH_bt$upr),RH_bt$RH_bt)
max
## plot fit
plot(RH_bt$RH_bt,pch = 19, xaxt='n',xlab="",ylab="",ylim=c(min,max),las=1)
axis(1,at=seq(1,68,5),labels=seq(1951,2018,5))## Draw the time axis
lines(x=seq(3,68,1),fitted(cont_RH_bt),col="black",lty="solid")
lines(fitted(full_RH_bt),col="red",lty="solid",lwd=2)
polygon(x=c(seq(4,68,1),rev(seq(4,68,1))),
        y=c(pred_RH_bt$lwr, rev(pred_RH_bt$upr)),
        col=adjustcolor("red",alpha.f = 0.2),border = NA)
mtext("Difference",cex=0.8,side=2,line = 3,las=0)
text(x=6,y=-24,labels = "Adjusted R-squared: 0.164", pos=4, col="black")
text(x=6,y=-31,labels = "Adjusted R-squared: 0.525", pos=4, col="red")
text(x=53,y=-31,labels = "P < 0.001", pos=4, col="red")
mtext("j) Percentage of BTs (by horsepower)",cex=0.8,side = 3,adj = 0, line=0.1)
### Fig. Sxc
pred_RC_bt_eez<-as.data.frame(na.omit(predict(full_RC_bt_eez,interval = "confidence")))
pred_RC_bt_eez
min<-min(min(pred_RC_bt_eez$lwr),RC_bt_eez$RC_bt_eez)
min
max<-max(max(pred_RC_bt_eez$upr),RC_bt_eez$RC_bt_eez)
max
## plot fit
plot(RC_bt_eez$RC_bt_eez,pch = 19, xaxt='n',xlab="",ylab="",ylim=c(min,max),las=1)
axis(1,at=seq(1,68,5),labels=seq(1951,2018,5))## Draw the time axis
lines(x=seq(1,68,1),fitted(cont_RC_bt_eez),col="black",lty="solid")
lines(fitted(opt_RC_bt),col="blue",lty="solid",lwd=2)
lines(fitted(full_RC_bt_eez),col="red",lty="solid",lwd=2)
polygon(x=c(seq(1,68,1),rev(seq(1,68,1))),
        y=c(pred_RC_bt_eez$lwr, rev(pred_RC_bt_eez$upr)),
        col=adjustcolor("red",alpha.f = 0.2),border = NA)
mtext("Difference",cex=0.8,side=2,line = 3,las=0)
text(x=4,y=-11,labels = "Adjusted R-squared: 0.021", pos=4, col="red")
text(x=53,y=-11,labels = "P = 0.385", pos=4, col="red")
mtext("k) % BTF catch from the claimed EEZ",cex=0.8,side = 3,adj = 0, line=0.1)
### Fig. Sxd
pred_FIB<-as.data.frame(na.omit(predict(full_FIB,interval = "confidence")))
pred_FIB
min<-min(min(pred_FIB$lwr),FIB$FIB)
min
max<-max(max(pred_FIB$upr),FIB$FIB)
max
## plot fit
plot(FIB$FIB,pch = 19, xaxt='n',xlab="",ylab="",ylim=c(min,max),las=1)
axis(1,at=seq(1,68,5),labels=seq(1951,2018,5))## Draw the time axis
lines(x=seq(3,64,1),fitted(cont_FIB),col="black",lty="solid")
lines(fitted(full_FIB),col="red",lty="solid",lwd=2)
polygon(x=c(seq(3,64,1),rev(seq(3,64,1))),
        y=c(pred_FIB$lwr, rev(pred_FIB$upr)),
        col=adjustcolor("red",alpha.f = 0.2),border = NA)
mtext("Difference",cex=0.8,side=2,line = 3,las=0)
text(x=1,y=-0.18,labels = "Adjusted R-squared: 0.083", pos=4, col="black")
text(x=1,y=-0.24,labels = "Adjusted R-squared: 0.317", pos=4, col="red")
text(x=52,y=-0.24,labels = "P < 0.05", pos=4, col="red")
mtext("l) Fishing-in-balance index",cex=0.8,side = 3,adj = 0, line=0.1)
dev.off()



############# Fig. 5 plot AHC and MDS plot
install.packages("dendextend")
install.packages("vegan")
#
library(dendextend)
library(vegan)

normalize <- function(x) {
   return ((x - min(x)) / (max(x) - min(x)))
}
## read data
setwd("C:/Users/xiong/Desktop")
mydata<-read.csv("sample.csv",header=TRUE,sep=',')
str(mydata)
## prepare data for cluster analysis
mydata2 <- na.omit(mydata)[,-1] # listwise deletion of missing
str(mydata2)
# scale data by range (0,1)
mydata2<-as.data.frame(lapply(mydata2, normalize))
summary(mydata2)
## calculate bray distance
row.names(mydata2)<-mydata$X
library(vegan)
dm<- vegdist(mydata2, method = "bray")
fit <- hclust(dm, method="ward.D") # link method "ward.D2"
fit
## plot the dendragram and nMDS in one figure
tiff("cluster-bray-ward.D2.tif", width = 11.5, height = 7.2, units = 'in', res = 300)
split.screen(rbind(c(0,0.3,0, 1), c(0.32, 1, 0, 1))) 
screen(1)
par(mar=c(1,1,1,8),oma=c(0,0,0,0))
fit
#plot(fit,labels = mydata$X,main="",las=2, cex=0.8) # display dendogram
#plot(as.dendrogram(fit),  xlab = "Height",
#    horiz = TRUE, col = mycol[g$group])
#rect.dendrogram(as.dendrogram(fit), 4, border=c("darkgreen","red","blue","purple"), horiz = TRUE)
groups <- cutree(fit, k=3) # cut tree into 5 clusters
groups$withinss
as.dendrogram(fit) %>%
   set("labels_col", value = c("red", "blue","purple"), k=3) %>%
   set("branches_k_color", value = c("red", "blue","purple"),k = 3) %>%
   plot(horiz=TRUE, axes=F)
text(x=1.1,y=8,labels = "Group 1",pos=3,col="purple",font=2,cex=1)
text(x=1.1,y=5.8,labels = "Group 2",pos=3,col="blue",font=2,cex=1)
text(x=1.1,y=2.5,labels = "Group 3",pos=3,col="red",font=2,cex=1)
mtext("a)",side=3,adj=0,cex=1)
#dend_height<-0
#for (i in 2:13) dend_height[i]<-fit$height[i-1]
#plot(13:1, dend_height, type = 'b', xlab="# of clusters", ylab="Dendrogram Height")
### NMDS scaling with correspondence plotting
library(vegan)
dm<- vegdist(mydata2, method = "bray") # distance matrix by bray or jaccard
dm
# calculates squared Mahalanobis  
set.seed(123)
scores<-metaMDS(dm,k=2,trymax = 500)
scores(scores)
scores$stress # stress = 0.12
group <- cutree(fit, k=3) # change to fit2, fit3...
group
mycol=c("purple","blue","red")
g<-as.data.frame(group)
str(g)
g
mygroup<-sort(unique(g$group))
# ANOSIM statistic
dm.ano<-anosim(dm,group)
dm.ano # ANOSIM statistic R = 0.745, p = 0.001
summary(dm.ano)
#plot(dm.ano, las=1)
# plot nMDS
screen(2)
par(mar=c(4,4,1,0.5))
plot(scores(scores),pch=c(17,17,17,17,18,19,18,18,17,19,18,19,17),xlim=c(-0.56,0.33), ylim=c(-0.42,0.4), col=mycol[g$group], cex=1.3,las=1) # ylim=c(-0.5,0.4),xlim=c(-0.5,0.5)
#abline(v=0,lty="dashed",col="grey")
#abline(h=0,lty="dashed",col="grey")
set.seed(123)
vf <- envfit(scores(scores), mydata2, perm = 999)
vf
plot(vf, len=0.1,cex=.8,col = "black") # c("blue","red","blue","darkgreen","red","purple","red","purple","purple","purple","purple")
ordihull(scores(scores),draw = "line", lty = "dotted", g$group, col=mycol,show.groups=mygroup)
text(scores(scores)-0.025, labels=mydata$X, cex=.8,col = mycol[g$group]) 
text(x=.28,y=.4,labels = "Stress: 0.12")
mtext("b)",side=3,adj=0,cex=1)
dev.off()
################# Fig. 6 comparisons among policy factors, types/levels
## read data
setwd("C:/Users/xiong/Desktop/Postdoc/dylm")
mydata<-read.csv("policy_2D.csv",header=TRUE,sep=',')
mydata2<-read.csv("policy_2D2.csv",header=TRUE,sep=',')
head(mydata2)
## points with labels
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(reshape2)
library(grid)

mydata3<-melt(mydata2[,c(1,6:8)], id ="Category")
str(mydata3)
summary(mydata3)
head(mydata3)
mydata3$Category <- factor(mydata3$Category, levels=unique(mydata3$Category))
levels(mydata3$Category)<-c("IA","OP","BP","IC","OC","LE","FS","ML","CL")
tiff("Fig. 6 performance indices2.tiff",width = 11.3, height = 7.5, units = 'in', res = 300)
theme_set(theme_bw())
p1<-ggplot(mydata, aes(x= IBI, y= MDI))+
   geom_point()
p11<-p1+geom_text_repel(data=mydata, aes(label=Policy),cex=3)+labs(x = "Influence breadth index (%)", y = "Mean dominance index (%)")

p2<-ggplot(mydata2, aes(x= Examined*100, y= Influence*100))+
   geom_point()+
   coord_cartesian(xlim = c(0,70), ylim = c(0,60))+
   geom_abline(slope = 1,intercept = 0, col= "red")
p22<-p2+geom_text_repel(data=mydata2, aes(label=Category),cex=3)+labs(x = "Policy contribution index (%)", y = "Influence contribution index (%)")


p3<-ggplot(mydata2, aes(x= Used*100, y= Usefulness*100))+
   geom_point()+
   coord_cartesian(xlim = c(25,100), ylim = c(25,100))+
   geom_abline(slope = 1,intercept = 0, col= "red")
p33<-p3+geom_text_repel(data=mydata2, aes(label=Category),cex=3)+labs(x = "Potential usefulness index (%)", y = "Significant usefulness index (%)")

p4<-ggplot(mydata3, aes(fill=variable, y=value, x=Category))+
   geom_bar(position="stack", stat="identity")
p44<-p4+scale_fill_manual("Policy", labels = c("Examined", "Unimplemented", "Unexaminable"), 
                          values = c("Number.of.examined.policies" = "royalblue1", 
                                     "Number.of.unimplemented.policies"="tomato",
                                     "Number.of.unexaminable"="grey"))+
   labs(x = "Policy type or level", y = "Number")

figure <- ggarrange(p11, p22, p33, p44,
                    labels = c("a)", "b)", "c)", "d)"),
                    ncol = 2, nrow = 2, vjust = c(1,1,0.8,0.8),
                    font.label = list(size = 12, face = "plain", color ="black"))
figure


dev.off()

