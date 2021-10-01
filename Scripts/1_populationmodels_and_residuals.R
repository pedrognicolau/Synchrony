# Population Models -------

library(INLA)
stackeddata <- readRDS("Data/stacked_data.rds")

# expanded stacked data with extra variables
stackeddata2 <- readRDS("~/OneDrive - UiT Office 365/Vole Synchrony/data/stacked_data2.rds")

# growth rates instead
GR_data0 <- readRDS("Vole Synchrony/data/log_centered_INLACR_growthrates.rds")
## Hansen 1999 Models =======
### Spring seasonal model ####
St_seasonal <- inla(S_t ~ 
                      Ft_1R1 + Ft_1R2 + Ft_1R3 +
                      St_1R1 + St_1R2 + St_1R3 +
                      Ft_2R1 + Ft_2R2 + Ft_2R3 +
                      St_2R1 + St_2R2 + St_2R3 
                    # + f(time, mod="rw2", scale.model = TRUE,
                    #            hyper = list(theta = list(prior="pc.prec", param=c(u=1,0.01))))
                    #, family="gaussian"
                    ,
                    control.predictor = list(compute = TRUE),
                    control.compute = list(config = TRUE),
                    data=stackeddata)

### Fall seasonal model ####
Ft_seasonal <- inla(F_t ~ StR1 + StR2 + StR3 +
                      Ft_1R1 + Ft_1R2 + Ft_1R3 +
                      St_1R1 + St_1R2 + St_1R3 +
                      Ft_2R1 + Ft_2R2 + Ft_2R3 
                    # + f(time, mod="rw2", scale.model = TRUE,
                    #            hyper = list(theta = list(prior="pc.prec", param=c(u=1,0.01))))
                    # , family="gaussian"
                    ,
                    control.predictor = list(compute = TRUE),
                    control.compute = list(config = TRUE),
                    data=stackeddata)

### Spring Yearly model ####
St_year <- inla(S_t ~ 
                  St_1R1 + St_1R2 + St_1R3 +
                  St_2R1 + St_2R2 + St_2R3 
                # + f(time, mod="rw2", scale.model = TRUE,
                #            hyper = list(theta = list(prior="pc.prec", param=c(u=1,0.01))))
                # , family="gaussian"
                ,
                control.predictor = list(compute = TRUE),
                control.compute = list(config = TRUE),
                data=stackeddata)

### Fall Yearly model ####
Ft_year <- inla(F_t ~ 
                  Ft_1R1 + Ft_1R2 + Ft_1R3 +
                  Ft_2R1 + Ft_2R2 + Ft_2R3 
                # + f(time, mod="rw2", scale.model = TRUE,
                #            hyper = list(theta = list(prior="pc.prec", param=c(u=1,0.01))))
                # , family="gaussian"
                ,
                control.predictor = list(compute = TRUE),
                control.compute = list(config = TRUE),
                data=stackeddata)

### List of All models ####
popmodels <- list(St_seasonal=St_seasonal,St_year=St_year,Ft_seasonal=Ft_seasonal,Ft_year=Ft_year)



## Obtain Residuals ####

source("Scripts/0_Important_functions.R")

# obtain residuals using Bayesian sample
St_seas_res <- matrix(residual_samples(stackeddata$S_t, St_seasonal, ns=200, mean = TRUE),ncol=19)
Ft_seas_res <- matrix(residual_samples(stackeddata$F_t, Ft_seasonal, ns=200, mean = TRUE),ncol=19)
St_year_res <- matrix(residual_samples(stackeddata$S_t, St_year, ns=200, mean = TRUE),ncol=19)
Ft_year_res <- matrix(residual_samples(stackeddata$F_t, St_year, ns=200, mean = TRUE),ncol=19)

popres <- list(St_seasonal=St_seas_res,Ft_seasonal=Ft_seas_res,St_year=St_year_res,Ft_year=Ft_year_res)

saveRDS(popres,"Data/popmodels_residuals.rds")

# combmat <- combn(19,2)
# library(dplyr)
# for(i in 1:ncol(combmat))
#     {
#     newdat <- data.frame(st_res_mat[,combmat[,i]])
#     nd2 <- arrange(newdat,X1)
#     plot(nd2, type="b", pch=19)
#     lines(smooth.spline(nd2), pch=19, col=2, lty=2)
# }
# 
# difmat <- cor(st_res_mat)-cor(st_res_mat, method = "spearman")
# par(mfrow=c(1,1))
# library(fields)
# image.plot(difmat)


