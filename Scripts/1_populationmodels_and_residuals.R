# Population Models -------

library(INLA)

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
                    ,data=stackeddata)
### Fall seasonal model ####
Ft_seasonal <- inla(F_t ~ StR1 + StR2 + StR3 +
                      Ft_1R1 + Ft_1R2 + Ft_1R3 +
                      St_1R1 + St_1R2 + St_1R3 +
                      Ft_2R1 + Ft_2R2 + Ft_2R3 
                    # + f(time, mod="rw2", scale.model = TRUE,
                    #            hyper = list(theta = list(prior="pc.prec", param=c(u=1,0.01))))
                    # , family="gaussian"
                    ,data=stackeddata)
summary(stacked_FallR3)

### Spring Yearly model ####
St_year <- inla(S_t ~ 
                  St_1R1 + St_1R2 + St_1R3 +
                  St_2R1 + St_2R2 + St_2R3 
                # + f(time, mod="rw2", scale.model = TRUE,
                #            hyper = list(theta = list(prior="pc.prec", param=c(u=1,0.01))))
                # , family="gaussian"
                ,data=stackeddata)
summary(stacked_ySpringR3)

### Fall Yearly model ####
Ft_year <- inla(F_t ~ 
                  Ft_1R1 + Ft_1R2 + Ft_1R3 +
                  Ft_2R1 + Ft_2R2 + Ft_2R3 
                # + f(time, mod="rw2", scale.model = TRUE,
                #            hyper = list(theta = list(prior="pc.prec", param=c(u=1,0.01))))
                # , family="gaussian"
                ,data=stackeddata)

### List of All models ####
popmodels <- list(St_seasonal=St_seasonal,St_year=St_year,Ft_seasonal=Ft_seasonal,Ft_year=Ft_year)



## Obtain Residuals ####