# library(devtools)
# install_github("andrewzm/STRbook")

library("sp")
library("spacetime")
library("ggplot2")
library("dplyr")
library("gstat")
library("RColorBrewer")
library("STRbook")
library("tidyr")


data("STObj3", package = "STRbook")
STObj4 <- STObj3[, "1993-07-01::1993-07-31"]
vv <- variogram(object = z ~ 1 + lat, # fixed effect component
                data = STObj4, # July data
                width = 80, # spatial bin (80 km)
                cutoff = 1000, # consider pts < 1000 km apart
                tlags = 0.01:6.01) # 0 days to 6 days



sepVgm <- vgmST(stModel = "separable",
                space = vgm(10, "Exp", 400, nugget = 0.1),
                time = vgm(10, "Exp", 1, nugget = 0.1),
                sill = 20)
sepVgm <- fit.StVariogram(vv, sepVgm)



metricVgm <- vgmST(stModel = "metric",
                   joint = vgm(100, "Exp", 400, nugget = 0.1),
                   sill = 10,
                   stAni = 100)
metricVgm <- fit.StVariogram(vv, metricVgm)

metricMSE <- attr(metricVgm, "optim")$value
sepMSE <- attr(sepVgm, "optim")$value


plot(vv, list(sepVgm, metricVgm),main = "Semi-variance")


spat_pred_grid <- expand.grid(
  lon = seq(-100, -80, length = 20),
  lat = seq(32, 46, length = 20)) %>%
  SpatialPoints(proj4string = CRS(proj4string(STObj3)))
gridded(spat_pred_grid) <- TRUE

temp_pred_grid <- as.Date("1993-07-01") + seq(3, 28, length = 6)

DE_pred <- STF(sp = spat_pred_grid, # spatial part
               time = temp_pred_grid) # temporal part


STObj5 <- as(STObj4[, -14], "STIDF") # convert to STIDF
STObj5 <- subset(STObj5, !is.na(STObj5$z)) # remove missing data

pred_kriged <- krigeST(z ~ 1 + lat, # latitude trend
                       data = STObj5, # data set w/o 14 July
                       newdata = DE_pred, # prediction grid
                       modelList = sepVgm, # semivariogram
                       computeVar = TRUE) # compute variances

color_pal <- rev(colorRampPalette(brewer.pal(11, "Spectral"))(16))

stplot(pred_kriged,
       main = "Predictions (degrees Fahrenheit)",
       layout = c(3, 2),
       col.regions = color_pal)

pred_kriged$se <- sqrt(pred_kriged$var1.var)
stplot(pred_kriged[, , "se"],
       main = "Prediction std. errors (degrees Fahrenheit)",
       layout = c(3, 2),
       col.regions = color_pal)



#Mab 4.2 Spatio Temporal Basis Functions with FRK

library("FRK")

data("STObj3", package = "STRbook") # load STObj3
STObj4 <- STObj3[, "1993-07-01::1993-07-31"] # subset time
STObj5 <- as(STObj4[, -14], "STIDF") # omit t = 14
STObj5 <- subset(STObj5, !is.na(STObj5$z)) # remove NAs

BAUs <- auto_BAUs(manifold = STplane(), # ST field on the plane
                  type = "grid", # gridded (not "hex")
                  data = STObj5, # data
                  cellsize = c(1, 0.75, 1), # BAU cell size
                  convex = -0.12, # hull extension
                  tunit = "days") # time unit is "days"

plot(as(BAUs[, 1], "SpatialPixels")) # plot pixel BAUs
plot(SpatialPoints(STObj5),
     add = TRUE, col = "red") # plot data points




BAUs_hex <- auto_BAUs(manifold = STplane(), # model on the plane
                      type = "hex", # hex (not "grid")
                      data = STObj5, # data
                      cellsize = c(1, 0.75, 1), # BAU cell size
                      nonconvex_hull = FALSE, # convex hull
                      tunit = "days") # time unit is "days"


plot(as(BAUs_hex[, 1], "SpatialPolygons"))


G_spatial <- auto_basis(manifold = plane(), # fns on plane
                        data = as(STObj5, "Spatial"), # project
                        nres = 2, # 2 res.
                        type = "bisquare", # bisquare.
                        regular = 0) # irregular


t_grid <- matrix(seq(1, 31, length = 20))

G_temporal <- local_basis(manifold = real_line(), # fns on R1
                          type = "bisquare", # bisquare
                          loc = t_grid, # centroids
                          scale = rep(2, 20)) # aperture par.


G <- TensorP(G_spatial,  G_temporal) # take the tensor product
BAUs$fs = 1
STObj5$std <- sqrt(0.049)
f <- z ~ lat + 1

S <- FRK(f = f, # formula
         data = list(STObj5), # (list of) data
         basis = G, # basis functions
         BAUs = BAUs, # BAUs
         n_EM = 3, # max. no. of EM iterations
         tol = 0.01) # tol. on change in log-likelihood
grid_BAUs <- predict(S)




data("NOAA_df_1990", package = "STRbook") # load NOAA data
NOAA_sub <- filter(NOAA_df_1990, # filter data to only
                   year == 1993 & # contain July 1993
                     month == 7 &
                     proc == "Tmax") # and just max. temp.
NOAA_sub_for_STdata <- NOAA_sub %>%
  transmute(ID = as.character(id),
            obs = z,
            date = date)

covars <- dplyr::select(NOAA_sub, id, lat, lon) %>%
  unique() %>%
  dplyr::rename(ID = id) # createSTdata expects "ID"



library("INLA")
library("dplyr")
library("tidyr")
library("ggplot2")
library("STRbook")
data("MOcarolinawren_long", package = "STRbook")

coords <- unique(MOcarolinawren_long[c("loc.ID", "lon", "lat")])
boundary <- inla.nonconvex.hull(as.matrix(coords[, 2:3]))

MOmesh <- inla.mesh.2d(boundary = boundary,
                       max.edge = c(0.8, 1.2), # max. edge length
                       cutoff = 0.1) # min. edge length


plot(MOmesh, asp = 1, main = "")
lines(coords[c("lon", "lat")], col = "red", type = "p")


spde <- inla.spde2.pcmatern(mesh = MOmesh,
                            alpha = 2,
                            prior.range = c(1, 0.01),
                            prior.sigma = c(4, 0.01))


n_years <- length(unique(MOcarolinawren_long$t))
n_spatial <- MOmesh$n
s_index <- inla.spde.make.index(name = "spatial.field",
                                n.spde = n_spatial,
                                n.group = n_years)


coords.allyear <- MOcarolinawren_long[c("lon", "lat")] %>%
  as.matrix()
PHI <- inla.spde.make.A(mesh = MOmesh,
                        loc = coords.allyear,
                        group = MOcarolinawren_long$t,
                        n.group = n_years)



## First stack: Estimation
n_data <- nrow(MOcarolinawren_long)
stack_est <- inla.stack(
  data = list(cnt = MOcarolinawren_long$cnt),
  A = list(PHI, 1),
  effects = list(s_index,
                 list(Intercept = rep(1, n_data))),
  tag = "est")

df_pred <- data.frame(lon = rep(MOmesh$loc[,1], n_years),
                      lat = rep(MOmesh$loc[,2], n_years),
                      t = rep(1:n_years, each = MOmesh$n))
n_pred <- nrow(df_pred)
PHI_pred <- Diagonal(n = n_pred)

## Second stack: Prediction
stack_pred <- inla.stack(
  data = list(cnt = NA),
  A = list(PHI_pred, 1),
  effects = list(s_index,
                 list(Intercept = rep(1, n_pred))),
  tag = "pred")

stack <- inla.stack(stack_est, stack_pred)

## PC prior on rho
rho_hyper <- list(theta=list(prior = 'pccor1',
                             param = c(0, 0.9)))
## Formula
formula <- cnt ~ -1 + Intercept +
  f(spatial.field,
    model = spde,
    group = spatial.field.group,
    control.group = list(model = "ar1",
                         hyper = rho_hyper))


output <- inla(formula,
               data = inla.stack.data(stack, spde = spde),
               family = "nbinomial",
               control.predictor = list(A = inla.stack.A(stack),
                                        compute = TRUE))

data("INLA_output", package = "STRbook")


output.field <- inla.spde2.result(inla = output,
                                  name = "spatial.field",
                                  spde = spde,
                                  do.transf = TRUE)
## plot p(beta0 | Z)
plot(output$marginals.fix$Intercept,
     type = 'l',
     xlab = expression(beta[0]),
     ylab = expression(p(beta[0]*"|"*Z)))
## plot p(rho | Z)

plot(output$marginals.hyperpar$`GroupRho for spatial.field`,
     type = 'l',
     xlab = expression(rho),
     ylab = expression(p(rho*"|"*Z)))

plot(output.field$marginals.variance.nominal[[1]],
     type = 'l',
     xlab = expression(sigma^2),
     ylab = expression(p(sigma^2*"|"*Z)))
## plot p(range | Z)
plot(output.field$marginals.range.nominal[[1]],
     type = 'l',
     xlab = expression(l),
     ylab = expression(p(l*"|"*Z)))

index_pred <- inla.stack.index(stack, "pred")$data
lp_mean <- output$summary.fitted.values$mean[index_pred]
lp_sd <- output$summary.fitted.values$sd[index_pred]

grid_locs <- expand.grid(
  lon = seq(min(MOcarolinawren_long$lon) - 0.2,
            max(MOcarolinawren_long$lon) + 0.2,
            length.out = 80),
  lat = seq(min(MOcarolinawren_long$lat) - 0.2,
            max(MOcarolinawren_long$lat) + 0.2,
            length.out = 80))

proj.grid <- inla.mesh.projector(MOmesh,
                                 xlim = c(min(MOcarolinawren_long$lon) - 0.2,
                                          max(MOcarolinawren_long$lon) + 0.2),
                                 ylim = c(min(MOcarolinawren_long$lat) - 0.2,
                                          max(MOcarolinawren_long$lat) + 0.2),
                                 dims = c(80, 80))

pred <- sd <- NULL
for(i in 1:n_years) {
  ii <- (i-1)*MOmesh$n + 1
  jj <- i*MOmesh$n
  pred[[i]] <- cbind(grid_locs,
                     z = c(inla.mesh.project(proj.grid,
                                             lp_mean[ii:jj])),
                     t = i)
  sd[[i]] <- cbind(grid_locs,
                   z = c(inla.mesh.project(proj.grid,
                                           lp_sd[ii:jj])),
                   t = i)
}


proj.grid <- inla.mesh.projector(MOmesh,
                                 xlim = c(min(MOcarolinawren_long$lon) - 0.2,
                                          max(MOcarolinawren_long$lon) + 0.2),
                                 ylim = c(min(MOcarolinawren_long$lat) - 0.2,
                                          max(MOcarolinawren_long$lat) + 0.2),
                                 dims = c(80, 80))

pred <- sd <- NULL
for(i in 1:n_years) {
  ii <- (i-1)*MOmesh$n + 1
  jj <- i*MOmesh$n
  pred[[i]] <- cbind(grid_locs,
                     z = c(inla.mesh.project(proj.grid,
                                             lp_mean[ii:jj])),
                     t = i)
  sd[[i]] <- cbind(grid_locs,
                   z = c(inla.mesh.project(proj.grid,
                                           lp_sd[ii:jj])),
                   t = i)
}

pred <- do.call("rbind", pred) %>% filter(!is.na(z))
sd <- do.call("rbind", sd) %>% filter(!is.na(z))





font_add_google("Aptos")

library(legocolors)








