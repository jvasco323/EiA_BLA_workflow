
library(splines)
library(Metrics)
library(ggplot2)
library(reshape2)

rm(list=ls())

raw_data <- read.csv('D:/# Jvasco/# Portfolio/Curriculum/6. CIMMYT Scientist/IDT Excellence in Agronomy/1-yield-gap-decomposition/external-consultant/tomas tenreiro/bla markdown final/BLA_Fermont_et.al.csv')
raw_data <- subset(raw_data, yld_t0 < 30)
names(raw_data)

# remover NA

# ------------------------------------------------------------------------------
# introduction


# ------------------------------------------------------------------------------
# 1) documentar funcao para identificacao dos pontos de fronteira

give_me_blp <- function(df, xvar){
   # As an example, we select the X_variable
   df <- df[,c("yld_t0", xvar)]
   # We modify column names for generic use of this script
   colnames(df)[1] = "Y"
   colnames(df)[2] = "X"
   df <- subset(df, X>0, select=c(Y, X)) # define here the X_variable of interest
   # 'NULL' values are excluded to avoid data transformation problems and calculation failures 
   # We correct NA values for both the Y and X 
   df$Y[is.na(df$Y)] <- mean(df$Y, na.rm=T)
   df$X[is.na(df$X)] <- mean(df$X, na.rm=T)
   # Fitting a BLA according to Fermont et al. (2009)
   x_subset_0.1 <- subset(df, X <= quantile(X, 0.1))
   x_subset_0.2 <- subset(df, X > quantile(X, 0.1) & X <= quantile(X, 0.2))
   x_subset_0.3 <- subset(df, X > quantile(X, 0.2) & X <= quantile(X, 0.3))
   x_subset_0.4 <- subset(df, X > quantile(X, 0.3) & X <= quantile(X, 0.4))
   x_subset_0.5 <- subset(df, X > quantile(X, 0.4) & X <= quantile(X, 0.5))
   x_subset_0.6 <- subset(df, X > quantile(X, 0.5) & X <= quantile(X, 0.6))
   x_subset_0.7 <- subset(df, X > quantile(X, 0.6) & X <= quantile(X, 0.7))
   x_subset_0.8 <- subset(df, X > quantile(X, 0.7) & X <= quantile(X, 0.8))
   x_subset_0.9 <- subset(df, X > quantile(X, 0.8) & X <= quantile(X, 0.9))
   x_subset_1.0 <- subset(df, X > quantile(X, 0.9) & X <= quantile(X, 1.0))
   # Select ymax to define boundary points subset
   BLP_subset_0.0 <- subset(x_subset_0.1, X == min(X))
   BLP_subset_0.1 <- subset(x_subset_0.1, Y == max(Y))
   BLP_subset_0.2 <- subset(x_subset_0.2, Y == max(Y))
   BLP_subset_0.3 <- subset(x_subset_0.3, Y == max(Y))
   BLP_subset_0.4 <- subset(x_subset_0.4, Y == max(Y))
   BLP_subset_0.5 <- subset(x_subset_0.5, Y == max(Y))
   BLP_subset_0.6 <- subset(x_subset_0.6, Y == max(Y))
   BLP_subset_0.7 <- subset(x_subset_0.7, Y == max(Y))
   BLP_subset_0.8 <- subset(x_subset_0.8, Y == max(Y))
   BLP_subset_0.9 <- subset(x_subset_0.9, Y == max(Y))
   BLP_subset_1.0 <- subset(x_subset_1.0, Y == max(Y))
   # Other option is to select the 0.95th quantile y-values to define the boundary points subset
   # BLP_subset_0.0 <- subset(x_subset_0.1, X == min(X)) 
   # BLP_subset_0.1 <- subset(x_subset_0.1, Y > quantile(Y, 0.95))
   # BLP_subset_0.2 <- subset(x_subset_0.2, Y > quantile(Y, 0.95))
   # BLP_subset_0.3 <- subset(x_subset_0.3, Y > quantile(Y, 0.95))
   # BLP_subset_0.4 <- subset(x_subset_0.4, Y > quantile(Y, 0.95))
   # BLP_subset_0.5 <- subset(x_subset_0.5, Y > quantile(Y, 0.95))
   # BLP_subset_0.6 <- subset(x_subset_0.6, Y > quantile(Y, 0.95))
   # BLP_subset_0.7 <- subset(x_subset_0.7, Y > quantile(Y, 0.95))
   # BLP_subset_0.8 <- subset(x_subset_0.8, Y > quantile(Y, 0.95))
   # BLP_subset_0.9 <- subset(x_subset_0.9, Y > quantile(Y, 0.95))
   # BLP_subset_1.0 <- subset(x_subset_1.0, Y > quantile(Y, 0.98))
   # bind all subsets together
   BLP_df <- rbind(BLP_subset_0.0, BLP_subset_0.1, BLP_subset_0.2, BLP_subset_0.3,
                   BLP_subset_0.4, BLP_subset_0.5, BLP_subset_0.6,
                   BLP_subset_0.7, BLP_subset_0.8, BLP_subset_0.9,
                   BLP_subset_1.0)
   return(BLP_df)
 }

# ------------------------------------------------------------------------------
# 2) estimar boundary lines + plotting

#variaveis <- c("pH_soil", "SOC", "totN_soil", "P_soil", "K_soil", "CBB_T0_9",
               #"CGM_T0_9", "RF_tot", "RF0_6M", "RF_.d", "clay", "days_to_harvest")

variaveis <- c("pH_soil", "totN_soil", "P_soil", "clay")


blp_new <- c()
raw_data_new <- c()

for(i in unique(variaveis)){
  
  print(i)
  
  # selecionar raw_data
  raw_data_subset <- raw_data[,c("year", "Trial", "Site", "Farm.no.", "yld_t0", i)]
  colnames(raw_data_subset)[5] = "Y"
  colnames(raw_data_subset)[6] = "X"
  raw_data_subset$variable <- i
  
  # 1) estimar pontos de fronteira
  data.bla <- give_me_blp(raw_data, i)
  data.bla$variable <- i
  
  # 2) estimar regressao
  model <- glm(Y ~ ns(X, df = 2), data = data.bla)
  # print(model)
  
  # 3) predict y_max
  # boundary points only
  data.bla$y_pred <- predict(model, newdata = data.bla)
  # raw data
  raw_data_subset$y_pred <- predict(model, newdata = raw_data_subset)

  # ERROR assessment
  print(rmse(data.bla$Y, data.bla$y_pred))
  
  # for plotting
  # (para melhorar, meter plots fora do loop)
  plot_lims <- range(raw_data_subset$X, na.rm=T)
  plot.grid <- seq(from=plot_lims[1], to=plot_lims[2])
  pred <- predict(model, newdata=list(X=plot.grid), se=T)
  spline.d <- as.data.frame(spline(data.bla$X, data.bla$y_pred))
  xmin <- min(raw_data_subset$X)
  xmax <- max(raw_data_subset$X)
  ymax <- max(raw_data_subset$Y)
  g1 <- ggplot(data.bla, aes(x=X, y=y_pred)) + 
    geom_point(size=6, alpha=0.6, color="blue") +
    geom_point(data=data.bla, aes(x=X, y=Y), size=3, alpha=0.9, color="red") +
    geom_point(data=raw_data_subset, aes(x=X, y=Y)) +
    coord_cartesian(xlim=c(xmin-2, xmax), ylim=c(0, ymax+5)) + 
    labs(title="BOUNDARY LINE ", subtitle="Spline method - Shatar & Mcbratney. (2004)", y="YIELD_tha",
         x="X", caption="[BL from 10 groups]") + 
    theme_bw() + 
    theme(text = element_text(size=15)) +
    geom_line(data = spline.d, aes(x = x, y = y), alpha=0.6, color="blue", size=2)
  print(g1)
  
  # final) juntar tudo
  blp_new <- rbind(blp_new, data.bla)
  raw_data_new <- rbind(raw_data_new, raw_data_subset)
  
  }

# ------------------------------------------------------------------------------
# 3) yield gap decomposition

# reshape df
raw_data_ygd <- dcast(raw_data_new, year + Trial + Site + Farm.no. + Y ~ variable, value.var='y_pred') 
raw_data_ygd$y_pred_min = apply(raw_data_ygd[6:9], 1, min) 
raw_data_ygd$factor_limitante = names(raw_data_ygd[6:9])[apply(raw_data_ygd[6:9], 1, which.min)]

# decompose yg
raw_data_ygd$IYG <- max(raw_data_ygd$Y, na.rm=T) - raw_data_ygd$y_pred_min
raw_data_ygd$UYG <- raw_data_ygd$y_pred_min - raw_data_ygd$Y


plot(raw_data_ygd$y_pred_min, raw_data_ygd$Y, xlim=c(0, 35), ylim=c(0, 35))
abline(a=0, b=1)
abline(a=0, b=0.5)


ggplot(raw_data_ygd, aes(x=Y, y=raw_data_ygd$y_pred_min, ymax = raw_data_ygd$y_pred_min+IYG)) + 
  geom_point(size=2.5, alpha=0.6, color="blue") +
  geom_abline() +
  geom_pointrange(aes(ymin = y_pred_min, ymax = y_pred_min+IYG), data = raw_data_ygd, width = 0.2, size=0.1, linetype='dashed', col = 'red')+
  geom_pointrange(aes(ymin = y_pred_min-UYG, ymax = y_pred_min), data = raw_data_ygd, width = 0.2, size=0.1, linetype='dashed', col = 'green')+
  # Coords and axis limits
  geom_hline(yintercept=max(raw_data_ygd$Y), size=1, linetype='dashed', col = 'red')+
  geom_hline(yintercept=0, linetype='dotted', col = 'red')+
  geom_vline(xintercept=0, linetype='dotted', col = 'red')+
  coord_cartesian(xlim=c(0, max(raw_data_ygd$Y)), ylim=c(0, max(raw_data_ygd$Y))) + 
  labs(title="Yield gap decomposition", subtitle="Following the method of Wairegi et al. (2010)", y="Predicted crop yield Y [t/ha]",
       x="Observed crop yield X [t/ha]", caption="") + theme_bw() + theme(text = element_text(size=15))

# ------------------------------------------------------------------------------
# recommendations / limitations


