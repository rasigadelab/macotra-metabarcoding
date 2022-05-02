library(data.table)
library(readxl)
library(dplyr)


d <- fread("../data/metadata_macotra_batch1_wSAabund.tsv", sep="\t", header = TRUE)
# Note : 8 carriers, 9 non carriers

names(d)

d$days_after_start_of_treatment
d$sample_CFU
d$SA_abund
d$carriage
d$participant_id


d[, carrier := (carriage == 1)]
d[, staphCFU := log10(sample_CFU + 1)]
d[, days := days_after_start_of_treatment]
d[, patient := factor(participant_id)]
d[, SA_prop := d$SA_abund/100]

d <- as.data.table(d %>% mutate(month = round(days/30.417, digit=1)))

cols <- c(rgb(0,0,0.8,0.15), rgb(0.8, 0, 0, 0.15))
cols_med <- c(rgb(0,0,0.8,0.4), rgb(0.8, 0, 0, 0.4))
cols_dark <- c(rgb(0,0,0.8,0.9), rgb(0.8, 0, 0, 0.9))


# Mean and CI of staphylococcus aureus log CFU

##plot(0, type = "n", main=expression( "Dynamics of " *italic(S.aureus)* " CFU in time"), xlim = range(d$month), xlab = "month after decolonization",  ylim = range(d$staphCFU), ylab=expression( "Log 10 of " *italic(S.~aureus)* " CFU in time"))

par(mar=c(5,6,4,1)+.1)
plot(0, type = "n", xlim = range(d$month), xlab = "Months after decolonization",  ylim = range(d$staphCFU), ylab=expression( "Log 10 of " *italic(S.~aureus)* " CFU in time"), cex.lab=1.3, cex.axis=1.3)
mtext("A", side = 3, adj = 0, line = 2, cex = 2, font= 2)

for(patient_current in levels(d$patient)) {
  dsub <- d[patient == patient_current]
  lines(dsub$month, dsub$staphCFU, col = cols_med[dsub$carrier[1] + 1], lty = 2)
}

legend("topright", bty = "n", title = "", fill = cols_med, legend = c("Non-carrier", "Carrier"), cex = 1.3)


cfumed <- d[, .(cfu_mean = mean(staphCFU), cfu_cilo = mean(staphCFU) - 1.96*sd(staphCFU)/sqrt(8), cfu_cihi = mean(staphCFU) + 1.96*sd(staphCFU)/sqrt(8),
                cfu_median = median(staphCFU), cfu_q25 = quantile(staphCFU, 0.25), cfu_q75 = quantile(staphCFU, 0.75)), by = .(month, carrier)  ]

cfumed[cfu_cilo < 0, cfu_cilo := 0]

x <- cfumed[carrier == TRUE]
lines(x$month, x$cfu_mean, col = cols_dark[2], lwd=2)
polygon(c(x$month, rev(x$month), x$month[1]), c(x$cfu_cilo, rev(x$cfu_cihi), x$cfu_cilo[1]), border = NA, col = cols[2])

cfumed <- d[, .(cfu_mean = mean(staphCFU), cfu_cilo = mean(staphCFU) - 1.96*sd(staphCFU)/sqrt(9), cfu_cihi = mean(staphCFU) + 1.96*sd(staphCFU)/sqrt(9),
                cfu_median = median(staphCFU), cfu_q25 = quantile(staphCFU, 0.25), cfu_q75 = quantile(staphCFU, 0.75)), by = .(month, carrier)  ]
cfumed[cfu_cilo < 0, cfu_cilo := 0]

x <- cfumed[carrier == FALSE]
lines(x$month, x$cfu_mean, col = cols_dark[1], lwd=2)
polygon(c(x$month, rev(x$month), x$month[1]), c(x$cfu_cilo, rev(x$cfu_cihi), x$cfu_cilo[1]), border = NA, col = cols[1])

svg("../outputs/SA_CFU_dynamics.svg", width = 866, height = 553)
dev.off()


# Median and IQR of staphylococcus reads abundance
# 
# png("../outputs/SA_abund_in_time_plot.png", width = 866, height = 553)
# plot(0, type = "n", xlim = range(d$days), xlab = "Days after decolonization", ylim = range(d$SA_abund), ylab ="Percentage of S. aureus reads")
# 
# for(patient_current in levels(d$patient)) {
#   dsub <- d[patient == patient_current]
#   lines(dsub$days, dsub$SA_abund, col = cols[dsub$carrier[1] + 1], lty = 2)
# }
# 
# legend("topleft", bty = "n", title = "", fill = cols, legend = c("Non-carrier", "Carrier"))
# 
# SAmed <- d[, .(SAabund_mean = mean(SA_abund), SA_cilo = as.numeric(mean(SA_abund) - 1.96*sd(SA_abund)/sqrt(8) ), SA_cihi = as.numeric(mean(SA_abund) + 1.96*sd(SA_abund)/sqrt(8))), by = .(days, carrier)  ]
# SAmed[SA_cilo < 0, SA_cilo := 0]
# 
# x <- SAmed[carrier == TRUE]
# lines(x$days, x$SAabund_mean, col = cols_dark[2], lwd=2)
# polygon(c(x$days, rev(x$days), x$days[1]), c(x$SA_cilo, rev(x$SA_cihi), x$SA_cilo[1]), border = NA, col = cols[2])
# 
# SAmed <- d[, .(SAabund_mean = mean(SA_abund), SA_cilo = as.numeric(mean(SA_abund) - 1.96*sd(SA_abund)/sqrt(9) ), SA_cihi = as.numeric(mean(SA_abund) + 1.96*sd(SA_abund)/sqrt(9))), by = .(days, carrier)  ]
# SAmed[SA_cilo < 0, SA_cilo := 0]
# 
# x <- SAmed[carrier == FALSE]
# lines(x$days, x$SAabund_mean, col = cols_dark[1], lwd=2)
# polygon(c(x$days, rev(x$days), x$days[1]), c(x$SA_cilo, rev(x$SA_cihi), x$SA_cilo[1]), border = NA, col = cols[1])
# 
# dev.off()
# 

# Median and IQR of staphylococcus reads abundance + log odds ratio on confidence intervals

df <- select(d, c(carrier, month, patient, SA_prop))
min_prop <- min(df[,SA_prop][which(df[,SA_prop]>0)])

# NOTE : Remove "0.00" values from proportions to not bias log means with -Inf values. Take the minimal value found and divide by e^1
df$SA_prop[df$SA_prop==0.00] <- min_prop/exp(1)
df$SAprop_logodds <- df[, log(SA_prop / (1 - SA_prop))]


#df$SAprop_logodds_sd <- df[,sd(log(SA_prop / (1 - SA_prop)))] # NAN
#df$SAprop_logodds_mean_exp <- exp(df2$SAprop_logodds)/ (1+exp(mean(df2$SAprop_logodds)))

# df2 <- df[carrier == TRUE, .(SAprop_logodds_mean = mean(SAprop_logodds),  
#                              SAprop_logodds_lower = mean(SAprop_logodds) - 1.96*sd(SAprop_logodds)/sqrt(8), 
#                              SAprop_logodds_upper = mean(SAprop_logodds) + 1.96*sd(SAprop_logodds)/sqrt(8)), by = .(days, carrier) ]
# 
# df2 <- rbind(df2, df[carrier == FALSE, .(SAprop_logodds_mean = mean(SAprop_logodds),  
#                              SAprop_logodds_lower = mean(SAprop_logodds) - 1.96*sd(SAprop_logodds)/sqrt(9), 
#                              SAprop_logodds_upper = mean(SAprop_logodds) + 1.96*sd(SAprop_logodds)/sqrt(9)), by = .(days, carrier) ])

df2 <- df[, .(SAprop_logodds_mean = mean(SAprop_logodds),  
                             SAprop_logodds_lower = mean(SAprop_logodds) - 1.96*sd(SAprop_logodds)/sqrt(8), 
                             SAprop_logodds_upper = mean(SAprop_logodds) + 1.96*sd(SAprop_logodds)/sqrt(8)), by = .(month, carrier) ]


df2$SAprop_mean <- exp(df2$SAprop_logodds_mean)/ (1 + exp(df2$SAprop_logodds_mean))
df2$SAprop_lower <- exp(df2$SAprop_logodds_lower)/ (1 + exp(df2$SAprop_logodds_lower))
df2$SAprop_upper <-exp(df2$SAprop_logodds_upper)/ (1 + exp(df2$SAprop_logodds_upper))

  

par(mar=c(5,6,4,1)+.1)
plot(0, type = "n", xlim = range(df2$month), xlab = "Months after decolonization", ylim = c(0, 1.1), ylab =expression("" *italic(S.~aureus)* " reads proportion"), cex.lab=1.3, cex.axis=1.3)
mtext("B", side = 3, adj = 0, line = 2, cex = 2, font= 2)
legend("topright", bty = "n", title = "", fill = cols_med, legend = c("Non-carrier", "Carrier"), cex = 1.3)

for(patient_current in levels(d$patient)) {
  dsub <- d[patient == patient_current]
  lines(dsub$month, dsub$SA_prop, col = cols_med[dsub$carrier[1] + 1], lty = 2)
}

x <- df2[carrier == TRUE]
lines(x$month, x$SAprop_mean, col = cols_dark[2], lwd=2)
polygon(c(x$month, rev(x$month), x$month[1]), c(x$SAprop_lower, rev(x$SAprop_upper), x$SAprop_lower[1]), border = NA, col = cols[2])

x <- df2[carrier == FALSE]
lines(x$month, x$SAprop_mean, col = cols_dark[1], lwd=2)
polygon(c(x$month, rev(x$month), x$month[1]), c(x$SAprop_lower, rev(x$SAprop_upper), x$SAprop_lower[1]), border = NA, col = cols[1])

png("../outputs/SA_abund_dynamics", width = 866, height = 553)
dev.off()


# EXAMPLE OF LOGODDS 95%CI
# 
# n <- 10
# x <- runif(n)*0.1
# 
# boxplot(x)
# 
# # Computing 95%CI on logodds
# x_logodds <- log(x / (1 - x))
# 
# hist(x_logodds)
# 
# # Confidence intervals
# xlo_mean <- mean(x_logodds)
# xlo_sd <- sd(x_logodds)
# 
# xlo_lower <- xlo_mean - 1.96*xlo_sd/sqrt(n)
# xlo_upper <- xlo_mean + 1.96*xlo_sd/sqrt(n)
# 
# x_mean <- exp(xlo_mean)  /(1 + exp(xlo_mean))
# x_lower <- exp(xlo_lower)  /(1 + exp(xlo_lower))
# x_upper <- exp(xlo_upper)  /(1 + exp(xlo_upper))
# 


