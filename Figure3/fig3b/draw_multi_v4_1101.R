# !/usr/bin/env Rscript
#
##############################
# Name:
# Author: Boyang JI
# Todo:
# Version:
# Reference
#
##############################

# -----------------------------
# load related library
# -----------------------------

# data import, tidy
library(tidyverse)
library(stringr)

# plotting
library(ggsci)
library(ggpubr)

# fonts
library(extrafont)
fonts()

library(RColorBrewer)
display.brewer.all()
library("extrafont")
loadfonts()
# -----------------------------
# load data
# -----------------------------
color_4M<-brewer.pal(n=8,"Dark2")[1:5]
color_12M<-c(brewer.pal(n=8,"Dark2")[1:5],c("#BC80BD", "#E6AB02", "#80B1D3"))
exp_file <- "experimentdata_20181031.txt"
prd_file<-"prediction_20181031.txt"

exp_data <- read.table(file = exp_file, sep = "\t", header = T, stringsAsFactors = F)
exp_data <- exp_data[, 1:3]
prd_data <- read.table(file = prd_file, sep = "\t", header = T, stringsAsFactors = F)
# prd_data <- read.table(file = prd_file, sep = ",", header = T, stringsAsFactors = F)


all_data <- merge(exp_data, prd_data, by = c("Type", "Species"))
all_data$Species <- factor(all_data$Species, labels = c(
  "Bfr", "Bth", "Bad", "Bbv",
  "Blg", "Ehal", "Fpr", "Rint"
))
all_data$Species <- factor(all_data$Species, levels = c(
  "Bth", "Bfr", "Blg", "Bbv",
  "Bad", "Ehal", "Fpr", "Rint"
))

## Add an alpha value to a colour
add.alpha <- function(col, alpha=1) {
  if (missing(col)) {
    stop("Please provide a vector of colours.")
  }
  apply(
    sapply(col, col2rgb) / 255, 2,
    function(x)
      rgb(x[1], x[2], x[3], alpha = alpha)
  )
}

col_8 <- brewer.pal(8, "Set1")
col_8 <- add.alpha(col_8, alpha = 0.9)
col_8<-c("#e95e66","#669eca","#78c47a","#b07ab6","#f6a168","#f6f278","#bc8162","#f2a4c7")
col_8<-c("#b73276", "#228fa0", "#325780", "#3a8f62", "#b48939","#8f5360", "#754577", "#a7624f")
# draw 4M data
M4_data <- all_data %>% filter(Type == "4 Months") %>% filter(prediction != 0)
# M4_others <- data.frame(Type="4 Months", Species="Others", experiment=1-sum(M4_data[,3]), prediction=1-sum(M4_data[,4]))
# M4_data_new <- rbind(M4_data, M4_others)
M4_data_new <- with(M4_data, M4_data[order(Species), ])

M4_sum_experiment <- round(sum(M4_data$experiment) * 100, 2)
M4_sum_experiment_label <- paste(M4_sum_experiment, "%", sep = " ")
M4_sum_prediction <- round(sum(M4_data$prediction) * 100, 2)
M4_sum_prediction_label <- paste(M4_sum_prediction, "%", sep = " ")

opar <- par(no.readonly = T)
pdf("4M_v2_1101.pdf", width = 11 / 2.54, height = 5 / 2.54, pointsize = 7)
par(mfrow = c(1, 2), mar = c(0, 0, 0, 0))

pie(
  M4_data_new[, 3],
  labels = M4_data_new$Species,
  col = c(col_8[1:5]),
  border="white",
  clockwise = T,
  init.angle = 0,cex=1.5
)
par(new = TRUE)
pie(1, radius = 0.4, col = "white", border = "white", lwd = 2, labels = "")
# text(0, 0, labels = M4_sum_experiment_label, cex = 1, font = 2)
text(0, 0, labels = "BF_exp", cex = 1.5)

pie(
  M4_data_new[, 4],
  labels = M4_data_new$Species,
  col = c(col_8[1:5]),
  border="white",
  clockwise = T,
  init.angle = 0,cex=1.5
)
par(new = TRUE)
pie(1, radius = 0.4, col = "white", border = "white", lwd = 2, labels = "")
# text(0, 0, labels = M4_sum_prediction_label, cex = 1, font = 2)
# text(0, 0, labels = "BF_pre", cex = 1, font =list(family="Arial", face=1))
text(0, 0, labels = "BF_pre", cex = 1.5)


dev.off()

# draw 4M data
M12_data <- all_data %>% filter(Type == "12 Months") %>% filter(prediction != 0)
# M12_others <- data.frame(Type="12 Months", Species="Others", experiment=1-sum(M12_data[,3]), prediction=1-sum(M12_data[,4]))
# M12_data_new <- rbind(M12_data, M12_others)
M12_data_new <- with(M12_data, M12_data[order(Species), ])

M12_sum_experiment <- round(sum(M12_data$experiment) * 100, 2)
M12_sum_experiment_label <- paste(M12_sum_experiment, "%", sep = " ")
M12_sum_prediction <- round(sum(M12_data$prediction) * 100, 2)
M12_sum_prediction_label <- paste(M12_sum_prediction, "%", sep = " ")


opar <- par(no.readonly = T)
pdf("12M_v2_1101.pdf", width = 11 / 2.54, height = 5 / 2.54, pointsize = 7)
par(mfrow = c(1, 2), mar = c(0, 0, 0, 0))

pie(
  M12_data_new[, 3],
  labels = M12_data_new$Species,
  col = c(col_8, "gray80"),
  border="white",
  clockwise = T,cex=1.5
)
par(new = TRUE)
pie(1, radius = 0.4, col = "white", border = "white", lwd = 2, labels = "")
# text(0, 0, labels = M12_sum_experiment_label, cex = 1, font = 2)
text(0, 0, labels = "SF_exp", cex = 1.5, font = 1)

pie(
  M12_data_new[, 4],
  labels = M12_data_new$Species,
  col = c(col_8, "gray80"),
  border="white",
  clockwise = T,cex=1.5
)
par(new = TRUE)
pie(1, radius = 0.4, col = "white", border = "white", lwd = 2, labels = "")
# text(0, 0, labels = M12_sum_experiment_label, cex = 1, font = 2)
text(0, 0, labels = "SF_pre", cex = 1.5, font =1)


dev.off()

