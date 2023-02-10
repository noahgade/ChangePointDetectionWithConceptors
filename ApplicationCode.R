# install.packages("conceptorCP_1.0.0.tar.gz", repos = NULL)
rm(list = ls())
library(conceptorCP)
library(ecp)
library(hdbinseg)
library(gridExtra)
library(tidyverse)

set.seed(102722)

# Delta band LFP data, C Varela & MA Wilson (2019, 2020)
LFP <- dplyr::as_tibble(read.csv("LFP_Section1_DeltaBand.csv", header = F))
colnames(LFP) <- c("Time", "HC", "PFC", "THAL")
Sleep <- dplyr::as_tibble(read.csv("SleepTimes_Section1.csv", header = F))
colnames(Sleep) <- c("SleepTimeStart", "SleepTimeEnd")

change_point1 <- as.numeric(Sleep[1,2])
change_point2 <- as.numeric(Sleep[2,1])
change_point3 <- as.numeric(Sleep[3,1])

TimeGrid <- matrix(c(change_point1 - 90, change_point2 - 60, change_point3 - 70, 
                     change_point1 + 10, change_point2 + 40, change_point3 + 30), ncol = 2)

example_data1 <- LFP[(LFP$Time >= TimeGrid[1, 1]) & (LFP$Time < TimeGrid[1, 2]),]
example_data2 <- LFP[(LFP$Time >= TimeGrid[2, 1]) & (LFP$Time < TimeGrid[2, 2]),]
example_data3 <- LFP[(LFP$Time >= TimeGrid[3, 1]) & (LFP$Time < TimeGrid[3, 2]),]

washoutL <- 40  # 10 seconds for washout (data frequency 4Hz)
trainL <- 120   # 30 seconds for training (data frequency 4Hz)

CCP_Example1 <- ccp(example_data1[,2:4], washoutL = washoutL, trainL = trainL, tol = 0.04, plot.it = F)
CCP_Example2 <- ccp(example_data2[,2:4], washoutL = washoutL, trainL = trainL, tol = 0.04, plot.it = F)
CCP_Example3 <- ccp(example_data3[,2:4], washoutL = washoutL, trainL = trainL, tol = 0.04, plot.it = F)

c(CCP_Example1$estimate, CCP_Example1$statistic, CCP_Example1$MBBquant, example_data1$Time[CCP_Example1$estimate])
c(CCP_Example2$estimate, CCP_Example2$statistic, CCP_Example2$MBBquant, example_data2$Time[CCP_Example2$estimate])
c(CCP_Example3$estimate, CCP_Example3$statistic, CCP_Example3$MBBquant, example_data3$Time[CCP_Example3$estimate])

ED_Example1 <- e.divisive(example_data1[(washoutL + trainL + 1):nrow(example_data1), 2:4])
ED_Example2 <- e.divisive(example_data2[(washoutL + trainL + 1):nrow(example_data2), 2:4])
ED_Example3 <- e.divisive(example_data3[(washoutL + trainL + 1):nrow(example_data3), 2:4])

SBS2_Example1 <- sbs.alg(t(example_data1[(washoutL + trainL + 1):nrow(example_data1), 2:4]), cp.type = 2)
SBS2_Example2 <- sbs.alg(t(example_data2[(washoutL + trainL + 1):nrow(example_data2), 2:4]), cp.type = 2)
SBS2_Example3 <- sbs.alg(t(example_data3[(washoutL + trainL + 1):nrow(example_data3), 2:4]), cp.type = 2)

KCP_Example1 <- kcpa(as.matrix(example_data1[(washoutL + trainL + 1):nrow(example_data1), 2:4]), L = 1, C = 2)
KCP_Example2 <- kcpa(as.matrix(example_data2[(washoutL + trainL + 1):nrow(example_data2), 2:4]), L = 1, C = 2)
KCP_Example3 <- kcpa(as.matrix(example_data3[(washoutL + trainL + 1):nrow(example_data3), 2:4]), L = 1, C = 2)

Example1PlotData <- dplyr::tibble(Time = rep(example_data1$Time, 3), 
                                  Values = c(example_data1$HC, example_data1$PFC, example_data1$THAL),
                                  Region = c(rep("HC", nrow(example_data1)), rep("PFC", nrow(example_data1)), rep("THAL", nrow(example_data1))))
Example2PlotData <- dplyr::tibble(Time = rep(example_data2$Time, 3), 
                                  Values = c(example_data2$HC, example_data2$PFC, example_data2$THAL),
                                  Region = c(rep("HC", nrow(example_data2)), rep("PFC", nrow(example_data2)), rep("THAL", nrow(example_data2))))
Example3PlotData <- dplyr::tibble(Time = rep(example_data3$Time, 3), 
                                  Values = c(example_data3$HC, example_data3$PFC, example_data3$THAL),
                                  Region = c(rep("HC", nrow(example_data3)), rep("PFC", nrow(example_data3)), rep("THAL", nrow(example_data3))))

Example1Plot <- ggplot(data = Example1PlotData, aes(x = Time, y = Values)) + 
  facet_wrap(~Region, ncol = 1, strip.position = "left", scales = "free_y") +
  theme_bw() + 
  ylab(NULL) +
  geom_line() + 
  scale_x_continuous("", expand = c(0, 0), breaks = seq(0, 2000, 20)) + 
  scale_color_manual(name = "", values = c("CCP" = "blue", "EDiv" = "firebrick2", "SBS2" = "darkgoldenrod1", "KCP" = "green3", "Actual" = "gray50")) + 
  scale_linetype_manual(name = "", values = c("CCP" = "longdash", "EDiv" = "dashed", "SBS2" = "twodash", "KCP" = "dotdash", "Actual" = "solid")) +
  geom_vline(aes(xintercept = change_point1, color = "Actual", linetype = "Actual"), size = 1) + 
  geom_vline(aes(xintercept = example_data1$Time[(washoutL + trainL + 1):nrow(example_data1)][ifelse(!is.na(ED_Example1$order.found[3]), ED_Example1$order.found[3], nrow(example_data1) - washoutL - trainL - 1)], color = "EDiv", linetype = "EDiv"), size = 1) + 
  geom_vline(aes(xintercept = example_data1$Time[(washoutL + trainL + 1):nrow(example_data1)][ifelse(SBS2_Example1$tree[[1]][3,1] > 0, SBS2_Example1$tree[[1]][3,1], nrow(example_data1) - washoutL - trainL - 1)] ,color = "SBS2", linetype = "SBS2"), size = 1) + 
  geom_vline(aes(xintercept = example_data1$Time[(washoutL + trainL + 1):nrow(example_data1)][ifelse(KCP_Example1[2] <= nrow(example_data1) - washoutL - trainL, KCP_Example1[2], nrow(example_data1) - washoutL - trainL)], color = "KCP", linetype = "KCP"), size = 1) +
  geom_vline(aes(xintercept = example_data1$Time[ifelse(CCP_Example1$MBBquant <= 0.05, CCP_Example1$estimate, nrow(example_data1))], color = "CCP", linetype = "CCP"), size = 1) + 
  guides(linetype = guide_legend(override.aes = list(size = 1))) +
  theme(strip.placement = "outside", strip.background = element_blank(), 
        strip.text = element_text(size = 14), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), legend.text = element_text(size = 16), 
        legend.position = "none", axis.text.x = element_text(size = 12), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.title.x = element_text(size = 16), panel.spacing = unit(0, "mm"))

Example2Plot <- ggplot(data = Example2PlotData, aes(x = Time, y = Values)) + 
  facet_wrap(~Region, ncol = 1, strip.position = "left", scales = "free_y") +
  theme_bw() + 
  ylab(NULL) +
  geom_line() + 
  scale_x_continuous("", expand = c(0, 0), breaks = seq(0, 2000, 20)) + 
  scale_color_manual(name = "", values = c("CCP" = "blue", "EDiv" = "firebrick2", "SBS2" = "darkgoldenrod1", "KCP" = "green3", "Actual" = "gray50")) + 
  scale_linetype_manual(name = "", values = c("CCP" = "longdash", "EDiv" = "dashed", "SBS2" = "twodash", "KCP" = "dotdash", "Actual" = "solid")) +
  geom_vline(aes(xintercept = change_point2, color = "Actual", linetype = "Actual"), size = 1) + 
  geom_vline(aes(xintercept = example_data2$Time[(washoutL + trainL + 1):nrow(example_data2)][ifelse(!is.na(ED_Example2$order.found[3]), ED_Example2$order.found[3], nrow(example_data2) - washoutL - trainL - 1)], color = "EDiv", linetype = "EDiv"), size = 1) + 
  geom_vline(aes(xintercept = example_data2$Time[(washoutL + trainL + 1):nrow(example_data2)][ifelse(SBS2_Example2$tree[[1]][3,1] > 0, SBS2_Example2$tree[[1]][3,1], nrow(example_data2) - washoutL - trainL - 1)] ,color = "SBS2", linetype = "SBS2"), size = 1) + 
  geom_vline(aes(xintercept = example_data2$Time[(washoutL + trainL + 1):nrow(example_data2)][ifelse(KCP_Example2[2] <= nrow(example_data2) - washoutL - trainL , KCP_Example2[2], nrow(example_data2) - washoutL - trainL - 1)], color = "KCP", linetype = "KCP"), size = 1) +
  geom_vline(aes(xintercept = example_data2$Time[ifelse(CCP_Example2$MBBquant <= 0.05, CCP_Example2$estimate, nrow(example_data2))], color = "CCP", linetype = "CCP"), size = 1) + 
  guides(linetype = guide_legend(override.aes = list(size = 1))) +
  theme(strip.placement = "outside", strip.background = element_blank(), 
        strip.text = element_text(size = 14), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), legend.text = element_text(size = 16), 
        legend.position = "none", axis.text.x = element_text(size = 12), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.title.x = element_text(size = 16), panel.spacing = unit(0, "mm"))

Example3Plot <- ggplot(data = Example3PlotData, aes(x = Time, y = Values)) + 
  facet_wrap(~Region, ncol = 1, strip.position = "left", scales = "free_y") +
  theme_bw() + 
  ylab(NULL) +
  geom_line() + 
  scale_x_continuous("Experiment Time (s)", expand = c(0, 0), breaks = seq(0, 2000, 20)) + 
  scale_color_manual(name = "", values = c("CCP" = "blue", "EDiv" = "firebrick2", "SBS2" = "darkgoldenrod1", "KCP" = "green3", "Actual" = "gray50")) + 
  scale_linetype_manual(name = "", values = c("CCP" = "longdash", "EDiv" = "dashed", "SBS2" = "twodash", "KCP" = "dotdash", "Actual" = "solid")) +
  geom_vline(aes(xintercept = change_point3, color = "Actual", linetype = "Actual"), size = 1) + 
  geom_vline(aes(xintercept = example_data3$Time[(washoutL + trainL + 1):nrow(example_data3)][ifelse(!is.na(ED_Example3$order.found[3]), ED_Example3$order.found[3], nrow(example_data3) - washoutL - trainL - 1)], color = "EDiv", linetype = "EDiv"), size = 1) + 
  geom_vline(aes(xintercept = example_data3$Time[(washoutL + trainL + 1):nrow(example_data3)][ifelse(SBS2_Example3$tree[[1]][3,1] > 0, SBS2_Example3$tree[[1]][3,1], nrow(example_data3) - washoutL - trainL - 1)] ,color = "SBS2", linetype = "SBS2"), size = 1) + 
  geom_vline(aes(xintercept = example_data3$Time[(washoutL + trainL + 1):nrow(example_data3)][ifelse(KCP_Example3[2] <=  nrow(example_data3) - washoutL - trainL, KCP_Example3[2], nrow(example_data3) - washoutL - trainL - 1)], color = "KCP", linetype = "KCP"), size = 1) +
  geom_vline(aes(xintercept = example_data3$Time[ifelse(CCP_Example3$MBBquant <= 0.05, CCP_Example3$estimate, nrow(example_data3))], color = "CCP", linetype = "CCP"), size = 1) + 
  guides(linetype = guide_legend(override.aes = list(size = 1))) +
  theme(strip.placement = "outside", strip.background = element_blank(), 
        strip.text = element_text(size = 14), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), legend.text = element_text(size = 16), 
        legend.position = "bottom", axis.text.x = element_text(size = 12), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.title.x = element_text(size = 16), panel.spacing = unit(0, "mm"))

ExamplePlot <- arrangeGrob(Example1Plot, Example2Plot, Example3Plot, widths = unit(8.3, "in"), heights = unit(c(2.3, 2.3, 2.9), "in"), ncol = 1)
plot(ExamplePlot)

ggsave(plot = ExamplePlot, file = "ExamplePlot.pdf", width = 8.3, height = 7.5)

CCPPlot1 <- plotCP(CCP_Example1)
CCPPlot2 <- plotCP(CCP_Example2)
CCPPlot3 <- plotCP(CCP_Example3)

quartz(file = "Example1CCP.pdf", type = "pdf", width = 13, height = 9)
CCPPlot1
dev.off()

quartz(file = "Example2CCP.pdf", type = "pdf", width = 13, height = 9)
CCPPlot2
dev.off()

quartz(file = "Example3CCP.pdf", type = "pdf", width = 13, height = 9)
CCPPlot3
dev.off()

