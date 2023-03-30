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

# Adjusting text sizes for plotting function in package conceptorCP  
# for publishable output instead of R output
library(scales)
plotCP <- function(ccp_output, nbreaks = 10) {
  if(check_integer(nbreaks) == FALSE) {
    stop("Enter an integer number of break points for plot.")
  }
  
  if(nbreaks < 2 | nbreaks > 40) {
    warning("Plot may be improved by adjusting number of break points.")
  }
  
  Time <- Window <- Values <- Reference <- WindowLength <- RCDF <- NULL
  Angles <- dplyr::tibble(Angles = ccp_output$angles)
  L <- nrow(Angles)
  Angles <- Angles %>% dplyr::mutate(Time = seq(1, L), timeMin = Time - 0.5, timeMax = Time + 0.5)
  PWAngles <- Angles[(ccp_output$netParams$washoutL + ccp_output$netParams$trainL + 1):L,]
  AngleMin <- min(PWAngles) - (1 - min(PWAngles)) / 100
  PWAngles <- PWAngles %>% dplyr::mutate(PWRanks = rank(PWAngles$Angles) / nrow(PWAngles))
  EndPts <- floor(seq(0, nrow(PWAngles), length.out = nbreaks + 1)) + ccp_output$netParams$washoutL + ccp_output$netParams$trainL
  PWAngles <- PWAngles %>% dplyr::rowwise() %>% dplyr::mutate(Window = sum(EndPts < Time))
  PWAngles$Values <- paste("Value", unlist(sapply(diff(EndPts), seq, simplify = F)))
  
  plotM <- ggplot2::ggplot(PWAngles, ggplot2::aes_string(x = "Time", y = "Angles")) +
    ggplot2::geom_rect(data = PWAngles, ggplot2::aes_string(xmin = "timeMin", xmax = "timeMax", ymin = "AngleMin", ymax = 1, fill = "PWRanks")) +
    ggplot2::scale_fill_gradient2(name = "",
                                  low = "red", mid = "white", high = "blue", midpoint = 0.5,
                                  limits = c(0, 1), na.value = "white",
                                  breaks = c(0, 0.5, 1),
                                  labels = c("Away from \nConceptor \nSpace", "\nMiddle: Percentiles of Cosine Similarities \n\nBottom: Relative ECDF Difference", "Towards \nConceptor \nSpace")) +
    ggplot2::geom_point(ggplot2::aes_string(x = "Time", y = "Angles"), shape = 20, size = 1) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(name = "Time",
                                expand = c(0, 0),
                                breaks = EndPts + c(rep(1, nbreaks), 0)) +
    ggplot2::scale_y_continuous(name = paste("Cosine Similaritiy\n Smaller  \u2194  Larger", sep = ''),
                                labels = scales::label_number(accuracy = 0.1),
                                limits = c(AngleMin, 1),
                                breaks = 1,
                                expand = c(0, 0)) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(color = "white", size = 2),
                   axis.ticks.y = ggplot2::element_blank(),
                   legend.position = "bottom",
                   legend.text = ggplot2::element_text(size = 16),
                   axis.title.x = ggplot2::element_text(size = 16),
                   axis.text.x = ggplot2::element_text(size = 12),
                   axis.title.y = ggplot2::element_text(size = 16),
                   legend.justification = "top",
                   plot.margin = ggplot2::unit(c(0.2, 0.4, 0.2, 0.2), "cm"),
                   legend.key.height = ggplot2::unit(0.3, "cm"),
                   legend.key.width = ggplot2::unit(3, "cm"),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank())
  
  CDF <- dplyr::tibble(Reference = rep(sort(PWAngles$Angles), nbreaks), RCDF = rep(seq(nrow(PWAngles)) / nrow(PWAngles), nbreaks), RCDFmin = RCDF - 1 / nrow(PWAngles), Window = rep(seq(1, nbreaks), each = nrow(PWAngles)), WindowLength = rep(diff(EndPts), each = nrow(PWAngles)))
  WCDF <- dplyr::select(PWAngles, Angles, Window, Values) %>% tidyr::pivot_wider(names_from = Values, values_from = Angles)
  WCDF <- dplyr::left_join(CDF, WCDF, by = "Window")
  WCDF <- dplyr::transmute(WCDF, dplyr::across(6:(max(diff(EndPts) + 5)), function(X) X <= Reference))
  CDF$WCDF <- dplyr::rowwise(WCDF) %>% dplyr::transmute(WCDF = sum(dplyr::c_across(cols = dplyr::everything()), na.rm = T))
  CDF <- dplyr::mutate(CDF, WCDF = WCDF$WCDF / WindowLength, Shading = RCDF - WCDF)
  
  plotB <- ggplot2::ggplot(data = CDF) +
    ggplot2::geom_rect(ggplot2::aes_string(xmin = "RCDFmin", xmax = "RCDF", ymin = 0, ymax = 1, fill = "Shading")) +
    ggplot2::scale_fill_gradient2(name = "",
                                  low = "red",  mid = "white", high = "blue",
                                  limits = c(-1, 1), na.value = "white", midpoint = 0) +
    ggplot2::facet_wrap(~Window, nrow = 1) +
    ggplot2::geom_point(data = CDF, ggplot2::aes_string(x = "RCDF", y = "WCDF"), col = "black", shape = 20, size = 0.2) +
    ggplot2::scale_x_continuous(name = "ECDF (Full Time Series)", limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(name = "ECDF\n(Time Window)", limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::theme(strip.text.x = ggplot2::element_blank(),
                   strip.text.y = ggplot2::element_blank(),
                   strip.background = ggplot2::element_rect(color = "black", fill = "white"),
                   panel.spacing.y = ggplot2::unit(0.05, "cm"),
                   panel.spacing.x = ggplot2::unit(0, "cm"),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.key.width = ggplot2::unit(0.3, "cm"),
                   legend.key.height = ggplot2::unit(0.5, "cm"),
                   axis.text.y = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_text(size = 16),
                   legend.justification = "bottom",
                   axis.text.x = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_text(size = 16),
                   legend.position = "none",
                   plot.margin = ggplot2::unit(c(0.2, 0.4, 0.2, 0.2), "cm"),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank())
  
  SSeries <- dplyr::tibble(Time = seq(ccp_output$netParams$washoutL + ccp_output$netParams$trainL + 1, L), Stat = ccp_output$statSeries)
  upper.limit <- max(SSeries$Stat, stats::quantile(ccp_output$MBBnull, 0.99)[[1]]) + 0.02
  
  plotT <- ggplot2::ggplot(SSeries) +
    ggplot2::geom_line(ggplot2::aes_string(x = "Time", y = "Stat")) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = stats::quantile(ccp_output$MBBnull, 0.90)[[1]], linetype = "Upper 10% \nMBB Null Dist."), color = "red", key_glyph = "cust") +
    ggplot2::geom_hline(ggplot2::aes(yintercept = stats::quantile(ccp_output$MBBnull, 0.95)[[1]], linetype = "Upper 5% \nMBB Null Dist."), color = "red", key_glyph = "cust") +
    ggplot2::geom_hline(ggplot2::aes(yintercept = stats::quantile(ccp_output$MBBnull, 0.99)[[1]], linetype = "Upper 1% \nMBB Null Dist."), color = "red", key_glyph = "cust") +
    ggplot2::geom_vline(ggplot2::aes(xintercept = ccp_output$estimate, linetype = "Most Likely \nChange Point"), color = "blue", key_glyph = "cust") +
    ggplot2::scale_linetype_manual("", values = c("Most Likely \nChange Point" = "solid", "Upper 10% \nMBB Null Dist." = "dotted", "Upper 5% \nMBB Null Dist." = "dashed", "Upper 1% \nMBB Null Dist." = "longdash"),
                                   guide = ggplot2::guide_legend(override.aes = list(colour = c("blue", "red", "red", "red")))) +
    ggplot2::scale_x_continuous(name = "Time", limits = c(ccp_output$netParams$washoutL + ccp_output$netParams$trainL + 1, L),
                                expand = c(0, 0), EndPts + c(rep(1, nbreaks), 0)) +
    ggplot2::scale_y_continuous(name = "Statistic", limits = c(0, upper.limit), labels = label_number(accuracy = 0.1),
                                expand = c(0, 0), breaks = seq(0, upper.limit, 0.1)) + ggplot2::theme_bw() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 12),
                   legend.position = "top",
                   legend.justification = "right",
                   legend.text = ggplot2::element_text(size = 16),
                   axis.title.x = ggplot2::element_text(size = 16),
                   axis.text.x = ggplot2::element_text(size = 12),
                   axis.title.y = ggplot2::element_text(size = 16),
                   plot.margin = ggplot2::unit(c(0.2, 0.4, 0.2, 0.2), "cm"),
                   plot.caption = ggplot2::element_text(size = 16),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank())
  
  plotS <- cowplot::axis_canvas(plotT, axis = "y", coord_flip = TRUE) +
    ggplot2::geom_histogram(ggplot2::aes(x = ccp_output$MBBnull), bins = 24, color = "black", fill = rgb(1, 0, 0, 0.2)) +
    ggplot2::coord_flip()
  suppressWarnings({
    plot1 <- cowplot::insert_yaxis_grob(plotT, plotS, grid::unit(0.1, "null"), position = "right")
    plot2 <- cowplot::plot_grid(plotM + ggplot2::theme(legend.position = "none"), plotB, ncol = 1, align = "hv", axis = "lr", rel_heights = c(1.1, 1))
    plot3 <- cowplot::insert_yaxis_grob(plot2, grid::nullGrob(), grid::unit(0.1, "null"), position = "right")
    plot4 <- cowplot::plot_grid(plot1, plot3, ncol = 1)
    plot5 <- cowplot::plot_grid(plot4, cowplot::get_legend(plotM), rel_heights = c(1, 0.15), ncol = 1, align = "hv", axis = "r")})
  return(plot5)
}

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

