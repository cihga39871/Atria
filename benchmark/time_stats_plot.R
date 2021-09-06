#!/usr/bin/env Rscript

if (is.na(packageDescription("plotly")[1])) install.packages("plotly")
library(plotly, quietly = T, warn.conflicts = F)
if (is.na(packageDescription("argparse")[1])) install.packages("argparse")
library(argparse, quietly = T)

parser <- ArgumentParser(description='Plot time stats (speed vs threads/CPU)')
parser$add_argument('-i', '--input', dest='input', metavar='FILE', type='character',
                    required=TRUE, nargs='+',
                    help='[REQUIRED] input time stats tables generated from time_stats.jl (1st=ungz, 2nd=gz')
parser$add_argument('-o', '--output', dest='out', metavar='PLOT', type='character',
                    default="time_stats_plot.html",
                    help='output html heatmap file (default: time_stats_plot.html)')

args <- parser$parse_args()

if (FALSE){
    setwd("~/analysis/atria-benchmark/simulate")
    args <- parser$parse_args(c("-i", "stats.time_benchmark3.df.txt", "stats.time_benchmark_gz3.df.txt", "-o", "time_stats_plot2.html"))
}

wrapper <- function(input_path, show_legend, is_gz){
    
    input <- read.delim(input_path, header=TRUE, as.is=TRUE)
    
    input$Trimmer <- ""
    input$Trimmer[grepl("atria", input$Command)] <- "Atria (consensus)"
    input$Trimmer[grepl("atria --no-consensus", input$Command)] <- "Atria"
    input$Trimmer[grepl("AdapterRemoval", input$Command)] <- "AdapterRemoval"
    input$Trimmer[grepl("skewer", input$Command)] <- "Skewer"
    input$Trimmer[grepl("trim_galore", input$Command)] <- "Trim Galore"
    input$Trimmer[grepl("trimmomatic", input$Command)] <- "Trimmomatic"
    input$Trimmer[grepl("ktrim", input$Command)] <- "Ktrim"
    input$Trimmer[grepl("atropos", input$Command)] <- "Atropos"
    input$Trimmer[grepl("fastp", input$Command)] <- "Fastp"
    input$Trimmer[grepl("SeqPurge", input$Command)] <- "SeqPurge"
    
    input_labels <- c("Atria (consensus)",
                      "Atria",
                      "AdapterRemoval",
                      "Skewer",
                      "Trim Galore",
                      "Trimmomatic",
                      "Ktrim",
                      "Atropos",
                      "Fastp",
                      "SeqPurge"
                      )
    
    input$Trimmer <- factor(input$Trimmer, input_labels)
    
    # input$Command <- NULL
    
    input_value = input
    for (j in 1:ncol(input)){
        for (i in 1:nrow(input)){
            input_value[i,j] <- sub(" ±.*", "", input[i,j])
        }
        if (!any(is.na(as.numeric(input_value[,j])))) {
            if (class(input_value[,j]) != "factor") {
                input_value[,j] <- as.numeric(input_value[,j])    
            }
        }
    }
    
    input_sd = input
    for (j in 1:ncol(input)){
        for (i in 1:nrow(input)){
            input_sd[i,j] <- sub(".*± ", "", input[i,j])
        }
        if (!any(is.na(as.numeric(input_sd[,j])))) {
            if (class(input_value[,j]) != "factor") {
                input_sd[,j] <- as.numeric(input_sd[,j])
                if (all(input_sd[,j] == input_value[,j])) {
                    input_sd[,j] <- 0
                }
            }
        }
    }
    
    input_high = input_value
    for (j in 1:ncol(input)){
        for (i in 1:nrow(input)){
            if (is.numeric(input_sd[i,j])){
                input_high[i,j] = input_value[i,j] + input_sd[i,j]
            }
        }
    }
    input_low = input_value
    for (j in 1:ncol(input)){
        for (i in 1:nrow(input)){
            if (is.numeric(input_sd[i,j])){
                input_low[i,j] = input_value[i,j] - input_sd[i,j]
            }
        }
    }
    
    
    if (sum(input_sd$Speed..M.Bases.s.) == 0){
        speed_error_y_array = NULL
        efficiency_error_y_array = NULL
    } else {
        speed_error_y_array <- input_sd$Speed..M.Bases.s.
        efficiency_error_y_array <- input_sd$Efficiency..M.Bases.s.CPU.
    }
    # writeLines(as.character(show_legend))
    if (is_gz){
        gz_title = " for Compressed Files"
    } else {
        gz_title = ""
    }
    x_tick_vals = unique(input_value$Threads)
    
    fig_speed <- plot_ly(x=input_value$Threads, y=input_value$Speed..M.Bases.s., color=input_value$Trimmer, legendgroup=input_value$Trimmer, error_y = list(array=speed_error_y_array), showlegend=show_legend) %>%
        add_lines(line = list(shape = "spline" )) %>%
        add_markers(showlegend = FALSE) %>%
        layout(
            xaxis = list(
                title = "Threads Assigned",
                tickvals = x_tick_vals
            ), yaxis = list(
                title = paste0("Speed", gz_title, "\n(M Bases / Second)")
            ))
    fig_speed
    
    fig_efficiency <- plot_ly(x=input_value$Threads, y=input_value$Efficiency..M.Bases.s.CPU., color=input_value$Trimmer, legendgroup=input_value$Trimmer, error_y = list(array=efficiency_error_y_array), showlegend=FALSE) %>%
        add_lines(line = list(shape = "spline")) %>%
        add_markers(showlegend = FALSE) %>%
        layout(
            xaxis = list(
                title = "Threads Assigned",
                tickvals = x_tick_vals
            ), yaxis = list(
                title = paste0("Efficiency", gz_title, "\n(M Bases / Second / CPU Usage)")
            ))
    fig_efficiency
    
    fig_speed_vs_realCPU <- plot_ly(x=input_value$CPU, y=input_value$Speed..M.Bases.s., color=input_value$Trimmer, legendgroup=input_value$Trimmer, error_y = list(array=efficiency_error_y_array), showlegend=FALSE) %>%
        add_trace(line = list(shape = "spline")) %>%
        add_markers(showlegend = FALSE) %>%
        layout(
            xaxis = list(
                title = "Real Average CPU Usage",
                tickvals = x_tick_vals
            ), yaxis = list(
                title = paste0("Speed", gz_title, "\n(M Bases / Second)")
            ))
    fig_speed_vs_realCPU
    
    return(list(
        input_value = input_value,
        input_sd = input_sd,
        input_high = input_high,
        input_low = input_low,
        fig_speed = fig_speed,
        fig_efficiency = fig_efficiency,
        fig_speed_vs_realCPU = fig_speed_vs_realCPU
    ))
}

stat_1 <- wrapper(args$input[1], T, F)
stat_2 <- wrapper(args$input[2], T, T)


p <- subplot(stat_1$fig_speed %>% layout(legend = list(orientation='h', 
                                                       y=1.3, 
                                                       bgcolor=rgb(0,0,0,0))), 
        stat_1$fig_speed_vs_realCPU,
        stat_2$fig_speed, 
        stat_2$fig_speed_vs_realCPU,
        nrows=2, shareX = T, shareY = T)

writeLines(sprintf("Output plot: %s", args$out))
htmlwidgets::saveWidget(as_widget(p), args$out)

plogx <- subplot(stat_1$fig_speed %>% layout(xaxis = list(type='log'),
                                             legend = list(orientation='h', 
                                                       y=1.3, 
                                                       bgcolor=rgb(0,0,0,0))), 
             stat_1$fig_speed_vs_realCPU %>% layout(xaxis = list(type='log')),
             stat_2$fig_speed %>% layout(xaxis = list(type='log')), 
             stat_2$fig_speed_vs_realCPU %>% layout(xaxis = list(type='log')),
             nrows=2, shareX = T, shareY = T)
plogx

outlogx = sub(".html$", ".logx.html", args$out)
writeLines(sprintf("Output logx plot: %s", outlogx))
htmlwidgets::saveWidget(as_widget(plogx), outlogx)
