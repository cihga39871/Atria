function julia_wrapper_rscript(rcode::String, args::Vector{String})
    try
        run(`Rscript --version`)
    catch e
        @error "Rscript not found in PATH. Please install R and export Rscript to PATH." _module=nothing _group=nothing _id=nothing _file=nothing
        return
    end
    tmpfile = tempname()::String
    tmpio = open(tmpfile, "w+")
    write(tmpio, rcode)
    close(tmpio)
    try
        run(`Rscript $tmpfile $args`)
    catch e
        # @error "Rscript exited with error." _module=nothing _group=nothing _id=nothing _file=nothing
    finally
        rm(tmpfile)
    end
    return
end

statplot_code = raw"""
#!/usr/bin/env Rscript

if (is.na(packageDescription("argparse")[1])) install.packages("argparse")
library(argparse, quietly = T)
if (is.na(packageDescription("plotly")[1])) install.packages("plotly")
library(plotly, quietly = T, warn.conflicts = F)
if (is.na(packageDescription("ggsci")[1])) install.packages("ggsci")
library("ggsci")
if (is.na(packageDescription("tidyverse")[1])) install.packages("tidyverse")
library("tidyverse")
if (is.na(packageDescription("stringr")[1])) install.packages("stringr")
library(stringr)

parser <- ArgumentParser(prog='atria statplot', description='Trimming performance plots with stat.tsv files')
parser$add_argument('-i', '--input', metavar='STAT.TSV', type='character',
                    nargs='+', default="auto",
                    help='input stat.tsv files; if auto, search files ended with r12.stat.tsv')
parser$add_argument('-o', '--outpref', metavar='PREF', type='character',
                    default="trimmer_accuracy",
                    help='Prefix of output files (default: trimmer_accuracy)')
parser$add_argument('-l', '--legend', metavar="DIR|DIR2|BASE|.",
                    default='DIR',
                    help = 'legend name: DIR = dir name, BASE = base name, DIR2 = second dir name (default: DIR)')
parser$add_argument('-a', '--adapter-length', metavar="INT",
                    help = 'Main adapter length for stats')
parser$add_argument('-F', '--no-format', action='store_true',
                    help = 'Stop formatting legend names (keep ^[A-Za-z]*, first uppercase)')

args <- parser$parse_args()

if (args$input == "auto") {
    # search files ended with r12.stat.tsv
    args$input <- try(system("find * -maxdepth 2 -type f | grep r12.stat.tsv", intern = TRUE))
    if (length(args$input) == 0) {
        stop("-i auto failed: Cannot find files ends with r12.stat.tsv under current working directory and 2-depth subdirectories.")
    }
    writeLines("-i auto: successfully found files:")
    writeLines(args$input)
}

dt = data.frame()
args$legend = toupper(args$legend)
args$input = sort(args$input)
for (i in args$input) {
    dt_i = read.delim(i)
    if (args$legend == 'DIR'){
        dt_i_file = dirname(i)
    } else if (args$legend == 'DIR2'){
        dt_i_file = str_split(i, pattern = '/')[[1]][2]
    } else if (args$legend == 'BASE') {
        dt_i_file = sub(".stat.tsv$", "", basename(i))
    } else {
        dt_i_file = sub(".stat.tsv$", "", i)
    }

    # format
    if (!args$no_format) {
        dt_i_file = paste(toupper(substr(dt_i_file, 1, 1)), substr(dt_i_file, 2, nchar(dt_i_file)), sep="")

        dt_i_file = str_extract(dt_i_file, "^([A-Za-z ]*)")
    }

    # find adapter length if indicated in input path name
    adapter_length = as.integer(str_extract(str_extract(i, "adapter_length_([0-9]+)"), "\\d+"))

    dt_i$file = dt_i_file
    dt_i$adapter_length = adapter_length
    dt_i$deviation_undertrim <- -dt_i$deviation_undertrim
    dt_i <- dt_i[order(dt_i$insert_size),]
    dt = rbind(dt, dt_i)
}

unique_files = unique(dt$file)
file_count = length(unique_files)

dt$has_adapter <- dt$seq_length > dt$insert_size

##### Full length adapter (if -a specified, else all adapter data). Stats of all error profiles.

args$adapter_length <- as.integer(args$adapter_length)

if (is.integer(args$adapter_length) && length(as.integer(args$adapter_length))==1 ) {
    dt_main <- filter(dt, adapter_length == args$adapter_length)
} else {
    dt_main <- dt
}

group_dt <- dt_main %>%
    group_by(file, has_adapter)

dt_summary <- group_dt %>% summarise(
    `Accurate Trim` = mean(precision),
    # accuracy_1bp = mean(rate_precision_in1) - mean(precision),
    `Over Trim = 1 bp` = mean(rate_overtrim) - mean(rate_overtrim_gt1),
    `Under Trim = 1 bp` =  mean(rate_undertrim) - mean(rate_undertrim_gt1),
    `Over Trim > 1 bp` = mean(rate_overtrim_gt1),
    `Under Trim > 1 bp` = mean(rate_undertrim_gt1)
)

dt_summary_long <- pivot_longer(dt_summary, `Accurate Trim`:`Under Trim > 1 bp`)

dt_summary_long$name <- factor(dt_summary_long$name, ordered = TRUE, levels = c("Accurate Trim", "Over Trim = 1 bp", "Under Trim = 1 bp", "Over Trim > 1 bp", "Under Trim > 1 bp") %>% rev)
rate_colors <- pal_locuszoom("default", alpha = 0.7)(5)


p11 <- plot_ly(dt_summary_long %>% filter(has_adapter==T),
               type='bar', orientation = 'h',
               colors=rate_colors,
               y = ~file, x = ~value, color=~name,
               legendgroup = "Statistics", showlegend=T) %>%
    layout(xaxis = list(title = 'Accumulate Rate (With Adapter)', tickformat=".1%"),
           yaxis = list(title = '', autorange='reversed'),
           barmode = 'stack')

p12 <- plot_ly(dt_summary_long %>% filter(has_adapter==F),
               type='bar', orientation = 'h',
               colors=rate_colors,
               y = ~file, x = ~value, color=~name,
               legendgroup = "Statistics", showlegend=F) %>%
    layout(xaxis = list(title = 'Accumulate Rate (No Adapter)', tickformat=".1%"),
           yaxis = list(title = '', autorange='reversed', side="left"),
           barmode = 'stack')

##### line plot: x = error rate, y = accuracy

group_dt <- dt_main %>% group_by(file, has_adapter, error_rate)

dt_summary <- group_dt %>% summarise(
    `Accurate Trim` = mean(precision),
    accuracy_1bp = mean(rate_precision_in1)
)

p21 <- plot_ly(dt_summary %>% filter(has_adapter==T)) %>%
    add_trace(
               mode='lines',
               y = ~`Accurate Trim`, x = ~error_rate, color=~file,
               legendgroup = "Trimmers", showlegend=T) %>%
    layout(yaxis = list(title = 'Accurate Trim Rate', tickformat=".1%"),
           xaxis = list(title = 'Error Rate (With Adapter)'))

p22 <- plot_ly(dt_summary %>% filter(has_adapter==F)) %>%
    add_trace(
        mode='lines',
        y = ~`Accurate Trim`, x = ~error_rate, color=~file,
        legendgroup = "Trimmers", showlegend=F) %>%
    layout(yaxis = list(title = 'Accurate Trim Rate', side="left", tickformat=".1%"),
           xaxis = list(title = 'Error Rate (No Adapter)'))

p31 <- plot_ly(dt_summary %>% filter(has_adapter==T)) %>%
    add_trace(
        mode='lines',
        y = ~accuracy_1bp, x = ~error_rate, color=~file,
        legendgroup = "Trimmers", showlegend=F) %>%
    layout(yaxis = list(title = 'Accurate Trim Rate (Allow 1 Bp Error)', tickformat=".1%"),
           xaxis = list(title = 'Error Rate (With Adapter)'))

p32 <- plot_ly(dt_summary %>% filter(has_adapter==F)) %>%
    add_trace(
        mode='lines',
        y = ~accuracy_1bp, x = ~error_rate, color=~file,
        legendgroup = "Trimmers", showlegend=F) %>%
    layout(yaxis = list(title = 'Accurate Trim Rate (Allow 1 Bp Error)', side="left", tickformat=".1%"),
           xaxis = list(title = 'Error Rate (No Adapter)'))

##### line plot: x = adapter length, y = accuracy

group_dt <- dt %>% group_by(file, has_adapter, adapter_length)

dt_summary <- group_dt %>% summarise(
    `Accurate Trim` = mean(precision),
    accuracy_1bp = mean(rate_precision_in1)
)

adapter_lengths <- sort(unique(dt_summary$adapter_length))
dtick <- if (length(adapter_lengths)>1) (adapter_lengths[2] - adapter_lengths[1]) else 4
tick0 <- min(adapter_lengths)

p41 <- plot_ly(dt_summary %>% filter(has_adapter==T)) %>%
    add_trace(
        mode='lines',
        y = ~`Accurate Trim`, x = ~adapter_length, color=~file,
        legendgroup = "Trimmers", showlegend=F) %>%
    layout(yaxis = list(title = 'Accurate Trim Rate', tickformat=".1%"),
           xaxis = list(title = 'Library Adapter Length (With Adapter)', dtick = dtick, tick0 = tick0))

p42 <- plot_ly(dt_summary %>% filter(has_adapter==F)) %>%
    add_trace(
        mode='lines',
        y = ~`Accurate Trim`, x = ~adapter_length, color=~file,
        legendgroup = "Trimmers", showlegend=F) %>%
    layout(yaxis = list(title = 'Accurate Trim Rate', side="left", tickformat=".1%"),
           xaxis = list(title = 'Library Adapter Length (No Adapter)', dtick = dtick, tick0 = tick0))

p51 <- plot_ly(dt_summary %>% filter(has_adapter==T)) %>%
    add_trace(
        mode='lines',
        y = ~accuracy_1bp, x = ~adapter_length, color=~file,
        legendgroup = "Trimmers", showlegend=F) %>%
    layout(yaxis = list(title = 'Accurate Trim Rate (Allow 1 Bp Error)', tickformat=".1%"),
           xaxis = list(title = 'Library Adapter Length (With Adapter)', dtick = dtick, tick0 = tick0))

p52 <- plot_ly(dt_summary %>% filter(has_adapter==F)) %>%
    add_trace(
        mode='lines',
        y = ~accuracy_1bp, x = ~adapter_length, color=~file,
        legendgroup = "Trimmers", showlegend=F) %>%
    layout(yaxis = list(title = 'Accurate Trim Rate (Allow 1 Bp Error)', side="left", tickformat=".1%"),
           xaxis = list(title = 'Library Adapter Length (No Adapter)', dtick = dtick, tick0 = tick0))

p_all<- subplot(p11, p12, p21, p22, p41, p42,
                nrows = 3, shareX = F, shareY = T) %>% layout(showlegend=T)

outname = paste0(args$outpref, '.plots', ".html")
htmlwidgets::saveWidget(as_widget(p_all), outname)
writeLines(paste("atria statplot: write to", outname))

"""
