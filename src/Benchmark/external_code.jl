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

parser <- ArgumentParser(prog='atria peReadStatsPlot', description='Trimming performance plots with stat.tsv files')
parser$add_argument('-i', '--input', metavar='STAT.TSV', type='character',
                    required=TRUE, nargs='+',
                    help='[REQUIRED] input stat.tsv files')
parser$add_argument('-o', '--outpref', metavar='PREF', type='character',
                    default="trimmer_performance",
                    help='Prefix of output files (default: trimmer_performance)')
parser$add_argument('-t', '--title',
                    default="", help = "suffix to figure title, eg. [R1] (default: nothing)")
parser$add_argument('-l', '--legend', metavar="DIR|BASE|.",
                    default='DIR',
                    help = 'legend name: DIR = dir name, BASE = base name, others = name (default: DIR)')
parser$add_argument('-H', '--highlight', metavar="STRING",
                    default="LIAN!GUN@JIAN#PAN%BA*SADKJWIDJSAPJDSLKJWQWASD",
                    help = 'a name in legend to be highlight')

args <- parser$parse_args()

seq_length_vlines <- function(p, x, y0 = 0, y1 = 1, showlegend = TRUE) {
    for (xi in x) {
        linename <- paste('Read Length:', xi)
        p <- add_trace(
            p, type = 'scatter', mode = 'lines',
            x = c(xi, xi), y = c(y0, y1), line=list(color=rgb(0,0,0), dash = '2px', width = 0.5),
            name = linename, hoverinfo = "name", hoverlabel = list(namelength=50),
            legendgroup = "Read Length", showlegend = showlegend
        )
    }
    p
}

legend_placeholder <- function(p) {
    add_trace(
        p, type = 'scatter', mode = 'markers', x = 0, y = 0, name = ' ',
        marker=list(color=rgb(0,0,0,0))
    )
}

yrange <- function(vec) {
    rg = range(vec)
    delta = (rg[2] - rg[1]) * 0.05
    c(rg[1]-delta, rg[2]+delta)
}
markers = c("101", "105", "124", "104")

dt = data.frame()
args$legend = toupper(args$legend)
args$input = sort(args$input)
for (i in args$input) {
    dt_i = read.delim(i)
    if (args$legend == 'DIR'){
        dt_i_file = dirname(i)
    } else if (args$legend == 'BASE') {
        dt_i_file = sub(".stat.tsv$", "", basename(i))
    } else {
        dt_i_file = sub(".stat.tsv$", "", i)
    }
    dt_i$file = dt_i_file
    dt_i$deviation_undertrim <- -dt_i$deviation_undertrim
    dt_i <- dt_i[order(dt_i$insert_size),]
    dt = rbind(dt, dt_i)
}

unique_files = unique(dt$file)
file_count = length(unique_files)
file_num_to_highlight = NULL
dt$imarker = 0
dt$idash = 0
for (i in 1:file_count) {
    ind_file <- dt$file == unique_files[i]
    dt$imarker[ind_file] <- i %% 5
    dt$idash[ind_file] <- i %% 3

    if (args$highlight == unique_files[i]){
        file_num_to_highlight = i
    }
}

# highlight
pal_colors = pal_d3("category20")(file_count)
if (!is.null(file_num_to_highlight)) {
    pal_colors_replace = pal_colors[file_num_to_highlight]
    pal_colors[file_num_to_highlight] <- pal_colors[1]
    pal_colors[1] <- pal_colors_replace
}

seq_lenths = unique(dt$seq_length)

######## Precision
p1 <- plot_ly(dt) %>%
    add_trace(
        type = 'scatter', x = ~insert_size, y = ~precision, color = ~file, frame = ~error_rate, mode = 'lines+markers',
        legendgroup = ~file, showlegend = TRUE, colors=pal_colors,
        hoverlabel = list(namelength = -1), line=list(width = 1.5, opacity = 0.82, dash = ~idash),
        marker = list(symbol = ~imarker, symbols=markers, size = 6, opacity = 0.82)
    ) %>%
    legend_placeholder() %>%
    seq_length_vlines(seq_lenths, y0 = 0, y1 = 1, showlegend = TRUE) %>%
    layout(#title = paste('Precision', args$title),
        xaxis = list(title = 'Original Insert Size'),
        yaxis = list(title = 'Precision', range=yrange(dt$precision))
    )

######## Rate Overtrim
p2 <- plot_ly(dt) %>%
    add_trace(
        type = 'scatter', x = ~insert_size, y = ~rate_overtrim, color = ~file, frame = ~error_rate, mode = 'lines+markers',
        legendgroup = ~file, showlegend = FALSE, colors=pal_colors,
        hoverlabel = list(namelength = -1), line=list(width = 1.5, opacity = 0.82, dash = ~idash),
        marker = list(symbol = ~imarker, symbols=markers, size = 6, opacity = 0.82)
    ) %>%
    seq_length_vlines(seq_lenths, y0 = 0, y1 = 1, showlegend = FALSE) %>%
    layout(#title = paste('Over-trim Rate', args$title),
        xaxis = list(title = 'Original Insert Size'),
        yaxis = list(title = 'Over-trim Rate', range=yrange(dt$rate_overtrim))
    )

######## Rate Undertrim
p3 <- plot_ly(dt) %>%
    add_trace(
        type = 'scatter', x = ~insert_size, y = ~rate_undertrim, color = ~file, frame = ~error_rate, mode = 'lines+markers',
        legendgroup = ~file, showlegend = FALSE, colors=pal_colors,
        hoverlabel = list(namelength = -1), line=list(width = 1.5, opacity = 0.82, dash = ~idash),
        marker = list(symbol = ~imarker, symbols=markers, size = 6, opacity = 0.82)
    ) %>%
    seq_length_vlines(seq_lenths, y0 = 0, y1 = 1, showlegend = FALSE) %>%
    layout(#title = paste('Under-trim Rate', args$title),
        xaxis = list(title = 'Original Insert Size'),
        yaxis = list(title = 'Under-trim Rate', range=yrange(dt$rate_undertrim))
    )

######## Deviation
p4 <- plot_ly(dt) %>%
    add_trace(
        type = 'scatter', x = ~insert_size, y = ~deviation, color = ~file, frame = ~error_rate, mode = 'lines+markers',
        legendgroup = ~file, showlegend = TRUE, colors=pal_colors,
        hoverlabel = list(namelength = -1), line=list(width = 1.5, opacity = 0.82, dash = ~idash),
        marker = list(symbol = ~imarker, symbols=markers, size = 6, opacity = 0.82)
    ) %>%
    legend_placeholder() %>%
    seq_length_vlines(seq_lenths, y0 = 0, y1 = max(dt$deviation), showlegend = TRUE) %>%
    layout(
        xaxis = list(title = 'Original Insert Size'),
        yaxis = list(title = 'Median Deviation', range=yrange(dt$deviation))
    )

######## Deviation Over-trim
p5 <- plot_ly(dt) %>%
    add_trace(
        type = 'scatter', x = ~insert_size, y = ~deviation_overtrim, color = ~file, frame = ~error_rate, mode = 'lines+markers',
        legendgroup = ~file, showlegend = FALSE, colors=pal_colors,
        hoverlabel = list(namelength = -1), line=list(width = 1.5, opacity = 0.82, dash = ~idash),
        marker = list(symbol = ~imarker, symbols=markers, size = 6, opacity = 0.82)
    ) %>%
    seq_length_vlines(seq_lenths, y0 = 0, y1 = max(dt$deviation_overtrim), showlegend = FALSE) %>%
    layout(
        xaxis = list(title = 'Original Insert Size'),
        yaxis = list(title = 'Median Deviation\n(Over-trim Portion)', range=yrange(dt$deviation_overtrim))
    )

######## Deviation Under-trim
p6 <- plot_ly(dt) %>%
    add_trace(
        type = 'scatter', x = ~insert_size, y = ~deviation_undertrim, color = ~file, frame = ~error_rate, mode = 'lines+markers',
        legendgroup = ~file, showlegend = FALSE, colors=pal_colors,
        hoverlabel = list(namelength = -1), line=list(width = 1.5, opacity = 0.82, dash = ~idash),
        marker = list(symbol = ~imarker, symbols=markers, size = 6, opacity = 0.82)
    ) %>%
    seq_length_vlines(seq_lenths, y0 = 0, y1 = max(dt$deviation_undertrim), showlegend = FALSE) %>%
    layout(
        xaxis = list(title = 'Original Insert Size'),
        yaxis = list(title = 'Median Deviation\n(Under-trim Portion)', range=yrange(dt$deviation_undertrim))
    )

############## integrate subplots

p123 <- subplot(p1, p2, p3,
                nrows = 3, shareX = T, shareY = T) %>%
    animation_opts(frame = 500, transition = 0) %>%
    animation_slider(currentvalue = list(prefix = "Error Rate = ")) %>%
    animation_button(label = "Error Rate", hide = TRUE) %>%
    layout(title = paste('Preciseness Benchmark', args$title))

p456 <- subplot(p4, p5, p6,
                nrows = 3, shareX = T, shareY = T) %>%
    animation_opts(frame = 500, transition = 0) %>%
    animation_slider(currentvalue = list(prefix = "Error Rate = ")) %>%
    animation_button(label = "Error Rate", hide = TRUE) %>%
    layout(title = paste('Median Deviation Benchmark', args$title))

outname123 = paste0(args$outpref, '.', 'preciseness', args$title, ".html")
outname456 = paste0(args$outpref, '.', 'deviation', args$title, ".html")

htmlwidgets::saveWidget(as_widget(p123), outname123)
htmlwidgets::saveWidget(as_widget(p456), outname456)

writeLines(paste("peReadSimulatorStatsPlot: Output:", outname123))
writeLines(paste("peReadSimulatorStatsPlot: Output:", outname456))

"""
