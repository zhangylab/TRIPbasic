#arg1: basis loop file
#arg2: qvalue threshold, default: 0.1
setwd("ov_analysis")
require(tidyverse)

args = commandArgs(trailingOnly = T)
loopbasisfn = args[1]
qval_threshold = as.numeric(args[2])

an1=read.table("an1_sorted_ov.bed",sep="\t")
an2=read.table("an2_sorted_ov.bed",sep="\t")
an1=an1[order(an1[,4]),]
an2=an2[order(an2[,4]),]
al=cbind(an1,an2)
al$ov=paste(al[,5],al[,11],sep="|")


loops=al[al$ov!="|",]
#db_loops=al[grepl("[0-9]+\\|[0-9]+",al$ov),]

loopori = read.table(loopbasisfn)
colnames(loopori)[c(9,15,25:26)]=c("isPeak1","isPeak2","pval", "qval")
allinfotab = merge(loopori, al, by.x=1:6, by.y=c(1:3, 7:9))

df=allinfotab[order(allinfotab$pval),]
df$type1 = df$ov!="|"
df$type2 = grepl("[0-9]+\\|[0-9]+",df$ov)
df$type3 = df$isPeak1 & df$isPeak2
df$type4 = df$type3 | df$type2
df$type5 = (df$isPeak1 | grepl("[0-9]+\\|",df$ov)) & (df$isPeak2 | grepl("\\|[0-9]+",df$ov))

fnbasis_matches = str_match(loopbasisfn, ".*FitHiChIP_loops\\/(.*)\\/([0-9]*)\\/")
fnbasis = paste0(fnbasis_matches[2], "_res", fnbasis_matches[3], "_")

determine_best_pval_cutoff = function()
{
    pvals = exp(-seq(-log(0.05), -log(0.000000001),by=0.25))

    test_FDR = function(pval_thre)
    {
        #message(".")
        df$type6 = ((df$isPeak1 | grepl("[0-9]+\\|",df$ov)) & (df$isPeak2 | grepl("\\|[0-9]+",df$ov))) | df$pval < pval_thre
        df$FP = cumsum(df$type6 == FALSE)
        df$TP = cumsum(df$type6 == TRUE)
        df$FDR = df$FP / (df$FP + df$TP)
        sum(df$FDR<0.05)/nrow(df)
    }

    fracs = sapply(pvals, test_FDR)

    x = seq(length(fracs))
    y = fracs
    d1 <- diff(y) / diff(x) # first derivative
    d2 <- diff(d1) / diff(x[-1]) # second derivative

    library(shiny)

    ui <- basicPage(
    plotOutput("plot1", click = "plot_click"),  
    plotOutput("plot2"),
    actionButton("myBtn", "Selection Done!"),
    verbatimTextOutput("info")
    )

    server <- function(input, output) {
        dataInput <- reactive({
            input$plot_click$x
        })

        output$plot1 <- renderPlot({        
            plot(abs(d2), ylab="Absolute second order difference", xlab="index");
            
        })

        output$plot2 <- renderPlot({        
            plot(fracs, ylab="Fraction of loops with qvalue < 0.05", xlab="index")
            abline(v=dataInput(),lty="dotted")
        })

    output$info <- renderText({
        paste0("x=", input$plot_click$x, "\ny=", input$plot_click$y)})

    observe({
        if(input$myBtn > 0){
            stopApp(floor(input$plot_click$x)+2)
        }})
    }

    pval_index = runApp(list(ui=ui,server=server))
    pval_cutoff = pvals[pval_index]

    require(ggplot2)
    require(patchwork)
    df2 = data.frame(x=1:length(pvals), y1=fracs, y2=c(NA, NA, abs(d2)))
    p1 = ggplot(df2) + geom_point(aes(x, y1)) + 
        geom_point(x=pval_index,y=fracs[pval_index],col="red")+
        ylab("Fraction of loops with q-value < 0.05")+xlab("Index")+
            theme(
        axis.text = element_text(family="Helvetica", size = 11, color="black"), 
        axis.title.x = element_text(family="Helvetica", size = 11, color="black"), 
        axis.title.y = element_text(family="Helvetica", size = 11, color="black"), 
        axis.text.x=element_text(size=11,hjust=1),
        #axis.text = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(fill=NA, colour = "black"),
        strip.text = element_text(size = 12, color = "black"), strip.background = element_blank(),
        text=element_text(family="Helvetica", size=11))

    p2 = ggplot(df2) + geom_point(aes(x, y2)) + 
        geom_point(x=pval_index,y=d2[pval_index-2],col="red")+
        ylab("Absolute second-order difference")+xlab("Index")+
            theme(
        axis.text = element_text(family="Helvetica", size = 11, color="black"), 
        axis.title.x = element_text(family="Helvetica", size = 11, color="black"), 
        axis.title.y = element_text(family="Helvetica", size = 11, color="black"), 
        axis.text.x=element_text(size=11,hjust=1),
        #axis.text = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(fill=NA, colour = "black"),
        strip.text = element_text(size = 12, color = "black"), strip.background = element_blank(),
        text=element_text(family="Helvetica", size=11))
    p1 / p2
    ggsave(paste0("../", fnbasis, "selection.tiff"), width=6, height = 6)
    pval_cutoff
}

pval_cutoff = determine_best_pval_cutoff()

df$type6 = ((df$isPeak1 | grepl("[0-9]+\\|",df$ov)) & (df$isPeak2 | grepl("\\|[0-9]+",df$ov))) | df$pval < pval_cutoff

df$FP = cumsum(df$type6 == FALSE)
df$TP = cumsum(df$type6 == TRUE)

df$FDR = df$FP / (df$FP + df$TP)

stats_txt = ""
stats_txt = paste0(stats_txt, "#original loops: ", nrow(al), "\n")
stats_txt = paste0(stats_txt, "#original loops without any feature overlaps: ", sum(al$ov=="|"), " (", sum(al$ov=="|")/nrow(al), ")", "\n")
stats_txt = paste0(stats_txt, "pvalue threshold: ", pval_cutoff, "\n")
stats_txt = paste0(stats_txt, "FDR threshold: ", qval_threshold, "\n")
stats_txt = paste0(stats_txt, "#filtered loops by FDR: ", sum(df$FDR<qval_threshold), " (", sum(df$FDR<qval_threshold)/nrow(df), ")", "\n")
fdrsummay = quantile(df$FDR, seq(0,1,0.1))
fdrsummay = paste(names(fdrsummay), fdrsummay, sep = ":", collapse = ", ")

stats_txt = paste0(stats_txt, "#filtered loops by ENCODE/1Dpeaks overlap: ", sum(df$type5), " (", sum(df$type5)/nrow(df), ")", "\n")
stats_txt = paste0(stats_txt, "#filtered loops by ENCODE/1Dpeaks overlap or significant loops: ", sum(df$type6), " (", sum(df$type6)/nrow(df), ")", "\n")
stats_txt = paste0(stats_txt, "FDR summary: ", fdrsummay, "\n")


# > sum(df$FDR<0.1)
# [1] 72270
# > sum(df$FDR<0.05)
# [1] 32354
# > sum(df$FDR<0.2)
# [1] 209845
# > sum(df$type2)
# [1] 113879
setwd("..")
writeLines(stats_txt, paste0(fnbasis, "stats.txt"))

fn = basename(loopbasisfn)
fn = str_match(fn, "(.*)\\.")[2]
fn = paste0(fn, "_filtered_newq", qval_threshold,".bedpe")
write.table(df[df$FDR<qval_threshold, c(1:6, 7, 9, 15, 25, 26)], fn, sep="\t", col.names = F, row.names = F, quote=F)

fn = basename(loopbasisfn)
fn = str_match(fn, "(.*)\\.")[2]
fn = paste0(fn, "_filtered_by_ENCODE_1Dpeaks_ov", ".bedpe")
write.table(df[df$type5,], fn, sep="\t", col.names = F, row.names = F, quote=F)

fn = basename(loopbasisfn)
fn = str_match(fn, "(.*)\\.")[2]
fn = paste0(fn, "_filtered_by_ENCODE_1Dpeaks_ov_or_sig_loops", ".bedpe")
write.table(df[df$type6,], fn, sep="\t", col.names = F, row.names = F, quote=F)

fn = basename(loopbasisfn)
fn = str_match(fn, "(.*)\\.")[2]
fn = paste0(fn, "_full_information.txt")

df_annot = df[,c(1:6, 9, 15, 25, 26, 33:42)]
colnames(df_annot)[1:6] = c("chr1", "start1", "end1", "chr2", "start2", "end2")
write.table(df_annot, fn, sep="\t", col.names = T, row.names = F, quote=F)