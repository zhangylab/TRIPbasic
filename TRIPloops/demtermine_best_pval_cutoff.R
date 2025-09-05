args = commandArgs(trailingOnly = T)
df_fn = args[1]
#df_fn="/STORE3/TRIP/24-10-15_H3K27ac/results/FitHiChIP_loops/TRIP_H3K27ac_pool_967M/500/FitHiChIP_Peak2ALL_b500_L1000_U2000000/P2PBckgr_0/Coverage_Bias/FitHiC_BiasCorr/filtered_loops/TRIP_H3K27ac_pool_967M_500.interactions_FitHiC_dmin3000_p0.05_full_information.txt"

df = read.table(df_fn, header=T)
#pvals = exp(-seq(-log(0.05), -log(0.000000001),by=-log(0.8)))
#pvals = exp(-seq(-log(0.05), -log(0.000000001),by=0.1))
#pvals = exp(-seq(-log(0.05), -log(0.000000001),by=0.05))
pvals = exp(-seq(-log(0.05), -log(0.000000001),by=0.25))
test_FDR = function(pval_thre)
{
    message(".")
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
#plot(abs(d2));abline(v=18,lty="dotted")
#plot(fracs);abline(v=18-1,lty="dotted")
#pvals[17]

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
        plot(abs(d2), ylab="Absolute 2nd order differences", xlab="index");
        
    })

    output$plot2 <- renderPlot({        
        plot(fracs, ylab="Fractions of loops wiht qvalue < 0.05", xlab="index")
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
pvals[pval_index]
#shinyApp(ui, server)
#stopApp()

require(ggplot2)
require(patchwork)
df2 = data.frame(x=1:length(pvals), y1=fracs, y2=c(NA, NA, abs(d2)))
p1 = ggplot(df2) + geom_point(aes(x, y1)) + 
    geom_point(x=pval_index,y=fracs[pval_index],col="red")+
    ylab("Fractions of loops wiht q-value < 0.05")+xlab("index")+
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
    ylab("Absolute 2nd order differences")+xlab("index")+
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
