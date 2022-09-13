# This script contains all of the code use to plot/visualize data in the mansuscript.


PlottingDir=paste0(getwd(), "/Figures")

LoadPackages = function() {
    library(ggplot2) # Used for plotting
    library(scales) # Used to transform axes into scientific notation
    library(ggpmisc) # Used for Sup. Fig 4 to add R2 to plot
    library(cowplot) # Used to stitch figures together
    library(data.table) # Used to aggregate individual plots
    #library(ggpubr)
    #require(gridExtra)
    #require(cowplot)
}

AddAxesScales = function(p) {
    plot_out = p +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) 
    return(plot_out)
}

AddColor = function(p, reversed) {
    if (reversed) {
        plot_out = p +
        scale_fill_manual(values = c("red","#7FBD32"), labels=c( 'Passengers', 'Drivers')) +
        scale_colour_manual(values = c("red","#7FBD32"), labels=c('Passengers', 'Drivers'))  
    } else {
        plot_out = p +
        scale_fill_manual(values = c("#7FBD32","red"), labels=c('Drivers', 'Passengers')) +
        scale_colour_manual(values = c("#7FBD32","red"), labels=c('Drivers', 'Passengers'))  
    }
    return(plot_out)
}

AddTheme = function(p) {
    plot_out = p +
    theme_classic() + 
    theme(legend.title=element_blank(), plot.title = element_text(hjust = 0.5, size=9, 
        face = "bold"), text = element_text(size=9)) 
    return(plot_out)
}

AdddNdscv = function(df, p) {
    df = subset(df, df$dNdS == 'dNdScv')
    df$high = NULL; df$low=NULL # Don't plot confidence intervals for dNdScv
    plot_out = p +
    geom_line(data=subset(df, df$Group == 'Passengers'), aes(x=Bin, y=true, color=Group)) +
    geom_line(data=subset(df, df$Group == 'Drivers_Bailey'), aes(x=Bin, y=true, color=Group)) 
    return(plot_out)
}


################
### Figure 2 ###
################
    
PlotFig2A = function(df) {
    # dN/dS Plot
    dNdSPlot = ggplot(data=df, aes(y = true, x=Bin, fill=Group, color=Group)) + 
            geom_point(size=.8) + geom_line(size=0.65) + # Add line/point
            labs(x="Total Number of Substitutions", y="dN/dS") + # Label axes
            geom_hline(aes(yintercept = 1)) + # Add neutral dN/dS line
            scale_y_continuous(limits= c(0.25,9), breaks = c(0.25,0.5,1,2,4,8), trans = log_trans()) +
            geom_ribbon(aes(ymin=low, ymax=high, fill=Group),  alpha=0.4, linetype=0)  ## Add CI
    dNdSPlot = AddTheme(dNdSPlot)
    dNdSPlot = AddColor(dNdSPlot, FALSE)
    dNdSPlot = AddAxesScales(dNdSPlot) + theme(legend.position='none') +   theme(strip.background = element_blank()) + facet_wrap(~dNdS) 
    # Counts Plot
    CountsTable = subset(df, df$Group == 'Passengers')
    CountsTable = unique(CountsTable[c('xbin','NumPatientsInBin','dNdS')])
    CountsTable$Bin = data.frame(do.call('rbind', strsplit(as.character(CountsTable$xbin),'-',fixed=TRUE)))$X1
    Counts = ggplot(data=CountsTable, aes(y=as.numeric(as.character(NumPatientsInBin)), 
        x=as.numeric(as.character(Bin)))) + geom_bar(position = "dodge", stat = "identity") + 
        scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
        theme_classic() + labs(y="Number of Tumors", x="") + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
        ggtitle('All SNVs')
    Counts = AddTheme(Counts)
    Counts = AddAxesScales(Counts) + facet_wrap(~dNdS) +  theme(strip.background = element_blank()) 
    PlotOut = plot_grid(Counts, dNdSPlot, ncol = 1, rel_heights = c(0.35, 0.66))
    #ggsave(paste0(PlottingDir, "Fig2A.pdf"))
    return(PlotOut)
}


PlotFig2B = function(df) {
    df[c('low','true','high','Bin')] = sapply(df[c('low','true','high','Bin')], as.numeric)
    PolyPlot = ggplot(data=df, aes(y = true, x=Bin, fill=Group, color=Group)) + 
            geom_point(size=.8) + geom_line(size=0.65) + # Add line/point
            labs(y="Fraction of Pathogenic Mutations", x="Total Number of Substitutions") + # Label axes
            geom_hline(aes(yintercept = 0.606300035)) + ### global COSMIC polyphen scores
            geom_ribbon(aes(ymin=low, ymax=high, fill=Group),  alpha=0.4, linetype=0) ## Add CI
    PolyPlot = AddColor(PolyPlot, TRUE)  
    PolyPlot = AddTheme(PolyPlot)
    PolyPlot = AddAxesScales(PolyPlot) + theme(legend.position='none')
    # Counts plot
    CountsTable = subset(df, df$Group == 'False')
    CountsTable = unique(CountsTable[c('xbin','NumPatientsInBin')])
    CountsTable$Bin = data.frame(do.call('rbind', strsplit(as.character(CountsTable$xbin),'-',fixed=TRUE)))$X1
    Counts = ggplot(data=CountsTable, aes(y=as.numeric(as.character(NumPatientsInBin)), x=as.numeric(as.character(Bin)))) + 
        geom_bar(position = "dodge", stat = "identity") + 
        scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
        theme_minimal() + labs(y="Number of Tumors", x="") + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
    Counts = AddTheme(Counts)
    Counts = AddAxesScales(Counts) + ggtitle('Polyphen2')
    PlotOut = plot_grid(Counts, PolyPlot, ncol = 1, rel_heights = c(0.35, 0.66))
    return(PlotOut)
}


PlotFig2C = function(df) {
    dEdIPlot = ggplot(data=df, aes(y = breakpointFreq_mean, x=as.numeric(as.character(Bin)), fill=Track, color=Track)) + 
            geom_point(size=.8) + geom_line(size=0.65) + # Add line/point
            labs(x="Total Number of CNAs", y="dE/dI") +
            geom_hline(aes(yintercept = 1)) + # Add neutral dE/dI line
            scale_y_continuous(limits= c(0.5,5), breaks = c(0.5,1,2,4), trans = log_trans()) +
            geom_ribbon(aes(ymin=breakpointFreq_low, ymax=breakpointFreq_high, fill=Track),  alpha=0.4, linetype=0) + # Add CI
            facet_wrap(~ Length_Category) + theme(legend.position='none')
    dEdIPlot = AddColor(dEdIPlot, FALSE)  
    dEdIPlot = AddTheme(dEdIPlot) 
    dEdIPlot = AddAxesScales(dEdIPlot) + theme(strip.background = element_blank()) + theme(legend.position='none')
    # Counts plot
    CountsTable = subset(df, df$Group == 'Passengers')
    CountsTable = unique(CountsTable[c('xbin','NumberOfTumorsInBin','Length_Category')])
    CountsTable$Bin = data.frame(do.call('rbind', strsplit(as.character(CountsTable$xbin),'-',fixed=TRUE)))$X1
    Counts = ggplot(data=CountsTable, aes(y=as.numeric(as.character(NumberOfTumorsInBin)), 
        x=as.numeric(as.character(Bin)))) + geom_bar(position = "dodge", stat = "identity") + 
        scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
        theme_classic() + labs(y="Number of Tumors", x="") + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
        facet_wrap(~ Length_Category) + ggtitle('CNAs')
    Counts = AddTheme(Counts)
    Counts = AddAxesScales(Counts) + theme(strip.background = element_blank()) 
    PlotOut = plot_grid(Counts, dEdIPlot, ncol = 1, rel_heights = c(0.35, 0.66)) 
    return(PlotOut)
}

PlotFig2D = function(df) {
    dnds = subset(df, (df$Group != 'Drivers') &  (df$Group != 'Passengers') )
    dNdSPlot = ggplot(data=dnds, aes(y = true, x=as.numeric(as.character(Bin)), fill=Group, color=Group)) +
        geom_point(size=.8) + geom_hline(aes(yintercept = 1)) + labs(x="Total Number of Substitutions", y="dN/dS") + 
        theme(text = element_text(size=12)) + theme(legend.position="bottom", legend.box="horizontal") + 
        geom_ribbon(aes(ymin=low, ymax=high, fill=Group),  alpha=0.4, linetype=0) + ## adds CI
        scale_colour_manual(breaks = c('Clonal_Drivers', 'Subclonal_Drivers','Clonal_Passengers', 'Subclonal_Passengers'),
                        values = c("#375F1B",  "#ABEF7B", "#54000b", "#ff576d"),
                        labels = c('Clonal Drivers','Subclonal Drivers','Clonal Passengers',  'Subclonal Passengers')) +
        geom_line(size=0.65) + scale_y_continuous(limits= c(0.19,16), breaks = c(0.25,0.5,1,2,4,8,16), trans = log_trans()) + 
        scale_fill_manual(breaks = c('Clonal_Drivers', 'Subclonal_Drivers','Clonal_Passengers', 'Subclonal_Passengers'),
                        values = c("#375F1B",  "#ABEF7B", "#54000b", "#ff576d"),
                        labels = c('Clonal Drivers', 'Subclonal Drivers','Clonal Passengers', 'Subclonal Passengers')) +
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) 
    dNdSPlot = AddTheme(dNdSPlot)+ theme(legend.position='none') + facet_wrap(~dNdS) +  theme(strip.background = element_blank()) 
    #dNdSPlot = AddAxesScales(dNdSPlot) + theme(legend.position='none')
    # Counts plot
    CountsTable = subset(df, df$Group == 'Passengers')
    CountsTable = unique(CountsTable[c('xbin','NumPatientsInBin','dNdS')])
    CountsTable$Bin = data.frame(do.call('rbind', strsplit(as.character(CountsTable$xbin),'-',fixed=TRUE)))$X1
    Counts = ggplot(data=CountsTable, aes(y=as.numeric(as.character(NumPatientsInBin)), 
        x=as.numeric(as.character(Bin)))) + geom_bar(position = "dodge", stat = "identity") + 
        scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
        theme_classic() + labs(y="Number of Tumors", x="") + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
    Counts = AddTheme(Counts) + facet_wrap(~dNdS) +  theme(strip.background = element_blank()) 
    Counts = AddAxesScales(Counts) + ggtitle('Clonal vs. Subclonal SNVs') 
    PlotOut = plot_grid(Counts, dNdSPlot, ncol = 1, rel_heights = c(0.35, 0.66))
    #ggsave(paste0(PlottingDir, "Fig2D.pdf"))
    return(PlotOut)
}


PlotFig2E = function(df) {
    PlotOut = ggplot(data=df, aes( x=as.character(Label), y= as.numeric(as.character(true)), fill=as.factor(MutationRateGroup))) + 
        geom_bar(stat='identity', position=position_dodge()) + 
        scale_fill_manual(values = c("black","grey"), labels=c('High Mutation Load', 'Low Mutation Load')) + 
        facet_wrap(~ Group, ncol = 4, scales='free_x') + coord_flip() + 
        geom_hline(aes(yintercept = 1), color="#434445", linetype=2) + 
        labs(x='', y='dN/dS') + theme_classic() + 
        theme(legend.position="bottom", legend.box="horizontal") +
        theme(text = element_text(size=9), legend.text=element_text(size=9), legend.title=element_blank()) 
    #ggsave(paste0(PlottingDir, "Fig2E_dNdSBroadCancerTypesDriversAndPassengers.pdf"))
    #print(paste0('Plot saved to: ', PlottingDir))
    return(PlotOut)
}  
      

PlotFig2F1 = function(df) {
    PlotOut = ggplot(data=df, aes( x=as.character(CancerType), y= as.numeric(as.character(true)), fill=as.factor(MutationRateGroup))) + 
        geom_bar(stat='identity', position=position_dodge()) + 
        scale_fill_manual(values = c("black","grey"), labels=c('High Mutation Rates', 'Low Mutation Rates')) + 
        facet_wrap(~ Group, ncol = 4, scales='free_x') + coord_flip() + 
        geom_hline(aes(yintercept = 1), color="#434445", linetype=2) + 
        labs(x='', y='dN/dS') + 
        theme(legend.position="bottom", legend.box="horizontal", legend.title=element_blank()) +
        theme_classic() + theme(text = element_text(size=9)) 
    #ggsave(paste0(PlottingDir, "Fig2F-1_dNdSSpecificCancerTypesDriversAndPassengers.pdf"))
    return(PlotOut)
}


PlotFig2F2 = function(df) {
    PlotOut = ggplot(data=df, aes( x=as.character(CancerType), y= as.numeric(as.character(breakpointFreq_mean)), fill=as.factor(MutationRateGroup))) + 
        geom_bar(stat='identity', position=position_dodge()) + 
        scale_fill_manual(values = c("black","grey"), labels=c('High Mutation Rates', 'Low Mutation Rates')) + 
        facet_wrap(~ Group, ncol = 4, scales='free_x') + coord_flip() + 
        geom_hline(aes(yintercept = 1), color="#434445", linetype=2) + 
        labs(x='', y='dE/dI') + 
        theme(legend.position="bottom", legend.box="horizontal", legend.title=element_blank()) +
        theme_classic() + theme(text = element_text(size=9)) + theme(legend.text=element_text(size=9))
    #ggsave(paste0(PlottingDir, "Fig2F-2_dEdISpecificCancerTypesDriversAndPassengers.pdf"))
    return(PlotOut)
}


PlotFig2G = function(df) {
    PlotOut = ggplot(data=df, aes( x=as.numeric(as.character(Bin)), y= as.numeric(as.character(true)), fill=Group)) + 
        geom_point(size=0.8) + geom_line(size=0.65) +
        scale_colour_manual(values = c("blue","orange","grey","purple"), 
            breaks= c('Chaperonins','HSP90','All','Proteasome'), labels=c('Chaperonins','HSP90','All','Proteasome')) +  
        labs(x="Total Number of Substitutions", y="Mean Expression of Gene Set\n(Relative to an Average Tumor)") + 
        geom_ribbon(aes(ymin=low, ymax=high, fill=Group),  alpha=0.4, linetype=0) + 
        scale_fill_manual(values = c("grey","purple", "blue","orange"), breaks=c('All','Proteasome','Chaperonins','HSP90'), 
            labels=c('All','Proteasome','Chaperonins','HSP90')) +
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10",math_format(10^.x))) +
        theme_classic() 
    PlotOut = AddTheme(PlotOut)
    PlotOut = AddAxesScales(PlotOut) +  theme(legend.position="bottom", legend.direction = "horizontal") +
         theme(legend.text=element_text(size=9))
    #ggsave(paste0(PlottingDir, "Fig2G_Expression.pdf"))
    return(PlotOut)
}

GetFigure2TopLegend = function(df) {
    dNdSPlot = ggplot(data=df, aes(y = true, x=as.numeric(as.character(Bin)), fill=Group, color=Group)) +
        geom_point(size=.8) + geom_hline(aes(yintercept = 1)) + labs(x="Total Number of Substitutions", y="dN/dS") + 
        theme(text = element_text(size=12)) +
        geom_ribbon(aes(ymin=low, ymax=high, fill=Group),  alpha=0.4, linetype=0) + ## adds CI
        scale_colour_manual(breaks = c('Clonal_Drivers', 'Drivers','Subclonal_Drivers','Clonal_Passengers','Passengers', 'Subclonal_Passengers'),
                        values = c("#375F1B", "#7FBD32", "#ABEF7B", "#54000b","#FE1601",  "#ff576d"),
                        labels = c('Clonal Drivers', 'All Drivers','Subclonal Drivers','Clonal Passengers', 'All Passengers', 'Subclonal Passengers')) +
        geom_line(size=0.65) + scale_y_continuous(limits= c(0.19,16), breaks = c(0.25,0.5,1,2,4,6,8,16), trans = log_trans()) + 
        scale_fill_manual(breaks = c('Clonal_Drivers', 'Drivers','Subclonal_Drivers','Clonal_Passengers','Passengers', 'Subclonal_Passengers'),
                        values = c("#375F1B", "#7FBD32", "#ABEF7B", "#54000b","#FE1601",  "#ff576d"),
                        labels = c('Clonal Drivers', 'All Drivers','Subclonal Drivers','Clonal Passengers', 'All Passengers', 'Subclonal Passengers')) 
    dNdSPlot = AddTheme(dNdSPlot)
    dNdSPlot = AddAxesScales(dNdSPlot) + theme(legend.position="bottom", legend.direction = "horizontal") +
        guides(color = guide_legend(nrow = 1)) + guides(fill = guide_legend(nrow = 1)) + theme(legend.text=element_text(size=9))
    return(dNdSPlot)
}

############################
### Supplemental Figures ###
############################

PlotSupFig2A = function(dnds) {
    # Aggregate all replicates and take the mean
    dNdS_Table = aggregate(dnds$neutral_dNdS, by=list(signature=dnds$SignatureNumber, # dN/dS before simulation of neg. selection
            nSites = dnds$ProportionOfSitesWithEachContext, truedNdS = dnds$True_Ka_Ks), FUN=mean) 
    Corrected_dNdS = aggregate(dnds$Ka_pKa_Ks_pKs, by=list(signature=dnds$SignatureNumber, # dN/dS after correcting for bias
            nSites = dnds$ProportionOfSitesWithEachContext, truedNdS = dnds$True_Ka_Ks), FUN=mean) 
    dNdS_Table$Corrected_dNdS = Corrected_dNdS$x
    PlotOut = ggplot(data=dNdS_Table, aes(x=truedNdS, y=Corrected_dNdS, color=x, shape=as.character(nSites))) +   geom_point(size=1.8)  +
      facet_wrap(~signature) + geom_abline(slope=1, intercept=0) + 
      labs(y = expression(frac(dn^(observed)/dn^(permuted),ds^(observed)/ds^(permuted)))) + 
      labs( x=c('True dN/dS ( 1 - Probability of Non-Syn Mut Removed)')) + 
      labs(color = "dN/dS Before Simulating\nPurifying Selection" ) + 
      labs(shape = "Number of Sites Simulated")  + 
      theme(legend.text=element_text(size=9)) + theme(legend.title=element_text(size=9)) + theme_classic()
    return(PlotOut)
}

PlotSupFig2B = function(df) {
    df[c('low','true','high','Bin')] = sapply(df[c('low','true','high','Bin')], as.numeric)
    df$BinID = '1'
    PlotOut = ggplot(data=df, aes(y = true, x= Bin, group=Type, fill=Type, color=Type, linetype=BinID)) + 
            geom_point(size=.8) + geom_line(size=0.65) + 
            scale_colour_manual(values = c("blue","red")) +  
            labs(x="Total Number of Substitutions", y="dN/dS") + # Label axes
            geom_hline(aes(yintercept = 1)) + # Add neutral dN/dS line
            geom_ribbon(aes(ymin=low, ymax=high, fill=Type),  alpha=0.4, linetype=0) + ## Add CI
            scale_fill_manual(values = c("blue", "red")) +  
            theme_classic() + theme(legend.position = "none") +
            scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))
    return(PlotOut)

}

PlotSupFig2C = function(df) {
    print(head(df))
    PlotOut = ggplot(data=df, aes(y = true, x=as.numeric(as.character(Bin)), fill=Group, color=Group, linetype=dNdS)) + 
            geom_point(size=.8) + geom_line(size=0.65) + # Add line/point
            scale_colour_manual(values = c("red"), labels=c('Passengers')) +  
            labs(x="Total Number of Substitutions", y="dN/dS") + # Label axes
            geom_hline(aes(yintercept = 1)) + # Add neutral dN/dS line
            geom_ribbon(aes(ymin=low, ymax=high, fill=Group),  alpha=0.4, linetype=0) + ## Add CI
            theme_classic() + theme(legend.position = "none") +
            scale_y_continuous(limits= c(0.25,9), breaks = c(0.25,0.5,1,2,4,6), trans = log_trans())
    PlotOut = AddAxesScales(PlotOut)
    return(PlotOut)
}


PlotSupFig3 = function(df) {
    PlotOut = ggplot(data=df, aes(y = true, x=as.numeric(as.character(Bin)), fill=Group, color=Group)) + 
            geom_point(size=.8) + geom_line(size=0.65) + # Add line/point
            labs(x="Total Number of Substitutions", y="dN/dS") + # Label axes
            geom_hline(aes(yintercept = 1)) + # Add neutral dN/dS line
            scale_y_continuous(limits= c(0.25,14), breaks = c(0.25,0.5,1,2,4,6,12), trans = log_trans()) +
            geom_ribbon(aes(ymin=low, ymax=high, fill=Group),  alpha=0.4, linetype=0) + ## Add CI
            scale_colour_manual(values = c("#7FBD32","red","blue","purple"), 
                breaks= c('Drivers_Bailey','Passengers','Oncogene','Tumor_Suppressor'), labels=c('Drivers (Bailey et.al 2018)','Passengers','Oncogenes','Tumor Suppressors')) +
            scale_fill_manual(values = c("#7FBD32","red","blue","purple"), breaks=c('Drivers_Bailey','Passengers','Oncogene','Tumor_Suppressor'), 
                labels = c('Drivers (Bailey et.al 2018)','Passengers','Oncogenes','Tumor Suppressors'))  +
            theme_classic() + theme(legend.title=element_blank()) +
            scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
            facet_wrap(~dNdS)
    ggsave(paste0(PlottingDir, "FigS3_dNdSOncogenesAndTumorSuppressors.pdf"), width=6, height=4, units="in")
    print(paste0('Plot saved to: ', PlottingDir))
    return(PlotOut)

}

PlotSupFig4 = function(df, MAF) { 
    PlotOut = ggplot(data=df, aes(y=value, x=BeginBin, fill=variable)) + 
        geom_bar(stat='identity',  position=position_dodge())  +   
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)), limits = c(1,1e4)) +
        theme_classic() + theme(legend.title=element_blank()) +
        scale_fill_manual(values= c('blue','red'),breaks=c('Fraction_Ka','Fraction_Ks'), labels=c('Non-Synonymous','Synonymous')) + 
        labs(x='Total Number of Substitutions', y='Fraction Overlap With\n1000 Genomes Polymorphisms')  +
        theme(text = element_text(size=9), plot.title = element_text(size=9)) 
    if (MAF  < 0.05) {
        PlotOut = PlotOut + ggtitle(paste0("Low-Frequency Polymorphisms (MAF > ", MAF, "%)")) + theme(legend.position = 'none')
    } else if (MAF == 0.05) {
        PlotOut = PlotOut + ggtitle(paste0("Common Polymorphisms (MAF > ", MAF, "%)")) + theme(legend.position = 'none') 
    } else if (MAF == 1) {
        PlotOut = PlotOut + ggtitle('All Polymorphisms')
    }
    print(paste0('Plot saved to: ', PlottingDir))
    return(PlotOut)
}



PlotSupFig5 = function(df) {
    PlotOut = ggplot(data=df, aes(y = true, x=as.numeric(as.character(Bin)), fill=Group, color=Group)) + 
            geom_point(size=.8) + geom_line(size=0.65) + # Add line/point
            labs(x="Total Number of Substitutions", y="dN/dS") + # Label axes
            geom_hline(aes(yintercept = 1)) + # Add neutral dN/dS line
            scale_y_continuous(limits= c(0.25,13), breaks = c(0.25,0.5,1,2,4,6,12), trans = log_trans()) +
            geom_ribbon(aes(ymin=low, ymax=high, fill=Group),  alpha=0.4, linetype=0) + ## Add CI
            scale_colour_manual(values = c("#CBE432","#4AC4EE","#7FBD32","#448D76","red"), 
                breaks=c('Drivers_Bailey','SpecificDriverGenes','Drivers_COSMIC','Drivers_Intogen','Passengers'),
                labels=c('Pan-cancer Drivers (Bailey.et.al 2018)','Cancer-Specific Drivers','Pan-cancer COSMIC Drivers', 'Pan-cancer Intogen Drivers', 'Passengers')) +  
            scale_fill_manual(values = c("#CBE432","#4AC4EE","#7FBD32","#448D76","red"), 
                breaks=c('Drivers_Bailey','SpecificDriverGenes','Drivers_COSMIC','Drivers_Intogen','Passengers'),
                labels=c('Pan-cancer Drivers (Bailey.et.al 2018)','Cancer-Specific Drivers','Pan-cancer COSMIC Drivers', 'Pan-cancer Intogen Drivers', 'Passengers')) +  
            theme_classic() + theme(legend.title=element_blank()) +
            scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
            facet_wrap(~dNdS)
    ggsave(paste0(PlottingDir, "FigS5_dNdSCancerSpecificDrivers.pdf"), width=7, height=4, units="in")
    print(paste0('Plot saved to: ', PlottingDir))
    return(PlotOut)

}

PlotSupFig6A = function(df) {
    PlotOut = ggplot(data=df, aes(y=as.numeric(as.character(purity)), x=as.numeric(as.character(mutRate)))) + 
        geom_point(size=1) + 
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
        theme_classic() + labs(x='Total Number of Substitutions', y='Tumor Purity') + 
        stat_poly_eq(formula = y ~ x, 
                    aes(label = paste( ..rr.label.., sep = "~~~")), label.y = 0.1, label.x = 0.9,
                    parse = TRUE) +
        geom_smooth(method='lm', formula= y~x) + theme( text = element_text(size=9))
    return(PlotOut)
}


PlotSupFig6B = function(df) {
    df$xbin3 = factor(df$xbin2, levels= c('1-3','3-10','10-10000'))
    print(unique(df$xbin3))
    PlotOut = ggplot(data=df, aes(y = as.numeric(as.character(purity)), x=xbin3)) + 
            geom_boxplot(position= "dodge2") + theme_classic() +
            labs(y='Tumor Purity', x='Total Number of Substitutions') +
            scale_colour_manual(values = c('blue')) + theme(legend.position = 'none', text = element_text(size=9))
    return(PlotOut)
}


PlotSupFig6C = function(df) {
    PlotOut = ggplot(data=df, aes(y = as.numeric(as.character(true)), x=as.numeric(as.character(Bin)), fill=Group, color=Group)) + 
            geom_point(size=.8) + geom_line(size=0.65) + # Add line/point
            labs(x="Total Number of Substitutions", y="dN/dS") + # Label axes
            geom_hline(aes(yintercept = 1)) + # Add neutral dN/dS line
            scale_y_continuous(limits= c(0.1,18), breaks = c(0.25,0.5,1,2,4,6, 12), trans = log_trans()) +
            geom_ribbon(aes(ymin=low, ymax=high, fill=Group),  alpha=0.4, linetype=0) + ## Add CI
            facet_wrap(~PurityThreshold) 
    PlotOut = AddTheme(PlotOut)
    PlotOut = AddColor(PlotOut, FALSE)
    PlotOut = AddAxesScales(PlotOut) + theme(legend.position='none')
    return(PlotOut)
}


PlotSupFig6D = function(df) {
    PlotOut = ggplot(data=df, aes(y=as.numeric(as.character(SummedMutsInBins)), x=as.numeric(as.character(Bin)))) + 
        geom_bar(position = "dodge", stat = "identity") + 
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
        scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
        theme_classic()  +  theme(legend.position = "none") + 
        labs(y="Number of Mutations", x="Total Number of Substitutions") +
        theme(axis.text.x = element_text(size=8), axis.text.y = element_text(size=8), text = element_text(size=8),
             legend.title=element_blank(), strip.text = element_text(size = 8), strip.background = element_blank()) +
        theme(legend.position = "none", text = element_text(size=8))+
        facet_wrap(~ PurityThreshold)
    return(PlotOut)
}

PlotSupFig7A = function(df) {
    #df = subset(df, (df$Group != "All") & (df$xbin != '10000-30000'))
    PlotOut = ggplot(data=df, aes(y = true, x=Bin, fill=Group, color=Group)) + 
            geom_point(size=.8) + geom_line(size=0.65) + # Add line/point
            labs(x="Total Number of Substitutions", y="dN/dS") + # Label axes
            geom_hline(aes(yintercept = 1)) + # Add neutral dN/dS line
            scale_y_continuous(limits= c(0.25,9), breaks = c(0.25,0.5,1,2,4,6), trans = log_trans()) +
            geom_ribbon(aes(ymin=low, ymax=high, fill=Group),  alpha=0.4, linetype=0) + ## Add CI
            theme_classic() + theme(legend.position='none') +
            scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))  +
            scale_fill_manual(breaks = c('All','Passengers','Drivers_Bailey'), values = c("grey",'red',"#7FBD32"), labels=c('All','Passengers','Drivers')) +
            scale_colour_manual(breaks = c('All','Passengers','Drivers_Bailey'), values = c("grey",'red',"#7FBD32"), labels=c('All','Passengers','Drivers'))   +
            facet_wrap(~dNdS)
    return(PlotOut)
}

PlotSupFig7B = function(df) {
    PlotOut = ggplot(data= df, aes(y=true, x=as.numeric(as.character(Bin)), color=Group, fill=Group)) + 
        scale_colour_manual(breaks = c('All','Passengers','Drivers_Bailey'), values = c("grey",'red',"#7FBD32"), 
            labels=c( 'All','Passengers','Drivers')) +  
        geom_line(size=0.65) + scale_y_continuous(limits= c(0.25,6), breaks = c(0.25,0.5,1,2,4,6), trans = log_trans()) +
        theme_classic() + labs(x="Total Number of Substitutions", y="dN/dS") + 
        theme(legend.title = element_blank()) + geom_hline(aes(yintercept = 1)) + ## neutral dN/dS line
        geom_ribbon(aes(ymin=low, ymax=high, fill=Group),  alpha=0.4, linetype=0) + ## adds
        scale_fill_manual(breaks = c('All','Passengers','Drivers_Bailey'), values = c("grey",'red',"#7FBD32"), 
            labels=c('All','Passengers','Drivers')) 
    return(PlotOut)
}



PlotSupFig8 = function(df) {
    dEdIPlot = ggplot(data=df, aes(y = mean, x=as.numeric(as.character(Bin)), fill=Group, color=Group)) + 
            geom_point(size=.8) + geom_line(size=0.65) + # Add line/point
            labs(x="Total Number of CNAs", y="dE/dI") +
            geom_hline(aes(yintercept = 1)) + # Add neutral dE/dI line
            scale_y_continuous(limits= c(0.5,2), breaks = c(0.5,1,2), trans = log_trans()) +
            geom_ribbon(aes(ymin=low, ymax=high, fill=Group),  alpha=0.4, linetype=0) + # Add CI
            facet_wrap(~ metric) 
    dEdIPlot = AddTheme(dEdIPlot) 
    dEdIPlot = AddAxesScales(dEdIPlot)
    ggsave(paste0(PlottingDir, "FigS8_dEdIOfNullPermutations.pdf"), width=6, height=4, units="in")
    return(dEdIPlot)
}

PlotSupFig9 = function(df) {
    dEdIPlot = ggplot(data=df, aes(y = fractionalOverlap_mean, x=as.numeric(as.character(Bin)), fill=Track, color=Track)) + 
            geom_point(size=.8) + geom_line(size=0.65) + # Add line/point
            labs(x="Total Number of CNAs", y="dE/dI") +
            geom_hline(aes(yintercept = 1)) + # Add neutral dE/dI line
            scale_y_continuous(limits= c(0.5,5), breaks = c(0.5,1,2,4), trans = log_trans()) +
            geom_ribbon(aes(ymin=fractionalOverlap_low, ymax=fractionalOverlap_high, fill=Track),  alpha=0.4, linetype=0) + # Add CI
            facet_wrap(~ Length_Category) 
    dEdIPlot = AddColor(dEdIPlot, FALSE)  
    dEdIPlot = AddTheme(dEdIPlot) + ggtitle('Fractional Overlap of CNAs')
    dEdIPlot = AddAxesScales(dEdIPlot) + theme(strip.background = element_blank())
    ggsave(paste0(PlottingDir, "FigS9_dEdIFractionalOverlap.pdf"), width=6, height=4, units="in")
    return(dEdIPlot)
}


PlotSupFig10 = function(df) {
    df = subset(df, df$Group != 'All')
    dNdSPlot = ggplot(data=df, aes(y = true, x=as.numeric(as.character(Bin)), fill=Group, color=Group)) +
        geom_point(size=.8) + geom_hline(aes(yintercept = 1)) + labs(x="Total Number of Substitutions", y="dN/dS") + 
        geom_ribbon(aes(ymin=low, ymax=high, fill=Group),  alpha=0.4, linetype=0) + ## adds CI
        scale_colour_manual(breaks = c('Clonal_Drivers', 'Drivers','Subclonal_Drivers','Clonal_Passengers','Passengers', 'Subclonal_Passengers'),
            values = c("#375F1B", "#7FBD32", "#ABEF7B", "#A70F01","#FE1601", "#FC8C82"),
            labels = c('Clonal Drivers', 'All Drivers','Subclonal Drivers','Clonal Passengers', 'All Passengers', 'Subclonal Passengers')) +
        geom_line(size=0.65) + scale_y_continuous(limits= c(0.1,32), breaks = c(0.1,0.25,0.5,1,2,4,6,8,16,32), trans = log_trans()) + 
        scale_fill_manual(breaks = c('Clonal_Drivers', 'Drivers','Subclonal_Drivers','Clonal_Passengers','Passengers', 'Subclonal_Passengers'),
            values = c("#375F1B", "#7FBD32", "#ABEF7B", "#A70F01","#FE1601", "#FC8C82"),
            labels = c('Clonal Drivers', 'All Drivers','Subclonal Drivers','Clonal Passengers', 'All Passengers', 'Subclonal Passengers')) +
        facet_wrap(~AF) + ggtitle('Subclonal SNVs by AF Threshold') +
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
        theme_classic() + theme(legend.title=element_blank(), plot.title = element_text(hjust = 0.5, size=9, 
        face = "bold"), text = element_text(size=9), legend.position = 'bottom', legend.box = "horizontal") +
        guides(fill = guide_legend(nrow = 1)) + guides(color = guide_legend(nrow = 1))
    ggsave(paste0(PlottingDir, "FigS10_dNdSClonalSubclonalWithDifferentThresholds.pdf"), width=7, height=7, units="in")
    return(dNdSPlot)
}


PlotSupFig11 = function(df) {
    PlotOut = ggplot(data=df, aes(fill=MutationRateGroup, x=Group, y=as.numeric(as.character(true)))) +
        geom_bar(stat='identity', position=position_dodge()) + geom_hline(aes(yintercept = 1), linetype=2) +
        scale_fill_manual(values = c("black","grey"), labels=c('Elevated Mutation Load', 'Low Mutation Load'))  +
        labs(x='', y='dN/dS') + # geom_hline(aes(yintercept = 1), linetype=2) + 
        coord_flip( ylim=c(0,1.5)) + theme_bw() + theme(legend.position  = "bottom", legend.title=element_blank(), legend.text=element_text(size=8))
    ggsave(paste0(PlottingDir, "FigS11_dNdSInFunctionalGeneSets.pdf"), width=5, height=4, units="in")
    return(PlotOut) 

}

PlotSupFig12A = function(df) {
    PlotOut = ggplot(data=df, aes(y=true, x=as.numeric(as.character(Bin)), color=Group)) + 
        geom_point(size=.8) + geom_line(size=0.65) +
        scale_y_continuous(limits= c(0.1,52), breaks = c(0.25,0.5,1,2,4,8,16,32), trans = log_trans()) +
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
        theme_classic() + facet_wrap(~ CancerType) + theme(legend.title=element_blank(), legend.position = 'bottom') + 
        labs(x="Total Number of Substitutions", y="dN/dS") + 
        geom_hline(aes(yintercept = 1)) + ## neutral dN/dS line
        geom_ribbon(aes(ymin=low, ymax=high, fill=Group),  alpha=0.4, linetype=0) + ## Add CI
        ggtitle('Whole-Exome SNVs') 
    PlotOut = AddColor(PlotOut, FALSE)  
    return(PlotOut)
}

PlotSupFig12B = function(df) {
    df$CancerTypeWithN = paste0(df$CancerType, '(n=', df$TotalPatientsInAllBins, ')')
    df$ID = 1
    df$MutationRateGroup = gsub("Low Mutation Rate", 'Low Burden', df$MutationRateGroup)
    df$MutationRateGroup = gsub("High Mutation Rate", 'High Burden', df$MutationRateGroup)
    df$true = format(round(df$true, 2), nsmall = 2)
    PlotOut = ggplot(data=df, aes(y=CancerTypeWithN, x=ID, fill=log10(as.numeric(as.character(true))))) +
        geom_tile() + geom_text(aes(label = true), size=3, parse=TRUE) + #scale_fill_manual(values = c("white")) +
        scale_fill_gradient2(low = '#E92E31', mid='white', high='#5CD843', midpoint = 0 , name='dN/dS') + 
        facet_wrap(Group ~ MutationRateGroup, ncol = 4) + theme_classic() +
        theme(axis.text.y = element_text(size=10),strip.text.x = element_text(size = 10),
            strip.text = element_text(size = 10), strip.background = element_blank(),  axis.ticks.x=element_blank(),
            axis.text.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
    return(PlotOut)
}

PlotSupFig13A = function(df) {
     PlotOut = ggplot(data=df, aes(y=fractionalOverlap_mean, x=as.numeric(as.character(Bin)), color=Group)) + 
        geom_point(size=.8) + geom_line(size=0.65) +
         scale_y_continuous(limits= c(0.4, 15), breaks = c(0.5,0.75,1,2,4,8,12), trans = log_trans()) +
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
        theme_classic() + facet_wrap(~ type) + theme(legend.title=element_blank(), legend.position = 'bottom') + 
        labs(x="Total Number of Substitutions", y="dE/dI") + 
        geom_hline(aes(yintercept = 1)) + ## neutral dN/dS line
        geom_ribbon(aes(ymin=fractionalOverlap_low, ymax=fractionalOverlap_high, fill=Group),  alpha=0.4, linetype=0) + ## Add CI
        ggtitle('Fractional Overlap') 
    PlotOut = AddColor(PlotOut, FALSE)  + theme(legend.position='none',strip.background = element_blank())
    return(PlotOut)
}

PlotSupFig13B = function(df) {
     PlotOut = ggplot(data=df, aes(y=fractionalOverlap_mean, x=as.numeric(as.character(Bin)), color=Group)) + 
        geom_point(size=.8) + geom_line(size=0.65) +
         scale_y_continuous(limits= c(0.4, 15), breaks = c(0.5,1,2,4,8), trans = log_trans()) +
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
        theme_classic() + facet_wrap(~  broad_type) + theme(legend.title=element_blank(), legend.position = 'bottom') + 
        labs(x="Total Number of Substitutions", y="dE/dI") + 
        geom_hline(aes(yintercept = 1)) + ## neutral dN/dS line
        geom_ribbon(aes(ymin=fractionalOverlap_low, ymax=fractionalOverlap_high, fill=Group),  alpha=0.4, linetype=0)  ## Add CI
    PlotOut = AddColor(PlotOut, FALSE)   + theme(strip.background = element_blank())
    return(PlotOut)
}

PlotSupFig13C = function(df) {
     PlotOut = ggplot(data=df, aes(y=breakpointFreq_mean, x=as.numeric(as.character(Bin)), color=Group)) + 
        geom_point(size=.8) + geom_line(size=0.65) +
        scale_y_continuous(limits= c(0.4, 15), breaks = c(0.5,0.75,1,2,4,8,12), trans = log_trans()) +
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
        theme_classic() + facet_wrap(~ type) + theme(legend.title=element_blank(), legend.position = 'bottom') + 
        labs(x="Total Number of Substitutions", y="dE/dI") + 
        geom_hline(aes(yintercept = 1)) + ## neutral dN/dS line
        geom_ribbon(aes(ymin=breakpointFreq_low, ymax=breakpointFreq_high, fill=Group),  alpha=0.4, linetype=0) + ## Add CI
        ggtitle('Breakpoint Frequency') + theme(legend.position='none', strip.background = element_blank())
    PlotOut = AddColor(PlotOut, FALSE)  
    return(PlotOut)
}

PlotSupFig13D = function(df) {
     PlotOut = ggplot(data=df, aes(y=breakpointFreq_mean, x=as.numeric(as.character(Bin)), color=Group)) + 
        geom_point(size=.8) + geom_line(size=0.65) +
        scale_y_continuous(limits= c(0.4, 15), breaks = c(0.5,1,2,4,8), trans = log_trans()) +
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
        theme_classic() + facet_wrap(~  broad_type) + theme(legend.title=element_blank(), legend.position = 'bottom') + 
        labs(x="Total Number of Substitutions", y="dE/dI") + 
        geom_hline(aes(yintercept = 1)) + ## neutral dN/dS line
        geom_ribbon(aes(ymin=breakpointFreq_low, ymax=breakpointFreq_high, fill=Group),  alpha=0.4, linetype=0)  ## Add CI
    PlotOut = AddColor(PlotOut, FALSE)  + theme(strip.background = element_blank())
    return(PlotOut)
}

PlotSupFig14A = function(df) {
    PlotOut = ggplot(data=df, aes( x=as.numeric(as.character(Bin)), y= as.numeric(as.character(true)), fill=Group)) + 
        geom_point(size=0.8) + geom_line(size=0.65) +
        scale_colour_manual(values = c("blue","orange","grey","purple"), 
            breaks= c('Chaperonins','HSP90','All','Proteasome'), labels=c('Chaperonins','HSP90','All','Proteasome')) +  
        labs(x="Total Number of Substitutions", y="Mean Expression of Gene Set\n(Relative to an Average Tumor)") + 
        geom_ribbon(aes(ymin=low, ymax=high, fill=Group),  alpha=0.4, linetype=0) + 
        scale_fill_manual(values = c("blue","orange","grey","purple"), breaks=c('Chaperonins','HSP90','All','Proteasome'), 
            labels=c('Chaperonins','HSP90','All','Proteasome')) +
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10",math_format(10^.x))) +
        theme_classic() 
    return(PlotOut)

}

PlotSupFig14B = function(df) {
    # unweighted = subset(unweighted, unweighted$r_value > 0)
    #print(subset(df, as.character(df$Proteasome == "True"))
    df$r2 = as.numeric(as.character(df$r_value))
    PlotOut = ggplot(data=df, aes(y=r2, x=reorder(GENE_NAME, r2), color='grey')) + 
        geom_bar(position = "dodge", stat = "identity", color='grey') +
        geom_segment(data = subset(df, df$Chaperonins == "True"), aes(x=GENE_NAME, xend=GENE_NAME, y=r2 + 0.1 , yend=r2), size=0.5,color='orange', arrow = arrow(length = unit(0.5, "cm"))) +
        geom_segment(data = subset(df, df$HSP90 == "True"), aes(x=GENE_NAME, xend=GENE_NAME, y=r2 + 0.1 , yend=r2), size=0.5, color='blue', arrow = arrow(length = unit(0.5, "cm"))) +
        geom_segment(data = subset(df, df$Proteasome == "True"), aes(x=GENE_NAME, xend=GENE_NAME, y=r2 + 0.1 , yend=r2), size=0.5 ,color='purple', arrow = arrow(length = unit(0.5, "cm"))) +
        theme_classic()  +
        geom_hline(aes(yintercept = 0.75), color='black', linetype='dashed') + ## global all mutations dN/dS line
        geom_hline(aes(yintercept = 0.5), color='black', linetype='dashed') + ## global all mutations dN/dS line
        geom_hline(aes(yintercept = 0.25), color='black', linetype='dashed') + ## global all mutations dN/dS line
        geom_hline(aes(yintercept = -0.75), color='black', linetype='dashed') + ## global all mutations dN/dS line
        geom_hline(aes(yintercept = -0.5), color='black', linetype='dashed') + ## global all mutations dN/dS line
        geom_hline(aes(yintercept = -0.25), color='black', linetype='dashed') + ## global all mutations dN/dS line
        #theme(legend.position = "none") + 
        #theme(text = element_text(size=12)) 
        labs(y="Correlation (R2) of Gene\nExpression With Mutational Burden", x="Gene") +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
        #theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size=10, face = "bold"), text = element_text(size=10)) 
        scale_color_identity(name = "Gene Sets", guide = "legend",  
                                breaks = c("red", "purple","blue"), labels = c("Chaperonins", "HSP90", "Proteasome"))
        return(PlotOut) 
}

PlotSupFig14C = function(df) {
    PlotOut = ggplot(data=df, aes(x=null)) + 
      geom_histogram() + 
      theme_classic() +
      theme(text = element_text(size=9), axis.text.y = element_text(size=9), strip.text.x = element_text(size = 9),legend.title=element_blank(), strip.text = element_text(size = 9)) +
      geom_vline(xintercept=0.66, color='red') + 
      scale_x_continuous(limits= c(-1, 1), breaks = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1)) +
      labs(x="Median Correlation Coefficient (r)\n of Randomly Sampled Sets of Genes (n=28) ", y="Counts") 
    return(PlotOut)

}

PlotSupFig14D = function(df) {
    df$ID = 1
    df$true = as.numeric(as.character(df$High)) - as.numeric(as.character(df$Low)) 
    df$true = format(round(df$true, 2), nsmall = 2)
    PlotOut = ggplot(data=df, aes(y=CancerType, x=ID, fill=(as.numeric(as.character(true))))) +
        geom_tile() + geom_text(aes(label = true), size=3, parse=TRUE) + #scale_fill_manual(values = c("white")) +
        scale_fill_gradient2(low = '#E92E31', mid='white', high='#5CD843', midpoint = -0.02 , name='') + 
        facet_wrap(~Group, ncol = 4) + theme_classic() +
        theme(axis.text.y = element_text(size=8),strip.text.x = element_text(size = 8),
            strip.text = element_text(size = 8), strip.background = element_blank(),  axis.ticks.x=element_blank(),
            axis.text.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), 
            plot.title = element_text(hjust = 0.5, size=12, face = "bold"), text = element_text(size=9), legend.position = 'bottom') +
        ggtitle('CNAs') + theme(legend.position = 'bottom')
    return(PlotOut)
} 

PlotSupFig14E = function(df) {
    df$ID = 1
    df$true = as.numeric(as.character(df$High)) - as.numeric(as.character(df$Low)) 
    df$true = format(round(df$true, 2), nsmall = 2)
    PlotOut = ggplot(data=df, aes(y=CancerType, x=ID, fill=(as.numeric(as.character(true))))) +
        geom_tile() + geom_text(aes(label = true), size=3, parse=TRUE) + #scale_fill_manual(values = c("white")) +
        scale_fill_gradient2(low = '#E92E31', mid='white', high='#5CD843', midpoint = -0.02 , name='') + 
        facet_wrap(~Group, ncol = 4) + theme_classic() +
        theme(axis.text.y = element_text(size=8),strip.text.x = element_text(size = 8),
            strip.text = element_text(size = 8), strip.background = element_blank(),  axis.ticks.x=element_blank(),
            axis.text.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), 
            plot.title = element_text(hjust = 0.5, size=12, face = "bold"), text = element_text(size=9), legend.position = 'bottom') + 
        ggtitle('SNVs')
    return(PlotOut)
} 




PlotSupFig21 = function(df) {
    dNdSPlot = ggplot(data=df, aes(y = true, x=Bin, fill=Group, color=Group)) + 
            geom_point(size=.8) + geom_line(size=0.65) + # Add line/point
            labs(x="Total Number of Substitutions", y="dN/dS") + # Label axes
            geom_hline(aes(yintercept = 1)) + # Add neutral dN/dS line
            scale_y_continuous(limits= c(0.25,9), breaks = c(0.25,0.5,1,2,4,8), trans = log_trans()) +
            geom_ribbon(aes(ymin=low, ymax=high, fill=Group),  alpha=0.4, linetype=0)
    dNdSPlot = AddTheme(dNdSPlot)
    dNdSPlot = AddColor(dNdSPlot, FALSE)
    dNdSPlot = AddAxesScales(dNdSPlot) + theme(legend.position='none') + theme(strip.background = element_blank()) + facet_wrap(Dataset~dNdS) 
    ggsave(paste0(PlottingDir, "FigS21_TCGA_MutationCalls.pdf"), width=5, height=5, unit='in')
}



PlotSupFig22A = function(df) {
    df = unique(df[c('Group','xbin','NumMutsInBin')])
    df$Bin = data.frame(do.call('rbind', strsplit(as.character(df$xbin),'-',fixed=TRUE)))$X1
    Counts = ggplot(data=df, aes(y=as.numeric(as.character(NumMutsInBin)), x=as.numeric(as.character(Bin)), fill=Group)) + 
        geom_bar(position = "dodge", stat = "identity") + 
        scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
        theme_minimal() + labs(y="Number of Mutations", x="Total Number of Substitutions")
    Counts = AddColor(Counts, FALSE)
    Counts = AddTheme(Counts) + theme(legend.position = 'none')
    Counts = AddAxesScales(Counts) + ggtitle('All SNVs')
    return(Counts)

}


PlotSupFig22B = function(df) {
    df = unique(df[c('Group','xbin','NumMutsInBin')])
    df$Bin = data.frame(do.call('rbind', strsplit(as.character(df$xbin),'-',fixed=TRUE)))$X1
    Counts = ggplot(data=df, aes(y=as.numeric(as.character(NumMutsInBin)), x=as.numeric(as.character(Bin)), fill=Group)) + 
        geom_bar(position = "dodge", stat = "identity") + 
        scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
        theme_minimal() + labs(y="Number of Mutations", x="Total Number of Substitutions") 
    Counts = AddTheme(Counts) + theme(legend.position = 'none')
    Counts = AddColor(Counts, TRUE)
    Counts = AddAxesScales(Counts) + ggtitle('PolyPhen2')
    return(Counts)

}


PlotSupFig22C = function(df) {
    df = unique(df[c('Group','xbin','NumMutsInBin','Length_Category')])
    df$Bin = data.frame(do.call('rbind', strsplit(as.character(df$xbin),'-',fixed=TRUE)))$X1
    Counts = ggplot(data=df, aes(y=as.numeric(as.character(NumMutsInBin)), x=as.numeric(as.character(Bin)), fill=Group)) + 
        geom_bar(position = "dodge", stat = "identity") + 
        scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
        theme_minimal() + labs(y="Number of Mutations", x="Total Number of CNAs") 
    Counts = AddTheme(Counts) + facet_wrap(~Length_Category) + theme(legend.position = 'none')
    Counts = AddColor(Counts, FALSE)
    Counts = AddAxesScales(Counts) + ggtitle('CNAs')
    return(Counts)
}



PlotSupFig22D = function(df) {
    df = unique(df[c('Group','xbin','NumMutsInBin')])
    df$Bin = data.frame(do.call('rbind', strsplit(as.character(df$xbin),'-',fixed=TRUE)))$X1
    Counts = ggplot(data=df, aes(y=as.numeric(as.character(NumMutsInBin)), x=as.numeric(as.character(Bin)), fill=Group)) + 
        geom_bar(position = "dodge", stat = "identity") + 
        scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
        theme_minimal() + labs(y="Number of Mutations", x="Total Number of Substitutions") +
        scale_fill_manual(breaks = c('Clonal_Drivers', 'Drivers','Subclonal_Drivers','Clonal_Passengers','Passengers', 'Subclonal_Passengers'),
                        values = c("#375F1B", "#7FBD32", "#ABEF7B", "#A70F01","#FE1601", "#FC8C82"),
                        labels = c('Clonal Drivers', 'All Drivers','Subclonal Drivers','Clonal Passengers', 'All Passengers', 'Subclonal Passengers')) 
    Counts = AddTheme(Counts) + ggtitle('Clonal vs. Subclonal SNVs') 
    Counts = AddAxesScales(Counts) + theme(legend.position = 'bottom', legend.box = "horizontal") + guides(fill = guide_legend(nrow = 1)) 
    return(Counts)

}


PlotReviewerResponse_1 = function(df) {
dNdSPlot = ggplot(data=df, aes(y = true, x=Bin, fill=Group, color=Group)) + 
            geom_point(size=.8) + geom_line(size=0.65) + # Add line/point
            labs(x="Total Number of Substitutions", y="dN/dS") + # Label axes
            geom_hline(aes(yintercept = 1)) + # Add neutral dN/dS line
            scale_y_continuous(limits= c(0.25,9), breaks = c(0.25,0.5,1,2,4,8), trans = log_trans()) +
            geom_ribbon(aes(ymin=low, ymax=high, fill=Group),  alpha=0.4, linetype=0) + ## Add CI +
            ggtitle('TCGA and ICGC')
    dNdSPlot = AddTheme(dNdSPlot)
    dNdSPlot = AddColor(dNdSPlot, FALSE)
    dNdSPlot = AddAxesScales(dNdSPlot) + theme(legend.position='none') +   theme(strip.background = element_blank()) + facet_wrap(~dNdS) 
    ggsave(paste0(PlottingDir, "ReviewerResponse_TCGAandICGC_PostFiltering.pdf"), width=5, height=4, unit='in')
}



PlotReviewerResponse_3 = function(df) {
    dNdSPlot = ggplot(data=df, aes(y = true, x=Bin, fill=Group, color=Group)) + 
            geom_point(size=.8) + geom_line(size=0.65) + # Add line/point
            labs(x="Total Number of Mutations", y="dN/dS") + # Label axes
            geom_hline(aes(yintercept = 1)) + # Add neutral dN/dS line
            scale_y_continuous(limits= c(0.25,12), breaks = c(0.25,0.5,1,2,4,8), trans = log_trans()) +
            geom_ribbon(aes(ymin=low, ymax=high, fill=Group),  alpha=0.4, linetype=0) 
    dNdSPlot = AddTheme(dNdSPlot)
    dNdSPlot = AddColor(dNdSPlot, FALSE)
    dNdSPlot = AddAxesScales(dNdSPlot) + theme(legend.position='none') + theme(strip.background = element_blank())  
    ggsave(paste0(PlottingDir, "ReviewerResponse_TCGAandTMB.pdf"), width=4, height=4, unit='in')
}

PlotCombinedPanels = function(FigNum, ...) {
    OutDirectory = "/labs/ccurtis2/tilk/scripts/hri/Figures/"
    AllDataFiles = list(...) # List of all dataframes 
    if (FigNum == 'Fig2') {   
        Out = ggdraw() + 
            draw_plot(PlotFig2A(AllDataFiles[[1]]), x = 0,    y = 0.49, width = 0.3, height = 0.49) + 
            draw_plot(PlotFig2B(AllDataFiles[[2]]), x = 0.3, y = 0.49, width = 0.165, height = 0.49) +
            draw_plot(PlotFig2C(AllDataFiles[[3]]) + theme(legend.position = "none"), x= 0.46, y = 0.49, width = 0.27, height = 0.50) +
            draw_plot(PlotFig2D(AllDataFiles[[4]]) + theme(legend.position = "none"), x = 0.72, y = 0.49, width = 0.28, height = 0.5) +
            draw_plot(get_legend(GetFigure2TopLegend(AllDataFiles[[4]])), x= 0.38, y = -0.04, width = 0.28, height = 1) +
            draw_plot(PlotFig2E(AllDataFiles[[5]]) + theme(legend.position = "none"), x = 0, y = 0.05, width = 0.3, height = 0.39) +
            draw_plot(PlotFig2F1(AllDataFiles[[6]]) + theme(legend.position = "none"),  x = 0.3, y = 0.05,   width = 0.22, height = 0.39) + 
            draw_plot(PlotFig2F2(AllDataFiles[[7]]) + theme(legend.position = "none"), x = 0.52, y = 0.05,   width = 0.23, height = 0.39) +
            draw_plot(PlotFig2G(AllDataFiles[[8]]) + theme(legend.position = "none"),  x = 0.756, y = 0.05,   width = 0.249, height = 0.39) +
            draw_plot(get_legend(PlotFig2G(AllDataFiles[[8]])), x = 0.68, y = -0.475, width = 0.249, height = 1) +
            draw_plot(get_legend(PlotFig2E(AllDataFiles[[5]])),  x = 0.3, y = -0.475, width = 0.249, height = 1) +
            draw_plot_label(label = c("A", "B", "C","D","E","F","G"), size = 12, x = c(0, 0.25, 0.5, 0.75,0,0.3,0.75), y = c(1,1,1,1,0.435,0.435,0.435))
        ggsave(plot=Out, paste0(OutDirectory, "Fig2.pdf"), width=10, height=8, units="in")
    } else if (FigNum == 'FigS2') {
        Out = ggdraw() +
            draw_plot(PlotSupFig2A(AllDataFiles[[1]]), x = 0, y = 0.3, width = 1, height = 0.7) + 
            draw_plot(PlotSupFig2B(AllDataFiles[[2]]), x= 0 , y= 0, width= 0.5 , height = 0.3) +
            draw_plot(PlotSupFig2C(AllDataFiles[[3]]), x= 0.5 , y= 0, width= 0.5, height = 0.3) +
            draw_plot_label(label = c("A","B","C"), size = 12, x = c(0, 0, 0.5), y = c(1,0.3,0.3)) 
        ggsave(plot=Out, paste0(OutDirectory, "FigS2_dNdSCorrectsForBiases.pdf"))
    } else if (FigNum == 'FigS4') {
        Out = plot_grid( PlotSupFig4(AllDataFiles[[1]], 1) +  theme(legend.position = 'none') , PlotSupFig4(AllDataFiles[[2]], 0.05), 
            PlotSupFig4(AllDataFiles[[3]], 0.01), PlotSupFig4(AllDataFiles[[4]], 0.005), labels = c('A', 'B', 'C', 'D'))
        Legend = get_legend(PlotSupFig4(AllDataFiles[[1]], 1) + theme(legend.position="bottom"))
        Combined = plot_grid( Out, Legend, ncol = 1, rel_heights = c(0.9, 0.08))
        ggsave(plot=Combined, paste0(OutDirectory, "FigS4_OverlapWithCommonPolymorphisms.pdf"), width=8, height=8, units='in')
    } else if (FigNum == 'FigS6') {    
        Panels = ggdraw() +
            draw_plot(PlotSupFig6A(AllDataFiles[[1]]), x = 0, y = 0.7, width = .5, height = .3) + 
            draw_plot(PlotSupFig6B(AllDataFiles[[2]]), x = 0.5, y = 0.7, width = 0.5, height = 0.3) + 
            draw_plot(PlotSupFig6C(AllDataFiles[[3]]), x = 0, y = 0, width = 0.5, height = 0.7) + 
            draw_plot(PlotSupFig6D(AllDataFiles[[4]]), x = 0.5, y = 0, width = 0.5, height = 0.7) + 
            draw_plot_label(label = c("A", "B", "C","D"), size = 12, x = c(0, 0.5, 0, 0.5), y = c(1, 1, 0.7,0.7))
        ggsave(plot=Panels, paste0(OutDirectory, "FigS6_dNdSPatternsPersistAcrossPurityThresholds.pdf"))
    } else if (FigNum == 'FigS7') {
        Out = ggdraw() +
                draw_plot(PlotSupFig7A(AllDataFiles[[1]]), x = 0, y = 0, width = .5, height = 1 ) +
                draw_plot(PlotSupFig7B(AllDataFiles[[2]]), x = 0.5, y = 0, width = .5, height = 1 ) + 
                draw_plot_label(label = c("A", "B"), size = 12, x = c(0, 0.5), y = c(1, 1))
        ggsave(plot=Out, paste0(OutDirectory, "FigS7_ComparisonOfdNdSToMartincorenaEtAl.pdf"), width=7, height=4, units="in")
    } else if (FigNum == 'FigS12') {
        Out = ggdraw() +
            draw_plot(PlotSupFig12A(AllDataFiles[[1]]), x = 0, y = 0, width = .45, height = 1) + 
            draw_plot(PlotSupFig12B(AllDataFiles[[2]]), x = 0.45, y = 0, width = 0.55, height = 1) +
            draw_plot_label(label = c("A", "B"), size = 12, x = c(0, 0.45 ), y = c(1, 1))
        ggsave(plot=Out, paste0(OutDirectory, "FigS12_dNdSAcrossCancerTypes.pdf"), width=11, height=7, units="in")
    } else if (FigNum == 'FigS13') {
        Out = plot_grid(PlotSupFig13A(AllDataFiles[[1]]),
                        PlotSupFig13C(AllDataFiles[[3]]), 
                        PlotSupFig13B(AllDataFiles[[2]]),
                        PlotSupFig13D(AllDataFiles[[4]]), labels = c('A', 'B', 'C', 'D'))
        ggsave(plot=Out, paste0(OutDirectory, "FigS13_dEdIAcrossCancerGroups.pdf"), width=8, height=8, units='in')
    } else if (FigNum == 'FigS14') {
        Out = ggdraw() +
            draw_plot(PlotSupFig14A(AllDataFiles[[1]]), x = 0, y = 0.64, width = .35, height = 0.35) + 
            draw_plot(PlotSupFig14B(AllDataFiles[[2]]), x= 0.35 , y= 0.64, width= 0.32 , height = 0.35) +
            draw_plot(PlotSupFig14C(AllDataFiles[[3]]), x= 0.68 , y= 0.64, width= 0.3, height = 0.35) +
            draw_plot(PlotSupFig14D(AllDataFiles[[4]]), x = 0, y = 0, width = 0.5, height = 0.625) + 
            draw_plot(PlotSupFig14E(AllDataFiles[[5]]), x = 0.5, y = 0, width = 0.5, height = 0.625) +
            draw_plot_label(label = c("A","B","C", "D","E"), size = 12, x = c(0,0.34,0.68, 0, 0.5), y = c(1,1,1,0.625,0.625))  
        ggsave(plot=Out, paste0(OutDirectory, "FigS14_UpregulationOfHeatshockPathways.pdf"), width=10, height=8, units="in")
    } else if (FigNum == 'FigS22') {
        Out = plot_grid(PlotSupFig22A(AllDataFiles[[1]]), 
                PlotSupFig22B(AllDataFiles[[2]]),
                PlotSupFig22C(AllDataFiles[[3]]),
                PlotSupFig22D(AllDataFiles[[4]]) + theme(legend.position = 'none'), ncol=2, labels=c('A','B','C','D'))
        Legend = get_legend(PlotSupFig22D(AllDataFiles[[4]]))
        Combined = plot_grid( Out, Legend, ncol = 1, rel_heights = c(0.9, 0.08))
        ggsave(plot=Combined, paste0(OutDirectory, "FigS22_CountsOfMutationsForFig2.pdf"), width=7, height=7, units="in")
        print(paste0('Plot saved to: ', OutDirectory))
    }
}

