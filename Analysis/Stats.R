
library(dplyr)

GetStats = function(Category) {
    DataDir = "/labs/ccurtis2/tilk/scripts/hri/Analysis/DataTables/"
    if (Category == 'FigS2B_Pvalue') {
        Samp = read.table(paste0(DataDir, "S2B_dNdSSubsampledFromHighMutationRateTCGA1000Reps"), sep=',', header=T)
        AvgdNdSForHighTMBInSameProportions = mean(subset(Samp, Samp$type == '1-3')$X0)
        print(paste0(c('Mean dN/dS for randomly sampled mutations from high TMB tumors in the same proportions as lowest bin Fig2A is: '), 
          as.character(AvgdNdSForHighTMBInSameProportions)))
        Obs = read.table(paste0(DataDir, "dNdSDriverAndPassengersTCGA"), sep='\t', header=T)
        ObsdNdSInLowestBin = unique(subset(Obs, (Obs$Group == 'Passengers') & (Obs$xbin == '1-3'))$true)
        print(paste0(c('Observed dN/dS for Fig2A is '), ObsdNdSInLowestBin))
        NumSimulationWithSimilarLevelsOfSelection = nrow(subset(Samp, as.numeric(as.character(Samp$X0)) < as.numeric(ObsdNdSInLowestBin)))
        print('Using a t-test, the observed mean in the lowest bin is significantly different from randomly sampled mutations in the same proportions.')
        print(t.test(subset(Samp, Samp$type == '1-3')$X0, mu = as.numeric(ObsdNdSInLowestBin)))
        print('Number of randomly sampled sets of mutations in lowest bin w/ similar patterns of selection')
        nrow(subset(Samp, (Samp$type == '1-3') & (Samp$X0 < 0.57)))/nrow(subset(Samp, (Samp$type == '1-3'))) * 100
    } else if (Category == 'Fig2C_Pvalue') {
        # Obs = na.omit(read.table(paste0(DataDir, "dEdITCGA"), sep='\t', header=T))
        # Obs=subset(Obs, Obs$Track == 'Drivers')
        # t.test(subset(Obs, Obs$Length_Category == "<100Kb")$breakpointFreq_mean, 
        #         subset(Obs, Obs$Length_Category == ">100Kb")$breakpointFreq_mean, paired=TRUE)
        # print(wilcox.test(as.numeric(subset(Obs, Obs$Length_Category == "<100Kb")$breakpointFreq_mean), 
        #                 as.numeric(subset(Obs, Obs$Length_Category == ">100Kb")$breakpointFreq_mean), paired = TRUE)$p.value)     
    } else if (Category == 'Fig2G_CNV_Weighted_R2') {
        df = read.table(paste0(DataDir, "CNVExpressionTCGA"), sep='\t', header=T)
        df$Bin = as.numeric(as.character(data.frame(do.call('rbind', strsplit(as.character(df$xbin),'-',fixed=TRUE)))$X1)) # Split bin column to get x axis values
        print('Weighted R2 of HSP90 gene set by mutational burden CNV bins is: ')
        print(summary(lm(true ~ log10(Bin), weights=(NumPatientsInBin), data=subset(df, df$geneGroup == "HSP90")))$r.squared)
        print('Weighted R2 of Chaperonins gene set by mutational burden CNV bins is: ')       
        summary(lm(true ~ log10(Bin), weights=(NumPatientsInBin), data=subset(df, df$geneGroup == "Chaperonins")))$r.squared
        print('Weighted R2 of Proteasome gene set by mutational burden CNV bins is: ')
        summary(lm(true ~ log10(Bin), weights=(NumPatientsInBin), data=subset(df, df$geneGroup == "Proteasome")))$r.squared 
    } else if (Category == 'Fig2G_SNV_Weighted_R2') {
        df = read.table(paste0(DataDir, "SNVExpressionTCGA"), sep='\t', header=T)
        df$Bin = as.numeric(as.character(data.frame(do.call('rbind', strsplit(as.character(df$xbin),'-',fixed=TRUE)))$X1)) # Split bin column to get x axis values
        print('Weighted R2 of HSP90 gene set by mutational burden SNV bins is: ')
        print(summary(lm(true ~ log10(Bin), weights=(NumPatientsInBin), data=subset(df, df$geneGroup == "HSP90")))$r.squared)
        print('Weighted R2 of Chaperonins gene set by mutational burden SNV bins is: ')       
        summary(lm(true ~ log10(Bin), weights=(NumPatientsInBin), data=subset(df, df$geneGroup == "Chaperonins")))$r.squared
        print('Weighted R2 of Proteasome gene set by mutational burden SNV bins is: ')
        summary(lm(true ~ log10(Bin), weights=(NumPatientsInBin), data=subset(df, df$geneGroup == "Proteasome")))$r.squared 
    } else if (Categ == 'FigS13_Pvalue') { # Wilcox rank test for gene sets
        df = na.omit(read.table(paste0(DataDir, "dNdSFunctionalGeneSetsTCGA"), sep='\t', header=T))
        df = subset(df, df$true != 'Inf')
        LowBins = data.frame(df %>% group_by(Group)  %>%  summarise(Bin = min(Bin)) %>% mutate(BinGroup='low'))
        HighBins = data.frame(df %>% group_by(Group)  %>%  summarise(Bin = max(Bin)) %>% mutate(BinGroup='high'))
        combined = rbind(merge(LowBins, df[c('xbin','true','Bin','Group')], by=c('Bin', 'Group')),
                      merge(HighBins, df[c('xbin','true','Bin','Group')], by=c('Bin', 'Group')))
        print('A Wilcoxon signed rank test between all gene groups in low vs high mutational burden groups yields a p-value of: ')
        print(wilcox.test(as.numeric(subset(combined, combined$BinGroup == "low")$true), 
                    as.numeric(subset(combined, combined$BinGroup == "high")$true), paired = TRUE)$p.value)
    } else if (Categ == 'FigS14_Pvalue') {
        df = na.omit(read.table(paste0(DataDir, "dNdSCancerTypeSpecificDriverAndPassengersTCGA"), sep='\t', header=T))
        df = subset(df, (df$true != 'Inf') & (df$true != 0))
        LowBins = data.frame(df %>% group_by(Group, Subtype)  %>%  summarise(Bin = min(Bin)) %>% mutate(BinGroup='low'))
        HighBins = data.frame(df %>% group_by(Group, Subtype)  %>%  summarise(Bin = max(Bin)) %>% mutate(BinGroup='high'))
        combined = rbind(merge(LowBins, df[c('xbin','true','Bin','Group','Subtype')], by=c('Bin', 'Group','Subtype')),
                      merge(HighBins, df[c('xbin','true','Bin','Group','Subtype')], by=c('Bin', 'Group', 'Subtype')))
         print('A Wilcoxon signed rank test between all cancer types for passengers in low vs high mutational burden groups yields a p-value of: ')
         print(wilcox.test(as.numeric(subset(combined, (combined$BinGroup == "low") & (combined$Group == 'Passengers'))$true), 
                           as.numeric(subset(combined, (combined$BinGroup == "high") & (combined$Group == 'Passengers'))$true), paired = TRUE)$p.value)
         print('A Wilcoxon signed rank test between all cancer types for drivers in low vs high mutational burden groups yields a p-value of: ')
         print(wilcox.test(as.numeric(subset(combined, (combined$BinGroup == "low") & (combined$Group == 'Drivers'))$true), 
                           as.numeric(subset(combined, (combined$BinGroup == "high") & (combined$Group == 'Drivers'))$true), paired = TRUE)$p.value)
    } else if (Categ == 'Fig2F2_Pvalue') {
        df = na.omit(read.table(paste0(DataDir, "dEdICancerTypeSpecificDriverAndPassengersTCGA"), sep='\t', header=T))
        df = df[df$type %in% c('UCEC','PRAD','LUAD','LGG','HSNC','BRCA'),]
        df$Bin = as.numeric(as.character(do.call(rbind,strsplit(as.character(df$xbin), "-"))[,1]))
        df$true = df$breakpointFreq_mean
        df$high = df$breakpointFreq_high
        df$low = df$breakpointFreq_low
        df = subset(df, (df$true != 'Inf') & (df$true != 0) & (df$high != 0) & (df$low != 0 ))
        LowBins = data.frame(df %>% group_by(Track, type)  %>%  summarise(Bin = min(Bin)) %>% mutate(BinGroup='low'))
        HighBins = data.frame(df %>% group_by(Track, type)  %>%  summarise(Bin = max(Bin)) %>% mutate(BinGroup='high'))
        combined = rbind(merge(LowBins, df[c('xbin','true','Bin','Track','type')], by=c('Bin', 'Track','type')),
                      merge(HighBins, df[c('xbin','true','Bin','Track','type')], by=c('Bin', 'Track', 'type')))
        print('A two-sided t-test between all cancer types whether dE/dI in low mutational burden bins is different from 1')
        t.test(as.numeric(subset(combined, (combined$BinGroup == "low") & (combined$Track == 'Passengers'))$true), mu = 1, alternative = "two.sided")
        print('A two-sided t-test between all cancer types whether dE/dI in low mutational burden bins is different from 1')
        t.test(as.numeric(subset(combined, (combined$BinGroup == "low") & (combined$Track == 'Drivers'))$true), mu = 1, alternative = "two.sided")
    } else if (Categ == 'FigS22B-C') {
        genes = read.table('/labs/ccurtis2/tilk/scripts/hri-data/Annotations/GeneSets/ProteinExpressionGeneSetUsed.txt', sep='\t', header=T)
        df = na.omit(read.table(paste0(DataDir, "S22B_R2MutationalBinByExpressionForAllGenes"), sep=',', header=T))
        df$Proteasome = df$GENE_NAME %in% subset(genes, genes$Group == 'Proteasome')$GeneName
        df$HSP90 = df$GENE_NAME %in% subset(genes, genes$Group == 'HSP90')$GeneName
        df$Chaperonins = df$GENE_NAME %in% subset(genes, genes$Group == 'Chaperonins')$GeneName
        ### Calculates number of genes in each gene set
        chaperonins_num= nrow(subset(df, df$Chaperonins == TRUE))
        proteasome_num = nrow(subset(df, df$Proteasome == TRUE))
        hsp90_num = nrow(subset(df, df$HSP90 == TRUE))
        CombinedNumGenes = chaperonins_num  + proteasome_num + hsp90_num 
        ### Calculates the median R2 of our chosen gene set
        MedianForGeneSet = median(c(subset(df, df$Chaperonins == TRUE)$r_value,
                                    subset(df, df$HSP90 == TRUE)$r_value,
                                    subset(df, df$Proteasome == TRUE)$r_value))
        ### Median r2 values for all genes is 0.66
        print(paste0('Median r2 values for all genes is: ', MedianForGeneSet))
    } else if (Categ == 'FigS22D') {
        df = read.table(paste0(DataDir, "S22D_CNVExpressionByCancerTypeTCGA"), sep='\t', header=T)
        df$Bin = as.numeric(as.character(data.frame(do.call('rbind', strsplit(as.character(df$xbin),'-',fixed=TRUE)))$X1)) # Split bin column to get x axis values
        LowBins = data.frame(df %>% group_by(geneGroup, type)  %>%  summarise(Bin = min(Bin)) %>% mutate(BinGroup='low'))
        HighBins = data.frame(df %>% group_by(geneGroup, type)  %>%  summarise(Bin = max(Bin)) %>% mutate(BinGroup='high'))
        combined = rbind(merge(LowBins, df[c('xbin','true','Bin','geneGroup','type')], by=c('Bin', 'geneGroup','type')),
                      merge(HighBins, df[c('xbin','true','Bin','geneGroup','type')], by=c('Bin', 'geneGroup', 'type')))
        print('A Wilcoxon signed rank test between HSP90 in low vs high mutational burden groups yields a p-value of: ')
        print(wilcox.test(
            as.numeric(combined[(combined$BinGroup == 'low') & (combined$geneGroup == 'HSP90'),]$true),
            as.numeric(combined[(combined$BinGroup == 'high') & (combined$geneGroup == 'HSP90'),]$true), paired = FALSE)$p.value)
        print('A Wilcoxon signed rank test between Chaperonins in low vs high mutational burden groups yields a p-value of: ') 
        print(wilcox.test(
            as.numeric(combined[(combined$BinGroup == 'low') & (combined$geneGroup == 'Chaperonins'),]$true),
            as.numeric(combined[(combined$BinGroup == 'high') & (combined$geneGroup == 'Chaperonins'),]$true), paired = FALSE)$p.value)  
        print('A Wilcoxon signed rank test between Proteasome in low vs high mutational burden groups yields a p-value of: ') 
        print(wilcox.test(
            as.numeric(combined[(combined$BinGroup == 'low') & (combined$geneGroup == 'Proteasome'),]$true),
            as.numeric(combined[(combined$BinGroup == 'high') & (combined$geneGroup == 'Proteasome'),]$true), paired = FALSE)$p.value)

    } else if (Categ == 'FigS22E') {
        df = read.table(paste0(DataDir, "S22E_SNVExpressionByCancerTypeTCGA"), sep='\t', header=T)
        df$Bin = as.numeric(as.character(data.frame(do.call('rbind', strsplit(as.character(df$xbin),'-',fixed=TRUE)))$X1)) # Split bin column to get x axis values
        LowBins = data.frame(df %>% group_by(geneGroup, type)  %>%  summarise(Bin = min(Bin)) %>% mutate(BinGroup='low'))
        HighBins = data.frame(df %>% group_by(geneGroup, type)  %>%  summarise(Bin = max(Bin)) %>% mutate(BinGroup='high'))
        combined = rbind(merge(LowBins, df[c('xbin','true','Bin','geneGroup','type')], by=c('Bin', 'geneGroup','type')),
                      merge(HighBins, df[c('xbin','true','Bin','geneGroup','type')], by=c('Bin', 'geneGroup', 'type')))
    
        print('A Wilcoxon signed rank test between HSP90 in low vs high mutational burden groups yields a p-value of: ')
        print(wilcox.test(
            as.numeric(combined[(combined$BinGroup == 'low') & (combined$geneGroup == 'HSP90'),]$true),
            as.numeric(combined[(combined$BinGroup == 'high') & (combined$geneGroup == 'HSP90'),]$true), paired = TRUE)$p.value)
        print('A Wilcoxon signed rank test between Chaperonins in low vs high mutational burden groups yields a p-value of: ')
        
        print(wilcox.test(
            as.numeric(combined[(combined$BinGroup == 'low') & (combined$geneGroup == 'Chaperonins'),]$true),
            as.numeric(combined[(combined$BinGroup == 'high') & (combined$geneGroup == 'Chaperonins'),]$true), paired = TRUE)$p.value)
        print('A Wilcoxon signed rank test between Proteasome in low vs high mutational burden groups yields a p-value of: ')
        print(wilcox.test(
            as.numeric(combined[(combined$BinGroup == 'low') & (combined$geneGroup == 'Proteasome'),]$true),
            as.numeric(combined[(combined$BinGroup == 'high') & (combined$geneGroup == 'Proteasome'),]$true), paired = TRUE)$p.value)
        
        
    }
}
