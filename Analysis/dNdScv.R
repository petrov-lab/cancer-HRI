# This script contains the set of functions required to run dNdScv.


LoadPackages = function() {
    library('dndscv')
    library('plyr')
}


BuildGRCh38Ref = function() {
    # CDs was downloaded from Ensembl and the raw column names were re-arranged/re-named to fit dndscv's format as follows:
    # 'Gene stable ID', 'Gene name', 'Protein stable ID', 'Chromosome/scaffold name', 'Genomic coding start', 
    # 'Genomic coding end', 'CDS start', 'CDS end', 'CDS Length', 'Strand'
    RefDir=paste0(getwd(), "/Analysis/dNdScvRef/")
    path_cds_table = paste0(RefDir, "GRCH38_CDS_Table.txt")
    path_genome_fasta = paste0(RefDir, "Homo_sapiens.GRCh38.dna.primary_assembly.fa")
    buildref(cdsfile=path_cds_table, genomefile=path_genome_fasta, outfile = paste0(RefDir, "GRCh38_HomoSapiens_refcds.rda"), 
        excludechrs="MT")
}

CalculatedNdScvHg38 = function(df) {
    ### Calculates global dNdScv of all mutations (not dN/dS per gene)
    RefDir=paste0(getwd(), "/Analysis/dNdScvRef/")
    tryCatch(
        expr = {
            GlobaldNdS = dndscv(df, max_coding_muts_per_sample = Inf,  
            max_muts_per_gene_per_sample = Inf, refdb = paste0(RefDir, "GRCh38_HomoSapiens_refcds.rda"))
            GlobaldNdSOut = as.data.frame(GlobaldNdS$globaldnds)
            GlobaldNdSOut = subset(GlobaldNdSOut, GlobaldNdSOut$name == 'wall')
            return(data.frame(xbin = df$xbin, GlobaldNdSOut))
        },
        error = function(e){ 
            print(e)
            return(data.frame(xbin = df$xbin, name= 'wall', mle=NA, cilow=NA, cihigh=NA))
        }
    )
}


ApplydNdScvToHg38 = function(df) {
    LoadPackages()
    dnds=unique(ddply(df, .(xbin), CalculatedNdScvHg38))
    dnds = subset(dnds, dnds$xbin != "")
    return(dnds)
}


CalculatedNdScvHg19 = function(df) {
    ### Calculates global dNdScv of all mutations (not dN/dS per gene)
    tryCatch(
        expr = {
            GlobaldNdS = dndscv(df, max_coding_muts_per_sample = Inf,  
            max_muts_per_gene_per_sample = Inf, refdb = "hg19")
            GlobaldNdSOut = as.data.frame(GlobaldNdS$globaldnds)
            GlobaldNdSOut = subset(GlobaldNdSOut, GlobaldNdSOut$name == 'wall')
            return(data.frame(xbin = df$xbin, GlobaldNdSOut))
        },
        error = function(e){ 
            print(e)
            return(data.frame(xbin = df$xbin, name= 'wall', mle=NA, cilow=NA, cihigh=NA))
        }
    )
}


ApplydNdScvToHg19 = function(df) {
    LoadPackages()
    dnds=unique(ddply(df, .(xbin), CalculatedNdScvHg19))
    dnds = subset(dnds, dnds$xbin != "")
    return(dnds)
}
