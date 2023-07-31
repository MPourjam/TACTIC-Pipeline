library(vegan) # check if needed
library(GUniFrac) # check if needed
library(ggplot2) # needed for rank abundance plot
library(iNEXT) # needed for rarefaction plot


#-----------------------
# read OTU table and normalize counts
#-----------------------
# Calculate the species richness in a sample
Species.richness <- function(x)
  {
    # Count only the OTUs that are present >0.5 normalized counts (normalization produces real values for counts)
    count=sum(x[x>0.5]^0)
    return(count)
  }

# Calculate the Effective species richness in each individual sample
Eff.Species.richness <- function(x)
  {
    # Count only the OTUs that are present more than the set proportion
    total=sum(x)
    count=sum(x[x/total>0.0025]^0)
    return(count)
  }


# Calculate the Normalized species richness in each individual sample
Norm.Species.richness <- function(x)
  {
    # Count only the OTUs that are present >0.5 normalized counts (normalization produces real values for counts)
    # Given a fixed Normalization reads depth
    total=sum(x)
    count=sum(x[1000*x/total>0.5]^0)
    return(count)
  }


# Calculate the Shannon diversity index
Shannon.entropy <- function(x)
  {
    total=sum(x)
    se=-sum(x[x>0]/total*log(x[x>0]/total))
    return(se)
  }

# Calculate the effective number of species for Shannon
Shannon.effective <- function(x)
  {
    total=sum(x)
    se=round(exp(-sum(x[x>0]/total*log(x[x>0]/total))),digits =2)
    return(se)
  }

# Calculate the Simpson diversity index
Simpson.concentration <- function(x)
  {
    total=sum(x)
    si=sum((x[x>0]/total)^2)
    return(si)
  }

# Calculate the effective number of species for Simpson
Simpson.effective <- function(x)
{
    total=sum(x)
    si=round(1/sum((x[x>0]/total)^2),digits =2)
    return(si)
}


otu_table <- read.table ('ZOTUs-table.final.tab',
                            check.names = FALSE,
                            header = TRUE,
                            dec = ".",
                            sep = "\t",
                            row.names = 1,
                            comment.char = "")

# Delete column with taxonomy information in dataframe
otu_table$taxonomy <- NULL

# Calculate the minimum sum of all columns/samples
min_sum <- min(colSums(otu_table))

# Divide each value by the sum of the sample and multiply by the minimal sample sum
norm_otu_table <- t(min_sum * t(otu_table) / colSums(otu_table))

# Order and transpose OTU-table
my_otu_table <- norm_otu_table
my_otu_table <-data.frame(t(my_otu_table))

richness <- apply(my_otu_table, 1, Species.richness)
eff_richness <- apply(my_otu_table, 1, Eff.Species.richness)
shannon <- apply(my_otu_table, 1, Shannon.entropy)
shannon_effective <- apply(my_otu_table,1,Shannon.effective)
simpson <- apply(my_otu_table,1, Simpson.concentration)
simpson_effective <-apply(my_otu_table, 1, Simpson.effective)


# write to file, so django can read them and add them to the database
richness_line <- paste('richness', richness, sep='\t')
eff_richness_line <- paste('eff_richness', eff_richness, sep='\t')
shannon_line <- paste('shannon', shannon, sep='\t')
shannon_eff_line <- paste('shannon_effective', shannon_effective, sep='\t')
simpson_line <- paste('simpson', simpson, sep='\t')
simpson_eff_line <- paste('simpson_effective', simpson_effective, sep='\t')

output_file <- paste('div_metrics.tab', sep='')
write(richness_line, output_file,append=TRUE)
write(eff_richness_line, output_file,append=TRUE)
write(shannon_line, output_file,append=TRUE)
write(shannon_eff_line, output_file,append=TRUE)
write(simpson_line, output_file,append=TRUE)
write(simpson_eff_line, output_file,append=TRUE)

if (TRUE){
#---------------------------------------
# RANK ABUNDANCE PLOT
#---------------------------------------
order_norm_counts <- sort(my_otu_table, decreasing=TRUE)
abund_data <- data.frame(counts = c(t(order_norm_counts)), positions=1:length(order_norm_counts))

g1 <- ggplot(abund_data, aes(x=positions,y=counts)) + geom_line(color='red') + geom_point()
g2 <- g1 + xlab('Abundance Rank') + ylab('Relative Abundance')
g3 <- g2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave('rank_abundance_plot.png', g3)

#---------------------------------------
# RAREFACTION CURVE PLOT
#---------------------------------------
inext_out <- iNEXT(otu_table[,1], q=c(0), datatype="abundance")
df <- fortify(inext_out, type=1)
df.point <- df[which(df$method=="observed"),]
df.line <- df[which(df$method!="observed"),]
df.line$method <- factor(df.line$method, c("interpolated",  "extrapolated"),c("interpolation",  "extrapolation"))

g4 <- ggplot(df, aes(x=x, y=y, colour=site)) + geom_point(aes(shape=site), size=5, data=df.point) 
g5 <- g4 + geom_line(aes(linetype=method), lwd=1.5, data=df.line)
g6 <- g5 + geom_ribbon(aes(ymin=y.lwr, ymax=y.upr,fill=site, colour=NULL), alpha=0.2) 
g7 <- g6 + labs(x="Number of reads" , y="Number of OTUs") 
g8 <- g7 + theme(legend.position = "bottom", legend.title=element_blank(),text=element_text(size=18),
                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave('rarafaction_curve.png', g8)
}