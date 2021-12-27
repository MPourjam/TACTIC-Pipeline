install.packages(c('igraph', 'fastmatch', 'quadprog'))
install.packages('vegan')
install.packages('GUniFrac')
install.packages('ggplot2')
install.packages('devtools')
library(devtools)
install.packages('iNEXT')

###############################
### EXTRA PACKAGES FOR RHEA ###
###############################
# Normalization, Alpha, Taxonomic require no further installations
# for the Beta step
install.packages('ade4')
install.packages("https://cran.r-project.org/src/contrib/Archive/phangorn/phangorn_2.5.5.tar.gz", repos=NULL, type='source', dependencies=T)
install.packages('cluster')
install.packages('fpc')
install.packages('compare')
# for the Serial steps
install.packages('plotrix')
install.packages('PerformanceAnalytics')
install.packages('reshape')
install.packages('gridExtra')
install.packages('grid')
install.packages('ggrepel')
install.packages('gtable')
install.packages('Matrix')
install.packages('cowplot')
install.packages('muStat')
install.packages('https://cran.r-project.org/src/contrib/Archive/latticeExtra/latticeExtra_0.6-28.tar.gz', repos=NULL, type='source', dependencies=T)
# install.packages('https://cran.r-project.org/src/contrib/Archive/Hmisc/Hmisc_4.3-0.tar.gz', repos=NULL, type='source', dependencies=T)
install.packages('Hmisc')

install.packages('corrplot')
