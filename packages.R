packages <- c('ggplot2', 'nlme', 'mgcv', 'TSA', 'forecast', 'devtools', 'ggfortify',
 'reshape2', 'stats', 'dlm', 'grid', 'gridExtra','scales')
lapply(packages, require, character.only = TRUE)
install_github('sinhrks/ggfortify')
