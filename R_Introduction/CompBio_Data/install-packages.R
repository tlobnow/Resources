install.packages(c("BiocManager", "bookdown"))

pack.intro <- c('BiocStyle', 'Biostrings', 'tidyverse', 'janitor' ,'GEOquery', 'biomaRt', 'cgdsr', 'ggdendro')
pack.intro <- c('ggdendro')
BiocManager::install(pack.intro, update = F)

#pack.amit <- c('GEOquery', 'pheatmap', 'oligo')
#install(pack.amit, update = F)

pack.schulz <- c('grid', 'org.Mm.eg.db', 'mogene10sttranscriptcluster.db', 'pheatmap', 'biomaRt')
install(pack.schulz, update = F)

pack.riemer <- c('tidyverse', 'pheatmap','DESeq2', 'goseq', 'geneLenDataBase', 'GO.db', 'org.Mm.eg.db')
BiocManager::install(pack.riemer, update = F)

#pack.bentele <- c('BSgenome.Ecoli.NCBI.20080805')
#install(pack.bentele, update = F)

pack.brandt <- c("tidyverse", "cowplot", "viridis", "scales", "sp")
BiocManager::install(pack.brandt, update = F)
