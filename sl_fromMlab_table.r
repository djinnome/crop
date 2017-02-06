##jmd
##10.15.12
##sl_fromMlab_table.r

source('farm_header.r')

### non-iso's ###########################################################################################################
##read
sl <- read.delim('matlab_data/nonIsozymeSynthLethalsPwys.tsv', as.is=TRUE)
#add pyr-1;uc-5
sl <- rbind(sl, c('NCU06532.5', 'NCU07334.5', '', 'PWY0-162,PWY-5686', ''))

## create 2 more columns, with symbol if available, else ncu
sl.names <- sl.add.names(sl, ncu.gene.map)
sl <- cbind(sl[,1:2], sl.names, sl[,-(1:2)])

write.table(sl, 'matlab_data/nonIsozymeSynthLethals_table.tsv', row.names=FALSE, quote=FALSE)

### isozyme SLs #########################################################################################################
isos <- read.table('matlab_data/isozymeSLs.txt', header=TRUE)

isos <- isos[!duplicated(isos),]
isos <- isos[order(isos[,1], isos[,2]),]

iso.names <- sl.add.names(isos, ncu.gene.map)
isos <- cbind(Gene1=isos[,1], Gene2=isos[,2], iso.names)

write.table(isos, 'matlab_data/isozymeSynthLethals_table.tsv', row.names=FALSE, quote=FALSE)
