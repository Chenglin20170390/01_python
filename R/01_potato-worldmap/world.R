install.packages("reshape")
install.packages("rworldmap")
install.packages("rworldxtra")
install.packages('sp')
require(reshape)
require(rworldmap)
require(rworldxtra)
rice <- read.table("accession_coordinates.v2.txt", head=T) #读入今天要用到的数据（请下载附件数据）


##更改颜色 potato world map 

op <- palette(c('green','yellow','orange','red'))
pdf('potato=world.pdf',width=8,height=7)
mapBubbles(rice,nameZSize="num",nameX="Longitude",nameY="Latitude",mapRegion='world',nameZs =c('outgroup','wild','landrace','candolleanum'),nameZColour="ID",colourPalette=op,symbolSize=0.3,lwd=0.5,barOrient='vert',oceanCol="skyblue",landCol="white",addColourLegend = TRUE, colourLegendPos='bottomleft',legendPos = "bottomright",addLegend = FALSE,xlim=c(-130,-80),ylim=c(-56,66))
dev.off()