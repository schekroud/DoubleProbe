# # # Double Probe - convert E1 and E2 data into whole csv files
fileList1 <- list.files(path= "/Users/user/Desktop/E1/csv/", pattern = ".csv")
dataFiles1 <- list(NULL)
count <- 1
for(i in fileList1){
  path <- paste0("/Users/user/Desktop/E1/csv/",i)
  dataFiles1[[count]] <- read.csv(path, header=TRUE,as.is=TRUE)
  count <- count + 1
}
df1 <- do.call("rbind",dataFiles1)
df1 <- df1[order(df1$subID),]

fileList2  <- list.files(path="/Users/user/Desktop/E2/csv/", pattern=".csv")
dataFiles2 <- list(NULL)
count <- 1
for(i in fileList2){
  path <- paste0("/Users/user/Desktop/E2/csv/",i)
  dataFiles2[[count]] <- read.csv(path, header=TRUE, as.is=TRUE)
  count <- count + 1
}
df2 <- do.call("rbind", dataFiles2)
df2 <- df2[order(df2$subID),]

fname <- "/Users/user/Desktop/DoubleProbe3/data/E1.csv"
write.csv(df1, file=fname, sep=',', eol='\n', dec='.', col.names=TRUE)

fname <- "/Users/user/Desktop/DoubleProbe3/data/E2.csv"
write.csv(df2, file=fname, sep=',', eol='\n', dec='.', col.names=TRUE)