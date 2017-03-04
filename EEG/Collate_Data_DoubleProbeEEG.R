# # # Double Probe EEG - convert all individual subject data into one csv file

fileList <- list.files(path= "/Users/user/Desktop/Experiments/Nick/DoubleProbe/EEG/data/csv/", pattern = ".csv")
dataFiles <- list(NULL)
count <- 1
for(i in fileList){
  path <- paste0("/Users/user/Desktop/Experiments/Nick/DoubleProbe/EEG/data/csv/",i)
  dataFiles[[count]] <- read.csv(path, header=TRUE,as.is=TRUE)
  count <- count + 1
}
df <- do.call("rbind",dataFiles)
df <- df[order(df$subID),]

fname <- "/Users/user/Desktop/Experiments/Nick/DoubleProbe/EEG/data/allEEG.csv"
write.csv(df, file=fname, sep=',', eol='\n', dec='.', col.names=TRUE)

