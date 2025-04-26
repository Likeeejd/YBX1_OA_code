#### GSE117999 ########
library(stringr)
library(AnnoProbe)
library(GEOquery)
library(limma)
gse="GSE117999"
geoChina(gse)

load("GSE117999_eSet.Rdata")
pd<-pData(gset[[1]])

raw_dir="/Users/zz/my_projects/实验数据/GEO+MR/data/GSE117999/GSE117999_RAW"
raw_datas=paste0(raw_dir,"/",dir(raw_dir))
raw_datas

raw_order <- str_extract(raw_datas, "GSM\\d+")

pd <- pd[match(raw_order, rownames(pd)), ]

if(any(is.na(match(raw_order, rownames(pd))))) {
  warning("Some GSM IDs in raw_datas don't match rownames(pd)")
}

raw_order=str_extract(raw_datas,"GSM\\d*")
pd=pd[match(raw_order,rownames(pd)),]

group_list<-ifelse(stringr::str_detect(pd$title,"HealthyControl"),"arthroscopic partial meniscectomy","osteoarthritis")
group_list<-factor(group_list,levels=c("arthroscopic partial meniscectomy","osteoarthritis"))

sample_file <- read.delim(gzfile(raw_datas[1]), nrows = 10, comment.char = "*")
colnames(sample_file)

#writeLines(readLines(gzfile(raw_datas[1]), n = 20))

sample_data <- read.delim(gzfile(raw_datas[1]), comment.char = "*", nrows = 10)
colnames(sample_data) 

for (file_path in raw_datas) {
  con <- gzfile(file_path, "r")
  lines <- readLines(con)
  close(con)
  if (length(lines) > 9) {
    lines <- lines[-(1:9)]
  } else {
    warning(paste("文件", file_path, "行数不足9行"))
  }
  
  # 写回文件
  con <- gzfile(file_path, "w")
  writeLines(lines, con)
  close(con)
  
  cat("已处理文件:", file_path, "\n")
}

for (file in raw_datas) {
  sample_data <- read.delim(gzfile(file), comment.char = "*", nrows = 5)
  cat("File:", file, "\n")
  print(colnames(sample_data))
  cat("\n")
}

test_file <- read.delim(gzfile(raw_datas[1]), sep = "\t", check.names = FALSE)
head(test_file[, c("rMedianSignal", "rBGMedianSignal", "rIsWellAboveBG")])


library(limma)

data_list <- lapply(raw_datas, function(file) {
  dat <- read.delim(gzfile(file), sep = "\t", check.names = FALSE)
  list(
    R = dat$rMedianSignal,        
    Rb = dat$rBGMedianSignal,      
    well_above = dat$rIsWellAboveBG 
  )
})

x <- new("RGList")
x$R <- do.call(cbind, lapply(data_list, `[[`, "R"))     
x$Rb <- do.call(cbind, lapply(data_list, `[[`, "Rb")) 
x$other <- do.call(cbind, lapply(data_list, `[[`, "well_above")) 
x$genes <- data.frame(ProbeName = test_file$ProbeName)   
#x$other <- list(rIsWellAboveBG = do.call(cbind, lapply(data_list, `[[`, "well_above")))

str(x)  

x$other

#names(x$other) 
library(limma)

x <- backgroundCorrect(x, method = "normexp")

x <- normalizeBetweenArrays(x, method = "quantile")

Control<-y$genes$ControlType==1L;table(Control)

NoSymbol<-is.na(y$genes$GeneName);table(NoSymbol)

IsExpr<-rowSums(y$other$gIsWellAboveBG>0)>=16;table(IsExpr)

Isdup<-duplicated(y$genes$GeneName);table(Isdup)

yfilt<-y[!Control&!NoSymbol&IsExpr&!Isdup,]
dim(yfilt)

exp=yfilt@.Data[[1]]
boxplot(exp)

colnames(exp)=str_extract(colnames(exp),"GSM\\d*")
exp[1:2,1:2]

anno=yfilt$genes
nrow(anno);nrow(exp)
rownames(exp)=rownames(anno)
ids=unique(anno$GeneName)
exp=exp[!duplicated(anno$GeneName),]
rownames(exp)=ids
exp[1:4,1:4]
