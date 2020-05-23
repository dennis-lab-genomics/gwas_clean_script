#This script will examine the numeric data in a file, 
#look for instances when the values are less than the minimum double precision value R can handle (as defined by the .Machine$double.xmin value)
#replace them with that minimum value, and then create a separate file that lists which rows had conversions and the original row values


args<-commandArgs(trailingOnly = TRUE)

target.path<-args[1]
output.path<-args[2]
sep.char<-args[3]
has.header<-as.logical(args[4])
target.cols <-args[5]

print(args)

if (1==2)
{
  #Rscript clean_file_for_R_small_values.R /data/davis_lab/dennisj/biomarkers/data/GWAS_summary_stats/lipids/jointGwasMc_LDL.txt ~/temp/jointGwasMc_LDL_minRVals.txt $'\t' FALSE 9,10
  
  
target.path<-"/data/davis_lab/dennisj/biomarkers/data/GWAS_summary_stats/lipids/jointGwasMc_LDL.txt"
output.path<-"/data/davis_lab/dennisj/biomarkers/data/GWAS_summary_stats/lipids/jointGwasMc_LDL_minRVals.txt"

target.path<-"X:\\dennisj\\biomarkers\\data\\GWAS_summary_stats\\lipids\\jointGwasMc_LDL.txt"
output.path<-"Z:\\temp\\jointGwasMc_LDL_minRVals.txt"

target.path<-"Z:\\temp\\jointGwasMc_LDL.txt"
output.path<-"Z:\\temp\\jointGwasMc_LDL_minRVals.txt"


target.path<-"Z:\\temp\\small_num_test.txt.txt"
output.path<-"Z:\\temp\\small_num_test_minRVals.txt"

target.path<-"~/temp/small_num_test.txt.txt"
output.path<-"~/temp/small_num_test_minRVals.txt"


sep.char<-"\t"
has.header<-FALSE
has.header<-TRUE
target.cols <-"9,10"
target.cols <-"2,3"
}


target.dat<-read.table(target.path, sep=sep.char, header=has.header, stringsAsFactors = FALSE, quote="", colClasses = "character")

orig.target.dat.dim<-dim(target.dat)
#orig.target.dat<-target.dat

target.col.split<-as.numeric(unlist(strsplit(target.cols, ",")))

min.R.val <- .Machine$double.xmin

mismatch.row<-c()

t=1
for (t in 1:length(target.col.split))
{
  this.col<-target.col.split[t]
  print(this.col)
  this.col.name<-colnames(target.dat)[this.col]
  #print("test1")
  target.dat$dup.col<-target.dat[,this.col]
  #print("test2")
  colnames(target.dat)[ncol(target.dat)]<-paste0(this.col.name, "_orig")
  
  target.dat[,this.col]<-as.numeric(as.character(target.dat[,this.col]))
  
  target.dat.na<-target.dat[is.na(target.dat[,this.col]),]
  
  target.dat.no.na<-target.dat[!is.na(target.dat[,this.col]),]
  
  target.dat.no.na.order<-target.dat.no.na[order(target.dat.no.na[,this.col]),]
  #head(target.dat.no.na.order)
  
  #now check to see if any values got converted to 0, and
  mismatch.target.dat<-target.dat.no.na.order[as.character(target.dat.no.na.order[,this.col])!=as.character(target.dat.no.na.order[,paste0(this.col.name, "_orig")]),]
  mismatch.target.dat.0<-mismatch.target.dat[mismatch.target.dat[,this.col]==0,]
  #check to make sure the original value was 0, or 0.000e+00 or something like that
  real.zero.rows<-grep("^0?\\.?0+e?\\+?0?0?$", mismatch.target.dat.0[,paste0(this.col.name, "_orig")])
  if (length(real.zero.rows)>0)
  {
    mismatch.target.dat.0<-mismatch.target.dat.0[-real.zero.rows,]
  }
  
  #now let's set the values of the probem 0 rows to the R minimum value  
  target.dat[as.numeric(rownames(mismatch.target.dat.0)), this.col]<-min.R.val
   
  #now that we have our problem rows identified, save them
  mismatch.row<-append(mismatch.row, as.numeric(rownames(mismatch.target.dat.0)))
}


#now print out the substituted file and the file showing the rows that original data in the substited file
write.table(target.dat[,1:orig.target.dat.dim[2]], sep=sep.char, file=output.path, col.names = has.header, row.names=FALSE, quote=FALSE)

#now print out a file show the rows that were substituted and row names
output.path.dir<-dirname(output.path)
output.base<-basename(output.path)
output.base.parts<-unlist(strsplit(output.base, "\\.")) 
output.subrows.base<-paste(paste(output.base.parts[1:length(output.base.parts)-1], collapse="."), "substituted_rows", output.base.parts[length(output.base.parts)],  sep=".")
write.table(target.dat[unique(mismatch.row),], sep=sep.char, file=file.path(output.path.dir,output.subrows.base), col.names = has.header, row.names=TRUE, quote=FALSE)
