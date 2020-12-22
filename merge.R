## library is reverse stranded
 files =list.files(path="star_align",recursive=T,pattern="ReadsPerGene.out.tab",full.names=T)
         for (file in files)
         {
          if(!exists("rc"))
          {rc=read.table(file,header=F,row.names=1,colClasses=c(NA,"NULL","NULL",NA),col.names=c("gene_id","NULL","NULL",file))}
          else if(exists("rc"))
          {
           temp_rc=read.table(file,header=F,row.names=1,colClasses=c(NA,"NULL","NULL",NA),col.names=c("gene_id","NULL","NULL",file))
           rc=cbind(rc,temp_rc)
           rm(temp_rc)
          }
         }
         head(rc)
         ncol(rc)
         sample_names=read.table("sample_names.txt")
         sample_names
         colnames(rc)=sample_names[,1]
         head(rc)
         write.csv(rc,"reads_count.csv")

