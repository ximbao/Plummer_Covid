# Covid RNAseq - human - full depth sequencing 
files <- dir("/home/ximbao/Documents/Covid/RNA_FullDepth/", pattern = "fastq")
path = "/home/ximbao/Documents/Covid/RNA_FullDepth/"
path.out = "/home/ximbao/Documents/Covid/RNA_FullDepth/aligned/"
for(i in 1:length(files)) {
  paste0("~/STAR/bin/Linux_x86_64/./STAR --genomeDir /home/ximbao/reference_files/ --runThreadN 10 --readFilesIn ",
         path, files[i], " --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ",path.out, files[i]) %>% system()
}



library(Rsubread)
bams = paste0("/home/ximbao/Documents/Covid/aligned/bam/", dir("/home/ximbao/Documents/Covid/aligned/bam/"))
bams <- bams[-1]

bams = paste0("/home/ximbao/Documents/Covid/newRNA2/aligned/bams/", dir("/home/ximbao/Documents/Covid/newRNA2/aligned/bams/", pattern = ".bam$"))

bams = paste0("/home/ximbao/Documents/Covid/RNA_FullDepth/aligned/bams/", dir("/home/ximbao/Documents/Covid/RNA_FullDepth/aligned/bams/", pattern = ".bam$"))

fc <- featureCounts(files = bams, annot.ext = "~/reference_files/gencode.v37.annotation.gtf",
                    isGTFAnnotationFile = T, nthreads = 20)

