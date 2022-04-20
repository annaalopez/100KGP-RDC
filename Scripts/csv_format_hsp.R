INFILES <- ("/home/alopezcarrasco/re_gecip/enhanced_interpretation/AnnaLopez/VEP/hsp/")
OUTFILES <- ("/home/alopezcarrasco/re_gecip/enhanced_interpretation/AnnaLopez/VEP/hsp/csv_format/")

file_name <- list.files(INFILES, pattern = ".txt")

read_files <- paste(INFILES, file_name, sep="") 
write_files <- paste(OUTFILES, paste(sub(".txt","", file_name),".csv"), sep="")

for (i in 1:length(read_files)) {
  csv_file <- (read.csv(file = read_files[i], header = TRUE, fill = FALSE, sep =""))
  write.csv(csv_file, file = write_files[i])
}



