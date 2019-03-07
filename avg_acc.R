source("/Users/aidanlakshman/Documents/github/wrightlab/syn_genome_rearrangement_v3_5.R")
#source("/Users/aidanlakshman/Documents/github/wrightlab/prev_alg_versions/syn_genome_rearrangement_v2_0.R")

file_dir <- "/Users/aidanlakshman/Documents/github/wrightlab/algorithm_results/tests/synteny_200/"
avg_out_file3 <- "/Users/aidanlakshman/Documents/github/wrightlab/algorithm_results/tests/avgACC3_200.txt"
min_out_file3 <- "/Users/aidanlakshman/Documents/github/wrightlab/algorithm_results/tests/minACC3_200.txt"
out_file2 <- "/Users/aidanlakshman/Documents/github/wrightlab/algorithm_results/tests/varACC2_200.txt"

ctr <- 1
avg_inv_3 <- vector(length=100)
avg_bi_3 <- vector(length=100)
avg_idt_3 <- vector(length=100)
min_inv_3 <- vector(length=100)
min_bi_3 <- vector(length=100)
min_idt_3 <- vector(length=100)
avg_inv_2 <- vector(length=100)
avg_bi_2 <- vector(length=100)
avg_idt_2 <- vector(length=100)

#cat("\n", file=min_out_file3)
#cat("\n", file=avg_out_file3)
#cat("\n", file=out_file2)

for (i in 101:2000){
  val <- toString(2*i)
  while(nchar(val) < 4)
    val <- paste("0", val, sep='')
  filename <- paste(file_dir, "Result", val, ".RData", sep="")
  
  load(filename)
  
  cat("\nResult", val, ".RData\n", sep='')
  inv_data <- output_lst[[400]]
  inv_syn <- output_lst[[399]]
  
  #cat("\tRunning version 3...\n")
  predictions <- GRarr_v3(inv_syn, verbose=F, mean=T, test_run=2)
  avg_inv_3[i] <- (predictions[1] - inv_data[1])
  avg_bi_3[i] <- (predictions[2] - inv_data[2])
  avg_idt_3[i] <- (predictions[3] - inv_data[3])
  
  cat("\n")
  
  predictions <- GRarr_v3(inv_syn, verbose=F, mean=F, test_run=2)
  min_inv_3[i] <- (predictions[1] - inv_data[1])
  min_bi_3[i] <- (predictions[2] - inv_data[2])
  min_idt_3[i] <- (predictions[3] - inv_data[3])
  
  cat("\n")
  #cat("\tRunning version 2...\n")
  predictions <- GRarr_v2(inv_syn, 10, F, c(), 2)
  avg_inv_2[i] <- (predictions[1] - inv_data[1])
  avg_bi_2[i] <- (predictions[2] - inv_data[2])
  avg_idt_2[i] <- (predictions[3] - inv_data[3])
  
  if ( i %% 100 == 0){
    cat(mean(avg_inv_3), sd(avg_inv_3), 
        mean(avg_bi_3), sd(avg_bi_3), 
        mean(avg_idt_3), sd(avg_idt_3), 
        "\n", file=avg_out_file3, sep=" ", append=TRUE)
    
    cat(mean(min_inv_3), sd(min_inv_3), 
        mean(min_bi_3), sd(min_bi_3), 
        mean(min_idt_3), sd(min_idt_3), 
        "\n", file=min_out_file3, sep=" ", append=TRUE)
    
    cat(mean(avg_inv_2), sd(avg_inv_2), 
        mean(avg_bi_2), sd(avg_bi_2), 
        mean(avg_idt_2), sd(avg_idt_2), 
        "\n", file=out_file2, sep=" ", append=TRUE)
    
    avg_inv_3 <- vector(length=100)
    avg_bi_3 <- vector(length=100)
    avg_idt_3 <- vector(length=100)
    min_inv_3 <- vector(length=100)
    min_bi_3 <- vector(length=100)
    min_idt_3 <- vector(length=100)
    avg_inv_2 <- vector(length=100)
    avg_bi_2 <- vector(length=100)
    avg_idt_2 <- vector(length=100)
  }
}