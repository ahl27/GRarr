##########
#Function to generate most likely rearrangement scenario
#Expects as input a synteny object
#Returns a list containing the rearrangement steps in order

##########
GRarr_v4 <- function(synt, num_runs=-1, verbose = FALSE, ultraverbose=FALSE, mean = FALSE, #main user parameters
                     opp_value = 0, min_unique_length = 10, #extra user parameters
                     actual=c(), test_run=0){ #internal parameters

  if(ultraverbose)
    verbose <- TRUE
  if(class(synt) != "Synteny")
    stop("Expected class of type 'Synteny'")
  
  
  # translocations <- 0
  # duplications <- 0
  # inversions <- 0
  # transpositions <- 0
  # insertdeletes <- 0
  
  num_genomes <- nrow(synt)
  rearrangements <- matrix(0, nrow=num_genomes, ncol=5)
  row_names <- vector(length=num_genomes)
  for (ind in 1:num_genomes){
    name <- paste("Genome", ind)
    row_names[ind] <- name
  }
  col_names <- c("Translocations", "Duplications", "Inversions", "Transpositions", "Insert/Deletes")
  rownames(rearrangements) <- row_names
  colnames(rearrangements) <- col_names
  rates <- vector("list", num_genomes)
  for ( i in 1:num_genomes ){
    for ( j in i:num_genomes ){
      if (j == i)
        next
      gen1 <- max(i, j)
      gen2 <- min(i, j)
      
      gen_info <- synt[[gen1, gen2]]
      
      num_chrom <- max(gen_info[,1],gen_info[,2]) 
      if(verbose)
        cat("\nComputing Scenario for Genomes ", i, " and ", j, "...", sep='')
      for ( chrom in 1:num_chrom ){
        rows <- gen_info[gen_info[,1] == chrom | gen_info[,2] == chrom,]
        if ( class(rows) != "matrix" ) #if only one row matches, it returns a vector
          rows <- matrix(rows, nrow=1)

        #eliminating rows we've already seen
        rows <- rows[rows[,1] >= chrom,]
        rows <- rows[rows[,2] >= chrom,]
        
        
        if (nrow(rows) == 0) #if no rows match, continue
          next()
        
        for (k in 1:nrow(rows)){
          if (rows[k, 1] != rows[k, 2]){
            #determine if translocation or duplication
            gen1_match <- gen_info[gen_info[,5]==rows[k,5] & gen_info[,7]==rows[k,7],]
            gen2_match <- gen_info[gen_info[,5]==rows[k,5] & gen_info[,7]==rows[k,7],]
            found <- F
            if(class(gen1_match) == "matrix"){
              rearrangements[gen1, 2] <- rearrangements[gen1, 2] + 1
              found <- T
            }
            if(class(gen2_match) == "matrix"){
              rearrangements[gen2, 2] <- rearrangements[gen2, 2] + 1
              found <- T
            }
            
            if (!found){
              rearrangements[gen1, 1] <- rearrangements[gen1, 1] + 1
              rearrangements[gen2, 1] <- rearrangements[gen2, 1] + 1
            }
              
          }
        }
        matching <- rows[rows[,1] == rows[,2],]
        if (class(matching) != "matrix"){
          matching <- matrix(matching, nrow=1)
        }
        if (nrow(matching) == 0){
          next()
        }
        if(verbose)
          cat("\n\tChromosome ", chrom, " of ", num_chrom, ":\n", sep='')  
        chrom_results <- rearrange_chromosome(matching, synt, verbose=verbose, ultraverbose=ultraverbose, mean=mean, chrom=chrom,
                                              opp_value = opp_value, min_unique_length=min_unique_length,
                                              first_gen = gen1, second_gen=gen2, test_run=2)
        
        rearrangements[gen1, 3] <- rearrangements[gen1, 3] + chrom_results[1]
        rearrangements[gen2, 3] <- rearrangements[gen2, 3] + chrom_results[1]
        rearrangements[gen1, 4] <- rearrangements[gen1, 4] + chrom_results[2]
        rearrangements[gen2, 4] <- rearrangements[gen2, 4] + chrom_results[2]
        rearrangements[gen1, 5] <- rearrangements[gen1, 5] + chrom_results[3]
        rearrangements[gen2, 5] <- rearrangements[gen2, 5] + chrom_results[3]
      }
    }
  }
  if(verbose)
    cat("\nComplete.\n")
  return(rearrangements)
}

rearrange_chromosome <- function(synblocks, syn, num_runs=-1, verbose = FALSE, ultraverbose=FALSE, mean = FALSE, #main user parameters
                                 opp_value = 0, min_unique_length = 10, #extra user parameters
                                 actual=c(), test_run=0, first_gen=1, second_gen=2, chrom=-1){ #internal parameters  
# ERROR CHECKING  --------------------------------------------------------------
  #set.seed(13)
  if(length(actual) > 5)
    warning("vector provided is of length greater than 5. Only the first 5 entries will be used.")
  else if(length(actual) < 5 && length(actual) != 0)
    stop("vector provided for actual values is of incorrect length. Expected length 0 or 5.")
  
  if ( second_gen == first_gen)
    stop("Error: Incorrect genome indices provided. Make sure indices are not the same.")
  if ( second_gen < first_gen ){
    temp <- second_gen
    second_gen <- first_gen
    first_gen <- temp
  }
  
  #case where both genomes are identical or have a single synteny block causes problems
  #this calculates the number of unique genome sequences and then returns before the main
  #function to avoid any potential errors
  if(nrow(syn[[second_gen,first_gen]]) == 1){
    len_1 <- syn[[1,1]]
    len_2 <- syn[[2,2]]
    start_1 <- synt[[2,1]][,5]
    start_2 <- synt[[2,1]][,6]
    end_1 <- synt[[2,1]][,7]
    end_2 <- synt[[2,1]][,8]
    
    ins_del_trans <- 0
    if(start_1 != 1)
      ins_del_trans <- 1
    if(start_2 != 1)
      ins_del_trans <- ins_del_trans + 1
    if(end_1 != len_1)
      ins_del_trans <- ins_del_trans + 1
    if(end_2 != len_2)
      ins_del_trans <- ins_del_trans + 1
    
    if(test_run == 2)
      return(c(0, 0, ins_del_trans))
    if(length(actual) != 5){
      cat("[0] Original: 1")
      cat(c("\nTotal Inversions: 0\n"))
      cat(c("Total Rearrangements: 0\n"))
      cat(c("Total Ins/Del/Transfer: ", ins_del_trans, "\n\n"))
    }
    else{
      cat("[0] Original: 1")
      cat(c("\nTotal Inversions: 0 (", actual[1], ")\n"))
      cat(c("Total Rearrangements: 0 (", actual[4], ")\n"))
      cat(c("Total Ins/Del/Transfer: ", ins_del_trans, " (", actual[3]+actual[5]+actual[2], ")\n\n"))
    }
    return()
  }
  
# END ERROR CHECKING------------------------------------------------------------
  
  # HELPER FUNCTIONS -------------------------------------------------------------
  
  ## Adds connection between two vertices of graph
  add_connection <- function(mat, i1, i2, col){
    mat[i1, col] <- i2
    mat[i2, col] <- i1
    return(mat)
  }
  
  ## Double cut and join operation
  DCJ <- function(graph, v1, v2){
    ##v1 and v2 are joined by a gray line
    v1b <- graph[v1, 1] #black connection of v1
    v2b <- graph[v2, 1] #black connection of v2
    
    graph <- add_connection(graph, v1, v2, 1)
    graph <- add_connection(graph, v1b, v2b, 1)
    
    return(graph)
  }
  
  ## Convert adjacency graph to genome (as vector)
  graph_to_genome <- function(graph){
    start <- nrow(graph)
    cur <- graph[start, 1]
    genome <- c()
    
    while (cur != (start - 1)){
      genome <- c(genome, (cur + (cur%%2)) / 2 * graph[cur,3])
      cur <- graph[cur + graph[cur,3], 1]
    }
    
    return(genome)
  }
  
  condense_genome <- function(genome){
    i <- 1
    while(i <= length(genome)){
      start <- i
      while(genome[i] + 1 == genome[i+1]){
        i <- i+1
        genome <- genome[-(i)]
      }
      if (start == i)
        i <- i+1
    }
    return(genome)
  }
  
  #faster and more memory efficient way to append to a list
  lappend <- function(lst, val, increase=10){
    
    #find first empty value
    first_empty <- length(lst) + 1
    for(i in 1:length(lst)){
      if(isEmpty(lst[[i]])){
        first_empty <- i
        break()
      }
    }
    
    #get current length of list
    len <- length(lst)
    tmp <- len
    
    #find new length to be able to add new value(s)
    #should overestimate slightly, this is intentional to reduce calls to this function
    while((tmp - first_empty) < length(val)){
      tmp <- tmp + increase
    }
    
    #if we increased the size of the list, copy the old list into the new one  
    if(tmp != len){
      new_lst <- vector("list", length = tmp)
      new_lst[1:len] <- lst
      lst <- new_lst
    }
    
    #add new value to list
    lst[first_empty:(first_empty+length(val)-1)] <- val 
    return(lst)
  }
  
  update_direc <- function(graph){
    start <- nrow(graph) #the last entry is the begin point for the mixed up genome
    cur <- graph[start, 1] #travel to the vertex connected by a black line
    
    #note here that the genome is represented as a sequence of permutations with a direction
    #ex. < -3 2 -4 5 1 > becomes < +3 -3 -2 +2 +4 -4 -5 +5 -1 +1 >, where genes are read from negative to positive
    #in the graph, instead of representing the ends of block n as -n and +n, we use 2n-1 and 2n
    #so the program represents the above as < 6 5 3 4 7 8 9 10 1 2 >
    #thus is a gene is even -> odd it's reversed, and if it's odd -> even it's in the correct orientation
    
    
    while (cur != (start-1)){
      if(cur%%2 == 1){ #if value is odd, gene is in the correct orientation (+)
        #we store a direction of 1 in both edges of the gene, and set cur to the other side of the block
        graph[cur, 3] <- 1 
        cur <- cur + 1
        graph[cur, 3] <- 1
      }
      else{ #if value is even, gene is in the reversed orientation (-)
        #same as above, except setting direction to -1
        graph[cur, 3] <- -1
        cur <- cur-1
        graph[cur, 3] <- -1
      }
      
      #follow the black path to the next gene in the mixed up genome
      cur <- graph[cur, 1]
    }
    
    
    #This may change the end caps' direction, which is not allowed
    #this fixes it
    graph[nrow(graph), 3] <- 1
    graph[nrow(graph)-1, 3] <- 1
    return(graph)
  }
  

  generate_rearrangements <- function(graph){
    rearrangements <- vector("list", length=nrow(graph))
    rearrangements[[1]] <- paste("Original: ", paste(graph_to_genome(graph), collapse=" "))

    block_count <- 0
    invert_count <- 0
    gen <- graph_to_genome(graph)
    t <- 0
    len <- length(gen) + 1
    gen <- c(0, gen, len)
    done <- FALSE
    while(!identical(as.integer(gen), 0:(length(gen)-1))){
      t <- t + 1
      if (t == (3*len)) {
        #The number of operations is at most double the number of blocks
        #(as we could sort any vector into place by transposing each block one at a time and then inverting it)
        #Thus if we ever get to triple the number of blocks, we can be sure it's in an infinite loop
        fname <- "~/ERRONEOUS_SYNT_OBJ.RData"
        save(synt, file=fname)
        err_msg <- paste("Infinite Loop, bug in code. Error-producing synteny object saved to ", fname, sep="")
        stop(err_msg)
      }
      
      gen <- graph_to_genome(graph)
      len <- length(gen) + 1
      gen <- c(0, gen, len)
      done <- TRUE
      order <- sample(1:len)
      for(i in order){
        sign1 <- sign(gen[i])
        sign2 <- sign(gen[i+1])
        if (sign1==0) sign1 <- 1
        
        #block interchange====
        for(p in c(1,-1)){
          if (done == FALSE)
            next()
          if (p == 1)
            subsec <- gen[i:len]
          else
            subsec <- gen[1:(i-1)]
          
          if(((gen[i]+p) %in% subsec) && (gen[i] + p) != gen[i+p]){
            #no block interchanges to connect inverted elements
            #getting length of interchange
            j <- i+p
            len_block <- 0
            cont <- FALSE
            while(gen[j] != (gen[i]+p)){
              len_block <- len_block + 1
              if(gen[j] == gen[length(gen)]){
                cont <- TRUE
                break()
              }
              j <- j+p
            }
            if(cont) next()
            #ex. 0 1 4 3 2 5, block interchange required
            #first cut creates the circular intermediate
            if((gen[i] < 0 && p==1) || (gen[i] > 0 && p==-1))
              v1 <-  abs(gen[i])*2 - 1
            else
              v1 <- abs(gen[i])*2
            
            new_graph <- DCJ(graph, v1, graph[v1,2])
            circular_intermediate <- sort(gen[(!gen%in%(graph_to_genome(new_graph)))])
            #we now have 0 1 2 5, with a CI of (4 3)
            #take genome as a vector and sort it, adding 0 to the front
            cut_vertex <- -1
            
            for(vert in circular_intermediate){
              if(vert == 0 || vert==max(circular_intermediate))
                next()
              if((vert-1)%in%gen && !(vert-1)%in%circular_intermediate){
                if (vert < 0){
                  cut_vertex <- (abs(vert))*2
                  break()
                }
                else{
                  cut_vertex <- (vert)*2-1
                  break()
                }
              }
              if((vert + 1)%in%gen && !(vert+1)%in%circular_intermediate){
                if (vert < 0){
                  cut_vertex <- (abs(vert))*2 -1
                  break()
                }
                else{
                  cut_vertex <- (vert)*2 
                  break()
                }
              }
              
              if((vert-1) %in% gen && !(vert-1)%in%circular_intermediate){
                if(vert < 0){
                  cut_vertex <- abs(vert)*2
                  break()
                }
                else{
                  cut_vertex <- (vert)*2 - 1
                  break()
                }
              }
            }
            
            if(cut_vertex != -1){
              graph <- new_graph
              #then we just do a DCJ on the left vertex of this block and its connected gray vertex
              graph <- DCJ(graph, cut_vertex, graph[cut_vertex, 2])
              rearrangements[[block_count + invert_count +2]] <- paste("block interchange: ", paste(graph_to_genome(graph), collapse=" "), "{", len_block, "}")
              block_count <- block_count + 1
              done <- FALSE
              break()   
            }
            else {
              gen <- graph_to_genome(graph)
              gen <- c(0, gen, length(gen) + 1)
            }
          }
          #end block interchange====
          
          #inversion====
          if (!(p==-1 && i==1) && gen[i] < 0 && ((-1*(gen[i]+p)) %in% gen)){
            v1 <- -1
            if(gen[i]*p > 0)
              v1 <- (abs(gen[i]*2)) 
            else 
              v1 <- (abs(gen[i])*2-1)
            
            #ex. 0 1 -3 -2 4
            i1 <- i
            i2 <- match((-1*(gen[i]+p)), gen)
            len_block <- abs(i1 - i2)
            
            if(v1 != -1){
              graph <- DCJ(graph, v1, graph[v1, 2])
              graph <- update_direc(graph)
              rearrangements[[invert_count + block_count + 2]] <- paste("inversion: ", paste(graph_to_genome(graph), collapse=" "), "{", len_block, "}")
              invert_count <- invert_count + 1
              done <- FALSE
              break()
            }
          }
        }
        #end inversion====
      }
    }
    total_ops <- block_count + invert_count
    
    rearrangements <- rearrangements[-((total_ops+2):length(rearrangements))]
    return(list("rearrangements"=rearrangements, "graph"=graph, "invert_count"=invert_count, "block_count"=block_count))
  }

  
  # END HELPER FUNCTIONS----------------------------------------------------------
  
  
# FUNCTION BODY # 
  
  # SYNTENY OBJECT TO VECTOR OF PERMUTATIONS 

  # ==== Synteny Object to Blocks ====
  block_matrix <- matrix(nrow=nrow(synblocks), ncol=4) #setting up the matrix to eventually return
  
  for(rowindex in 1:nrow(synblocks)){
    row <- synblocks[rowindex,] #grab information on the block
    
    start1 <- row[5] #start of block on the first genome
    blocklength <- row[7] - row[5] #length of block
    
    start2 <- row[6]
    #get direction of the block
    if(row[3] == opp_value){ 
      dir <- 1
    }
    else{ 
      dir <- -1
    }
    
    block_matrix[rowindex, 1] <- start1
    block_matrix[rowindex, 2] <- start2
    block_matrix[rowindex, 3] <- blocklength
    block_matrix[rowindex, 4] <- dir
  }
  
  # ==== end synteny object to blocks ====

  # ==== Blocks to Permutation Matrix ====
  sorted_mat <- block_matrix
  
  sorted_mat <- sorted_mat[order(sorted_mat[,1]),] #sort based on the first genome
  indices <- 1:nrow(block_matrix) #simple vector from 1 to {n}, to order the blocks
  sorted_mat <- cbind(sorted_mat, indices) #attach the indices, now sorted_mat[,5] represents permutation order for genome1
  
  sorted_mat <- sorted_mat[order(sorted_mat[,2]),] #sort based on start coord on genome2
  sorted_mat <- cbind(sorted_mat, indices) #append indices, now sorted_mat[,6] represents permutation order for genome2
  sorted_mat[,6] <- sorted_mat[,4] * sorted_mat[,6] #multiply by direction to get pos/neg permutations for genome2
  
  sorted_mat <- sorted_mat[order(sorted_mat[,5]),] #sort by permutation order for genome1
  
  permutation_order <- cbind(sorted_mat[,5], sorted_mat[,6]) #now columns 5 and 6 represent permutation order for genomes 1 and 2, resp.
  genome <- permutation_order[,2]
  
  # ==== end blocks to matrix ====
  
  
  # ==== Finding genomic blocks not shared between genomes ====
  start_ind_gen1 <- cbind(sorted_mat[,1], sorted_mat[,3])
  start_ind_gen2 <- sorted_mat[,2:3]
  
  gen1 <- start_ind_gen1[sort.list(start_ind_gen1[,1]),]
  gen2 <- start_ind_gen2[sort.list(start_ind_gen2[,1]),]
  
  unique_1 <- c()
  unique_2 <- c()
  
  cur <- 1
  prev <- 1
  for (i in 1:nrow(gen1)){
    cur <- gen1[i,1]
    if(cur != prev && ((cur-prev) > min_unique_length)){
      unique_1 <- c(unique_1, prev, cur-prev)
    }
    prev <- cur + gen1[i,2]+1
  }
  
  cur <- 1
  prev <- 1
  for (i in 1:nrow(gen2)){
    cur <- gen2[i,1]
    if(cur != prev && ((cur-prev) > min_unique_length)){
      unique_2 <- c(unique_2, prev, cur-prev)
    }
    prev <- cur + gen2[i,2] + 1
  }
  # ==== end finding unique genomic blocks ====
  
  # ==== Simplifying vector of permutations ====
  i <- 1
  indices_to_remove <- c()
  while(i < length(genome)){
    start <- i
    while(i < length(genome) && (genome[i] + 1) == genome[i+1]){
      i <- i+1
      indices_to_remove <- c(indices_to_remove, i)
    }
    if (start == i)
      i <- i+1
  }
  
  ctr <- 0
  for(ind in indices_to_remove){
    genome <- genome[-(ind-ctr)]
    ctr <- ctr + 1
  }
  
  new <- 1:length(genome)
  temp <- sort(abs(genome))
  
  conv <- as.list(setNames(new, temp))
  
  new_genome <- c()
  for(entry in genome){
    if(entry < 0)
      mult <- -1
    else
      mult <- 1
    
    new_genome <- c(new_genome, conv[[toString(abs(entry))]]*mult)
  }
    
  genome <- new_genome
  num_blocks <- length(genome)
  # ==== end simplifying vector of permutation
  
  # ===== End synteny to vector of permutations =====
  
  # ===== Breakpoint Graph from Genome =====
  len <- 2 * length(genome)
  
  bpg <- matrix(nrow=len + 2, ncol = 3) #BreakPoint Graph
  
  #adding caps for the genomes
  #len + 1 is the end cap
  #len + 2 is the beginning cap
  
  #every iteration we create a black line between this vertex and the previous vertex in A
  prev_black <- len + 2 
  
  for (i in 1:length(genome)){
    val <- genome[i] #get permutation number
    
    if (val < 0){ #if negative, map permutation n to 2n and 2n-1 (in that order)
      v1 <- abs(val) * 2
      v2 <- (abs(val) * 2) - 1
      bpg[v1, 3] <- -1
      bpg[v2, 3] <- -1
      bpg <- add_connection(bpg, v2, v2 - 1, 2) #add gray line between left vertex of block and previous block's right vertex
      bpg <- add_connection(bpg, v1, v1 + 1, 2) #add gray line between right vertex of block and next block's left vertex
    }
    else { #if positive, map permutation n to 2n-1 and 2n (in that order)
      v1 <- (val * 2) - 1
      v2 <- val * 2
      bpg[v1, 3] <- 1
      bpg[v2, 3] <- 1
      bpg <- add_connection(bpg, v1, v1 - 1, 2) #add gray line between left vertex of block and previous block's right vertex
      bpg <- add_connection(bpg, v2, v2 + 1, 2) #add gray line between right vertex of block and next block's left vertex
    }
    
    bpg <- add_connection(bpg, v1, prev_black, 1) #add black line between current vertex and prev. vertex in A
    prev_black <- v2 #store right vertex to connect to the next one
  }
  #black line between last value of A and end cap of A
  bpg <- add_connection(bpg, prev_black, len+1, 1)
  
  #gray line between last value of B and end cap of B
  #and gray line between first value of B and start cap of B
  bpg <- add_connection(bpg, len, len + 1, 2)
  bpg <- add_connection(bpg, 1, len + 2, 2)
  # ===== finished creating breakpoint graph =====
  
  # ===== Rearrangement Scenarios =====
  if (num_runs == -1)
    num_runs <- length(graph_to_genome(bpg))
  
  if(verbose || ultraverbose)
    progress <- txtProgressBar(max = num_runs, char = "=", style=3)
  
  if(mean){
    invert_count <- 0
    block_count <- 0
    
    for(i in 1:num_runs){
      rearr <- generate_rearrangements(bpg)
      invert_count <- invert_count + rearr$invert_count
      block_count <- block_count + rearr$block_count
      if (verbose || ultraverbose)
        setTxtProgressBar(progress, i)
    }
    rearrangements <- rearr$rearrangements
    invert_count <- invert_count / num_runs
    block_count <- block_count / num_runs
  }
  else{
    block_count <- 3*length(graph_to_genome(bpg))
    invert_count <- block_count
    
    for (i in 1:num_runs){
      rearr <- generate_rearrangements(bpg)
      if((rearr$block_count + rearr$invert_count) < (block_count+invert_count)){
        block_count <- rearr$block_count
        invert_count <- rearr$invert_count
        rearrangements <- rearr$rearrangements
        graph <- rearr$rearrangements
      }
      if (verbose || ultraverbose)
        setTxtProgressBar(progress, i)
    }
  }

  # ===== end rearrangement scenarios =====
  
  # OUTPUT ---------------------------------------------------------------------
  if (ultraverbose){
    cat(c("\nUnique Genomic Sequences: (min. length ", min_unique_length, ")\n\n"))
    cat("\tGenome 1:\n")
    if(length(unique_1) != 0){
      for(strand in 1:(length(unique_1)/2)){
        st <- unique_1[strand*2-1]
        ln <- unique_1[strand*2]
        cat(c("\t[", toString(strand), "]", st, "(Length:", ln, ")\n"))
      }
    }
    else
      cat("\tNone\n")
    cat("\n\tGenome 2: \n")
    if(length(unique_2) != 0){
      for(strand in 1:(length(unique_2)/2)){
        st <- unique_2[strand*2-1]
        ln <- unique_2[strand*2]
        cat(c("\t[", strand, "]", st, " (Length:", ln, ")\n"))
      }
    }
    else
      cat("\tNone\n")
    if(!mean){
      cat("\nRearrangements:")
      for (key in 1:length(rearrangements)){
        cat(paste(paste("\n[",toString(key-1),"]", sep=""), rearrangements[[key]]))
      }
    }
    cat("\n")
    cat(c("\nUnique Genomic Sequences in Genome 1:", length(unique_1)/2, "\n"))
    cat(c("Unique Genomic Sequences in Genome 2:", length(unique_2)/2, "\n"))
  }
  
  if(test_run == 1){
    cat(c("\n", invert_count, "(", actual[1], ")\n"))
    cat(c(block_count, "(", actual[4], ")\n"))
    cat(c(length(unique_1)/2 + length(unique_2)/2, "(", actual[3]+actual[5]+actual[2], ")\n\n"))
  }
  else if(test_run == 2){
    return(c(invert_count, block_count, length(unique_1)/2 + length(unique_2)/2))
  }
  
  else if(length(actual) != 5){
    if(verbose)
      cat(c("\n[0]", rearrangements[[1]]))
    else
      cat("\nNumber of Blocks:", num_blocks, "\n")
    cat(c("\nTotal Inversions:", invert_count, "\n"))
    cat(c("Total Rearrangements:", block_count, "\n"))
    cat(c("Total Ins/Del/Transfer:", length(unique_1)/2 + length(unique_2)/2, "\n\n"))
  }
  else{
    if(ultraverbose)
      cat(c("\n[0]", rearrangements[[1]]))
    else
      cat("\nNumber of Blocks:", num_blocks, '\n')
    cat(c("\nTotal Inversions:", invert_count, "(", actual[1], ")\n"))
    cat(c("Total Rearrangements:", block_count, "(", actual[4], ")\n"))
    cat(c("Total Ins/Del/Transfer:", length(unique_1)/2 + length(unique_2)/2, "(", actual[3]+actual[5]+actual[2], ")\n\n"))
  }
  # END OUTPUT -----------------------------------------------------------------
  return(list("invert"=invert_count, "block"=block_count, "insdel"=(length(unique_1)/2 + length(unique_2)/2)))
}
