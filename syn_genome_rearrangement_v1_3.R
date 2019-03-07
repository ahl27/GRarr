##########
#Function to generate most likely rearrangement scenario
#Expects as input a synteny object
#Returns a list containing the rearrangement steps in order

#ERROR: Can't handle case where synteny object only has one block
##########
syn_genome_rearrangement <- function(synt, min_unique_length = 0, verbose = FALSE, actual=c()){
  
# ERROR CHECKING  --------------------------------------------------------------
  if(class(synt) != "Synteny")
    stop("Expected class of type 'Synteny'")
  if(length(actual) > 5)
    warning("vector provided is of length greater than 5. Only the first 5 entries will be used.")
  else if(length(actual) < 5 && length(actual) != 0)
    stop("vector provided for actual values is of incorrect length. Expected length 0 or 5.")
  
  #case where both genoems are identical or have a single synteny block causes problems
  #this calculates the number of unique genome sequences and then returns before the main
  #function to avoid any potentioal errors
  if(nrow(synt[[2,1]]) == 1){
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
  graph_to_genome <- function(graph, AA){
    if(!AA){
      start <- nrow(graph)
    }
    else{
      start <- nrow(graph) - 2
    }
    cur <- graph[start, 1]
    genome <- c()
    
    while (cur != (start - 1)){
      genome <- c(genome, (cur + (cur%%2)) / 2 * graph[cur,3])
      cur <- graph[cur + graph[cur,3], 1]
    }
    
    return(genome)
  }
  
  ##faster and more memory efficient way to append to a list
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
  
  #Finds inversions
  inversions <- function(graph, start=1, best=-1){
    ret_lst <- list()
    # ===== Inversions =====
    invert_count <- 0
    done <- FALSE
    better<- TRUE
    while (!done){ #inversions can create oriented gray edges, so we need to recheck all vertices if we make an inversion
      done <- TRUE
      j <- start-1
      if(invert_count > best && best != -1){
        better <- FALSE
        break()
      }
      while(j != (start + nrow(graph))){ #iterating over graph vertices looking for oriented gray lines
        i <- (j%%nrow(graph))+1
        j <- j + 1
        gen <- -1 * rev(graph_to_genome(graph, AA))
        if (identical(as.integer(gen), 1:length(gen))){
          next
        }
        if(AA && (i==nrow(graph) || i==(nrow(graph) - 1)))
          next
        if(graph[i,3] + graph[graph[i,2],3] == 0){ #if we find one, set done to FALSE and perform DCJ
          if(AA && (graph[i,2] == nrow(graph) || graph[i,2] == (nrow(graph)-1))){
            next
          }
          done <-FALSE
          graph <- DCJ(graph, i, graph[i,2]) #DCJ operation is performed on the two vertices incident to the gray line
          if(AA && graph[nrow(graph) - 2, 1] == (nrow(graph) - 3)){
            graph <- graph[-(nrow(graph) - 2),]
            graph <- graph[-(nrow(graph) - 2),]
            graph <- add_connection(graph, n+1, graph[n+1, 1], 1)
            graph <- add_connection(graph, n+1, graph[n+1, 2], 2)
            graph <- add_connection(graph, n+2, graph[n+2, 1], 1)
            graph <- add_connection(graph, n+2, graph[n+2, 2], 2)
            AA <- FALSE
          }
          
          # ---- Update Directionality of Graph ----
          start <- nrow(graph) #the last entry is the begin point for the mixed up genome
          if(AA)
            start <- nrow(graph) - 2
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
          
          
          #This algorithm may change the end cap's direction, which is not allowed
          #this fixes it
          graph[nrow(graph), 3] <- dir
          graph[nrow(graph)-1, 3] <- dir
          if (AA){
            graph[nrow(graph)-2, 3] <- dir
            graph[nrow(graph)-3, 3] <- dir
          }
          
          # ---- end updating directionality of graph ----
          ret_lst[[length(ret_lst)+1]] <- paste("inversion: ", paste(graph_to_genome(graph, AA), collapse=" "))
          
          invert_count <- invert_count + 1 #increment value (could also save what block(s) were inverted here)
          break()
        }
      }
    }
    
    return(list("better"=better, "graph"=graph, "count"=invert_count, "ops"= ret_lst))
  }
  
  # END HELPER FUNCTIONS----------------------------------------------------------
  
  
# FUNCTION BODY # 
  
  # SYNTENY OBJECT TO VECTOR OF PERMUTATIONS 
  
  # ==== Synteny Object to Blocks ====
  synblocks <- synt[[2,1]] #this grabs the information on the blocks
  block_matrix <- matrix(nrow=nrow(synblocks), ncol=4) #setting up the matrix to eventually return
  
  for(rowindex in 1:nrow(synblocks)){
    row <- synblocks[rowindex,] #grab information on the block
    
    start1 <- row[5] #start of block on the first genome
    blocklength <- row[7] - row[5] #length of block
    
    opp_value <- 0
    
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
  # ==== end simplifying vector of permutation
  
  # ===== End synteny to vector of permutations =====
  
  ## ===== Breakpoint Graph from Genome =====
  len <- 2 * length(genome)
  
  bpg <- matrix(nrow=len + 4, ncol = 3) #BreakPoint Graph
  
  #adding caps for the genomes
  #len + 1 is the end of genome A 
  #len + 2 is the beginning of genome A
  #len + 3 is the end of genome B
  #len + 4 is the beginning of genome B
  
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
  bpg <- add_connection(bpg, len, len + 3, 2)
  bpg <- add_connection(bpg, 1, len + 4, 2)
  
  bpg[len+1, 2] <- NA #removing an incorrect line generated by the last iteration of the above for loop
  
  # ===== finished creating breakpoint graph =====
  
  ## ===== Condensing End Caps of Graph =====
  graph <- bpg
  n <- len
  x <- graph[n+2, 1]
  val <- 1
  val_mod <- 1
  
  #finding cycle from the A begin cap
  while (x != (n + 2)){
    val <- val + val_mod
    val_mod <- val_mod * -1
    
    if (x == (n+1) || x == (n+3) || x == (n+4)){
      dir <- (x%%2) * -2 + 1
      end <- x
      x <-n+2
    }
    else{
      x <- graph[x, val]
    }
  }
  
  #identifying A cap with the B cap it leads to
  if (end != (n+1)){
    graph[n+2, 2] <- graph[end, 2]
    AA <- FALSE
  }
  else{
    graph <- add_connection(graph, n+1, n+2, 2)
    AA <- TRUE
  }
  
  #Finding cycle from the cap not identified with the A begin cap
  if(AA){
    dir <- 1
    graph <- add_connection(graph, n+3, n+4, 1)
  }
  else if (end == n+3){
    graph[n+1, 2] <- graph[n+4, 2]
  }
  else{
    graph[n+1, 2] <- graph[n+3, 2]
  }
  
  graph[n+1, 3] <- 1
  graph[n+2, 3] <- 1
  
  #identifying these two caps with each other
  if(!AA){
    graph <- graph[-(n+4),]
    graph <- graph[-(n+3),]
    
    #correcting some values to avoid out of bounds errors later
    graph <- add_connection(graph, n+1, graph[n+1, 1], 1)
    graph <- add_connection(graph, n+1, graph[n+1, 2], 2)
    graph <- add_connection(graph, n+2, graph[n+2, 1], 1)
    graph <- add_connection(graph, n+2, graph[n+2, 2], 2)
    graph[n+1, 3] <- dir
    graph[n+2, 3] <- dir
  }
  
  else{
    graph[n+3, 3] <- 1
    graph[n+4, 3] <- 1
  }
  # ===== finished condensing end caps =====
  
  # REARRANGEMENT OPERATIONS ---------------------------------------------------
  rearrangements <- list()
  rearrangements[[1]] <- paste("Original: ", paste(graph_to_genome(graph, AA), collapse=" "))
  
  # ===== Inversions =====
  cur <- inversions(graph)
  best <- cur$count
  for(j in 2:nrow(graph)){
    test <- inversions(graph, j, best)
    if(test$better){
      cur <- test
      best <- cur$count
    }
  }
  
  invert_count <- cur$count
  graph <- cur$graph
  rearrangements <- c(rearrangements, cur$ops)
  # ===== end of inversions =====
  
  # ===== Block Interchange =====
  block_count <- 0
  done <- FALSE
  
  ##Fixing AA and BB path 
  if ( AA && graph[nrow(graph)-3, 1] != graph[nrow(graph)-3, 2]){
    graph <- DCJ(graph, nrow(graph), graph[nrow(graph), 2])
    graph <- DCJ(graph, (nrow(graph) - 3), graph[nrow(graph)-3, 2])
    graph <- graph[-(nrow(graph) - 2),]
    graph <- graph[-(nrow(graph) - 2),]
    graph <- add_connection(graph, n+1, graph[n+1, 1], 1)
    graph <- add_connection(graph, n+1, graph[n+1, 2], 2)
    graph <- add_connection(graph, n+2, graph[n+2, 1], 1)
    graph <- add_connection(graph, n+2, graph[n+2, 2], 2)
    AA <- FALSE
    
    if(length(unique(graph[1:(nrow(graph)-2), 3])) == 1 && graph[1,3] == -1){
      graph[,3] <- abs(graph[,3])
    }
    
    rearrangements[[length(rearrangements)+1]] <- paste("block interchange: ", paste(graph_to_genome(graph, AA), collapse=" "))
    block_count <- block_count + 1
  }
  
  while (!done){
    done <- TRUE
    
    for (i in nrow(graph):1){
      #looking for two vertices that are out of order
      if(graph[i,1] != graph[i, 2]){ 
        done <- FALSE
        
        #If we find them, we perform a double cut and join
        graph <- DCJ(graph, i, graph[i, 2])
        
        #This double cut and join will create a circular intermediate (CI)
        #We can find this section of the genome and find a gray line attached to it
        #to add the CI back into the genome
        
        #convert genome into a vector and sort it, adding 0 to the front
        gen <- sort(graph_to_genome(graph, AA))
        gen <- c(0, gen)
        block_to_dcj <- 0
        for (i in 1:length(gen) - 1){
          #Now i represents what the next gene should be after gen[i]
          block_to_dcj <- i
          if (gen[i+1] != (i)){
            #once we find a missing gene, break
            break
          }
        }
        
        #then we just do a DCJ on the left vertex of this block and its connected gray vertex
        vertex <- 2*block_to_dcj - 1
        graph <- DCJ(graph, vertex, graph[vertex, 2])
        
        rearrangements[[length(rearrangements)+1]] <- paste("block interchange: ", paste(graph_to_genome(graph, AA), collapse=" "))
        block_count <- block_count + 1
        break
        
      }
    }#for i in nrow(graph)
    
  }#while(!done)
  
  # ===== end of block interchange =====
  
  # END REARRANGEMENT OPERATIONS -----------------------------------------------
  
  
  # OUTPUT ---------------------------------------------------------------------
  if (verbose){
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
    cat("\nRearrangements:")
    for (key in 1:length(rearrangements)){
      cat(paste(paste("\n[",toString(key-1),"]", sep=""), rearrangements[[key]]))
    }
    cat("\n")
    cat(c("\nUnique Genomic Sequences in Genome 1:", length(unique_1)/2, "\n"))
    cat(c("Unique Genomic Sequences in Genome 2:", length(unique_2)/2, "\n"))
  }
  
  if(length(actual) != 5){
    cat(c("\n[0]", rearrangements[[1]]))
    cat(c("\nTotal Inversions:", invert_count, "\n"))
    cat(c("Total Rearrangements:", block_count, "\n"))
    cat(c("Total Ins/Del/Transfer:", length(unique_1)/2 + length(unique_2)/2, "\n\n"))
  }
  else{
    cat(c("\n[0]", rearrangements[[1]]))
    cat(c("\nTotal Inversions:", invert_count, "(", actual[1], ")\n"))
    cat(c("Total Rearrangements:", block_count, "(", actual[4], ")\n"))
    cat(c("Total Ins/Del/Transfer:", length(unique_1)/2 + length(unique_2)/2, "(", actual[3]+actual[5]+actual[2], ")\n\n"))
  }
  # END OUTPUT -----------------------------------------------------------------
  #return(rearrangements)
}
