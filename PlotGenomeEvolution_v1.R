library(DECIPHER)

L <- 5e6
maxCount <- 200L # number of operations to perform
#probabilities of performing an operation on each of the genomes
probG1 <- 0.5 #0.6
probG2 <- 0.5 #0.2
probG3 <- 0 #0.2
probG4 <- 0

invert <- function(x, mean=12, sd=2) {
	repeat {
		l <- rlnorm(1, mean, sd)
		if (l < width(x)[1])
			break
	}
	s <- sample(width(x)[1] - l, 1)
	replaceAt(x,
		IRanges(s, s + l - 1),
		reverseComplement(subseq(x, s, s + l - 1)))
}

insert <- function(x, mean=10, sd=1) {
	l <- rlnorm(1, mean, sd)
	s <- sample(width(x)[1], 1)
	replaceAt(x,
		IRanges(s, s - 1),
		paste(sample(DNA_BASES, l, T), collapse=""))
}

delete <- function(x, mean=10.7, sd=1) {
	repeat {
		l <- rlnorm(1, mean, sd)
		if (l < width(x)[1])
			break
	}
	s <- sample(width(x)[1] - l, 1)
	replaceAt(x,
		IRanges(s, s + l - 1))
}

rearrange <- function(x, mean=12, sd=2) {
	repeat {
		l <- rlnorm(1, mean, sd)
		if (l < width(x)[1])
			break
	}
	s <- sample(width(x)[1] - l, 1)
	z <- subseq(x, s, s + l - 1)
	y <- replaceAt(x,
		IRanges(s, s + l - 1))
	s <- sample(width(y)[1], 1)
	replaceAt(y,
		IRanges(s, s - 1),
		z)
}

transfer <- function(x, y, mean=10, sd=1) {
	repeat {
		l <- rlnorm(1, mean, sd)
		if (l < width(x)[1] & l < width(y)[1])
			break
	}
	s1 <- sample(width(x)[1] - l, 1)
	s2 <- sample(width(y)[1] - l, 1)
	replaceAt(y,
		IRanges(s2, s2 + l - 1),
		subseq(x, s1, s1 + l - 1))
}

evo1 <- DNAStringSet(paste(sample(DNA_BASES, L, T), collapse=""))
evo4 <- evo3 <- evo2 <- evo1
ops <- c(0, 0, 0, 0, 0)
ops_1 <- c(0, 0, 0, 0, 0)
ops_2 <- c(0, 0, 0, 0, 0)
count <- 0L
first <- TRUE
repeat {
	
	count <- count + 1L
	if (count > maxCount)
		break
	
	s1 <- sample(5, 1)
	s2 <- sample(4, 1, prob=c(probG1, probG2, probG3, probG4)) #probabilities to change the three genomes (G1, G2, G3)
	ops[s1] <- ops[s1] + 1
	if (s2 == 1) ops_1[s1] <- ops_1[s1] + 1 else if(s2==2) ops_2[s1] <- ops_2[s1] + 1
	if (s1==1) {
		message("\ncount = ", count, " invert", sep="")
		if (s2==1) {
			evo1 <- invert(evo1)
		} else if (s2==2) {
			evo2 <- invert(evo2)
		} else if (s2==3) {
			evo3 <- invert(evo3)
		} else if (s2==4) {
		  evo4 <- invert(evo4)
		}
	} else if (s1==2) {
		message("\ncount = ", count, " insert", sep="")
		if (s2==1) {
			evo1 <- insert(evo1)
		} else if (s2==2) {
			evo2 <- insert(evo2)
		} else if (s2==3) {
			evo3 <- insert(evo3)
		} else if (s2==4) {
		  evo4 <- insert(evo4)
		}
	} else if (s1==3) {
		message("\ncount = ", count, " delete", sep="")
		if (s2==1) {
			evo1 <- delete(evo1)
		} else if (s2==2) {
			evo2 <- delete(evo2)
		} else if (s2==3) {
			evo3 <- delete(evo3)
		} else if (s2==4) {
		  evo4 <- delete(evo4)
		}
	} else if (s1==4) {
		message("\ncount = ", count, " rearrange", sep="")
		if (s2==1) {
			evo1 <- rearrange(evo1)
		} else if (s2==2) {
			evo2 <- rearrange(evo2)
		} else if (s2==3) {
			evo3 <- rearrange(evo3)
		} else if (s2==4) {
		  evo4 <- rearrange(evo4)
		}
	} else if (s1==5) {
		message("\ncount = ", count, " transfer", sep="")
		s3 <- sample(0:2, 1)
		message(s3)
		if (s2==1) {
			if (s3==0) {
				evo2 <- transfer(evo1, evo2) # evo1 to evo2
			} else if(s3==1){
				evo3 <- transfer(evo1, evo3) # evo1 to evo3
			} else {
			  evo4 <- transfer(evo1, evo4)
			}
		} else if (s2==2) {
		  if (s3==0) {
		    evo1 <- transfer(evo2, evo1) # evo1 to evo2
		  } else if(s3==1){
		    evo3 <- transfer(evo2, evo3) # evo1 to evo3
		  } else {
		    evo4 <- transfer(evo2, evo4)
		  }
		} else if (s2==3) {
		  if (s3==0) {
		    evo1 <- transfer(evo3, evo1) # evo1 to evo2
		  } else if(s3==1){
		    evo2 <- transfer(evo3, evo2) # evo1 to evo3
		  } else {
		    evo4 <- transfer(evo3, evo4)
		  }
		}
		else {
		  if (s3==0) {
		    evo1 <- transfer(evo4, evo1) # evo1 to evo2
		  } else if(s3==1){
		    evo2 <- transfer(evo4, evo2) # evo1 to evo3
		  } else {
		    evo3 <- transfer(evo4, evo3)
		  }
		}
	}
	message(s2)
}

tf <- tempfile()
dbConn <- dbConnect(SQLite(), tf)
Seqs2DB(evo1, "XStringSet", dbConn, "Genome 1", v=F)
Seqs2DB(evo2, "XStringSet", dbConn, "Genome 2", v=F)
Seqs2DB(evo3, "XStringSet", dbConn, "Genome 3", v=F)
Seqs2DB(evo4, "XStringSet", dbConn, "Genome 4", v=F)
syn <- FindSynteny(dbConn, v=T)
dbDisconnect(dbConn)
