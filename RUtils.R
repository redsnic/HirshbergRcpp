# INFORMAZIONI GENERALI

nucleotides.order <- "ATGCSWRYKMBVHDN"

DNAfull.score.matrix <- matrix(
  c( 5,  -4,  -4,  -4,  -4,   1,   1,  -4,  -4,   1,  -4,  -1,  -1,  -1,  -2,
     -4,   5,  -4,  -4,  -4,   1,  -4,   1,   1,  -4,  -1,  -4,  -1,  -1,  -2,
     -4,  -4,   5,  -4,   1,  -4,   1,  -4,   1,  -4,  -1,  -1,  -4,  -1,  -2,
     -4,  -4,  -4,   5,   1,  -4,  -4,   1,  -4,   1,  -1,  -1,  -1,  -4,  -2,
     -4,  -4,   1,   1,  -1,  -4,  -2,  -2,  -2,  -2,  -1,  -1,  -3,  -3,  -1,
     1,   1,  -4,  -4,  -4,  -1,  -2,  -2,  -2,  -2,  -3,  -3,  -1,  -1,  -1,
     1,  -4,   1,  -4,  -2,  -2,  -1,  -4,  -2,  -2,  -3,  -1,  -3,  -1,  -1,
     -4,   1,  -4,   1,  -2,  -2,  -4,  -1,  -2,  -2,  -1,  -3,  -1,  -3,  -1,
     -4,   1,   1,  -4,  -2,  -2,  -2,  -2,  -1,  -4,  -1,  -3,  -3,  -1,  -1,
     1,  -4,  -4,   1,  -2,  -2,  -2,  -2,  -4,  -1,  -3,  -1 , -1,  -3,  -1,
     -4,  -1,  -1,  -1,  -1,  -3,  -3,  -1,  -1,  -3,  -1,  -2,  -2,  -2,  -1,
     -1,  -4,  -1,  -1,  -1,  -3,  -1,  -3,  -3,  -1,  -2,  -1,  -2,  -2,  -1,
     -1,  -1,  -4,  -1,  -3,  -1,  -3,  -1,  -3,  -1,  -2,  -2,  -1,  -2,  -1,
     -1,  -1,  -1,  -4,  -3,  -1,  -1,  -3,  -1,  -3,  -2,  -2,  -2,  -1,  -1,
     -2,  -2,  -2,  -2,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1),
  ncol=15, byrow = TRUE)

DNA.base.matrix <- matrix(
  c(
    1,-1,-1,-1,
    -1,1,-1,-1,
    -1,-1,1,-1,
    -1,-1,-1,1
  ), ncol = 4, byrow = TRUE)

# PROCEDURA DI TESTING

test.score <- function(s1,s2,mode){
  
  a_NWSW <- align_NWSW(s1, s2, mode,
                       nucleotides.order,
                       DNA.base.matrix,
                       -1)
  
  a_HIRSH <- align_Hirshberg(s1, s2, mode,
                             nucleotides.order,
                             DNA.base.matrix,
                             -1)
  
  a1 <- score_alignment(s1, s2, mode,
                        nucleotides.order,
                        DNA.base.matrix,
                        -1, a_NWSW)
  
  a2 <- score_alignment(s1, s2, mode,
                        nucleotides.order,
                        DNA.base.matrix,
                        -1, a_HIRSH)
  
  if(a1!=a2){
    print("Incongruenza individuata:")
    print_alignment(s1,s2,a_NWSW)
    print_alignment(s1,s2,a_HIRSH)
    print(a1)
    print(a2)
  }
  return(a1==a2)
} 


# PROFILING TEMPI E MEMORIA UTILIZZATA

profile_time <- function(mode, N, M=N, nrep=100){

  timeHirshberg <- rep(0,times=nrep)
  timeNWSW <- rep(0,times=nrep)
  
  pb <- txtProgressBar(min = 0, max = nrep, style = 3)
  
  for(rep in 1:nrep){ 
    setTxtProgressBar(pb, rep)
    
    s1 <- paste(collapse = "", sample(c("A","C","G","T"), sample(N,1,replace = TRUE), replace=TRUE))
    s2 <- paste(collapse = "", sample(c("A","C","G","T"), sample(M,1,replace = TRUE), replace=TRUE))
    
    # NWSW
    
    if(N<=20000){
      timeNWSW[rep] <- system.time({
        align_NWSW(s1,
                   s2,
                   mode,
                   nucleotides.order,
                   DNA.base.matrix, -1)
      }, gcFirst = FALSE)[3]
    }
    
    # HIRSHBERG
    
    timeHirshberg[rep] <- system.time({
      align_Hirshberg(s1,
                      s2,
                      mode,
                      nucleotides.order,
                      DNA.base.matrix, -1)
    }, gcFirst = FALSE)[3]
    
    
  }
  cat("\n")
  return(list(th=timeHirshberg,tnw=timeNWSW))
} 
