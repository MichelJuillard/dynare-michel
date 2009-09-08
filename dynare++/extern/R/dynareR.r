## $Id: dynareR.r 862 2006-08-04 17:34:56Z tamas $

## Copyright 2006, Tamas K Papp

dyn.load("dynareR.so")                  # FIXME: make it platform-independent

## FIXME hide auxiliary functions in a namespace

dynareR.indextensor <- function(ord, nume, nums) {
  nume*((nums^ord-1)/(nums-1))
}

dynareR.extracttensor <- function(tensor, ord, nume, nums) {
  aperm(array(tensor[dynareR.indextensor(ord,nume,nums)+(1:(nume*nums^ord))],
              c(nume,rep(nums,ord))),(ord+1):1)
}

dynareR.errormessages <- c("Sylvester exception",
                           "Dynare exception",
                           "OGU exception",
                           "Tensor library exception",
                           "K-order expansion library exception",
                           "Error matching names")

calldynare <- function(modeleq, endo, exo, parameters, expandorder,
                       parval, vcovmatrix, initval=rep(1,length(endo)),
                       numsteps=0, jnlfile="/dev/null") {
  ## check type of parameters
  local({
    is.charvector <- function(cv) { is.character(cv) && is.vector(cv) }
    stopifnot(is.charvector(modeleq) && is.charvector(endo) &&
              is.charvector(exo) && is.charvector(parameters) &&
              is.charvector(jnlfile))
  })
  stopifnot(is.numeric(expandorder) && is.vector(expandorder) &&
            (length(expandorder) == 1) && (expandorder >= 0))
  stopifnot(length(jnlfile) == 1)
  local({                               # variable names
    checkvarname <- function(v) {
      stopifnot(length(grep("[^a-zA-Z].*",v)) == 0) # look for strange chars 
    }
    checkvarname(endo)
    checkvarname(exo)
    checkvarname(parameters)
  })
  stopifnot(is.vector(parval) && is.numeric(parval))
  stopifnot(is.vector(initval) && is.numeric(initval))
  stopifnot(is.matrix(vcovmatrix) && is.numeric(vcovmatrix))
  stopifnot(is.numeric(numsteps) && is.vector(numsteps) &&
            (length(numsteps)==1))
  ## append semicolons to model equations if necessary
  modeleq <- sapply(modeleq, function(s) {
    if (length(grep("^.*; *$",s))==1)
      s
    else
      sprintf("%s;",s)
  })
  ## then concatenate into a single string
  modeleq <- paste(modeleq, collapse=" ")
  ## call dynareR
  nume <- length(endo)
  maxs <- length(endo)+length(exo)
  dr <- .C("dynareR",
           endo,as.integer(nume),
           exo,as.integer(length(exo)),
           parameters,as.integer(length(parameters)),
           modeleq,as.integer(expandorder),jnlfile,
           as.double(parval),as.double(vcovmatrix),
           as.double(initval),
           as.integer(numsteps),
           tensorbuffer=double(dynareR.indextensor(expandorder+1,nume,maxs)),
           numstate=integer(1), orderstate=integer(maxs),
           orderendo=integer(nume),
           orderexo=integer(length(exo)),
           newinitval=double(nume),
           error=integer(1),
           errormessage=character(1),
           kordcode=integer(1))
  ## check for errors
  kordcode <- 0
  if (dr$error == 0) {
    if (dr$error == 5) {
      list(kordcode=dr$kordcode - 251)  # magic dynare++ constant
    } else {
      ## return result
      with(dr, {
        nums <- numstate+length(exo)
        list(ss=dynareR.extracttensor(dr$tensorbuffer,0,nume,nums), # ss
             rule=sapply(1:expandorder,function (o) { # decision rule
               dynareR.extracttensor(dr$tensorbuffer,o,nume,nums)
             }),                            
             orderstate=orderstate[1:numstate], # state ordering
             orderendo=orderendo,           # endog. ordering
             orderexo=orderexo,             # exog. ordering
             newinitval=newinitval,         # new init values
             kordcode=0)
      })
    }
  } else {
    stop(sprintf("%s (\"%s\")",dynareR.errormessages[dr$error],
                 dr$errormessage))
  }
}
