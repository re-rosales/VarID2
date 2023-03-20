## VarID2 is part of the RaceID (>= v0.2.4) package.
## This script contains additional functions used in the script VarID2_Figure1.R, with the code to reproduce the initial analysis and benchmarking of VarID2 method
## contained in Figures 1 and S1 of Rosales-Alvarez et al., manuscript

#Compute the standard deviation of a NB distribution:
sigmaNB <- function(z, r) {sqrt(z + 1/r * z^2)}
#z represents the mean, and r represents the dispersion parameter of the NB distribution

#Compute 1D or 2D estimation:

fitNBtbTest <- function(z, gamma=2, x0=0, lower=0, upper=100, grad=TRUE){
  z <- round(z,0)
  rt <- fitGammaRt(apply(z,2,sum))
  
  w <- apply(z,1,function(z,gamma,x0,rt,grad){
    mu <- mean(z)
    if (!grad){
      maxfun <- function(k){ v <- postfntb(k[2],z,x0,gamma,k[1],rt); if ( is.na(v) ) v = 0; if ( v == Inf ) v = 2e16; if ( v ==  -Inf ) v = -2e16; v }
      if ( mu == 0 ){
        eps = NA
      }else{
        opt <- optim(c(5,.5*(lower+upper)),maxfun, lower = c(0,lower), upper = c(200,upper))#,  method = "L-BFGS-B")
        mu <- opt$par[1]
        eps <- opt$par[2]
      }
    }else{
      gf <- function(x) gradp(x,z,x0,gamma,mu,rt)
      
      gu <- gf(upper)
      gl <- gf(lower)
      
      if ( gu > 0 & gl >= 0 ){
        eps <- lower
      }else if ( gu <= 0 & gl < 0 ){
        eps <- upper
      }else if ( gu < 0 & gl > 0 ){
        eps <- NA
      }else if ( gu == 0 & gl == 0 ){
        eps <- NA
      }else{  
        if ( mu == 0 ){
          eps = NA
        }else{
          #eps <- .External2(stats:::C_zeroin2, gf, lower, upper, f.lower=gl, f.upper=gu, .Machine$double.eps^0.25, maxiter=1000)[1]
          eps <- eval(parse(text = ".External2(stats:::C_zeroin2, gf, lower, upper, f.lower=gl, f.upper=gu, .Machine$double.eps^0.25, maxiter=1000)[1]"))
        }
      }
    }
    return( c(mu,eps,rt) )
    
  },gamma=gamma,x0=x0,rt=rt,grad=grad)
  
  w <- as.data.frame( t(w) )
  colnames(w) <- c("mu","epsilon","rt")
  w$alphaG <- (1+w$rt*w$epsilon)/w$rt
  return(w)
}

#Internal function for fitNBtbTest() function
gradp <- function(eps,z,x0,gamma,mu,rt){
  r <- rt/(1+eps*rt)
  sum( ( digamma( r + z ) - digamma( r ) + log( r ) + 1 - log( r + mu ) - ( r + z )/(r + mu )  ) * r**2 )  + 2 * ( eps - x0 )/gamma**2 * 1/( 1 + ( eps - x0 )**2/gamma**2)
}