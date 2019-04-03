"power.exp" <- function(rm=5, df.num=1, df.denom=10, delta=1, sigma=1, 
                        alpha=0.05, print=FALSE)
{
  if (df.denom < 1e-03)
  {
    warning("cannot calculate power because denominator degrees of freedom are close to zero")
    powr <- NA
  } else
  {
    fcrit <- qf(1-alpha, df.num, df.denom)
    lambda <- rm / 2 * (delta / sigma)^2
    powr <- 1 - pf(fcrit, df.num, df.denom, lambda)
    if (print == TRUE)
    { print.data.frame(data.frame(rm, df.num, df.denom, alpha, delta, sigma, lambda, powr))}
    powr

	}
  return(powr) 
}

"power.diff.r" <- function(r=5, multiple=1, df.num=1, df.denom=expression((df.num+1)*(r-1)), 
                           delta=1, sigma=1, alpha=0.05, power=0.8, print=FALSE)
{
	rm <- r * multiple
	powrdiff <- abs(power-power.exp(rm=rm, df.num=df.num, df.denom=eval(df.denom), 
	                                delta=delta, sigma=sigma, alpha=alpha, print=print))
  return(powrdiff)
}

"no.reps" <- function(multiple=1, df.num=1, df.denom=expression((df.num+1)*(r-1)), 
                        delta=1, sigma=1, alpha=0.05, power=0.8, tol = 0.1, 
                        print=FALSE)
{
	rmini <- 2
	rmaxi <- 50
	local <- TRUE
#
# loop while minimum is local
	while (local)
	{
		local <- FALSE
#
# get a minimum
		minimum <- optimize(power.diff.r, interval=c(rmini,rmaxi), tol=tol, print=print,  
					 multiple=multiple, df.num=df.num, df.denom=df.denom, 
					 delta=delta, sigma=sigma, alpha=alpha, power=power)
		if (minimum$minimum < 2) #stop if need less than 2 reps
		{
		  r <- 2
		} else
		{
		  #
		  # check for local minimum				
		  if (minimum$objective > tol)
		  {
		    local <- TRUE
		    r <- minimum$minimum
		    rm <- r * multiple
		    if (power < power.exp(rm=rm, df.num=df.num, df.denom=eval(df.denom), 
		                          delta=delta, sigma=sigma))
		    {
		      rmaxi <- floor(minimum$minimum)
		      rmini <- rmaxi / 2
		    }
		    else
		    {
		      rmini <- ceiling(minimum$minimum)
		      rmaxi <- 2 * rmaxi
		    }
		  }
		}
		r <- ceiling(minimum$minimum)
	}
#
# compute integer number of pure reps and the corresponding power
	rm <- r * multiple
	power <- power.exp(rm=rm, df.num=df.num, df.denom=eval(df.denom), delta=delta, sigma=sigma)
	list(nreps=r, power=power)
}

"power.diff.delta" <- function(delta=1, rm=5, df.num=1, df.denom=10, 
                               sigma=1, alpha=0.05, power=0.8, print=FALSE)
{
  powrdiff <- abs(power-power.exp(rm=rm, df.num=df.num, df.denom=df.denom, 
                                  delta=delta, sigma=sigma, alpha=alpha, print=print))
  return(powrdiff)
}

"detect.diff" <- function(rm=5, df.num=1, df.denom=10, sigma=1, alpha=0.05, power=0.8, tol = 0.001, 
                          print=FALSE)
{
  deltamini <- 0.01
  deltamaxi <- 1
  local <- TRUE
  #
  # loop while minimum is local
  while (local)
  {
    local <- FALSE
    #
    # get a minimum
    minimum <- optimize(power.diff.delta, interval=c(deltamini,deltamaxi), tol=tol, print=print,  
                        rm=rm, df.num=df.num, df.denom=df.denom, 
                        sigma=sigma, alpha=alpha, power=power)
    #
    # check for local minimum  			
    if (minimum$objective > tol)
    {
      local <- TRUE
      delta <- minimum$minimum
      if (power < power.exp(rm=rm, df.num=df.num, df.denom=eval(df.denom), 
                            delta=delta, sigma=sigma))
      {
        deltamaxi <- floor(minimum$minimum)
        deltamini <- deltamaxi / 2
      }
      else
      {
        deltamini <- ceiling(minimum$minimum)
        deltamaxi <- 2 * deltamaxi
      }
    }
  }
  #
  # compute integer number of pure reps and the corresponding power
  delta <- minimum$minimum
  power <- power.exp(rm=rm, df.num=df.num, df.denom=df.denom, delta=delta, sigma=sigma, alpha=alpha)
  return(delta)
}
