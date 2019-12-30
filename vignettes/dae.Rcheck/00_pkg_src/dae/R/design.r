designLatinSqrSys <- function(order, start = NULL)
{
  #process start argument
  if (is.null(start))
    start <- 1:order
  if (length(unique(start))!= order | !all(start >0 & start <= order))
    stop("start does not consist of order unique values between one and order")
  
  #generate design
  des <- lapply(start, 
                function(k, order)
                {
                  if (k == 1)
                    row <- 1:order
                  else
                    row <- c(1:order)[c(k:order,1:(k-1))]
                  return(row)
                }, order = order)
  des <- unlist(des)
  return(des)
}
