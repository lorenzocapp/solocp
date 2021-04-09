#' Subset the change points computed either basad.cp or solo.cp
#' Corresponds to point 2-4 in Algorithm 1 in the manuscript 
#'
#' @param ratio marginal inclusion probabilities, as computed for example by the solocp_single function
#' @param del Delta parameter in the parameter, mininum spacing between change points to be considered nonconsecutive
#' (we set it equal to 5 by defauly)
#' 
#' @return subset of change points indexes
#'  
#' @export

subset_changepoints <- function(ratio,del=5){
  first <- which(ratio >=.5)
  n.c <- length(first)
  
  low <- first - del
  #id.low <- which(low[2:(n.c)]<= first[1:(n.c-1)])
  #id.up <- which(up[1:(n.c-1)]>= first[2:(n.c)]) you do not need both sides: it is symmetric!
  id.split <-which((low[2:(n.c)]<= first[1:(n.c-1)])==FALSE)
  
  i.low.bl <-first[c(1,id.split+1)]
  i.up.bl <- first[c(id.split,n.c)]
  
  
  #sweep through the block to pick a change point.
  change.points <- c()
  for (i in 1: length(i.low.bl)){
    f <- which(ratio[i.low.bl[i]:i.up.bl[i]]==max(ratio[i.low.bl[i]:i.up.bl[i]]))
    change.points <- c(change.points, i.low.bl[i]+ f[ceiling(length(f)/2)] -1)
  }
  #remove change.points right at the beginning and at the end
  change.points<-change.points[change.points>5 & change.points<(length(ratio)-5)]
  return(change.points)
}
