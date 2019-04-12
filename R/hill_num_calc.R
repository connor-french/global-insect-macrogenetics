## Get one hill number from a list of genetic distances. Original python code written by Isaac Overcast
hill_calc <- function(dists, order = 1) { 
  if (order == 0) {
    return(length(dists))
  }
  if (order == 1) {
    h1 = exp(entropy::entropy(dists))
    return(h1)
  }
  else {
    tot = sum(dists)
    proportions = dists/tot
    prop_order = proportions**order
    h2 = sum(prop_order)**(1/(1-order))
    return(h2)
  }
}