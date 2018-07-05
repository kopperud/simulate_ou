## Simulate OU process


library(slouch)
library(hashmap)
library(ape)

data(neocortex)
data(artiodactyla)


concat.factor <- function(...){
  as.factor(do.call(c, lapply(list(...), as.character)))
}


dW <- function(s){
  x <- rnorm(1, mean = 0, sd = s)
  return(x)
}



## dy = alpha * (theta - y) + s*dW
dy <- function(y, theta, alpha, s, infinitesimal){
  res <- alpha *(theta - y)*infinitesimal + dW(s)
  return(res)
}


simulate_edge <- function(edge, y, regime, depths, theta, alpha, s, infinitesimal){
  steps <- round(edge$edge.length / infinitesimal)
  theta <- theta[[regime]]
  
  timeseries <- list(rep(0, steps))
  timeseries[1] <- y + dy(y, theta, alpha, s, infinitesimal)
  for (i in seq(2, steps)){
    y <- y + dy(y, theta, alpha, s, infinitesimal)
    timeseries[i] <- y
  }
  timeseries <- as.numeric(timeseries)
  
  x = seq(0, edge$edge.length, length.out = steps) + depths[edge$ancestor] #treelength
  df <- list(timeseries = timeseries,
                   x = x,
                   regime = regime)
  
  return(df)
}

traverse <- function(node, y, phy, edges, env, depths, theta, alpha, s, infinitesimal){
  is_descendant = edges[,1] == node
  descendant_edges <- edges[is_descendant, ]
  ndescendants <- nrow(descendant_edges)
  
  y_terminal = list(rep(NA, ndescendants))
  for (i in seq(ndescendants)){
    edge <- descendant_edges[i, ]
    #edge_idx <- #simdata$keys()[i]
    edge_idx <- as.numeric(rownames(edge))  
    regime <- edge$regimes
    
    simulate <- simulate_edge(edge, y, regime, depths, theta, alpha, s, infinitesimal)
    y_terminal[[i]] <- tail(simulate$timeseries, n = 1)
    
    env$simdata[[edge_idx]] <- simulate
  }
  
  descendant_nodes <- descendant_edges$descendant
  is_internal <- descendant_nodes > length(phy$tip.label)
  
  if (length(descendant_nodes) > 0){
    mapply(traverse, 
                  descendant_nodes[is_internal],
                  y_terminal[is_internal],
                  MoreArgs = list(phy = phy,
                                  edges = edges, 
                                  env = env,
                                  depths = depths,
                                  theta = theta, 
                                  alpha = alpha, 
                                  s = s, 
                                  infinitesimal = infinitesimal))
    #res <- do.call(c, res)
    #return(res)
  }else{
    return(NULL)
  }
}



simulate_ou <- function(phy,
                        regimes,
                        infinitesimal,
                        Ya = 2.0,
                        theta = NULL,
                        alpha = 0.5,
                        s = 1.2){
  if(missing(infinitesimal)){
    infinitesimal <- max(node.depth.edgelength(phy))/1000
  }
  
  stopifnot(is.rooted(phy))
  
  edges <- data.frame(phy$edge, 
                      "edge.length" = phy$edge.length, 
                      "regimes" = regimes[phy$edge[,2]])
  
  colnames(edges)[1:2] <- c("ancestor", "descendant")
  nedges <- nrow(phy$edge)
  depths <- node.depth.edgelength(phy)
  treelength <- max(depths)
  

  
  rootnode <- length(phy$tip.label) + 1
  env <- environment()
  env$simdata <- list()
  
  traverse(rootnode, y = Ya, phy, edges, env, depths, theta, alpha, s, infinitesimal)
  
  return(env$simdata)
}

.addlines <- function(i, simdata, colmap){
  x <- simdata[[i]]
  lines(x$x, x$timeseries, col = colmap[[x$regime]], lwd = 1.3)
  
  #text(tail(x$x, n = 1), tail(x$timeseries, n = 1), i)
  #text(head(x$x, n = 1), head(x$timeseries, n = 1), i)
  return(NULL)
}

ou_plot <- function(phy, simdata, colmap){
  treelength <- max(node.depth.edgelength(phy))
  
  plot(0, type = "n", 
       xlim = range(do.call(c, lapply(simdata, function(e) range(e$x, na.rm = TRUE))), na.rm = TRUE), 
       ylim = range(do.call(c, lapply(simdata, function(e) range(e$timeseries, na.rm = TRUE))), na.rm = TRUE),
       xlab = "time",
       ylab = "y")
  for (i in seq(simdata)){
    .addlines(i, simdata, colmap)
  }
}

theta <- hashmap(c("MF", "Gr", "Br"), c(1, 8, 3))

regimes <- concat.factor(neocortex$diet, artiodactyla$node.label)

simdata <- lapply(1:50, function(e) simulate_ou(artiodactyla, 
                                          regimes, 
                                          theta = theta,
                                          alpha = 0.05,
                                          s = 0.01,
                                          infinitesimal = 0.001))



phy <- reorder.phylo(artiodactyla, order = "postorder")
par(mfrow = c(1,1))
plot(phy); edgelabels()
rTraitCont(phy, model = "OU",
           theta = c(1, 8, 3)[as.numeric(regimes[phy$edge[,2]])],
           alpha = 0.3,
           sigma = 0.5)








dat <- data.frame(artiodactyla$edge)
dat <- dat[(dat[,2] %in% 1:43),]
dat <- dat[order(dat$X2, rownames(dat)), ]

bar <- function(simdata){
  d <- simdata[as.numeric(rownames(dat))]
  
  y = sapply(d, function(x) tail(x$timeseries, n = 1))
  mv <- rnorm(length(y), 0, sd = 0.3*sd(y))
  predictions <- data.frame(species = artiodactyla$tip.label,
                            y = y + mv,
                            diet = sapply(d, function(x) x$regime))
  
  m1 <- slouch.fit(artiodactyla,
                   species = predictions$species,
                   response = predictions$y,
                   mv.response = (0.15)^2*mv^2,
                   fixed.fact = predictions$diet,
                   hillclimb = TRUE)
  return(m1)
}

m <- lapply(simdata, bar)



par(mfrow = c(2, 3))

# Optima 
coefs <- data.frame(t(sapply(m, function(e) e$beta_primary$coefficients[,1])))

hist(coefs$Br, main = "Browser", xlab = "Optimum", breaks = 15); abline(v = 3, col = "red", lwd = 3)
hist(coefs$Gr, main = "Grazer", xlab = "Optimum", breaks = 15); abline(v = 8, col = "red", lwd = 3)
hist(coefs$MF, main = "Mixed feeder", xlab = "Optimum", breaks = 15); abline(v = 1, col = "red", lwd = 3)

# OU parameters
evolpar <- data.frame(t(sapply(m, function(e) e$evolpar)))

hist(sapply(evolpar$hl, function(e) e), breaks = 15, main = "Half-life", xlab = "Time (myr)")
abline(v = log(2)/0.25, col = "red", lwd = 3)

hist(sapply(evolpar$vy, function(e) e), breaks = 15, xlim = c(-0.2, 1), main = "Stationary variance", xlab = "Trait squared")
abline(v = (0.01^2) / (2*0.25), col = "red", lwd = 3)




colmap <- hashmap(theta$keys(), c("orange", "black", "#2ca25f"))

par(mfrow = c(1,1))
#ou_plot(artiodactyla, simdata[[1]],  colmap)




