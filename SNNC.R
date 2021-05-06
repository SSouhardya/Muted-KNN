
distance<-function(x,y){
  return(sqrt(sum((x-y)**2)))
}


set_border_points<-function(D,k,alpha,border_times,distance_rank){
  #function to identify the border points

  di_values<-function(){	#gets the density influence values
    #N is the nearest neighbor matrix. But not only 0/1 but infact the ranks of the first k elements
    f<-function(vec){
      return(which(vec == sort(vec)[k+1])[1])
    }
    largest = apply(distance_rank,1,f) #the index of the kth nearest neighbor-vector for each data point
    sigma = diag(D[,largest])#the denominator of the gaussian kernel

    h<-function(vec){
      return(as.numeric(vec > 0 & vec <= sort(vec)[k+1]))
    }
  	R = apply(distance_rank,1,h) 	#RNN matrix with entries 0/1, note that the transpose has already been taken care of
    diag(R) = 0
    assign_di<-function(i){
      #a function to assign density influence values to the ith index
      return(sum((exp(-(D[i,]**2)/sigma))*R[i,]))
    }
    b = sapply(1:nrow(D), assign_di)  #density influence values
  	return(b)
    #return(b)
  }

  n = nrow(D)
  #t() is required as apply stacks as columns
  f<-function(vec){
    vec[-which(vec %in% 0:k)] = 0
    return(vec)
  }

  curr_indices = 1:n
  border_indices = NULL
  number_of_points = n
  for(iter in 1:border_times){
    #N = t(apply(distance_rank,1,f))
    #nearest neighbor matrix
    j = order(di_values())[1:ceiling(alpha*number_of_points)]
    b = curr_indices[j]
    number_of_points = number_of_points - length(b)
    border_indices = c(border_indices, b)
    distance_rank = distance_rank[-j,][,-j]
    curr_indices = setdiff(curr_indices,b)
    D = D[-j,][,-j]
  }
  return(border_indices)
}



form_graph<-function(X, border, k, D, N, mode, cutoff){
  #N is the nearest neighbor rank matrix as used earlier
  n = nrow(X)
  d = ncol(X)
  G = matrix(0,nrow = n, ncol = n)

  core = setdiff(1:n,border)
  #print(length(border))

  vec = rep(0,n)  #will be used in the function below
  f<-function(i){
    #i is the index of a core point
    vec[order(D[i,])[1:(min(N[i,border])+1)]] = 1
    return(c(vec, order(D[i,])[min(N[i,border])+1]))  #the last column contains the index of the first hit border point
  }
  Gcore = t(sapply(core,f))
  #print(dim(Gcore))
  hit_by_core = unique(Gcore[,(n+1)]) #contains 1 for those border points that have been hit by core
  G[core,] = Gcore[,-(n+1)]
  #diag(G) = 1
  #return(G)
  #step 1 done

  not_hit_by_core = setdiff(border, hit_by_core)

  m = length(hit_by_core)
  vec = rep(0,m)
  g<-function(i){ #returns the index of the hit_by_core point
    #i is the index of a not_hit_by_core point
    vec[order(N[i,hit_by_core])[1]]=1
    return(c(vec,hit_by_core[order(N[i,hit_by_core])[1]]))
  }

  Gcore = t(sapply(not_hit_by_core,g))  #reusing Gcore, as not needed anymore

  nearest_hit_by_core = Gcore[,(m+1)]
  G[hit_by_core,not_hit_by_core] = t(Gcore[,-(m+1)])


  if(mode == 1){
    vec = rep(0,length(not_hit_by_core))
    f<-function(i){
      #i is a memeber of nearest_hit_by_core
      index_core = order(D[i,])[2:(k+1)]
      index_totest = not_hit_by_core[which(G[i,not_hit_by_core] == 1)]
      index_total = c(index_totest,index_core)
      z = prcomp(X[index_total,])$x[,1:2]
      distances = mahalanobis(x = z[1:length(index_totest),], center = apply(z,2,mean), cov = cov(z))
      #distances = k*((z[1:length(index_totest)] - m)/sd)^2
      vec[which(G[i,not_hit_by_core] == 1)] = as.numeric((distances < qchisq(0.99,df = 2)))
      return(vec)
    }
      G[nearest_hit_by_core,not_hit_by_core] = t(sapply(nearest_hit_by_core,f))
  }

  if(mode == 2){
    if(is.null(cutoff)){
      g<-function(i){
          v = var(sort(D[i,])[2:(k+1)])*(k-1)
          m = mean(sort(D[i,])[2:(k+1)])
          index_totest = not_hit_by_core[which(G[i,not_hit_by_core] == 1)]
          return(((D[i,index_totest] - m)^2/v)*(k/(k+1)))
      }

      vec = sapply(nearest_hit_by_core,g)
      l = length(vec)
      var_ratio = NULL
      for(i in 1:l){
        var_ratio = c(var_ratio,vec[[i]])
      }
      plot(sort(var_ratio), pch = 19, type = 'l', ylab = 'DSS', lwd = 2)
      return(0)
    }
    vec = rep(0,length(not_hit_by_core))
    g<-function(i){
        v = var(sort(D[i,])[2:(k+1)])*(k-1)
        m = mean(sort(D[i,])[2:(k+1)])
        index_totest = not_hit_by_core[which(G[i,not_hit_by_core] == 1)]
        distances = ((D[i,index_totest] - m)^2/v)*(k/(k+1))
        vec[which(G[i,not_hit_by_core] == 1)] = as.numeric(distances < cutoff)
        return(vec)
    }
      G[nearest_hit_by_core,not_hit_by_core] = t(sapply(nearest_hit_by_core,g))
  }
  diag(G)=1
  return(G)
}




SNNC<-function(X,k=20,alpha=0.1,border_times=5, mode = 1, cutoff = NULL)
{
  #X is the data matrix
  #k1 is the number of nearest neighbors to calculate for obtaining the di values
  #k2 is for the revision step
  #alpha is the peeling percentage at each iteration
  require(igraph)
	n = nrow(X)
	D = as.matrix(dist(X))	#the distance matrix))
  give_rank<-function(vec){
    return(rank(vec,ties.method = 'random')-1)
  }
  distance_rank = t(apply(D,2,give_rank))  #each row contains rankings. Diagonals are 0
	mute = set_border_points(D,k,alpha,border_times,distance_rank)
	G<-form_graph(X,mute,k = k, D,distance_rank, mode, cutoff)
  if((mode == 2) & is.null(cutoff)){
    return()
  }
	gr = graph_from_adjacency_matrix(G,mode='undirected',diag=FALSE)
	l = components(gr)
	labels = l$membership
	return(labels)
}

