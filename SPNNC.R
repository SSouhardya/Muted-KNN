distance<-function(x,y)#Utility function to calculate distance between two vectors
{
	return(norm(as.matrix(x)-as.matrix(y),type='F'))
}

di_order<-function(X,k,D){	#ets the density influence values
	n = nrow(X)
	N = R = matrix(rep(0,n**2),nrow=n)
	#N and R are the NN and RNN matrices respectively
	largest=vector(length=n)
	sigma=vector(length=n)	#denominator of the gaussian kernel
	b=vector(length=n)	#density influence values

	#forming the nearest neighbour matrix
	for(i in 1:n){
		nghs=order(D[i,])[1:k]
		largest[i]=nghs[k]
		sigma[i]=distance(X[i,],X[largest[i],])
		for(j in 1:k){
			N[i,nghs[j]]=1
		}
	}
	#N is the nearest neighbour matrix
	R = t(N) 	#surprise!
	for(i in 1:n){
		b[i]=0
		for(j in 1:n){
			if(j!=i){
				b[i] = b[i] + exp((-(D[i,j]**2))/sigma[j])*R[i,j]
			}
		}
	}
	p=order(b)
	return(p)
}


NN<-function(X,mute,K2,D){
	n=nrow(X)	#number of data points
	N=matrix(rep(0,n**2),nrow=n)	#nearest neighbour matrix
	G = N 	#the graph adjacency matrix
	N0 = N
	Dist=vector(length=n)
	marked=vector(length=n)
	index = 0 #index for handling marked 
	for(i in 1:n){
		if(i %in% mute == FALSE){
			ord=order(D[i,])
			for(j in ord){
				N[i,j]=1
				if(j %in% mute){
					index = index + 1
					marked[index] = j
					break
				}
			}#end of for
		}
	}
	for(i in 1:n){
		for(j in 1:n){
			if(i %in% mute | j %in% mute){
				G[i,j] = as.numeric(N[i,j] | N[j,i])
			}
			else{
				G[i,j] = as.numeric(N[i,j] | N[j,i])
			}
		}
	}
	marked=unique(marked[1:index])
	unmarked=mute[which(mute %in% marked == FALSE)]
	lu = length(unmarked)
	lm = length(marked)
	Dist = vector(length=lm)
	for(i in 1:lu){
		for(j in 1:lm){
			Dist[j] = D[unmarked[i],marked[j]]
		}
		closest = which(Dist == min(Dist))[1]	#closest is the index of the closest marked boundary point from the current unmarked boundary point
		N0[unmarked[i],marked[closest]]=1
	}
	Dist_mat=D[marked,unmarked]
	for(i in 1:lm){
		accept = which(Dist_mat[i,] %in% Dist_mat[i,order(Dist_mat[i,])][1:K2] == TRUE)	#accept is the set of K2 nearest neighbors of current marked boundary point
		for(j in accept){
			N0[marked[i],unmarked[j]]=1
		}
	}
	for(i in 1:n){
		for(j in 1:n){
			G[i,j] = G[i,j] + N0[i,j] * N0[j,i]
		}
	}
	for(i in 1:n){
		G[i,i]=1
	}
	return(G)
}

SPNNC<-function(X,k1=40,k2=13,alpha=0.4)
{
	n = nrow(X)
	D = as.matrix(dist(X))	#the distance matrix
	bval<-di_order(X,k1,D)
	lim=ceiling(nrow(X)*alpha)
	mute=bval[1:lim]
	G<-NN(X,mute,k2,D)
	return(G)
}