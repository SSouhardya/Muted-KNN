distance<-function(x,y)#Utility function to calculate distance between two vectors
{
return(norm(as.matrix(x)-as.matrix(y),type='F'))
}

setval<-function(X,k){
	n=nrow(X)
	N=R=matrix(rep(0,n**2),nrow=n)
	#N and R are the NN and RNN matrices respectively
	dist=vector(length=n)
	largest=vector(length=n)
	b=vector(length=n)	#density influence values

	#forming the nearest neighbour matrix
	for(i in 1:n){
		z=X[i,]
		for(j in 1:n){
			dist[j]=distance(z,X[j,])
		}
		nghs=order(dist)[1:k]
		largest[i]=nghs[k]
		for(j in 1:k){
			N[i,nghs[j]]=1
		}
	}
	#N is the nearest neighbour matrix

	R = t(N) 	#surprise!
	sigma=vector(length=n)	#denominator of the gaussian kernel
	for(i in 1:n){
		sigma[i]=distance(X[i,],X[largest[i],])
	}
	for(i in 1:n){
		b[i]=0
		for(j in 1:n){
			if(j!=i){
				b[i] = b[i] + exp(-(distance(X[i,],X[j,])**2)/sigma[j])*R[i,j]
			}
		}
	}
	p=order(b)
	return(p)
}


opt.kNN<-function(X,mute,K2){
	n=nrow(X)	#number of data points
	N=matrix(rep(0,n**2),nrow=n)	#nearest neighbour matrix
	G = N 	#the graph adjacency matrix
	N1 = N
	N0 = N
	dist=vector(length=n)
	marked=NULL
	for(i in 1:n){
		if(i %in% mute == FALSE){
			for(j in 1:n){
				dist[j]=distance(X[i,],X[j,])
			}
			ord=order(dist)
			for(j in ord){
				N[i,j]=1
				if(j %in% mute){
					marked = c(marked,j)
					break
				}
			}#end of for
		}
	}
	for(i in 1:n){
		for(j in 1:n){
			if(i %in% mute | j %in% mute){
				G[i,j]=N[i,j]+N[j,i]
			}
			else{
				G[i,j]=as.numeric(N[i,j] | N[j,i])
			}
		}
	}
	marked=unique(marked)
	unmarked=mute[which(mute %in% marked == FALSE)]
	lu = length(unmarked)
	lm = length(marked)
	dist = vector(length=lm)
	for(i in 1:lu){
		for(j in 1:lm){
			dist[j] = distance(X[unmarked[i],],X[marked[j],])
		}
		help = which(dist == min(dist))[1]
		N0[unmarked[i],marked[help]]=1
	}
	dist.mat=matrix(0,nrow=lm,ncol=lu)
	for(i in 1:lm){
		for(j in 1:lu){
			dist.mat[i,j] = distance(X[marked[i],],X[unmarked[j],])
		}
	}
	for( i in 1:lm){
		dist=dist.mat[i,]
		send=which(dist %in% dist[order(dist)][1:K2] == TRUE)
		for(j in send){
			N0[marked[i],unmarked[j]]=1
		}
	}
	for(i in 1:n){
		for(j in 1:n){
			N1[i,j]=N0[i,j]*N0[j,i]
		}
	}
	G = G + N1
	for(i in 1:n){
		G[i,i]=1
	}
	return(G)
}

SPNNC<-function(X,k1=40,k2=13,alpha=0.4)
{
	bval<-setval(X,k1)
	lim=ceiling(nrow(X)*alpha)
	mute=bval[1:lim]
	G<-opt.kNN(X,mute,k2)
	return(G)
}
