DFS_graph<-function(X)
{
	#X is the companion matrix of the graph
	G<-X
	n<-dim(X)[1]
	flag=0
	curr=0
	change=0
	compo<-NULL
	compoL<-NULL
	stack<-NULL
	temp<-NULL
	repeat{
		i=1
		for(i in 1:n){
			if(max(G[i,])!=0){
				curr=i
				flag=1
				break
			}
			flag=0
		}
		if(flag==0){
			break
		}
		sub_poin=curr
		repeat{
			i=1
			for(i in 1:n){
				if(G[sub_poin,i]==1){
					stack<-c(stack,i)
					change=1
					G[sub_poin,i]=0
					sub_poin=i
					break	#breaks from for
				}
					change=0
			}#end of for
			if(change==0){
				temp<-c(temp,stack[length(stack)])
				stack<-stack[-length(stack)]	#stack's length decreases by 1
				sub_poin=stack[length(stack)]
			}
			if(length(stack)==0){
				temp<-unique(temp)
				compo<-c(compo,temp)
				compoL<-c(compoL,length(temp))
				temp<-NULL
				break	#breaks from inner repeat
			}
		}#end of inner repeat
	}#end of outer repeat
	toreturn<-list(compo,compoL)
	return(toreturn)
}

labels<-function(R,p){
	x<-p[[1]]
	y<-p[[2]]
	n=length(y)
	s<-NULL
	i=1
	label=rep(1,nrow(R))
	for(i in 2:n){
		low<-sum(y[1:(i-1)])+1
		high<-sum(y[1:i])
		s<-x[low:high]
		label[s]=i
	}
	return(label)
}

master.SPNNC<-function(X,G){
	p<-DFS_graph(G)
	lab = labels(X,p)
	return(lab)
}

