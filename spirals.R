rm(list=ls())
library(rgl)
library(mlbench)

#Calculo das verossimilhanças
pdfnvar<- function(x,m,K,n) ((1/(sqrt((2*pi)^n*(det(K)))))*exp(-0.5*(t(x-m) %*% (solve(K)) %*% (x-m))))

#FUNÇÃO K-MEANS
mykmeans <- function(x,k)
{
  
  #Parâmetros:
  #x: matriz de n dados de c dimensões
  #k: nº de clusters
  
  #Número de dimensões e dados
  c=ncol(x)
  n=nrow(x)
  
  #Inicializando os centróides aleatóriamente
  u=matrix(0,nrow=k,ncol=n)
  ix=sample(n)
  centroids=x[ix[1:k],]
  
  #Variável de parada, contagem de iterações e vetor de função de custo
  change=1
  iter=0
  J=c()
  
  while(change)
  {
    
    iter=iter+1
    
    #Cálculo da nova matriz de pertinência 
    u.old=u
    u=matrix(0,nrow=k,ncol=n)
    aux=c()
    Jaux=0
    for(i in 1:n)
    {
      for(j in 1:k)
      {
        aux[j]=t(centroids[j,]-x[i,])%*%(centroids[j,]-x[i,])
      }
      minimum.idx=which.min(aux)
      u[minimum.idx,i]=1
    }
    
    #Cálculo da função de custo
    dist=0
    for(j in 1:k)
    {
      idx=which(u[j,]==1)
      for(i in 1:(length(idx)))
      {
        dist=dist+t(centroids[j,]-x[idx[i],])%*%(centroids[j,]-x[idx[i],])
      }
    }
    J[iter]=dist
    
    #Cálculo dos centróides
    for(i in 1:k)
    {
      centroids[i,]=colMeans(x[which(u[i,]==1),])
    }
    
    #Teste do critério de parada: matriz de pertinêncial igual à antiga
    change=sum(1*(u!=u.old))
    
  }
  
  out=list(C=centroids,U=u,J=J,it=iter)
  return(out)
}

#Cálculo dos parâmetros das distribuições
mixofgaussians <- function(x,k)
{
  #Parâmetros:
  #x: lista de dados, com d classes de n dados de c dimensões
  #k: vetor de nº de gaussianas por classe

  #Número de classes
  d=length(x)

  outofkmeans=list()
  n=c()
  out=list()
  for(i in 1:d)
  {
    #Número de dados em cada classe
    n=nrow(x[[i]])
    
    #Chamada do cluster
    outofkmeans=mykmeans(x[[i]],k[i])
    U=outofkmeans[[2]]
    
    #Cálculo dos parâmetros
    Klist=list()
    mlist=list()
    plist=list()
    for(j in 1:k[i])
    {
      idx1=which(U[j,]=='1')
      Klist[[j]]=cov(x[[i]][idx1,])
      mlist[[j]]=as.matrix(colMeans(x[[i]][idx1,]))
      plist[[j]]=nrow(x[[i]][idx1,])/nrow(x[[i]])
    }
    out[[i]]=list(Klist,mlist,plist)
  }
  return(out)
}

#Cálculo da verossimilhança
systempdf <- function(x,w)
{
  #Parâmetros
  #x: lista de dados, com d classes de n dados de c dimensões
  #w: lista de parâmetros
  
  #Número de classes
  d=length(w)
  
  px<-list()
  for(i in 1:d)
  {
    #Número de dados da classe
    n=nrow(x)
    c=ncol(x)
    
    #Cálculo das verossimilhanças
    pxC=c()
    for(j in 1:n)
    {
      aux=0
      k=length(w[[i]][[i]])
      for(l in 1:k)
      {
        aux=aux+w[[i]][[3]][[l]]*pdfnvar(x[j,],w[[i]][[2]][[l]],
                                          w[[i]][[1]][[l]],c)
      }
      pxC[j]=aux
    }
  px[[i]]=pxC
  }
  return(px)
}

#Definição do problema
N<-1500
spirals<-mlbench.spirals(N,cycles=1.5,sd=0.05)

ic1<-which(spirals$classes=='1')
ic2<-which(spirals$classes=='2')

x<-unlist(spirals[[1]])
plot(x[ic1,1],x[ic1,2],col='magenta',xlim=c(-1.3,1.3),ylim=c(-1.3,1.3),
     xlab='x1',ylab='x2')
par(new=T)
plot(x[ic2,1],x[ic2,2],col='blue',xlim=c(-1.3,1.3),ylim=c(-1.3,1.3),
     xlab='',ylab='',main='Distribuição dos dados')

#Divisão das classes
x<-list(x[ic1,],x[ic2,])
k<-c(12,12)       
w<-mixofgaussians(x,k)

#spacemap
seqi<-seq(-1.3,1.3,0.1)
seqj<-seq(-1.3,1.3,0.1)
M1<-matrix(0,nrow=length(seqi),ncol=length(seqj))
M2<-matrix(0,nrow=length(seqi),ncol=length(seqj))
ci<-0
for(i in seqi)
{
  ci<-ci+1
  cj<-0
  for(j in seqj)
  {
    cj<-cj+1
    xt<-t(as.matrix(c(i,j)))
    aux<-systempdf(xt,w)
    M1[ci,cj]<-aux[[1]]
    M2[ci,cj]<-aux[[2]]
  }
}

par(new=T)
contour(seqi,seqj,M1,xlim=c(-1.3,1.3),ylim=c(-1.3,1.3),col='magenta')
par(new=T)
contour(seqi,seqj,M2,xlim=c(-1.3,1.3),ylim=c(-1.3,1.3),col='blue')

M<-1*((M1*nrow(x[[1]])/nrow(x[[1]]+x[[2]]))
      >=(M2*nrow(x[[2]])/nrow(x[[1]]+x[[2]])))

plot(x[[1]][,1],x[[1]][,2],col='magenta',xlim=c(-1.3,1.3),ylim=c(-1.3,1.3),
     xlab='',ylab='')
par(new=T)
plot(x[[2]][,1],x[[2]][,2],col='blue',xlim=c(-1.3,1.3),ylim=c(-1.3,1.3),
     xlab='x1',ylab='x2',main='Função de separação')
par(new=T)
contour(seqi,seqj,M,xlim=c(-1.3,1.3),ylim=c(-1.3,1.3),
        nlevels=1,drawlabels = FALSE)

persp3d(seqi,seqj,M1,zlim=c(0,3),col='blue')
persp3d(seqi,seqj,M2,zlim=c(0,3),col='magenta',add=T)
persp3d(seqi,seqj,M,zlim=c(0,3),col='grey',add=T)


# n<-4
# 
# N<-30
# g1<-matrix(rnorm(N),ncol = 2)*0.8 + matrix(c(2,2),N/2,2,byrow=T)
# g2<-matrix(rnorm(N),ncol = 2)*0.8 + matrix(c(2,4),N/2,2,byrow=T)
# g3<-matrix(rnorm(N),ncol = 2)*0.8 + matrix(c(4,4),N/2,2,byrow=T)
# g4<-matrix(rnorm(N),ncol = 2)*0.8 + matrix(c(4,2),N/2,2,byrow=T)
# 
# plot(g1[,1],g1[,2],col='blue',xlim=c(0,6),ylim=c(0,6))
# par(new=T)
# plot(g2[,1],g2[,2],col='red',xlim=c(0,6),ylim=c(0,6))
# par(new=T)
# plot(g3[,1],g3[,2],col='magenta',xlim=c(0,6),ylim=c(0,6))
# par(new=T)
# plot(g4[,1],g4[,2],col='black',xlim=c(0,6),ylim=c(0,6))
# 
# x<-rbind(g1,g2,g3,g4)
