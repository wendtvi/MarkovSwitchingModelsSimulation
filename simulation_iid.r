#Ex 1
#################################################
#####################CASO IID####################
#################################################
#Algoritmo EM para maximização de Log Verossimilhanças (caso iid)
set.seed(43985)
n_states=2 #define numero de estados
TT=150 #define quantas observações no tempo temos de yt
y_1=rnorm(TT,0.1,1.1) #gera um processo iid para estado 1
y_2=rnorm(TT,4.1,1.1)#gera um processo iid para estado 2
y=vector()
y_state=vector()
init_state=1 #define o estado inicial do processo
y[1]=y_1[1] #seta valor inicial de y como sendo do estado 1
y_state[1]=init_state #seta estado inicial
p11=1/4
p22=3/4
p12=1-p11
p21=1-p22
P=matrix(data = c(p11,p12,p21,p22),ncol = 2)
P #define matriz de transição do proceso

#simulando dados
for (t in 2:TT){
  if(y_state[t-1]==1) ind_state=rbinom(1,1,p11) #define estado 1 como sucesso em uma dist de bernoulli
  if(y_state[t-1]==2) ind_state=rbinom(1,1,p21) #define estado 1 como sucesso em uma dist de bernoulli
  
  if(ind_state==1){y[t]=y_1[t];y_state[t]=1} #gera dados y_t
  if(ind_state==0){y[t]=y_2[t];y_state[t]=2} #gera estados
}
#plota processo orginal yt e variável latente
par(mfrow=c(3,1))
plot(y,col = ifelse(y_state == 1,'red','blue'))
lines(y)
plot(y,col = ifelse(y_state == 1,'red','blue'))
plot(y_state,col = ifelse(y_state == 1,'red','blue'))
sum(y_state==1)/TT

#estimando parametros
#sabendo que existem 2 estados e sendo que var resposta para ambos vem de uma dist normal
#inicializando probabilidades de cada estado
epsolon_tt_hat=matrix(data = NA,ncol = n_states,nrow = TT)
epsolon_tt_hat[1,]=t(c(0.5,0.5))
ni=pnorm(y,)
L=rep(-Inf,n_states)
L1=vector()
L2=vector()
state_hat=vector()

#supondo que sabemos o vetor de parâmetros real e a distribuição da variável resposta em cada estado
p1_hat=dnorm(y,mean = .1,sd = 1.1)
p2_hat=dnorm(y,mean = 4.1,sd = 1.1)
epsolon_tt_hat=matrix(NA,ncol = n_states+1,nrow = TT)
epsolon_tt_hat[,1]=p1_hat
epsolon_tt_hat[,2]=p2_hat
epsolon_tt_hat[,3]=2
y_state_hat=vector(length = TT)
epsolon_tt_hat[epsolon_tt_hat[,1]>=epsolon_tt_hat[,2],3]=1
par(mfrow=c(3,1))
par(mfrow=c(3,1))
plot(y,col = ifelse(y_state == 1,'red','blue'),main = "Plot variável resposta y")
plot(y_state,col = ifelse(y_state == 1,'red','blue'),main = "Plot dos estados da var resposta y")
plot(epsolon_tt_hat[,3],col = ifelse(epsolon_tt_hat[,3] == 1,'red','blue'),main = "Plot dos estados estimados da var resposta y")
sum(epsolon_tt_hat[,3]==y_state)/TT

#supondo que sei apenas a distribuição da var resposta y (Inferência ótima)
init_mu1=-5
init_sd1=2
init_mu2=5
init_sd2=2
pi1=sum(dnorm(y,mean = init_mu1,sd = init_sd1))/TT
pi2=sum(dnorm(y,mean = init_mu2,sd = init_sd2))/TT
init_state=1
n_max_test=100
teta_hat=matrix(NA,ncol = 3*n_states+1,nrow = n_max_test)
teta_hat[1,]=c(init_mu1,init_mu2,init_sd1,init_sd2,pi1,pi2,init_state)
epsolon_tt_hat1=vector()
epsolon_tt_hat2=vector()
L=matrix(NA,nrow = n_max_test,ncol = n_states)
L[1,]=c(0,0)

#########################################
#################smoothed inference######
#########################################
for (k in 2:n_max_test){
  epsolon_tt_hat1=dnorm(y,mean = teta_hat[k-1,1],sd = teta_hat[k-1,3])
  epsolon_tt_hat2=dnorm(y,mean = teta_hat[k-1,2],sd = teta_hat[k-1,4])
  epsolon_tt_hat1=epsolon_tt_hat1/(epsolon_tt_hat1+epsolon_tt_hat2)
  epsolon_tt_hat2=epsolon_tt_hat2/(epsolon_tt_hat1+epsolon_tt_hat2)
  teta_hat[k,1]=sum(y*epsolon_tt_hat1)/sum(epsolon_tt_hat1)
  teta_hat[k,2]=sum(y*epsolon_tt_hat2)/sum(epsolon_tt_hat2)
  teta_hat[k,3]=sum(((y-teta_hat[k,1])^2)*epsolon_tt_hat1)/sum(epsolon_tt_hat1)
  teta_hat[k,4]=sum(((y-teta_hat[k,2])^2)*epsolon_tt_hat2)/sum(epsolon_tt_hat2)
  teta_hat[k,5]=sum(epsolon_tt_hat1)/sum(epsolon_tt_hat1+epsolon_tt_hat2)
  teta_hat[k,6]=sum(epsolon_tt_hat2)/sum(epsolon_tt_hat1+epsolon_tt_hat2)
  
  L[k,1]=-sum(log(teta_hat[k,5]/sqrt(2*pi*sqrt(teta_hat[k,3])))) + sum((y-teta_hat[k,1])^2)/(2*teta_hat[k,1]^2)
  L[k,2]=-sum(log(teta_hat[k,6]/sqrt(2*pi*sqrt(teta_hat[k,4])))) + sum((y-teta_hat[k,2])^2)/(2*teta_hat[k,2]^2)
}

p1_hat=dnorm(y,mean = -0.08810256,sd = sqrt(1.0406942561287246817869345250074))
p2_hat=dnorm(y,mean = 4.227197,sd = sqrt(0.275404972))
epsolon_tt_hat=matrix(NA,ncol = n_states+1,nrow = TT)
epsolon_tt_hat[,1]=p1_hat
epsolon_tt_hat[,2]=p2_hat
epsolon_tt_hat[,3]=2
y_state_hat=vector(length = TT)
epsolon_tt_hat[epsolon_tt_hat[,1]>=epsolon_tt_hat[,2],3]=1
par(mfrow=c(3,1))
plot(y,col = ifelse(y_state == 1,'red','blue'),main = "Plot variável resposta y")
plot(y_state,col = ifelse(y_state == 1,'red','blue'),main = "Plot dos estados da var resposta y")
plot(epsolon_tt_hat[,3],col = ifelse(epsolon_tt_hat[,3] == 1,'red','blue'),main = "Plot dos estados estimados da var resposta y")
sum(epsolon_tt_hat[,3]==y_state)/TT
count_11_hat=0
count_22_hat=0
for(l in 2:TT){
  if(epsolon_tt_hat[l,3]==epsolon_tt_hat[l-1,3] && epsolon_tt_hat[l,3]==1)count_11_hat=count_11_hat+1
  if(epsolon_tt_hat[l,3]==epsolon_tt_hat[l-1,3] && epsolon_tt_hat[l,3]==2)count_22_hat=count_22_hat+1
}
p11_hat=count_11_hat/sum(epsolon_tt_hat[,3]==1)
p22_hat=count_22_hat/sum(epsolon_tt_hat[,3]==2)
P_hat=matrix(data = c(p11_hat,1-p11_hat,1-p22_hat,p22_hat),ncol = 2)

