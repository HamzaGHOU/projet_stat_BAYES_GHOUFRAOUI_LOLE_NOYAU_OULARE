#librairie 
library(mvtnorm)
library(purrr)
library(stats)
library(matlib)

## on charge les données 
y <- structure(c(151, 145, 147, 155, 135, 159, 141, 159, 177, 134, 
                 160, 143, 154, 171, 163, 160, 142, 156, 157, 152, 154, 139, 146, 
                 157, 132, 160, 169, 157, 137, 153, 199, 199, 214, 200, 188, 210, 
                 189, 201, 236, 182, 208, 188, 200, 221, 216, 207, 187, 203, 212, 
                 203, 205, 190, 191, 211, 185, 207, 216, 205, 180, 200, 246, 249, 
                 263, 237, 230, 252, 231, 248, 285, 220, 261, 220, 244, 270, 242, 
                 248, 234, 243, 259, 246, 253, 225, 229, 250, 237, 257, 261, 248, 
                 219, 244, 283, 293, 312, 272, 280, 298, 275, 297, 350, 260, 313, 
                 273, 289, 326, 281, 288, 280, 283, 307, 286, 298, 267, 272, 285, 
                 286, 303, 295, 289, 258, 286, 320, 354, 328, 297, 323, 331, 305, 
                 338, 376, 296, 352, 314, 325, 358, 312, 324, 316, 317, 336, 321, 
                 334, 302, 302, 323, 331, 345, 333, 316, 291, 324), .Dim = c(30, 
                                                                             5))

x <- c(8.0, 15.0, 22.0, 29.0, 36.0)

##graph intro 

matplot(c(x), t(y), type = 'b', pch=1, ylab = 'Poids', xlab = 'Jour X_j', xaxp=c(8,36,4))

library(ggplot2)
library(tidyr)
data_long <- data.frame(x = rep(x, nrow(y)), y = c(t(y)), Rat = factor(rep(1:nrow(y), each = length(x))))

palette_couleurs <- rainbow(length(unique(data_long$Rat)))

ggplot(data_long, aes(x = x, y = y, group = Rat, color = Rat, linetype = Rat)) +
  geom_path() +
  geom_point(shape = 1) +
  labs(x = "Jour X_j", y = "Poids") 
scale_color_manual(values = palette_couleurs) +  # Utiliser une palette de couleurs
  scale_linetype_manual(values = rep("dashed", length(unique(data_long$Rat))))  # Utiliser des lignes pointillées

# = 1/R (matrice de cov , R=matrice de précision)
R_inv <- structure(c(0.005, 0, 0, 5), .Dim = c(2, 2)) 

#par init
sp = solve(R_inv)
beta <- structure(c(100,100,100,100,100,
                    100,100,100,100,100,
                    100,100,100,100,100,
                    100,100,100,100,100,
                    100,100,100,100,100,
                    100,100,100,100,100,
                    6,6,6,6,6,
                    6,6,6,6,6,
                    6,6,6,6,6,
                    6,6,6,6,6,
                    6,6,6,6,6,
                    6,6,6,6,6),
                  .Dim = c(30, 2))

mu0 <- c(0,0)
S0 <- matrix(c(0.0001,0,0,0.0001),nrow=2,ncol=2) # mat préc de mu_B

###MCMC -> ech de gibbs 

birats_gibbs <-  function(N,a0,b0,mu0,S0,R,beta,x,y){
  n = dim(y)[1]
  # To return
  mu_B = matrix(NaN,2*N,nrow=N,ncol=2,dimnames = list(c(),c("mu_B1","mu_B2"))) #dim = (Nchain x 2)
  tau_c = matrix(NaN,N,nrow=N,ncol=1,dimnames = list(c(),c("tau_c"))) #dim = (Nchain x 1)
  Beta_1 = matrix(NaN,nrow=30,ncol=N)
  Beta_2 = matrix(NaN,nrow=30,ncol=N)
  cors = matrix(NaN,nrow=1,ncol=N)
  
  # Initialisation
  mu_B[1,1:2] = c(0,0)
  tau_c[1] = 1
  Omega =  matrix(c(1,0,0,1),nrow = 2, ncol = 2) #precision matrix
  Beta = beta # dim = (30 x 2)
  
  Beta_1[,1] = Beta[,1]
  Beta_2[,1] = Beta[,2]
  
  cors[1] = 1
  
  for(t in 2:N){
    cat("t= ", t,"\n")
    
    # On update mu_B
    # Prec_post = diag(diag(n*Omega + S0)) # dim = (2x2)
    # Mu_post = solve(Prec_post)%*%(n*Omega%*%colMeans(Beta)) # + S0%*%matrix(mu0,nrow=2,ncol=1) = 0
    #
    # mu_B[t,1:2] = rmvnorm(1,mean=Mu_post,sigma = solve(Prec_post))
    
    # On update mu_B1
    var0 = 1/S0[1,1]
    
    s_star_1 = 1/(1/var0 + n*Omega[1,1])
    mu_star_1 = s_star_1 * (Omega[1,1]*sum(Beta[,1]) + Omega[1,2]*sum(Beta[,2]) -
                              n/2*mu_B[t-1,2]*(Omega[2,1] + Omega[1,2]))
    
    mu_B[t,1] = rnorm(1, mu_star_1, sqrt(s_star_1))
    
    # On update mu_B2
    
    s_star_2 = 1/(1/var0 + n*Omega[2,2])
    mu_star_2 = s_star_2 * (Omega[2,1]*sum(Beta[,1]) + Omega[2,2]*sum(Beta[,2]) -
                              n/2*mu_B[t,1]*(Omega[2,1] + Omega[1,2]))
    
    mu_B[t,2] = rnorm(1, mu_star_2, sqrt(s_star_2))
    
    # On update Omega
    somme = matrix(c(0,0,0,0),nrow=2,ncol=2)
    
    for (i in 1:30){
      Beta_i_moins_mu = matrix(Beta[i,1:2] - mu_B[t,1:2],nrow=2,ncol=1) #dim = (2,1)
      somme = somme + Beta_i_moins_mu%*%t(Beta_i_moins_mu)
    }
    
    Omega = rWishart(1,Sigma = solve(somme+solve(R)),df = n+2)[,,1]
    
    # Calcul de la corrélation entre Beta_1 et Beta_2
    SIGMA = solve(Omega)
    cors[t] = SIGMA[1,2]/sqrt(SIGMA[1,1]*SIGMA[2,2])
    
    # On update tau_c
    mu = matrix(NaN,30*5,nrow=30,ncol=5) #mu_ij de dim = (30,5)
    y_moins_mu_2 = matrix(NaN,30*5,nrow=30,ncol=5) # = (y_ij - mu_ij)^2
    
    for(i in 1:30){
      for (j in 1:5){
        mu[i,j] = Beta[i,1] + Beta[i,2]*x[j]
      }
    }
    y_moins_mu_2 = (y - mu)^2
    
    tau_c[t] = rgamma(1,shape = a0 + 5*n/2, rate = b0 + sum(y_moins_mu_2)/2) 
    
    #On update les betas
    for (i in 1:30){
      somme_Beta_A_i = matrix(c(0,0,0,0),nrow=2,ncol=2)
      somme_Beta_m_i = matrix(c(0,0),nrow=2,ncol=1)
      
      for (j in 1:5){
        Cj = matrix(c(1,x[j]),nrow=2,ncol=1)
        somme_Beta_A_i = somme_Beta_A_i + Cj%*%t(Cj)
        
        somme_Beta_m_i =  somme_Beta_m_i + y[i,j]*Cj
      }
      
      
      A = Omega + tau_c[t]*somme_Beta_A_i
      m = Omega%*%mu_B[t,1:2] + tau_c[t]*somme_Beta_m_i
      
      Beta[i,1:2] = rmvnorm(1,mean = solve(A)%*%m , sigma = solve(A))
    }
    Beta_1[,t] = Beta[,1]
    Beta_2[,t] = Beta[,2]
  }
  return(list(mu_B,tau_c,Beta_1,Beta_2,cors))
}

##On Lance la fonction avec 1100 itérations on stock les differents param 
res = birats_gibbs(N = 11000,a0=0.001,b0=0.001,mu0 = mu0,S0 = S0,R=R_inv,beta=beta,x=x,y=y)
mu_1 = res[[1]][,1] 
mu_2 = res[[1]][,2] 
tau_c = res[[2]]
Beta_1 = res[[3]] 
Beta_2 = res[[4]] 
cors = res[[5]][1000:length(res[[5]])] 

# Calcul des statistiques
moyennes <- c(mean(mu_1), mean(mu_2), mean(tau_c))
ecart_types <- c(sqrt(var(mu_1)), sqrt(var(mu_2)), sqrt(var(1/sqrt(tau_c))))

tableau <- data.frame(
  Paramètre = c("mu_1", "mu_2", "tau_c"),
  `Moyenne empirique` = moyennes,
  `Ecart-type empirique` = ecart_types
)

print(tableau)

# Affichage graphique des chaînes obtenues par le sampler : 
par(mfrow = c(3,2))
# Chaînes des Beta_1 pour chaque rat
plot(1000:11000,Beta_1[1,1000:11000],type="l",xlab='itérations',ylab='beta_1')
abline(h=106.6,col="red")
# Chaînes des Beta_2 pour chaque rat
plot(1000:11000,Beta_2[1,1000:11000],type="l",xlab='itérations',ylab='beta_2')
abline(h=6.185,col="red")
# Chaine pour mu_Beta_1
plot(1000:11000,mu_1[1000:11000],type="l",xlab='itérations',ylab='mu_beta_1')
abline(h=106.6,col="red")
abline(h=mean(mu_1),col="blue")
#chaine pour mu_Beta_2
plot(1000:11000,mu_2[1000:11000],type="l",xlab='itérations',ylab='mu_beta_2')
abline(h=6.185,col="red")
abline(h=mean(mu_2),col="blue")
# chaine pour tau_c
plot(1000:11000,tau_c[1000:11000],type="l",xlab='itérations',ylab='tau_c')
abline(h=mean(tau_c),col="blue")
abline(h=1/(6.136^2),col="red")

# Lien linéaire entre Beta_2 et Beta 1 , on utulise cors : 

# Intervalle de crédibilité de de niveau 95%  : 

len = c()
I = seq(0.001,0.05,0.001)
for (i in 1:length(B)){
  QI = quantile(cors[1000:length(cors)],probs = I[i])
  Q1mapI = quantile(cors[1000:length(cors)],probs = 0.95+I[i])
  len[i] = Q1mapI - QI
}

I_min = I[which.min(len)] # On récupère I tel que [q_I  ; Q_I+1-alpha] DE longueur minimale 
Q_I = quantile(cors[1000:length(cors)],probs = I_min)
Q_1mapI = quantile(cors[1000:length(cors)],probs = 0.95+I_min)

paste0("Intervalle de crédibilité de longeur minimale [",round(Q_I,3)," ; ",round(Q_1mapI,3),"]")

# Densité de la loi des corrélations + Intervalle de crédibilité pour la corrélation
plot(density(cors[1000:length(cors)]),type="l",main = "densité corrélation entre Beta_1 et Beta_2")
abline(v=c(round(Q_I,3),round(Q_1mapI,3)),col="blue",ylim = c(0,1))
legend("topright",legend = c("Densité corrélations","I-C"),
       pch=c("-","-"),col=c("black","blue"),bg = "white")

