rm(list = ls())
graphics.off()
setwd('~')
setwd('Documents/ENSTA_2021/STA203/Projet')


## Importation des données
load('cookie.app.RData')
load('cookie.val.RData')


#################################
## Analyse exploratoire
#################################

#########
## Q1
#########

#Lecture des données

xtrain <- cookie.app[-1]
head(xtrain)
ytrain <- cookie.app[1]
head(ytrain)

xtest <- cookie.val[-1]
head(xtrain)
ytest <- cookie.val[1]
head(ytrain)

#Etude uni-varié
boxplot(xtrain)


# Affichage des spectres
sapply( 1:nrow(xtrain), function(i){
                    matplot(1:700, t(xtrain[i,]), type = 'l',
                            main = paste("Spectre de l'individu: ", i," .")
                            , xlab="Fréquence dans le proche infra-rouge"
                            , ylab="Absorbance") })

#Corrélation entre les mesures aux différentes fréquences
require(corrplot)
corrplot(cor(xtrain),is.corr = F)
C <- cor(xtrain)

#########
## Q2
#########
#ACP
require(FactoMineR)
res.acp <- PCA(xtrain, ncp = 700,graph=F)

#Graphe des valeurs propres
barplot(res.acp$eig[,1] / sum(res.acp$eig[,1]) * 100, las = 2)
abline(h = 100 / 700, col = "red")
which(res.acp$eig[,3] > 100 - 1/7)
## Au nombre de 39 -> Prblm de dégénérescence du rang de X

# Tracer des nuages sur les six premiers axes principaux
par(mfrow = c(1, 2))

plot(res.acp, axes = c(1, 2), choix = "ind",
     graph.type="classic")

plot(res.acp, axes = c(1, 2), choix = "var",
     graph.type="classic")

plot(res.acp, axes = c(2, 3), choix = "ind",
     graph.type="classic")

plot(res.acp, axes = c(2, 3), choix = "var",
     graph.type="classic")

plot(res.acp, axes = c(3, 4), choix = "ind",
     graph.type="classic")

plot(res.acp, axes = c(3, 4), choix = "var",
     graph.type="classic")

plot(res.acp, axes = c(4, 5), choix = "ind",
     graph.type="classic")

plot(res.acp, axes = c(4, 5), choix = "var",
     graph.type="classic")

plot(res.acp, axes = c(5, 6), choix = "ind",
     graph.type="classic")

plot(res.acp, axes = c(5, 6), choix = "var",
     graph.type="classic")

## On retrouve bien que les variables sont représentées seulement dans les deux premiers axes principaux
## 



#########
## Q3
#########

reconstruct <- function(res, nr, Xm, Xsd){
  Fs <- as.matrix(res$ind$coord[,1:nr])
  us <- as.matrix(res$svd$V[,1:nr])
  rep <- Fs[,1]%*%t(us[,1])
  if(nr > 1){
    for( i in 2:nr){
      rep <- rep + Fs[,i]%*%t(us[,i])
    }
  }
  rep<- t(apply(rep, 1, function(x){x*Xsd+Xm}))
  return(rep)
}

#Verification
Xm <- apply(xtrain, 2, mean)
Xsd <- apply(xtrain, 2, function(x){sqrt((var(x)*(length(x)-1))/length(x))})

reconstruct(res.acp, 1, Xm, Xsd)
sum((reconstruct(res.acp, 39, Xm, Xsd) - xtrain)^2)/28000
sum(abs(reconstruct(res.acp, 39, Xm, Xsd) - xtrain))/28000

plot(1:700,reconstruct(res.acp, 39, Xm, Xsd)[1,], type = 'l')
for(i in 2:40){
  lines(1:700, reconstruct(res.acp, 39, Xm, Xsd)[i,], type = 'l')
}

#Représenter la reconstruction
par(mfrow = c(2, 3))
nr = c(1, 2, 3, 4, 5, 39)
for(i in nr){
  X <- reconstruct(res.acp, i, Xm, Xsd)
  RMSE <- round(sqrt(sum((X - xtrain)^2)/28000), 5)
  MAE <- round(sum(abs(X - xtrain))/28000, 5)
  plot(1:700,X[1,], type = 'l',
       main = paste("Reconstruction pour nr =", i, "; RMSE =", RMSE, "MAE =", MAE,"."),
       xlab = "Fréquence (pas en HZ)",
       ylab = "Amplitude")
  for(j in 2:nrow(X)){
    lines(1:700,X[j,], type = 'l')
  }
}

par(mfrow = c(1, 1))
nr = c(1, 2, 3, 4, 5, 39)
X <- reconstruct(res.acp, 1, Xm, Xsd)
X24 <- X[,24]
for(i in nr[-1]){
  X <- reconstruct(res.acp, i, Xm, Xsd)
  X24 <- cbind(X24, X[,24])
}
plot(X24[,1], col = 1)
for(j in 2:6){
  points(X24[,j], col = j)
}

#################################
## Régression pénalisée 
#################################
require(glmnet)
require(MASS)
#########
## Q1
#########
glmnet(x = xtrain, y = unlist(ytrain), grid=10^seq(6,-10,length=100))
#Le paramètre estimé tend vers 0


#########
## Q2
#########
Y <- as.matrix(ytrain)
X <- as.matrix(xtrain)
lm.ridge( Y ~ X, grid=10^seq(6,-10,length=100))


#################################
## Régression logistique pénalisée 
#################################

