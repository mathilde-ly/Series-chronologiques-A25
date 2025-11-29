
######################################################
###  Consommation d'électricité mensuelle aux USA  ###
######################################################

### Lecture des données

data<-read.table("data/Elec.dat")
data<-ts(data,,start=c(1983,1),frequency=12)


### Graphe de la série brute

autoplot(data,main="Consommation d'électricité",xlab="x",ylab="Consomm",col="blue")

autoplot(elec) +
  ggtitle("Consommation mensuelle d'électricité aux Etats-Unis (1983-2020)") +
  xlab("Année")+
  ylab("Consommation (en milliards de kWh)") +
  col='blue'

## Décomposition

fit1 <- decompose(data, type="additive")

autoplot(fit1) + xlab("Année") +ggtitle("Décomposition additive de la série Elec")

fit2 <- decompose(data, type="multiplicative")

autoplot(fit2) + xlab("Année") +ggtitle("Décomposition multiplicative de la série Elec")


### Homoscédasticité

BoxCox.lambda(data)
lambda<- BoxCox.lambda(data)

### Série des racines carrées

rdata <- sqrt(data)

cbind("Nombre de prescriptions" = data,
      "Racine carrée des prescriptions"=rdata) %>%
  autoplot(facets=TRUE) + xlab("Années") + ylab("")

### Corrélogrammes pour la série des racines carrées

rdata %>%  ggtsdisplay(lag=48)


### Corrélogrammes pour la série des racines carrées différenciée saisonnièrement (D=1)

rdata %>% diff(lag=12) %>%
  ggtsdisplay(xlab="Années", main="Série racine carrée des prescriptions avec D=1") ## commenter


### D=1 et d=1 ou seulement D=1 ?

kpss.test(diff(rdata,12))

### Corrélogrammes pour la série des racines carrées différenciée saisonnièrement et non saisonnièrement (D=1 et d=1)

acf2(diff(diff(rdata,12)))

### Modèle (i)

(fit1<-Arima(data,order=c(4,1,1),seasonal=list(order=c(0,1,2),period=12),fixed=c(NA,NA,0,NA,NA,NA,NA),lambda=lambda))

t_stat(fit1)

### Modèle (ii)

(fit2<-Arima(data,order=c(4,1,1),seasonal=list(order=c(0,1,2),period=12),fixed=c(NA,0,0,NA,NA,NA,NA),lambda=lambda))

t_stat(fit2)

autoplot(fit2)

residuals(fit2) %>% ggtsdisplay

t.test(residuals(fit2))

checkresiduals(fit2)

shapiro.test(fit2$resid)

### Modèle (iii)

(fit3<-Arima(data,order=c(4,1,1),seasonal=list(order=c(0,1,2),period=12),fixed=c(0,0,0,NA,NA,NA,NA),lambda=lambda))

t_stat(fit3)

autoplot(fit3)

residuals(fit3) %>% ggtsdisplay

t.test(residuals(fit3))

checkresiduals(fit3)

shapiro.test(fit3$resid)



###  Autre modèle trouvé après quelques essais

(fit<-Arima(data,order=c(1,1,1),seasonal=list(order=c(0,1,2),period=12),lambda=lambda))

t_stat(fit)

autoplot(fit)

residuals(fit) %>% ggtsdisplay

t.test(residuals(fit))

checkresiduals(fit)

shapiro.test(fit$resid)


### Prévisions sur deux ans avec les intervalles de prévision

h=12

fit3 %>% forecast(h=h) %>% autoplot + ylab("Nombre de prescriptions")+xlab("Années")


## ou en zoomant

fit3 %>% forecast(h=12) %>% autoplot(xlim=c(2015,2021),ylim=c(12,35)) + ylab("Nombre de prescriptions")+xlab("Années") ## zoom

### Prévisions par ets ###

h=12
fit4<-ets(data)
forecast(fit4,h=h) %>% autoplot(xlab="Années")
forecast(fit4,h=h) %>% autoplot(xlab="Années",xlim=c(2014,2020),ylim=c(12,38)) ## zoom



### Essai de robustesse sur deux ans.


train<-window(data,end=c(2016,6))

test<-window(data,start=c(2016,7))

h<-length(test)

(modele<-Arima(train,order=c(4,1,1),seasonal=list(order=c(0,1,2),period=12),fixed=c(0,0,0,NA,NA,NA,NA),lambda=lambda))

t_stat(modele)

autoplot(modele)

t.test(residuals(modele))

residuals(modele) %>% ggtsdisplay

checkresiduals(modele)

shapiro.test(modele$resid)




plot(forecast(modele,h=h))
lines(data)

plot(forecast(modele,h=h),xlim=c(2014.1,2018.6),ylim=c(12,25.5)) ## zoom
lines(data)



accuracy(modele)  ## mesures d'erreur


### Calcul des racines des polynômes autorégressif et moyenne mobile

p=4; P=0
q=1 ; Q=2
d=1 ; D=1
S=12

(prod.fit1<-Arima(train,order=c(4,1,1),seasonal=list(order=c(0,1,2),period=12),fixed=c(0,0,0,NA,NA,NA,NA),lambda=lambda))

# Racines de la partie "Autorégressive"
if(p+P>0) {alpha=c(1,rep(0,P*S+p))
if(P==0) {for (i in 1:p){alpha[i+1]=-prod.fit1$coef[i]}} else {for (i in 1:P) {alpha[S*i+1]=-prod.fit1$coef[p+q+i]}}
if(p>0) { for (i in 1:p) { alpha[1+i]=-prod.fit1$coef[i]
if(P>0) {for (j in 1:P) alpha[j*S+i+1]=prod.fit1$coef[i]*prod.fit1$coef[p+q+j]}}}
polyroot(alpha) #racines du polynôme autorégressif
abs(polyroot(alpha))  }  # modules des racines du polynôme autorégressif

# Racines de la partie "Moyenne Mobile"
if(q+Q>0) {beta=c(1,rep(0,Q*S+q))}
if(Q==0) {for (i in 1:q){beta[i+1]=prod.fit1$coef[p+i]}} else {for (i in 1:Q) {beta[S*i+1]=prod.fit1$coef[p+q+P+i]}}
if(q>0) { for (i in 1:q) { beta[1+i]=prod.fit1$coef[p+i]
if(Q>0) {for (j in 1:Q) beta[j*S+i+1]=prod.fit1$coef[p+i]*prod.fit1$coef[p+q+P+j]}}}
polyroot(beta)  #racines du polynôme moyenne mobile
abs(polyroot(beta))  #modules des racines du polynôme moyenne mobile

#### Robustesse par ets

fit2.modele<-ets(train)
plot(forecast(fit2.modele))
lines(data)

plot(forecast(fit2.modele),xlim=c(2015.1,2018.6),ylim=c(12,27)) ## zoom
lines(data)

accuracy(fit2.modele) ### le modèle ets semble jusqu'ici le meilleur !

