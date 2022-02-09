library(brms)  # for models
set.seed(123456)

database <- read.csv("https://raw.githubusercontent.com/ProfNascimento/WBA/main/database_entropy_sham.csv", sep=";")[,-1]
names(database)

database$Dose=as.factor(database$Dose)
database$Dose=ordered(database$Dose,c("0", "1", "2", "3"))
database$Side=as.factor(database$Side)
database$Groups=as.factor(database$Groups)

database$Groups=as.factor(database$Groups)
database$Groups=factor(database$Groups, 
                levels = c("baseline","AC","CC","SHAM-AC","SHAM-CC"))

prior1 <- prior(normal(0,10), class = b) + prior(cauchy(0,2), class = sd)
fit.mod30=brm(Entropy~-1+Groups*Dose+Time*Side+(1|ID), data=database,
            family=Gamma(link="log"),chains=4,thin=10,iter=5000,warmup=500,
            prior = prior1)

## Summary of the Bayesian adopted Model
print(fit.mod30,digits=3)

## Convergence Checking
plot(fit.mod30, plotfun = "trace")
