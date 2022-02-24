##########################
## WBA ENTROPY ANALYSIS ##
##########################
library(plyr)
library(ggplot2)
library(lme4)
library(lmerTest)

set <- read.csv("https://raw.githubusercontent.com/ProfNascimento/WBA/main/database_entropy_sham.csv", sep=";")[,-1]
names(set)

set$Montage=as.factor(set$Montage)
set$Dose=as.factor(set$Dose)
set$Dose=ordered(set$Dose,c("0", "1", "2", "3"))
set$Replica=as.factor(set$Replica)
set$Groups=as.factor(set$Groups)
set$Groups = relevel(set$Groups, "baseline")
set$Side=as.factor(set$Side)

## FIGURE 04
a=ggplot(data = set,aes(x = Side,fill=Side))+ggtitle("Force Plates Variation")+
  stat_summary_bin(aes(y = Entropy), 
                   fun.y = mean,
                   fun.ymin = function(x) mean(x), 
                   fun.ymax = function(x) mean(x) + sd(x), geom = "pointrange")+
  stat_summary_bin(aes(y = Entropy), fun.y = "mean", geom = "bar")+labs(x="Side")+labs(y="Entropy")

b=ggplot(data=set,aes(x=Time,y = Entropy,color=Side)) + geom_smooth()+labs(x="Trial Time (minute)")

cols=c("L_Cathode"="#0000cc", "L_Anode"="red", "L_Sham-CC"="#9999cc", "L_Sham-AC"="#FF9999",
       "R_Cathode"="#0000cc", "R_Anode"="red", "R_Sham-CC"="#9999cc", "R_Sham-AC"="#FF9999")
c=ggplot()+ggtitle("Dose-response per Montage")+
  geom_smooth(data = set[set$Groups=="AC" & set$Side=="L",],se=F,
              aes(x = Dose, y = Entropy, group = Groups, colour="L_Cathode"),method = "loess",linetype="twodash") +
  geom_smooth(data = set[set$Groups=="CC" & set$Side=="L",],se=F,
              aes(x = Dose, y = Entropy, group = Groups, colour="L_Anode"),method = "loess",linetype="twodash") +
  geom_smooth(data = set[set$Groups=="SHAM-AC" & set$Side=="L",],se=F,
              aes(x = Dose, y = Entropy, group = Groups, colour="L_Sham-AC"),method = "loess",linetype="twodash") +
  geom_smooth(data = set[set$Groups=="SHAM-CC" & set$Side=="L",],se=F,
              aes(x = Dose, y = Entropy, group = Groups, colour="L_Sham-CC"),method = "loess",linetype="twodash") +
  geom_smooth(data = set[set$Groups=="AC" & set$Side=="R",],se=F,
              aes(x = Dose, y = Entropy, group = Groups, colour="R_Cathode"),method = "loess") +
  geom_smooth(data = set[set$Groups=="CC" & set$Side=="R",],se=F,
              aes(x = Dose, y = Entropy, group = Groups, colour="R_Anode"),method = "loess") +
  geom_smooth(data = set[set$Groups=="SHAM-AC" & set$Side=="R",],se=F,
              aes(x = Dose, y = Entropy, group = Groups,colour="R_Sham-AC"),method = "loess") +
  geom_smooth(data = set[set$Groups=="SHAM-CC" & set$Side=="R",],se=F,
              aes(x = Dose, y = Entropy, group = Groups,colour="R_Sham-CC"),method = "loess") +
  scale_colour_manual(name="Legend",values=cols)+
  scale_linetype_manual(values=c("line","line","line","line","twodash","twodash","twodash","twodash"))+labs(x="Dose (mA)")+labs(y="Entropy")

cowplot::plot_grid(cowplot::plot_grid(a,b, nrow = 2, labels = c("A", "C")),
                   c, nrow = 1, labels = c("", "B"), rel_widths = c(1,3))

## FIGURE 05
set$Groups=factor(set$Groups,levels = c("baseline","AC","SHAM-AC","CC","SHAM-CC"))

bb=ggplot(set[set$Dose==0&set$Groups!="baseline",],aes(x=Trial,y=Entropy,color=Side))+
  geom_jitter()+geom_smooth()+facet_wrap(~Groups,ncol=4)+ylim(0,0.175)+
  theme(text=element_text(size=16))

aa=ggplot(set[set$Dose==0&set$Groups=="baseline",],aes(x=Trial,y=Entropy,color=Side))+ylim(0,0.185)+
  geom_violin()+theme(text=element_text(size=16),legend.position = "none")

cowplot::plot_grid(aa,bb, ncol = 2, labels = c("", ""), rel_widths = c(1,5))

## TABLE 1
ddply(set[set$Groups=="SHAM-AC" | set$Groups=="SHAM-CC",],
      .(Side,Dose,Groups),summarise,
      Median=quantile(Entropy,0.5, na.rm =T),
      MIN=min(Entropy, na.rm =T),
      MAX=max(Entropy, na.rm =T))

## TEST U (SHAM-CC vs SHAM-AC)
# RIGHT SIDE ONLY
wilcox.test(set[set$Dose==0 & set$Groups=="SHAM-AC" & set$Side=="R","Entropy"], 
            set[set$Dose==0 & set$Groups=="SHAM-CC" & set$Side=="R","Entropy"])$p.val

wilcox.test(set[set$Dose==1 & set$Groups=="SHAM-AC" & set$Side=="R","Entropy"], 
            set[set$Dose==1 & set$Groups=="SHAM-CC" & set$Side=="R","Entropy"])$p.val

wilcox.test(set[set$Dose==2 & set$Groups=="SHAM-AC" & set$Side=="R","Entropy"], 
            set[set$Dose==2 & set$Groups=="SHAM-CC" & set$Side=="R","Entropy"])$p.val

wilcox.test(set[set$Dose==3 & set$Groups=="SHAM-AC" & set$Side=="R","Entropy"], 
            set[set$Dose==3 & set$Groups=="SHAM-CC" & set$Side=="R","Entropy"])$p.val

# LEFT SIDE ONLY
wilcox.test(set[set$Dose==0 & set$Groups=="SHAM-AC" & set$Side=="L","Entropy"], 
            set[set$Dose==0 & set$Groups=="SHAM-CC" & set$Side=="L","Entropy"])$p.val

wilcox.test(set[set$Dose==1 & set$Groups=="SHAM-AC" & set$Side=="L","Entropy"], 
            set[set$Dose==1 & set$Groups=="SHAM-CC" & set$Side=="L","Entropy"])$p.val

wilcox.test(set[set$Dose==2 & set$Groups=="SHAM-AC" & set$Side=="L","Entropy"], 
            set[set$Dose==2 & set$Groups=="SHAM-CC" & set$Side=="L","Entropy"])$p.val

wilcox.test(set[set$Dose==3 & set$Groups=="SHAM-AC" & set$Side=="L","Entropy"], 
            set[set$Dose==3 & set$Groups=="SHAM-CC" & set$Side=="L","Entropy"])$p.val

## Gamma-Regression (GLMM)
set.seed(1990)
#set$TimeH=(set$Time/60) # UNIT per hour 
fit.m=glmer(Entropy ~ -1 + Groups*Dose + Time*Side + (1|ID), 
            data=set,family=Gamma(link="log"))
summary(fit.m)
AIC(fit.m);BIC(fit.m)

require(jtools)
plot_summs(fit.m,scale=TRUE,exp=TRUE)
cat_plot(fit.m, mod2=Dose, pred=Side,modx=Groups,
         colors=c("black","red","#FF9999","#0000cc","#9999cc"))

hist(residuals(fit.m),freq = F,breaks = seq(-2.2,2.2,0.2),
     main="Residual hist. - Mixed-Model Regression",xlab="Residuals")
car::qqPlot(residuals(fit.m),envelope=.99,cex=0.5)


## BAYESIAN VERSION
library(brms)  # for models

prior1 <- prior(normal(0,10), class = b) + prior(cauchy(0,2), class = sd)
fit.mod30=brm(Entropy~-1+Groups*Dose+Time*Side+(1|ID), data=set,
              family=Gamma(link="log"),chains=4,thin=10,iter=15000,warmup=5000, 
              prior = prior1,seed=123456,cores=getOption("mc.cores",4))

## Summary of the Bayesian adopted Model
print(fit.mod30,digits=3)

## Convergence Checking
plot(fit.mod30, plotfun = "trace")
