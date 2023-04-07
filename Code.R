library(emmeans)
library(survival)
library(rawr)
library(survminer)
library(ggplot2)
library(ggplotify)

data<-read.csv("scavenging.csv")
head(data)
data$carcass[data$carcass == "MOUSE"] <- "Mouse"
data$carcass[data$carcass == "RAT"] <- "Rat"
data$carcass[data$carcass == "SNAKE"] <- "Brown Tree Snake"

#GLM with binomial distribution- response is whether it was scavenged
#Dredging then fitting the top model
modelcomplete<-glm(Scaved~(season+habitat+carcass)^3, family = "binomial", data=data, na.action=na.fail)
modddredge<-dredge(modelcomplete)
fitted<-glm(Scaved~carcass+habitat+season+season:habitat+carcass:habitat, family = "binomial", data=data, na.action=na.fail)
summary(fitted)


#Parsing out the habitat carcass interaction"
emmeans(fitted, pairwise~habitat|carcass, adjust="Tukey", type='response')
predict(fitted, data.frame(carcass="Rat", season="Wet", habitat="Coastal"), se.fit=TRUE, type=c('response'))
predict(fitted, data.frame(carcass="Rat", season="Wet", habitat="Upland"), se.fit=TRUE, type=c('response'))
emmip(fitted, ~habitat|carcass, adjust="tukey", type="response", CIs=TRUE)

#Parsing out the season habitat interaction#
emmeans(fitted, pairwise~season|habitat)
predict(fitted, data.frame(carcass="Rat", season="Dry", habitat="Upland"), se.fit=TRUE, type=c('response'))
predict(fitted, data.frame(carcass="Rat", season="Wet", habitat="Upland"), se.fit=TRUE, type=c('response'))
emmip(fitted, ~season|habitat, adjust="tukey", type="response", CIs=TRUE)

#Figure
dat<-emmip(fitted, ~habitat, CIs=TRUE, plotit=FALSE)
emmip(data=dat, habitat~weight, at=mylist, CIs=TRUE)
ggplot(data=dat, aes(x=weight, y=yvar, color=habitat))+geom_line()+
  geom_ribbon(aes(ymax=UCL, ymin=LCL,fill=habitat, alpha=0.4))+theme(strip.background=element_rect(color="black", fill="darkslategray3"), panel.border = element_rect(colour='black', fill='NA'), 
                  panel.background = element_rect(fill='gray97'), panel.grid = element_line(colour = NA), panel.spacing = unit(1, "lines"),
                legend.key = element_rect(fill='white'), legend.position="bottom", legend.text = element_text(size=14),
                legend.title = element_text(size=14), axis.ticks.x = element_blank(),
           axis.text=element_text(size=14),
          axis.title.y = element_text(size=14, vjust=3.9), strip.text = element_text(size=14, face="bold"))+
 ylab("Estimated proportion scavenged")+ scale_y_continuous(limits = c(0,1))



#Species richness#
datarich<-data[!(data$Rich=="0"),] #Removing the ones that weren't scavenged
nrow(datarich)
mean(datarich$Rich)
xtabs(~Rich, data=datarich)
modelrich<-glm(Rich~(season+carcass+habitat)^2, family = "poisson", data=datarich, na.action=na.fail)
summary(modelrich)
dredge(modelrich)
#top model was the null, nothing further to explore

a1<-data.frame(dredge(modelcomplete))
b1<-data.frame(dredge(modelrich))
write.csv(bind_rows(a1, b1), "bigmod.csv")

#Binary GLm by species#
#For each: dredging then making predictions from the top model
nosnake<-data[!(data$carcass=="SNAKE"),] #removing the BTS carcasses because they weren't scavenged by BTS
nrow(nosnake)
snakefull<-glm(BTSscav~habitat+season+carcass, family="binomial", data=nosnake, na.action=na.fail)
summary(snakefull)
dredge(snakefull)

snakefit<-glm(BTSscav~habitat, family="binomial", data=nosnake, na.action=na.fail)
summary(snakefit)
predict(snakefit, data.frame(habitat="Coastal"), se.fit=TRUE, type=c('response'))
predict(snakefit, data.frame(habitat="Upland"), se.fit=TRUE, type=c('response'))

toadfull<-glm(Ctscav~habitat+season+carcass, family="binomial", data=data2, na.action=na.fail)
dredge(toadfull)
predict(toadfull, data.frame(habitat="Upland", season="Wet", carcass="RAT"), se.fit=TRUE, type=c('response'))
predict(toadfull, data.frame(habitat="Upland", season="Dry",carcass="RAT"), se.fit=TRUE, type=c('response'))
predict(toadfull, data.frame(habitat="Upland", season="Wet", carcass="RAT"), se.fit=TRUE, type=c('response'))
predict(toadfull, data.frame(habitat="Coastal", season="Wet",carcass="RAT"), se.fit=TRUE, type=c('response'))
emmip(toadfull, ~carcass)

MLfull<-glm(MLscav~habitat+season+carcass, family="binomial", data=data, na.action=na.fail)
dredge(MLfull)
MLfit<-glm(MLscav~habitat+season, family="binomial", data=data, na.action=na.fail)
predict(MLfit, data.frame(habitat="Upland", season="Dry"), se.fit=TRUE, type=c('response'))
predict(MLfit, data.frame(habitat="Upland", season="Wet"), se.fit=TRUE, type=c('response'))
predict(MLfit, data.frame(habitat="Coastal", season="Wet"), se.fit=TRUE, type=c('response'))
predict(MLfit, data.frame(habitat="Coastal", season="Dry"), se.fit=TRUE, type=c('response'))

hermfull<-glm(Hermitscav~habitat+season+carcass, family="binomial", data=data, na.action=na.fail)
hermfit<-glm(Hermitscav~habitat+carcass, family="binomial", data=data, na.action=na.fail)
dredge(hermfull)
summary(hermfull)
emmeans(hermfit, pairwise~carcass, type="response")
emmip(hermfit, ~carcass, type="response")
predict(hermfit, data.frame(habitat="Coastal", carcass="RAT"), se.fit=TRUE, type=c('response'))
predict(hermfit, data.frame(habitat="Upland", carcass="RAT"), se.fit=TRUE, type=c('response'))

sm<-data.frame(dredge(snakefull))
tm<-data.frame(dredge(toadfull))
mlm<-data.frame(dredge(MLfull))
hcm<-data.frame(dredge(hermfull))
write.csv((bind_rows(sm, tm, mlm, hcm)), "models.csv")

#decomposition survival analysis#

survdata<-read.csv("decomp.csv")
survdata$Days<-(survdata$Time/24)
head(survdata)
nrow(survdata)
tabl<-xtabs(~habitat+carcass+season+status, data=survdata)
ftable(tabl)
t<-survreg(Surv(Time, status)~(season+habitat+carcass+scale(weight)),
           data=survdata, dist='weibull', na.action=na.fail)
summary(t)
dredge(t)

fittedsurv<-survreg(Surv(Days, status)~carcass+season, data=survdata, dist='weibull', na.action=na.fail)
fittedseason<-survfit(Surv(Days, status)~season, data=survdata)
fittedcarc<-survfit(Surv(Days, status)~carcass, data=survdata)
survdiff_pairs(fittedcarc)
head(survdata)
custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      axis.title.y=element_text(angle=90, hjust=0.5, vjust=5),
      axis.title.x = element_text(size=12, vjust=-2),
      plot.margin=unit(c(0.5,0.5,0.5,1), "cm")
    )
}

cbbPalette <- c("#E69F00", "#56B4E9", "#009E73")

ggsurvplot(fittedseason, censor=FALSE, ggtheme = custom_theme(), palette =cbbPalette,
           legend="bottom", legend.title="Season", legend.labs=c("Dry", "Wet"), xlab= "Time (Days)", ylab="Proportion carcasses not decomposed")
ggsurvplot(fittedcarc, censor=FALSE, ggtheme = custom_theme(), palette =cbbPalette, 
           legend="bottom", legend.title="Carcass", legend.labs=c("Brown tree snake", "Mouse", "Rat"), xlab= "Time (Days)", ylab="Proportion carcasses not decomposed")

predict(fittedsurv, data.frame(carcass="RAT", season="WET", habitat="T"), se.fit=TRUE, type=c('quantile'))
predict(fittedsurv, data.frame(carcass="RAT", season="WET", habitat="U"), se.fit=TRUE, type=c('response'))
predict(fittedsurv, data.frame(carcass="RAT", season="DRY", habitat="U"), se.fit=TRUE, type=c('response'))

survdata$carcass[survdata$carcass == "MOUSE"] <- "Mouse"
survdata$carcass[survdata$carcass == "RAT"] <- "Rat"
survdata$carcass[survdata$carcass == "BTS"] <- "Brown tree snake"

#Histogram of days to decomposition
ggplot(aes(x=Days,fill=carcass), data=survdata)+geom_histogram(color="black", fill="grey80", binwidth = 0.75)+theme(strip.background=element_rect(color="black", fill="lightblue"), panel.border = element_rect(colour='black', fill='NA'), 
                                                       panel.background = element_rect(fill='white'), panel.grid = element_line(colour = NA), panel.spacing = unit(1, "lines"),
                                                       legend.key = element_rect(fill='white'), legend.position="none", 
                                                       axis.text=element_text(size=12), axis.title.x = element_text(size=14, vjust=-1),
                                                       axis.title.y = element_text(size=14, vjust=3), strip.text = element_text(size=12, face="bold"))+scale_color_manual(values=c('coral2', 'steelblue', 'black'))+
  theme(plot.margin=unit(c(0.5,0.5,0.5,1), "cm"))+xlab("Days to decomposition")+ylab("Number of carcasses")+facet_wrap(~carcass)

#kaplan-meier diagnostic plots for the weibull regression for survival analysis
survdata2<-read.csv("survival.csv")
survdata3<-survdata2[!(survdata2$TimeRight =="No"),] #removing ones with incorrect time stamps
s<-with(survdata3,Surv(Time, Status))
fKM <- survfit(s ~ 0,data=survdata3)
sWei <- survreg(s ~ as.factor(habitat)+as.factor(carcass),dist='weibull',data=survdata3)
plot(fKM, col = "white")
pct<-seq(0.01, 0.99, by=0.01)
lines(predict(sWei, newdata=list(habitat="T", carcass="Rat"),type="quantile",p=seq(.01,.99,by=.01)),seq(.99,.01,by=-.01), col="red")
lines(predict(sWei, newdata=list(habitat="T", carcass="Mouse"),type="quantile",p=seq(.01,.99,by=.01)),seq(.99,.01,by=-.01),col="blue")
lines(predict(sWei, newdata=list(habitat="T", carcass="Snake"),type="quantile",p=seq(.01,.99,by=.01)),seq(.99,.01,by=-.01),col="green")

lines(predict(sWei, newdata=list(habitat="U", carcass="Rat"),type="quantile",p=seq(.01,.99,by=.01)),seq(.99,.01,by=-.01),col="red", lty=2)
lines(predict(sWei, newdata=list(habitat="U", carcass="Mouse"),type="quantile",p=seq(.01,.99,by=.01)),seq(.99,.01,by=-.01),col="blue", lty=2)
lines(predict(sWei, newdata=list(habitat="U", carcass="Snake"),type="quantile",p=seq(.01,.99,by=.01)),seq(.99,.01,by=-.01),col="green", lty=2)

sWei2 <- survreg(s ~ as.factor(habitat)+as.factor(season),dist='weibull',data=survdata3)
lines(predict(sWei2, newdata=list(habitat="T", season="WET"),type="quantile",p=seq(.01,.99,by=.01)),seq(.99,.01,by=-.01),col="red")
lines(predict(sWei2, newdata=list(habitat="U", season="WET"),type="quantile",p=seq(.01,.99,by=.01)),seq(.99,.01,by=-.01),col="blue")
lines(predict(sWei2, newdata=list(habitat="T", season="DRY"),type="quantile",p=seq(.01,.99,by=.01)),seq(.99,.01,by=-.01),col="green")
lines(predict(sWei2, newdata=list(habitat="U", season="DRY"),type="quantile",p=seq(.01,.99,by=.01)),seq(.99,.01,by=-.01),col="black")

#survival analysis- days to carcass consumption

t2<-survreg(Surv(Time, Status)~((season+carcass+habitat)^3),
            data=survdata3, dist='weibull', na.action=na.fail) #weibull regression
summary(t2)
dredge(t2)
fittedsurv<-survreg(Surv(Time, Status)~carcass+season+habitat+season:habitat+carcass:habitat, data=survdata3, dist='weibull', na.action=na.fail)
summary(fittedsurv)

survdata3$int<-with(survdata3, interaction(carcass, habitat)) #posthoc comparisons for the interactions
sfit<-survfit(Surv(Time, Status)~int, data=survdata3)
survdiff_pairs(sfit, method='BH',  rho=0)

survdata3$int2<-with(survdata3, interaction(season, habitat))
sfit2<-survfit(Surv(Time, Status)~int2, data=survdata3)
combine_table(surv_table(sfit2))
survdiff_pairs(sfit2, method='BH', rho=0)


##Figures##
#time to carcass removal by carcass type
ggsurvplot_facet(sfit, facet.by="carcass", censor=FALSE, legend.title="Habitat", xlab="Time (days)", ylab= "Proportion not fully scavenged",
                 data=survdata3, legend.labs=c("Coastal", "Coastal", "Coastal", "Upland", "Upland", "Upland"),
                 panel.labs = list(c("Mouse", "Rat", "Snake")), short.panel.labs = TRUE)+theme_classic()+
  theme(legend.position="bottom", panel.spacing.x = unit(1, "lines"),
        strip.background = element_rect(color="black", size=0.5, fill="lightblue"), panel.border = element_rect(colour='black', fill='NA'), 
        panel.background = element_rect(fill='gray97'), panel.grid = element_line(colour = NA), panel.spacing = unit(1, "lines"),
        legend.key = element_rect(fill='white'), legend.title = element_text(size=12, face="bold"),
        axis.text=element_text(size=12), axis.title.x = element_text(size=14, vjust=-1),
        axis.title.y = element_text(size=14, vjust=2), strip.text = element_text(size=12, face="bold"))+
  theme(plot.margin=unit(c(0.5,0.5,0.8,0.8), "cm"))

#time to carcass removal by habitat
ggsurvplot_facet(sfit2, facet.by="habitat", censor=FALSE, legend.title="Habitat", xlab="Time (days)", ylab= "Proportion not fully scavenged",
                 data=survdata3, legend.labs=c("Dry", "Wet", "Dry", "Wet"),
                 panel.labs = list(habitat=c("Coastal", "Upland")), short.panel.labs = TRUE)+theme_classic()+
  theme(legend.position="bottom", panel.spacing.x = unit(1, "lines"),
        strip.background = element_rect(color="black", size=0.5, fill="lightblue"), panel.border = element_rect(colour='black', fill='NA'), 
        panel.background = element_rect(fill='gray97'), panel.grid = element_line(colour = NA), panel.spacing = unit(1, "lines"),
        legend.key = element_rect(fill='white'), legend.title = element_text(size=12, face="bold"),
        axis.text=element_text(size=12), axis.title.x = element_text(size=14, vjust=-1),
        axis.title.y = element_text(size=14, vjust=2), strip.text = element_text(size=12, face="bold"))+
  theme(plot.margin=unit(c(0.5,0.5,0.8,0.8), "cm"))

#Proportion scavenged by species and habitat then combining them to one figure
my_labeller<-as_labeller(function(x){
  return(paste0("", x))
})
a<-emmip(fitted, ~habitat|carcass, CIs=TRUE, CIarg=aes(lwd=1, color="black", linetype=1), type="response")+
  facet_wrap(~carcass, labeller = my_labeller)+
  theme(strip.background=element_rect(color="black", fill="lightblue"), panel.border = element_rect(colour='black', fill='NA'), 
        panel.background = element_rect(fill='gray97'), panel.grid = element_line(colour = NA), panel.spacing = unit(1, "lines"),
        legend.key = element_rect(fill='white'), legend.position="bottom", 
        axis.text=element_text(size=12), axis.title.x = element_blank(),
        axis.title.y = element_blank(), strip.text = element_text(size=12, face="bold"))+
  theme(plot.margin=unit(c(0.7,0.5,0.3,1), "cm"))+xlab("Habitat")+ylab("Proportion of carcasses scavenged")
b<-emmip(fitted, ~season|habitat, CIs=TRUE, CIarg=aes(lwd=1, color="black", linetype=1), type="response")+
  facet_wrap(~habitat, labeller = my_labeller)+
  theme(strip.background=element_rect(color="black", fill="lightblue"), panel.border = element_rect(colour='black', fill='NA'), 
        panel.background = element_rect(fill='gray97'), panel.grid = element_line(colour = NA), panel.spacing = unit(1, "lines"),
        legend.key = element_rect(fill='white'), legend.position="bottom", 
        axis.text=element_text(size=12), axis.title.x = element_blank(),
        axis.title.y = element_blank(), strip.text = element_text(size=12, face="bold"))+
  theme(plot.margin=unit(c(0.5,0.5,0.3,1), "cm"))+xlab("Season")+ylab("Proportion of carcasses scavenged")

annotate_figure((ggarrange(a,b, ncol=1, labels=c("A", "B"), hjust=c(-3,-3), vjust=c(1.6,1))), 
                left=text_grob("Proportion of carcasses scavenged", size = 14, rot=90, vjust = 1.3, hjust=0.4))







