#Full analyses of visual acuity estimates from behavioural & morphological assays
#Presented in Wright et al. (2023) "Quantifying visual acuity in Heliconius butterflies"

library(MASS)
library(lme4)
library(ggplot2)
library(stats)
library(effects)
library(splines)
library(xlsx)
library(car)
library(lattice)
library(grid)
library(gridExtra)
library(pbkrtest)
library(multcomp)
library(dplyr)

####BEHAVIOURAL VISUAL ACUITY####

#import data from csv
#csv file "opt.cyr"
opt.cyr <- read.csv("~/Documents/LMU/Optomotor_cyrbia/opt_cyr.csv")


#####Data Organization#####

#check the data structure
str(opt.cyr)
opt.cyr$weather <- as.factor(opt.cyr$weather)
opt.cyr$id <- as.factor(opt.cyr$id)
opt.cyr$species <- as.factor(opt.cyr$species)
opt.cyr$sex <- as.factor(opt.cyr$sex)
opt.cyr$origin <- as.factor(opt.cyr$origin)
opt.cyr$observer <- as.factor(opt.cyr$observer)
str(opt.cyr) #ok

#how many per sex
t.first <- opt.cyr[match(unique(opt.cyr$id), opt.cyr$id),]
t.first %>%
  group_by(sex) %>%
  summarize(count = n())
#23 females
#26 males

#####Mean Visual Acuity Estimates#####

######overall######
opt.cyr %>%
  group_by(species) %>%
  summarise_at(vars(threshold), list(name = mean, sd))

#mean = 0.491
#std = 0.189

#standard error of the overall mean = 0.027
0.189/(sqrt(49))

######males######
opt.cyr %>%
  group_by(sex) %>%
  summarise_at(vars(threshold), list(name = mean, sd))

#male mean = 0.547
#male std = 0.224

#male se = 0.04393001
0.224/ (sqrt(26))

######females######
opt.cyr %>%
  group_by(sex) %>%
  summarise_at(vars(threshold), list(name = mean, sd))

#female mean = 0.427
#female std = 0.114

#female se = 0.02377064
0.114/ (sqrt(23))


#####Analyses of Behavioural Visual Acuity#####

#what does the response variable look like?
hist(opt.cyr$threshold) #Gamma family for GLM seems appropriate. 

######set the random effect structure######
c1 <- glmer(threshold~sex + time + (1|weather) + (1|observer),
            data = opt.cyr, family = "Gamma",
            control=glmerControl(calc.derivs = F,optCtrl = list(maxfun=20000)))

c1.1 <- glmer(threshold~sex + time +  (1|observer),
              data = opt.cyr, family = "Gamma",
              control=glmerControl(calc.derivs = F,optCtrl = list(maxfun=20000)))

AIC(c1,c1.1) 
anova(c1, c1.1) #no difference, use only observer

c1.2 <- glmer(threshold~sex + time +  (1|weather),
              data = opt.cyr, family = "Gamma",
              control=glmerControl(calc.derivs = F,optCtrl = list(maxfun=20000)))

AIC(c1.1,c1.2)
anova(c1.1, c1.2)  #no difference, stick with only observer

#model without random effects
c1.3 <- glm(threshold~sex + time, data = opt.cyr)

AIC(c1.1,c1.3) #model with observer random effect is better

#check model assumptions
plot(c1.1) 
hist(resid(c1.1)) 

qqnorm(resid(c1.1))
qqline(resid(c1.1)) 


######check for overdispersion######
overdisp <- function(model) {
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  (rdf <- nrow(model@frame)-model.df)
  rp <- residuals(model)
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE,log.p=TRUE)
  c(chisq=Pearson.chisq,ratio=prat,p=exp(pval))
}

overdisp(c1.1) #model is not overdispersed


######simplify the fixed effects structure######
drop1(c1.1, test="Chisq") #time ns
c2 <- update(c1.1,.~.-time)
drop1(c2, test="Chisq") #sex significant

Anova(c2,test="Chisq") 
#sex p=0.02071

plot(allEffects(c2)) #plotted sex effects of the m.a.m
#plotting inverse, higher values are lower (x-axis inverted?)

qqnorm(resid(c2))
qqline(resid(c2)) 

Anova(c1,test="Chisq") 
#time p=0.4


#####**Behaviour Sex Effect Plot**#####

#theme settings for plots
require(grid)
pub<- theme_update(
  panel.grid.major=element_line(colour=NA),
  panel.grid.minor=element_line(colour=NA),
  panel.background = element_rect(colour = NA,fill=NA,linewidth = 1),
  panel.border = element_rect(linetype = "solid", colour = "gray47",fill=NA),
  axis.line.x = element_line(color="black"),
  axis.line.y = element_line(color="black"),
  axis.title.x=element_text(size=15,face="bold",hjust=0.5,vjust=0.5,angle=0),
  axis.title.y=element_text(size=15,face="bold",hjust=0.5,vjust=1,angle=90),
  axis.text.x=element_text(colour="black",angle=0,size=15),
  axis.text.y=element_text(colour="black",angle=0,size=15),
  axis.ticks=element_line(colour="black",linewidth=1))

#order sex as a factor 
opt.cyr$sex <- factor(opt.cyr$sex, c("M","F"))

#import additional data for males tested at >1.0 cpd (group two)
#csv file "male tests"

male_tests <- read.csv("~/Documents/LMU/Optomotor_cyrbia/male tests.csv")

str(male_tests)
male_tests$sex <- as.factor(male_tests$sex)

#restrict data only males responding at 1.0 cpd or higher

male_tests2 <- male_tests[which(male_tests$threshold >= 1),]

#dataset to calculate hypothetical mean if all males from group one had 1.4 cpd
#only males that responded at 1.0 cpd (4 males)

opt.cyr2 <- opt.cyr
opt.cyr2$threshold[opt.cyr2$threshold == 1.0] <- 1.4
opt.cyr2.1 <- opt.cyr2[which(opt.cyr2$sex == "M"),]


library(ggsignif)

behav <- ggplot(opt.cyr,aes(x=sex,y=threshold))+
  coord_cartesian(ylim=c(0.3, 1.52))+ 
  geom_point(data=opt.cyr, aes(x = sex, y = threshold, fill=sex),
             position = position_jitterdodge(0.45), size=1, shape=1, color="black", alpha=0.75)+
  geom_point(data=male_tests2, aes(x = sex, y = threshold, fill=sex),
             position = position_jitterdodge(0.45), size=1, shape=17, color="grey75")+
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar",width = 0.08, size=1, 
               color = "black")+
  stat_summary(fun = mean, geom = "point",size = 4)+
  stat_summary(data=opt.cyr2.1, fun.data = mean_cl_boot, geom = "errorbar",width = 0.08, size=1, 
               color = "grey75", position = position_nudge(0.11))+
  stat_summary(data=opt.cyr2.1,fun = mean, geom = "point",size = 4,
               position = position_nudge(0.11), color = "grey")+
  
  labs(x="",y="visual acuity (cpd)")+
  scale_y_continuous(breaks = c(0.3,0.5,0.7,0.9,1.1,1.3, 1.5))+
  scale_x_discrete(labels=c("male", "female"))+
  geom_signif(comparisons = list(c("M", "F")),map_signif_level=F,
              annotations = c("*"),y_position = c(1.15),
              color = "black", size = 0.5, textsize = 6, tip_length = 0.02,
              vjust=0.4)+
  theme(legend.position='none')+
  theme(axis.text.x=element_text(colour="black",angle=0,size=15, face="bold"))+
  theme(axis.text.y=element_text(colour="black",angle=0,size=14))+
  theme(plot.margin = margin(1,0.5,0.5,0.5, "cm"))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(text = element_text(family = "Times New Roman"))


####PERCEPTUAL ESTIMATES OF VISUAL ACUITY#######
#using AcuityView package

library(AcuityView)

#####males#####
#mean visual acuity = 0.547
1/0.547

#5cm
AcuityView(photo = NULL, distance = 5, realWidth = 8.8,
           eyeResolutionX = 1.828, eyeResolutionY = NULL, plot = T,
           output = "male5.jpg")

#25cm
AcuityView(photo = NULL, distance = 25, realWidth = 8.8,
           eyeResolutionX = 1.828, eyeResolutionY = NULL, plot = T,
           output = "male25.jpg")

#50cm
AcuityView(photo = NULL, distance = 50, realWidth = 8.8,
           eyeResolutionX = 1.828, eyeResolutionY = NULL, plot = T,
           output = "male50.jpg")

#max visual acuity = 1.4
1/1.4

#5cm
AcuityView(photo = NULL, distance = 5, realWidth = 8.8,
           eyeResolutionX = 0.7142, eyeResolutionY = NULL, plot = T,
           output = "male5.jpg")

#25cm
AcuityView(photo = NULL, distance = 25, realWidth = 8.8,
           eyeResolutionX = 0.7142, eyeResolutionY = NULL, plot = T,
           output = "male25.jpg")

#50cm
AcuityView(photo = NULL, distance = 50, realWidth = 8.8,
           eyeResolutionX = 0.7142, eyeResolutionY = NULL, plot = T,
           output = "male50.jpg")


#####females####
#mean visual acuity = 0.427 
1/0.427 

#5cm
AcuityView(photo = NULL, distance = 5, realWidth = 8.8,
           eyeResolutionX = 2.341, eyeResolutionY = NULL, plot = T,
           output = "female5.jpg")

#25cm
AcuityView(photo = NULL, distance = 25, realWidth = 8.8,
           eyeResolutionX = 2.341, eyeResolutionY = NULL, plot = T,
           output = "female25.jpg")

#50cm
AcuityView(photo = NULL, distance = 50, realWidth = 8.8,
           eyeResolutionX = 2.341, eyeResolutionY = NULL, plot = T,
           output = "female50.jpg")

#####bird predator####
#mean visual acuity =  22.05
(25.6 + 18.5) / 2
1/22.05

#5cm
AcuityView(photo = NULL, distance = 5, realWidth = 8.8,
           eyeResolutionX = 0.045, eyeResolutionY = NULL, plot = T,
           output = "bird5.jpg")

#25cm
AcuityView(photo = NULL, distance = 25, realWidth = 8.8,
           eyeResolutionX = 0.045, eyeResolutionY = NULL, plot = T,
           output = "bird25.jpg")

#50cm
AcuityView(photo = NULL, distance = 50, realWidth = 8.8,
           eyeResolutionX = 0.045, eyeResolutionY = NULL, plot = T,
           output = "bird50.jpg")


####EYE MORPHOLOGY####

#import data from csv
#csv file "eye_morph_cyr"
eye_morph <- read.csv("~/Documents/LMU/Optomotor_cyrbia/eye_morph_cyr.csv")

#####Data Organization####

#####Correlation Between Left vs. Right####
cor.test(eye_morph$L_count, eye_morph$R_count)

#new variables using only left side, or right side if left is missing
eye_morph$whole_count <- ifelse(is.na(eye_morph$L_count), eye_morph$R_count, eye_morph$L_count)
eye_morph$area <- ifelse(is.na(eye_morph$L_area), eye_morph$R_area, eye_morph$L_area)
eye_morph$tibia <- ifelse(is.na(eye_morph$L_tibia), eye_morph$R_tibia, eye_morph$L_tibia)


#####Mean ommatidia counts####

######males######
eye_morph %>%
  group_by(sex) %>%
  summarise_at(vars(whole_count), list(name = mean, sd))

#male mean = 14089
#male std = 879

#male se = 332.2308
879/(sqrt(7))


######females######
eye_morph %>%
  group_by(sex) %>%
  summarise_at(vars(whole_count), list(name = mean, sd))

#female mean = 13245
#female std = 493

#female se = 186.3365
493/(sqrt(7))


#####Analysis of ommatidia count#####

#using log10 values to account for allometric relationships

hist(log10(eye_morph$whole_count))

v5 <- lm((log10(whole_count)) ~ sex + (log10(tibia)), data=eye_morph)

plot(v5, which=1)
hist(resid(v5)) 

qqnorm(resid(v5))
qqline(resid(v5)) 

drop1(v5, test="Chisq") #tibia ns
v6 <- update(v5,.~.-log10(tibia))
drop1(v6, test="Chisq") #sex significant

Anova(v6, test="Chisq")
#sex p = 0.046

plot(allEffects(v6))

qqnorm(resid(v6))
qqline(resid(v6)) 

Anova(v5, test="Chisq")
#tibia p = 0.76


#####**Morphology Sex Effect Plot**#####
eye_morph$sex <- factor(eye_morph$sex, c("M","F"))

morph <- ggplot(eye_morph,aes(x=sex,y=whole_count))+
  coord_cartesian(ylim=c(12470, 15520))+ 
  geom_point(data=eye_morph, aes(x = sex, y = whole_count, fill=sex),
             position = position_jitterdodge(0.5), size=1, shape=1, color="black", alpha = 0.75)+
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar",width = 0.08, size=1, color = "black")+
  stat_summary(fun = mean, geom = "point",size = 4)+
  labs(x="",y="ommatidia count")+
  scale_y_continuous(breaks = c(12500, 13000, 13500, 14000, 14500, 15000, 15500))+
  scale_x_discrete(labels=c("male", "female"))+
  geom_signif(comparisons = list(c("M", "F")),map_signif_level=F,
              annotations = c("*"),y_position = c(15025),
              color = "black", size = 0.5, textsize = 6, tip_length = 0.02,
              vjust=0.4)+
  theme(legend.position='none')+
  theme(axis.text.x=element_text(colour="black",angle=0,size=15, face="bold"))+
  theme(axis.text.y=element_text(colour="black",angle=0,size=14))+
  theme(plot.margin = margin(1,0.5,0.5,0.5, "cm"))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(text = element_text(family = "Times New Roman"))


####FINAL COMBO PLOT OF BEHAVIOUR & MORPHOLOGY####
library(ggpubr)
ggarrange(behav, morph,
          vjust=1.25,
          ncol = 2, nrow = 1, widths = c(0.95, 1), heights = c(1,1),
          font.label = list(size = 15), align = "h")


