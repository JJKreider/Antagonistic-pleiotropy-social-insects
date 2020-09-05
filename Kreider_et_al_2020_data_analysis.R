# Publication: Antagonistic pleiotropy and the evolution of ageing in social insects
# Authors: Jan J. Kreider, Ido Pen, Boris H. Kramer
# data analysis R-script

library(ggplot2)
library(cowplot)
library(viridis)
library(grid)
library(dplyr)
library(ggpubr)
library(emmeans)
library(multcomp)
library(betareg)

# read all data
setwd("") # set working directory to folder with transformed data (see data transformation R-script)
setwd("mutAcc")
mutAccSur <- read.csv("survivalData.csv")
mutAccFec <- read.csv("fecundityData.csv")
mutAccSur$ext <- as.factor(mutAccSur$ext)
mutAccFec$ext <- as.factor(mutAccFec$ext)
mutAccFec2 <- aggregate(mutAccFec$fecundity, by=list(replicate=mutAccFec$replicate, ext=mutAccFec$ext), FUN=sum)
mutAccSur2 <- aggregate(mutAccSur$survival, by=list(replicate=mutAccSur$replicate, ext=mutAccSur$ext, caste=mutAccSur$caste), FUN=sum)

setwd('..')
setwd("WCWT")
WCWTSur <- read.csv("survivalData.csv")
WCWTFec <- read.csv("fecundityData.csv")
WCWTSur$ext <- as.factor(WCWTSur$ext)
WCWTFec$ext <- as.factor(WCWTFec$ext)
WCWTFec2 <- aggregate(WCWTFec$fecundity, by=list(replicate=WCWTFec$replicate, ext=WCWTFec$ext), FUN=sum)
WCWTSur2 <- aggregate(WCWTSur$survival, by=list(replicate=WCWTSur$replicate, ext=WCWTSur$ext, caste=WCWTSur$caste), FUN=sum)

setwd('..')
setwd("WCBT")
WCBTSur <- read.csv("survivalData.csv")
WCBTFec <- read.csv("fecundityData.csv")
WCBTSur$ext <- as.factor(WCBTSur$ext)
WCBTFec$ext <- as.factor(WCBTFec$ext)
WCBTFec2 <- aggregate(WCBTFec$fecundity, by=list(replicate=WCBTFec$replicate, ext=WCBTFec$ext), FUN=sum)
WCBTSur2 <- aggregate(WCBTSur$survival, by=list(replicate=WCBTSur$replicate, ext=WCBTSur$ext, caste=WCBTSur$caste), FUN=sum)

setwd('..')
setwd("BCWT")
BCWTSur <- read.csv("survivalData.csv")
BCWTFec <- read.csv("fecundityData.csv")
BCWTSur$ext <- as.factor(BCWTSur$ext)
BCWTFec$ext <- as.factor(BCWTFec$ext)
BCWTFec2 <- aggregate(BCWTFec$fecundity, by=list(replicate=BCWTFec$replicate, ext=BCWTFec$ext), FUN=sum)
BCWTSur2 <- aggregate(BCWTSur$survival, by=list(replicate=BCWTSur$replicate, ext=BCWTSur$ext, caste=BCWTSur$caste), FUN=sum)

setwd('..')
setwd("BCBT")
BCBTSur <- read.csv("survivalData.csv")
BCBTFec <- read.csv("fecundityData.csv")
BCBTSur$ext <- as.factor(BCBTSur$ext)
BCBTFec$ext <- as.factor(BCBTFec$ext)
BCBTFec2 <- aggregate(BCBTFec$fecundity, by=list(replicate=BCBTFec$replicate, ext=BCBTFec$ext), FUN=sum)
BCBTSur2 <- aggregate(BCBTSur$survival, by=list(replicate=BCBTSur$replicate, ext=BCBTSur$ext, caste=BCBTSur$caste), FUN=sum)

setwd("..")

# transform data
parameter <- c(rep("mutAcc", 80), rep("WCWT", 80), rep("WCBT", 80), rep("BCWT", 80), rep("BCBT", 80))
Sur <- rbind(mutAccSur2,WCWTSur2,WCBTSur2,BCWTSur2,BCBTSur2)
Sur2 <- cbind(Sur, parameter)
Sur2$parameter <- as.factor(Sur2$parameter) 
Sur2$parameter <- factor(Sur2$parameter, levels(Sur2$parameter)[c(3,5,4,2,1)])

parameter2 <- c(rep("mutAcc", 40), rep("WCWT", 40), rep("WCBT", 40), rep("BCWT", 40), rep("BCBT", 40))
Fec <- rbind(mutAccFec2,WCWTFec2,WCBTFec2,BCWTFec2,BCBTFec2)
Fec2 <- cbind(Fec, parameter2)
Fec2$parameter2 <- as.factor(Fec2$parameter2)
Fec2$parameter2 <- factor(Fec2$parameter2, levels(Fec2$parameter2)[c(3,5,4,2,1)])

# lifespan statistics
Sur2$x2 <- Sur2$x/20
Sur2Q <- Sur2[which(Sur2$caste=="queen"),]
Sur2W <- Sur2[which(Sur2$caste=="worker"),]

# compare queen lifespans with beta regression
model1 <- betareg(x2 ~ parameter * ext, Sur2Q)
summary(model1)
cld(emmeans(model1, ~parameter * ext),sort=FALSE, Letters=letters)

par(mfrow=c(2,3))
plot(model1)
plot(model1, which = 5, type = "deviance", sub.caption = "")
plot(model1, which = 1, type = "deviance", sub.caption = "")

# compare worker lifespans with beta regression
model2 <- betareg(x2 ~ parameter * ext, Sur2W)
summary(model2)
cld(emmeans(model2, ~parameter * ext),sort=FALSE, Letters=letters)

plot(model2)
plot(model2, which = 5, type = "deviance", sub.caption = "")
plot(model2, which = 1, type = "deviance", sub.caption = "")

# compare queen fecundities with beta regression
Fec2$x2 <- Fec2$x/100
model3 <- betareg(x2 ~ parameter2 * ext, Fec2)
summary(model3)
cld(emmeans(model3, ~parameter2 * ext),sort=FALSE, Letters=letters)

plot(model3)
plot(model3, which = 5, type = "deviance", sub.caption = "")
plot(model3, which = 1, type = "deviance", sub.caption = "")
par(mfrow=c(1,1))

# obtain means and standard deviations
Sur4 <- Sur2 %>%
  group_by(ext, caste) %>% 
  summarise(
    "mean" = mean(x),
    "sd" = sd(x)
  )

Sur5 <- Sur2 %>%
  group_by(caste, parameter) %>% 
  summarise(
    "mean" = mean(x),
    "sd" = sd(x)
  )

Fec3 <- Fec2 %>%
  group_by(parameter2) %>% 
  summarise(
    "mean" = mean(x),
    "sd" = sd(x)
  )

Sur6 <- Sur2[which(Sur2$parameter!="mutAcc"),]
Sur7 <- Sur6 %>%
  group_by(caste) %>% 
  summarise(
    "mean" = mean(x),
    "sd" = sd(x)
  )

Sur8 <- Sur6[which(Sur6$parameter!="WCWT" & Sur6$parameter!="WCBT"),]
Sur9 <- Sur8 %>%
  group_by(caste) %>% 
  summarise(
    "mean" = mean(x),
    "sd" = sd(x)
  )

Sur10 <- Sur6[which(Sur6$parameter!="BCWT" & Sur6$parameter!="BCBT"),]
Sur11 <- Sur10 %>%
  group_by(caste) %>% 
  summarise(
    "mean" = mean(x),
    "sd" = sd(x)
  )

mutAccSur3 <- mutAccSur[which(mutAccSur$ext=="0"),]

# Fig 2, survival curve
plot1 <- ggplot(mutAccSur3, aes(ageclass, survival, color = caste, fill = caste)) +
  stat_summary(fun = median, geom="line", size=2) +
  stat_summary(fun.data = "median_hilow", fun.args=(conf.int=1),
               colour = "transparent",
               geom = "ribbon", 
               alpha = 0.2) +
  scale_color_manual(values = c(c(plasma(20))[c(1,12)]),labels = c("Queen", "Worker")) +
  scale_fill_manual(values = c(c(plasma(20))[c(1,12)]),labels = c("Queen", "Worker")) +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
  scale_x_continuous(limits = c(0,20))+
  labs(x = "Age", y = "Survival") +
  theme_minimal()+
  theme(
    plot.title = element_blank(),
    axis.text=element_text(size=15,face="plain",color="black"),
    legend.text = element_text(colour="black", size=17, face="plain",angle=0),
    legend.background = element_rect(color = NA,fill="transparent", size = 0.2, linetype = "solid"),
    legend.title = element_blank(),
    legend.position = c(0.1,0.1), 
    axis.title = element_text(size = 21),
    axis.line = element_line(color="black", size = 1.5),
    panel.border = element_rect(colour = "darkgray", fill=NA, size=1),
    plot.background =element_blank(), 
    legend.key.width = unit(1.1,"cm"), 
    plot.margin=margin(t = 4, r = 4, b = 8, l = 8, unit = "pt"),
    panel.background = element_rect(fill = "white", colour = "grey50"),
    strip.text.x = element_text(size= 17, color = "black", face = "italic"),
    strip.text.y = element_text(size = 17, color = "black")
  )
plot1

pdf(file="survival_curve.pdf", onefile=T,  width = 8, height = 6) 
plot(plot1)
dev.off()

# Fig 3, boxplots with survival and fecundity
Sur2$parameterShort = factor(Sur2$parameter, levels=c("mutAcc", "WCWT",   "WCBT" ,  "BCWT",   "BCBT"),
                          labels = c("Baseline", "WCWT",   "WCBT" ,  "BCWT",   "BCBT"))                          
                          
plot2a <- ggplot(Sur2, aes(ext,x, color= caste))+
  stat_summary(fun.data=median_hilow, fun.args=(conf.int=1), geom="pointrange", size = 1) +
  geom_jitter(shape=16, position=position_jitter(width = 0.2, height = 0), alpha = 0.2)+
  scale_color_manual(name = "caste", labels = c("Queen", "Worker"), 
                     values=c(plasma(20))[c(1,12)]) +
  scale_fill_manual(name = "caste", labels = c("Queen", "Worker"), 
                    values=c(plasma(20))[c(1,12)]) +
  facet_grid(~parameterShort,labeller=label_wrap_gen(width = 10,multi_line = TRUE) )+
  theme_minimal()+
  labs(y="Lifespan")+
  scale_y_continuous(limits = c(0,20)) +
  theme(
    plot.title = element_blank(),
    axis.text=element_text(size=15,face="plain",color="black"),
    legend.position = c(0.9,0.5),
    legend.text = element_text(colour="black", size=17, face="plain",angle=0), 
    legend.background = element_rect(color = NA,fill="transparent", size = 0.2, linetype = "solid"), 
    legend.title = element_blank(),
    axis.title = element_text(size = 21),
    axis.line = element_line(color="black", size = 1.5),
    panel.border = element_rect(colour = "darkgray", fill=NA, size=1),
    plot.background =element_blank(), 
    legend.key.width = unit(1.1,"cm"), 
    panel.background = element_rect(fill = "white", colour = "grey50"),
    strip.text.x = element_text(size= 17, color = "black", face = "plain"),
    strip.text.y =element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x= element_blank(),
    axis.title.x = element_blank()
  )

ann_text <- data.frame(ext = rep(c("0","0.2"), 5), x = rep(15, 10),
                       parameterShort = factor(c(rep("Baseline", 2),rep("WCWT", 2),rep("WCBT", 2),rep("BCWT", 2),rep("BCBT", 2)),levels = c("Baseline","WCWT","WCBT", "BCWT", "BCBT")), 
                       caste = rep("queen", 10))
plot2a <- plot2a + geom_text(data = ann_text,label = c("a", "a", "bc", "c", "b", "c", "bc", "c", "ad", "d"), color = "black")

ann_text2 <- data.frame(ext = rep(c("0","0.2"), 5), x = rep(2.5, 10),
                       parameterShort = factor(c(rep("Baseline", 2),rep("WCWT", 2),rep("WCBT", 2),rep("BCWT", 2),rep("BCBT", 2)),levels = c("Baseline","WCWT","WCBT", "BCWT", "BCBT")), 
                       caste = rep("worker", 10))
plot2a <- plot2a + geom_text(data = ann_text2,label = c("a", "a", "b", "b", "a", "ab", "c", "c", "c", "c"), color = "black ")


plot2a
plot2b <- plot2a + theme(legend.position = "none")
 

# take legend 
leg1 <- get_legend(plot2a) 
leg2 <- get_legend(plot2a) %>% as_ggplot() %>% ggplotGrob()
as_ggplot(leg1)
grid.draw(leg2)

plot3 <- ggplot(Fec2, aes(ext,x, color= caste))+
  stat_summary(fun.data=median_hilow, fun.args=(conf.int=1),
               geom="pointrange", color = plasma(20)[1], size = 1) +
  geom_jitter(shape=16, position=position_jitter(width = 0.2, height = 0), color = plasma(20)[1], alpha = 0.2)+
  scale_y_continuous(limits=c(0, 100)) +
  scale_x_discrete(labels = c("0.0", "0.2", "0.0", "0.2","0.0", "0.2","0.0", "0.2","0.0", "0.2"))+
  labs(x="Extrinsic mortality",y="Fecundity") +
  facet_grid(~parameter2)+
  theme_minimal()+
  theme(
    plot.title = element_blank(),
    axis.text=element_text(size=15,face="plain",color="black"),
    legend.position = "null",#c(0.83,0.76),
    axis.title = element_text(size = 21),
    axis.line = element_line(color="black", size = 1.5),
    panel.border = element_rect(colour = "darkgray", fill=NA, size=1),
    plot.background =element_blank(), 
    legend.key.width = unit(1.1,"cm"), 
    panel.background = element_rect(fill = "white", colour = "grey50"),
    strip.text.x = element_blank(),
    strip.text.y =element_blank()
  )

ann_text <- data.frame(ext = rep(c("0","0.2"), 5), x = rep(62.5, 10),
                       parameter2 = factor(c(rep("mutAcc", 2),rep("WCWT", 2),rep("WCBT", 2),rep("BCWT", 2),rep("BCBT", 2)),levels = c("mutAcc","WCWT","WCBT", "BCWT", "BCBT")),
                       caste = rep("queen", 10))
plot3 <- plot3 + geom_text(data = ann_text,label = c("a", "a", "b", "bc", "b", "bc", "d", "d", "c", "bc"), color = "black")


fig3 <- align_plots(plot2b, plot3, align = "v")
combined_plots <- plot_grid(fig3[[1]], fig3[[2]], labels = c("A", "B"), ncol = 1, rel_heights = c(0.6,0.45), label_size = 20, vjust = c(1.2,-0.2))
vpleg <- viewport(width = 0.5, height = 0.6, x=0.69, y=0.2)
pushViewport(vpleg)
grid.draw(leg2)

pdf(file="lifespan_boxplot.pdf", onefile=T,  width = 8, height = 6) 
plot(combined_plots)
vpleg <- viewport(width = 0.5, height = 0.6, x=0.69,y=0.2)
pushViewport(vpleg)
grid.draw(leg2)
dev.off()

