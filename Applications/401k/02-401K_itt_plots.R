###########################################################################
# Application: 401k - ITT Plots
###########################################################################
# Startup - clear memory, load packages, and set parameters 
# Clear memory
rm(list = ls(all = T))
#----------------------------------------------------------------------------
# Set seed
set.seed(1234)
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# Install and load packages
library(reshape2)
library(ggplot2)
library(latex2exp)
library(ggpubr)
library(here)
#-----------------------------------------------------------------------------
# ggplot2 theme
theme_set(
  #theme_clean() + 
  theme_classic() +
    theme(plot.background = element_blank(),
          legend.background = element_rect(color = "white"),
          panel.grid.major.y = element_line(colour = "grey95"), 
          panel.grid.minor.y = element_blank(),
          panel.spacing = unit(10, "lines"))
)
#----------------------------------------------------------------------------
# Load ITT results
load(here("Applications/401k/Results/401k_itt.RData"))
#----------------------------------------------------------------------------

# quantile parameters and outcomes
tau = seq(0.1,0.9, 0.01)
tw = data401k$tw
nfa = data401k$net_tfa

# exponential
qte.ps.exp2.tw <- IPS::QTE(tw,
                           treat, 
                           baselineX, 
                           ps.exp$fitted.values, 
                           ps.exp$lin.rep,
                           tau, bw = "nrd0",
                           whs = NULL )
qte.ps.exp2.nfa <- IPS::QTE(nfa,
                            treat, 
                            baselineX, 
                            ps.exp$fitted.values, 
                            ps.exp$lin.rep,
                            tau, bw = "nrd0",
                            whs = NULL )

# proj
qte.ps.proj2.tw <- IPS::QTE(tw, 
                            treat, 
                            baselineX, 
                            ps.proj$fitted.values, 
                            ps.proj$lin.rep,
                            tau, bw = "nrd0",
                            whs = NULL )

qte.ps.proj2.nfa <- IPS::QTE(nfa, 
                             treat, 
                             baselineX, 
                             ps.proj$fitted.values, 
                             ps.proj$lin.rep,
                             tau, bw = "nrd0",
                             whs = NULL )
#CBPS Just identified
cbps.lin.rep <- inflc_cbps(treat, baselineX, 
                           ps.cbps$fitted.values, method = "exact")

qte.ps.cbps2.tw <- IPS::QTE(tw, 
                            treat, 
                            baselineX, 
                            ps.cbps$fitted.values, 
                            cbps.lin.rep,
                            tau, bw = "nrd0",
                            whs = NULL )

qte.ps.cbps2.nfa <- IPS::QTE(nfa, 
                             treat, 
                             baselineX, 
                             ps.cbps$fitted.values, 
                             cbps.lin.rep,
                             tau, bw = "nrd0",
                             whs = NULL )

#CBPS Overidentified
cbps.over.lin.rep <- inflc_cbps(treat, baselineX, 
                                ps.cbps.over$fitted.values, method = "over")

qte.ps.over.cbps2.tw <- IPS::QTE(tw, 
                                 treat, 
                                 baselineX, 
                                 ps.cbps.over$fitted.values, 
                                 cbps.over.lin.rep,
                                 tau, bw = "nrd0",
                                 whs = NULL )


qte.ps.over.cbps2.nfa <- IPS::QTE(nfa, 
                                  treat, 
                                  baselineX, 
                                  ps.cbps.over$fitted.values, 
                                  cbps.over.lin.rep,
                                  tau, bw = "nrd0",
                                  whs = NULL )

#MLE
glm.lin.rep <- inflc_glm(treat, baselineX, 
                         ps.glm$fitted.values)

qte.ps.glm2.tw <- IPS::QTE(tw, 
                           treat, 
                           baselineX, 
                           ps.glm$fitted.values, 
                           glm.lin.rep,
                           tau, bw = "nrd0",
                           whs = NULL )

qte.ps.glm2.nfa <- IPS::QTE(nfa, 
                            treat, 
                            baselineX, 
                            ps.glm$fitted.values, 
                            glm.lin.rep,
                            tau, bw = "nrd0",
                            whs = NULL )


dt.plot.tw = data.frame(tau,IPS_exp = qte.ps.exp2.tw$qte, 
                        IPS_proj =qte.ps.proj2.tw$qte,
                        CBPS_just =qte.ps.cbps2.tw$qte,
                        CBPS_over =qte.ps.over.cbps2.tw$qte,
                        ML =qte.ps.glm2.tw$qte)

dt.plot.gg <- melt(dt.plot.tw, id="tau")  # convert to long format

p1 <- ggplot(data=dt.plot.gg, aes(x=tau, y=value, colour=variable, linetype = variable)) +
  geom_line(size = 1) + 
  scale_x_continuous(breaks=seq(0.1,0.9,.1)) +
  expand_limits(y=0) +
  expand_limits(y=30000) +
  labs(x="Quantile index", y="Total Wealth (US$)") +
  ggtitle("QTE of 401(k) participation") 

p1



max.cbps = max( qte.ps.over.cbps2.tw$qte,  qte.ps.over.cbps2.nfa$qte)
max.diff = max( qte.ps.over.cbps2.tw$qte - qte.ps.proj2.tw$qte,
                qte.ps.over.cbps2.nfa$qte - qte.ps.proj2.nfa$qte)+500


p2.tw <- ggplot() + 
  geom_line(aes(y = qte.ps.over.cbps2.tw$qte, x = tau, linetype = "dashed"), size=.5)+ 
  scale_x_continuous(breaks=seq(0.1,0.9,.1)) +
  expand_limits(y=0) +
  expand_limits(y=max.cbps) +
  labs(x="Quantile index", y="US$") +
  ggtitle("QTE on total wealth") + 
  geom_line(aes(y = qte.ps.proj2.tw$qte, x = tau, linetype = "solid"), size=.5) +
  geom_line(aes(y = qte.ps.glm2.tw$qte, x = tau, linetype = "dotted"), size=.5) +
  theme(legend.title = element_blank()) +
  guides(linetype=guide_legend(keywidth = 2, keyheight = 1),
         colour=guide_legend(keywidth = 2, keyheight = 1)) +
  scale_linetype_manual(name="linetype", values=c("dashed", "solid", "dotted"),
                        labels = unname(TeX(c("$CBPS_{over}$   ", 
                                              "$IPS_{proj}$   ",
                                              "$MLE$"))))+
  theme(legend.text = element_text( size = 8),
        legend.spacing.x = unit(.2, 'cm'))+
  theme(legend.title = element_text(size=8))+
  theme(legend.title=element_blank())+
  theme(legend.position="bottom")+
  theme(plot.title = element_text(color="darkgray", face="bold", size=8)) + 
  theme(axis.title = element_text(color="black",  size=8))

p2.tw




p2.tw.diff <- ggplot() + 
  geom_line(aes(y = (qte.ps.over.cbps2.tw$qte - qte.ps.proj2.tw$qte), x = tau, 
                linetype = "dashed"), size=.5)+ 
  scale_x_continuous(breaks=seq(0.1,0.9,.1)) +
  expand_limits(y=0) +
  expand_limits(y=max.diff) +
  labs(x="Quantile index", y="US$") +
  ggtitle("Difference between QTE's based on the overidentified 
CBPS and based on the IPS with projection weigthing 
function. Outcome: Total wealth.") + 
  theme(legend.position = "none") +
  theme(legend.text = element_text( size = 8))+
  theme(legend.title = element_text(size=8))+
  theme(legend.title=element_blank())+
  theme(plot.title = element_text(color="darkgray", face="bold", size=8)) + 
  theme(axis.title = element_text(color="black",  size=8))

p2.tw.diff



p2.nfa <- ggplot() + 
  geom_line(aes(y = qte.ps.over.cbps2.nfa$qte, x = tau, linetype = "dashed"), size=.5)+ 
  scale_x_continuous(breaks=seq(0.1,0.9,.1)) +
  expand_limits(y=0) +
  expand_limits(y=max.cbps) +
  labs(x="Quantile index", y="US$") +
  ggtitle("QTE on net financial assets") + 
  geom_line(aes(y = qte.ps.proj2.nfa$qte, x = tau, linetype = "solid"), size=.5) +
  geom_line(aes(y = qte.ps.glm2.nfa$qte, x = tau, linetype = "dotted"), size=.5) +
  theme(legend.title = element_blank()) +
  guides(linetype=guide_legend(keywidth = 2, keyheight = 1),
         colour=guide_legend(keywidth = 2, keyheight = 1)) +
  scale_linetype_manual(name="linetype", values=c("dashed", "solid", "dotted"),
                        labels = unname(TeX(c("$CBPS_{over}$   ", 
                                              "$IPS_{proj}$   ",
                                              "$MLE$"))))+
  theme(legend.text = element_text( size = 8),
        legend.spacing.x = unit(.2, 'cm'))+
  theme(legend.title = element_text(size=8))+
  theme(legend.title=element_blank())+
  theme(legend.position="bottom")+
  theme(plot.title = element_text(color="darkgray", face="bold", size=8)) + 
  theme(axis.title = element_text(color="black",  size=8))

p2.nfa



p2.nfa.diff <- ggplot() + 
  geom_line(aes(y = (qte.ps.over.cbps2.nfa$qte - qte.ps.proj2.nfa$qte), x = tau, 
                linetype = "dashed"), size=.5)+ 
  scale_x_continuous(breaks=seq(0.1,0.9,.1)) +
  expand_limits(y=0) +
  expand_limits(y=max.diff) +
  labs(x="Quantile index", y="US$") +
  ggtitle("Difference between QTE's based on the overidentified 
CBPS and based on the IPS with projection weigthing 
function. Outcome: Net financial assets.") + 
  theme(legend.position = "none") +
  theme(legend.text = element_text( size = 8))+
  theme(legend.title = element_text(size=8))+
  theme(legend.title=element_blank())+
  theme(plot.title = element_text(color="darkgray", face="bold", size=8)) + 
  theme(axis.title = element_text(color="black",  size=8))

p2.nfa.diff



ggarrange(p2.nfa,  p2.tw , p2.nfa.diff, p2.tw.diff,
          ncol = 2, nrow = 2)

ggsave(here("Applications/401k/plots/plots-qte.pdf"),
       width = 8, 
       height = 6, 
       units = "in")


