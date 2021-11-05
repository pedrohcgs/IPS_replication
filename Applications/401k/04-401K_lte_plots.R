###########################################################################
# Application: 401k - LTE Plots
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
load(here("Applications/401k/Results/401k_lte.RData"))
#----------------------------------------------------------------------------

tau = seq(0.1,0.9, 0.01)
tw = data401k$tw
nfa = data401k$net_tfa

# exponential
lqte.ps.exp2.tw <- IPS::LQTE(tw,
                             z, treat, 
                             baselineX, 
                             ps.exp$fitted.values, 
                             ps.exp$lin.rep,
                             tau, bw = "nrd0",
                             whs = NULL )
lqte.ps.exp2.nfa <- IPS::LQTE(nfa,
                              z, treat, 
                              baselineX, 
                              ps.exp$fitted.values, 
                              ps.exp$lin.rep,
                              tau, bw = "nrd0",
                              whs = NULL )

# proj
lqte.ps.proj2.tw <- IPS::LQTE(tw, 
                              z,treat, 
                              baselineX, 
                              ps.proj$fitted.values, 
                              ps.proj$lin.rep,
                              tau, bw = "nrd0",
                              whs = NULL )

lqte.ps.proj2.nfa <- IPS::LQTE(nfa, 
                               z,treat, 
                               baselineX, 
                               ps.proj$fitted.values, 
                               ps.proj$lin.rep,
                               tau, bw = "nrd0",
                               whs = NULL )
#CBPS Just identified
cbps.lin.rep <- inflc_cbps(z, baselineX, 
                           ps.cbps$fitted.values, method = "exact")

lqte.ps.cbps2.tw <- IPS::LQTE(tw, 
                              z,treat, 
                              baselineX, 
                              ps.cbps$fitted.values, 
                              cbps.lin.rep,
                              tau, bw = "nrd0",
                              whs = NULL )

lqte.ps.cbps2.nfa <- IPS::LQTE(nfa, 
                               z,treat, 
                               baselineX, 
                               ps.cbps$fitted.values, 
                               cbps.lin.rep,
                               tau, bw = "nrd0",
                               whs = NULL )

#CBPS Overidentified
cbps.over.lin.rep <- inflc_cbps(z, baselineX, 
                                ps.cbps.over$fitted.values, method = "over")

lqte.ps.over.cbps2.tw <- IPS::LQTE(tw, 
                                   z,treat, 
                                   baselineX, 
                                   ps.cbps.over$fitted.values, 
                                   cbps.over.lin.rep,
                                   tau, bw = "nrd0",
                                   whs = NULL )


lqte.ps.over.cbps2.nfa <- IPS::LQTE(nfa, 
                                    z,treat, 
                                    baselineX, 
                                    ps.cbps.over$fitted.values, 
                                    cbps.over.lin.rep,
                                    tau, bw = "nrd0",
                                    whs = NULL )

#MLE
glm.lin.rep <- inflc_glm(z, baselineX, 
                         ps.glm$fitted.values)

lqte.ps.glm2.tw <- IPS::LQTE(tw, 
                             z,treat, 
                             baselineX, 
                             ps.glm$fitted.values, 
                             glm.lin.rep,
                             tau, bw = "nrd0",
                             whs = NULL )

lqte.ps.glm2.nfa <- IPS::LQTE(nfa, 
                              z,treat, 
                              baselineX, 
                              ps.glm$fitted.values, 
                              glm.lin.rep,
                              tau, bw = "nrd0",
                              whs = NULL )


dt.plot.tw = data.frame(tau,IPS_exp = lqte.ps.exp2.tw$lqte, 
                        IPS_proj =lqte.ps.proj2.tw$lqte,
                        CBPS_just =lqte.ps.cbps2.tw$lqte,
                        CBPS_over =lqte.ps.over.cbps2.tw$lqte,
                        ML =lqte.ps.glm2.tw$lqte)

dt.plot.gg <- melt(dt.plot.tw, id="tau")  # convert to long format

p1 <- ggplot(data=dt.plot.gg, aes(x=tau, y=value, colour=variable, linetype = variable)) +
  geom_line(size = 1) + 
  scale_x_continuous(breaks=seq(0.1,0.9,.1)) +
  expand_limits(y=0) +
  expand_limits(y=30000) +
  labs(x="Quantile index", y="Total Wealth (US$)") +
  ggtitle("LQTE of 401(k) participation") 

p1



max.cbps = max( lqte.ps.over.cbps2.tw$lqte,  lqte.ps.over.cbps2.nfa$lqte)
max.diff = max( lqte.ps.over.cbps2.tw$lqte - lqte.ps.proj2.tw$lqte,
                lqte.ps.over.cbps2.nfa$lqte - lqte.ps.proj2.nfa$lqte)+500


p2.tw <- ggplot() + 
  geom_line(aes(y = lqte.ps.over.cbps2.tw$lqte, x = tau, linetype = "dashed"), size=.5)+ 
  scale_x_continuous(breaks=seq(0.1,0.9,.1)) +
  expand_limits(y=0) +
  expand_limits(y=max.cbps) +
  labs(x="Quantile index", y="US$") +
  ggtitle("LQTE on total wealth") + 
  geom_line(aes(y = lqte.ps.proj2.tw$lqte, x = tau, linetype = "solid"), size=.5) +
  geom_line(aes(y = lqte.ps.glm2.tw$lqte, x = tau, linetype = "dotted"), size=.5) +
  theme(legend.title = element_blank()) +
  guides(linetype=guide_legend(keywidth = 2, keyheight = 1),
         colour=guide_legend(keywidth = 2, keyheight = 1)) +
  scale_linetype_manual(name="linetype", values=c("dashed", "solid", "dotted"),
                        labels = unname(TeX(c("$CBPS_{over}$   ", 
                                              "$LIPS_{proj}$   ",
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
  geom_line(aes(y = (lqte.ps.over.cbps2.tw$lqte - lqte.ps.proj2.tw$lqte), x = tau, 
                linetype = "dashed"), size=.5)+ 
  scale_x_continuous(breaks=seq(0.1,0.9,.1)) +
  expand_limits(y=0) +
  expand_limits(y=max.diff) +
  labs(x="Quantile index", y="US$") +
  ggtitle("Difference between LQTE's based on the overidentified 
CBPS and based on the LIPS with projection weigthing 
function. Outcome: Total wealth.") + 
  theme(legend.position = "none") +
  theme(legend.text = element_text( size = 8))+
  theme(legend.title = element_text(size=8))+
  theme(legend.title=element_blank())+
  theme(plot.title = element_text(color="darkgray", face="bold", size=8)) + 
  theme(axis.title = element_text(color="black",  size=8))

p2.tw.diff



p2.nfa <- ggplot() + 
  geom_line(aes(y = lqte.ps.over.cbps2.nfa$lqte, x = tau, linetype = "dashed"), size=.5)+ 
  scale_x_continuous(breaks=seq(0.1,0.9,.1)) +
  expand_limits(y=0) +
  expand_limits(y=max.cbps) +
  labs(x="Quantile index", y="US$") +
  ggtitle("LQTE on net financial assets") + 
  geom_line(aes(y = lqte.ps.proj2.nfa$lqte, x = tau, linetype = "solid"), size=.5) +
  geom_line(aes(y = lqte.ps.glm2.nfa$lqte, x = tau, linetype = "dotted"), size=.5) +
  theme(legend.title = element_blank()) +
  guides(linetype=guide_legend(keywidth = 2, keyheight = 1),
         colour=guide_legend(keywidth = 2, keyheight = 1)) +
  scale_linetype_manual(name="linetype", values=c("dashed", "solid", "dotted"),
                        labels = unname(TeX(c("$CBPS_{over}$   ", 
                                              "$LIPS_{proj}$   ",
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
  geom_line(aes(y = (lqte.ps.over.cbps2.nfa$lqte - lqte.ps.proj2.nfa$lqte), x = tau, 
                linetype = "dashed"), size=.5)+ 
  scale_x_continuous(breaks=seq(0.1,0.9,.1)) +
  expand_limits(y=0) +
  expand_limits(y=max.diff) +
  labs(x="Quantile index", y="US$") +
  ggtitle("Difference between LQTE's based on the overidentified 
CBPS and based on the LIPS with projection weigthing 
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

ggsave(here("Applications/401k/plots/plots-lqte.pdf"),
       width = 8, 
       height = 6, 
       units = "in")



