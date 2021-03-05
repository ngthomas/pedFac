### polygamous study - sampling fraction


###tasks to do:
#need to clean up full sib half sib for sequoia
#make sure sf>0, sequoia can tell half-sib

#correct final ROC rate when positive is not completed
#check out half-sib relationship for index plot for the tail end ..looks funny for pedFac



library(tidyverse)
library(xtable)
library(kableExtra)
library(gridExtra)

load("/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/fullsib/halfsib.sf.pedfac.186v")
load("/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/fullsib/halfsib.sf.COLONY.rerun")
load("/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/fullsib/halfsib.sf.sequoia.rerun")

Longo <- expand_grid(
  samp_frac_idx = 1:5,
  rep.indx = 1:5
) %>%
  mutate(
    sampl.frac = c(1, 0.75, 0.5, 0.25, 0)[samp_frac_idx],
    COLONY = map2(.x = samp_frac_idx, .y = rep.indx, .f = function(x, y) halfsib.sf.COLONY.rerun[[x]][[y]]),
    sequoia = map2(.x = samp_frac_idx, .y = rep.indx, .f = function(x, y) halfsib.sf.sequoia.rerun[[x]][[y]]),
    pedFac = map2(.x = samp_frac_idx, .y = rep.indx, .f = function(x, y) halfsib.sf.pedfac.187v[[x]][[y]])
  ) %>%
  pivot_longer(
    cols = COLONY:pedFac,
    names_to = "software",
    values_to = "output"
  )


# make a function to get runtimes
get_runtime_mins <- function(L) {
  L$runtime/60
}

# get the fraction n.match
fract_nmatch <- function(Tib) {
  mean(Tib$n.match)
}

## making parentage assignment
cleanOutParentTbl <- function(Tib) {
  if ("LLRpair" %in% colnames(Tib)) {
    Tib <- Tib %>% mutate(norm.LLRpair = ifelse(is.na(LLRpair),
                                                0,
                                                exp(LLRpair)),
                          prob=norm.LLRpair/max(norm.LLRpair)) }

  Tib %>%
    arrange(desc(prob)) %>%
    group_by(kid.id) %>% slice_head(n=1) %>%
    rename(num.correct.calls = n.match) %>%
    mutate(pa.id = as.numeric(pa.id),
           ma.id = as.numeric(ma.id),
           n.obs.parents = (pa.T!=-1) +(ma.T!=-1)) %>%
    ungroup() %>%
    arrange(desc(prob)) %>%
    mutate(rank = row_number())
}

# Longo %>%
#   group_by(sampl.frac, rep.indx, software) %>%
#   summarise(
#     mins = map_dbl(output, get_runtime_mins)
#   ) %>% View

## preparing for parentage analysis

compiled.parentage.tbl <- Longo %>% filter(sampl.frac != 0) %>%
  group_by(sampl.frac, rep.indx, software) %>%
  mutate(
    software = ifelse(software == "colony_no_sibprior", "COLONY", software),
    #seconds = map_dbl(output, "runtime"),
    parent_tbl_prob = case_when(
      software %in% c("pedFac") ~ map(output, "parent.prob.tbl"),
      software == "sequoia" ~ map(output, "parent.tbl"),
      TRUE ~ map(output, "parent.tbl")),
    parent_org_tbl_prob = map(parent_tbl_prob, cleanOutParentTbl)) %>%
  ungroup() %>%
  select(-parent_tbl_prob, -output) %>%
  unnest(cols=c(parent_org_tbl_prob))

parentage.stat <- compiled.parentage.tbl %>%
  ungroup() %>%
  mutate(n.miss.obs.parent = (pa.T != -1 & pa.T != pa.id) +
           (ma.T != -1 & ma.T != ma.id),
         n.miss.unobs.parent = (pa.T == -1 & pa.T != pa.id) +
           (ma.T == -1 & ma.T != ma.id)) %>%
  select(rep.indx, sampl.frac, software, rank,
         prob, num.correct.calls, n.obs.parents,
         n.miss.obs.parent, n.miss.unobs.parent)


prob.rank.plot.tbl <- parentage.stat %>% #filter(!software %in% c("pedFac_w_swap","pedFac_w_swap_sub")) %>%
  mutate(software = factor(software, levels=c("pedFac",  "COLONY", "sequoia"), labels = c("pedFac",  "COLONY", "sequoia")),
         sampl.frac = factor(sampl.frac, levels=c(1,0.75,0.5, 0.25), labels = c("sf:\n1.0","sf:\n0.75","sf:\n0.5","sf:\n0.25")))

prob.rank.err.line2.tbl <- prob.rank.plot.tbl %>%
  filter(num.correct.calls != 2) %>%
  mutate(x=rank, xend=rank, y= prob, yend = prob+0.1*(n.miss.obs.parent))
prob.rank.err.line3.tbl <- prob.rank.plot.tbl %>%
  filter(num.correct.calls != 2) %>%
  mutate(x=rank, xend=rank, y= prob, yend = prob-0.1*(n.miss.unobs.parent))

prob.err.cum.sum <- prob.rank.plot.tbl %>%
  arrange(rank) %>%
  group_by(software, sampl.frac, rep.indx) %>%
  mutate(y=cumsum(2-num.correct.calls)/(2*n()))

parentage.prob.ggplot <- list()

for (i in 1:5) {

  parentage.prob.plot <-
    ggplot(data=prob.rank.plot.tbl %>% filter(rep.indx==i)) +
    #geom_segment(data=prob.rank.nobs.line.tbl, aes(x=x, y=y, xend=xend, yend=yend), color =
    #              "grey", alpha=0.4) +
    geom_segment(data=prob.rank.err.line2.tbl %>% filter(rep.indx==i), aes(x=x, y=y, xend=xend, yend=yend), color ="#a45eff" ,alpha=0.8) +
    geom_segment(data=prob.rank.err.line3.tbl %>% filter(rep.indx==i), aes(x=x, y=y, xend=xend, yend=yend),color="#a6d7ff" ,alpha=0.9)+
    geom_line(aes(x=rank, y=prob), )+
    #geom_line(data=prob.err.cum.sum %>% filter(rep.indx==i), aes(x= rank, y= y),color="#5499C7", alpha=0.9)+
    #scale_color_manual(values = c("#5499C7", "#FF5733", "orange"), labels=c("n","d","d1"))+
    facet_grid(sampl.frac~software, switch="y")+
    theme_light()+
    scale_y_continuous("Prob", position = "right", breaks=seq(0,1,0.2), minor_breaks = seq(-0.1,1.1,0.1))+
    scale_x_continuous("Individual index")+
    theme(#strip.placement = "outside",
      strip.background.y = element_rect(fill= "dark grey"),
      strip.background.x = element_rect(fill= "dark grey"),
      strip.text.y.left = element_text(angle = 0,colour = "white"))
  #panel.grid.minor.y = element_blank())

  legend.parentage <- ggplot()+
    geom_text(aes(x=0.08, y=0.6, label= "Parent that is misassigned is:"),size=3.5, hjust="left", vjust="middle")+
    geom_segment(aes(x=0.55, y=0.0, xend=0.55, yend=1), color ="#a45eff" ,alpha=0.8, size=0.9)+
    geom_text(aes(x=0.58, y=0.6, label= "observed"),size=3.5, hjust="left", vjust="middle")+
    geom_segment(aes(x=0.75, y=0.0, xend=0.75, yend=1), color ="#a6d7ff" ,alpha=0.9, size=0.9)+
    geom_text(aes(x=0.78, y=0.6, label= "unobserved"),size=3.5, hjust="left", vjust="middle")+
    scale_y_continuous("", position = "right", breaks=FALSE,limits = c(0,1),minor_breaks = NULL,labels = NULL)+
    scale_x_continuous("", limits=c(0,1))+
    theme(#strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0,colour = "white"),
      panel.grid.minor.y = element_blank(),
      axis.line=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      legend.position="none",
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())

  ggsave(paste0("notes/fig/ch2_poly_parent_indxplot_",i,".pdf"), parentage.prob.plot,  device="pdf", height=5, width=9, units="in", dpi=500)


  parentage.prob.ggplot[[i]] <- grid.arrange(parentage.prob.plot, legend.parentage, ncol=1,heights=c(10.5,1),
                                             padding = unit(0, "line"))
  ggsave(paste0("notes/fig/ch2_poly_parent_indxplot_",i,".png"), parentage.prob.ggplot[[i]],  device="png", height=5, width=9, units="in", dpi=500)

}



## need to define the probability of parental assignment is the fraction of cases/sweep where we observe the parental pair assignmet for the recent generation cohort
## need a legend that the size and orientation of the line
# here we have a line of individuals sorted by a descreasing order of the probability value. We are showing one of five replicates in this parentage assignment study. The column are divided by the method (pedfac plain prior, COLONY - no sibship size prior but have true s.f estimate in param)
#blood orange line indicatesthe  number of observed parents of the indicative individual being missassigned while  (downward line) in gold indiciates the unobserved missasigned parent
#


## prepare for the ROC curve

roc.line.df <- parentage.stat %>%
  arrange(rank, desc(prob)) %>%
  group_by(software, rep.indx, sampl.frac) %>%
  mutate(n.true.case = cumsum(num.correct.calls == 2),
         n.false.case = cumsum(num.correct.calls != 2),
         tot.true.case = sum(n.true.case),
         tot.false.case = sum(n.false.case)) %>%
  group_by(software, rep.indx, sampl.frac, prob) %>%
  summarise(tpr = ifelse(tot.true.case[1]==0,0,
                         max(n.true.case)/tot.true.case[1]),
            fpr = ifelse(tot.false.case[1]==0,0,
                         max(n.false.case)/tot.false.case[1]),
            n.true.case = tot.true.case[1],
            n.false.case = tot.false.case[1])

roc.pt.fixed.1 <- roc.line.df %>%
  group_by(software, rep.indx, sampl.frac) %>%
  summarise(posterior = 0,
            tpr = 1,
            fpr = 1.0001,
            n.true.case=max(n.true.case),
            n.false.case=max(n.false.case))

roc.pt.fixed.2 <- roc.line.df %>%
  group_by(software, rep.indx, sampl.frac) %>%
  summarise(posterior = 1,
            tpr = -0,
            fpr = -0.00001,
            n.true.case=max(n.true.case[posterior==1]),
            n.false.case=max(n.false.case[posterior==1]))

roc.compiled.df <- bind_rows(roc.pt.fixed.2, roc.line.df, roc.pt.fixed.1)

Calculate.AUC <- function(prob, categ){
  match.score <- prob[categ==2] # num.correct.calls
  mismatch.score <- prob[categ!=2]
  if (length(match.score)==0) return(0)
  if (length(mismatch.score)==0) return(1)
  o <- outer(match.score, mismatch.score, "-")
  auc <- mean((o>0) + .5*(o==0))
}

parentage.stat.tbl <- parentage.stat %>%
  mutate(software = factor(software, levels=rev(c("pedFac",  "COLONY", "sequoia")), labels = rev(c("pedFac",  "COLONY", "sequoia"))),
         sampl.frac = factor(sampl.frac, levels=c(1,0.75,0.5, 0.25), labels = c("sf:\n1.0","sf:\n0.75","sf:\n0.5","sf:\n0.25"))) %>%
  arrange(rank, desc(prob)) %>%
  group_by(software, rep.indx, sampl.frac) %>%
  summarise (AUC = Calculate.AUC(prob, num.correct.calls),
             AR = mean(num.correct.calls/2)) %>%
  pivot_longer(
    cols = AUC:AR,
    names_to = "stat",
    values_to = "score"
  )

# AUC.label.tbl <- parentage.stat.tbl %>%
#   arrange(rank, desc(prob)) %>%
#   group_by(software, rep.indx, sampl.frac) %>%
#   mutate(x=0.8,
#          y=0.2*i+0.1,
#          auc.label = paste0("AUC: ",round(auc.score,2)))

## since most of the assignment is right on target, we skip to report the stat

parentage.boxstat.tbl <- parentage.stat.tbl %>%
  group_by(software, sampl.frac, stat) %>%
  summarise(min = min(score),
            max = max(score),
            mean = mean(score),
            sd = sd(score))

parentage.stat.tbl.full  <- parentage.stat.tbl %>%
  mutate(
    stat = case_when(
      stat=="AUC"~"Area Under the Curve",
      stat=="AR"~"Assignment Rate"
    )
  )

parentage.boxstat.tbl.full <- parentage.boxstat.tbl %>%
  mutate(
    stat = case_when(
      stat=="AUC"~"Area Under the Curve",
      stat=="AR"~"Assignment Rate"
    )
  )

g <- ggplot(data=parentage.stat.tbl.full) +
  geom_rect(data = subset(parentage.stat.tbl.full, sampl.frac %in% c("sf:\n1.0","sf:\n0.5")), aes(xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf), alpha = 0.01, fill="#AED6F1") +
  geom_point(aes(y=software,x=score), alpha=0.2, stat = "identity", size=3)+
  #geom_point(data= parentage.boxstat.tbl, aes(y=software,x=q.25), color="orange", shape=124, size=2)+
  #geom_point(data= parentage.boxstat.tbl, aes(y=software,x=q.75), color="orange", shape=124, size=2)+
  geom_point(data= parentage.boxstat.tbl.full, aes(y=software,x=mean), color="red", shape=124, size=4, stroke=4)+
  facet_grid(sampl.frac~stat,scales = "free_x", switch="y")+
  theme_light()+
  scale_y_discrete("", position = "right")+
  scale_color_discrete(guide=FALSE)+
  #scale_x_continuous("Individual index (ordered by desc. prob)")+
  theme(strip.placement = "outside",
        panel.border = element_rect( fill=NA, size=c(0.2)),
        strip.background.y = element_rect(fill= "dark grey"),
        strip.background.x = element_rect(fill= "dark grey"),
        strip.text.y.left = element_text(angle = 0,colour = "white", size=12),
        strip.text.x.top = element_text(size=11),
        axis.text.y.right = element_text(size=12, angle = 10),
        panel.grid.minor.y = element_blank())

ggsave("notes/fig/ch2_poly_parent_stat.pdf", g, device="pdf", height=5, width=9, units="in", dpi=500)

ggsave("notes/fig/ch2_poly_parent_stat.png", g,  device="png", height=5, width=9, units="in", dpi=500)



###create xtable in respect to graph
parentage.prep.xtable <- parentage.boxstat.tbl %>% ungroup() %>%
  mutate(
    software = factor(software, levels=c("pedFac","COLONY", "sequoia")),
    sampl.frac = factor(sampl.frac, labels=c(1.0,0.75,0.5, 0.25), levels = c("sf:\n1.0","sf:\n0.75","sf:\n0.5","sf:\n0.25")),
         min = round(min, 2),
         max = round(max, 2),
         mean = round(mean, 2),
         sd = round(sd, 2)) %>%
  pivot_wider(names_from = stat, names_glue = "{stat}_{.value}", values_from=min:sd) %>%
  arrange(sampl.frac, desc(software)) %>%
  group_by(sampl.frac, software) %>%
  summarise("AR_mean w/ sd" = sprintf("%.2f +/- %.2f", AR_mean, AR_sd),
            "AR_range" = sprintf("(%.2f,%.2f)", AR_min, AR_max),
            "AUC_mean w/ sd" = sprintf("%.2f +/- %.2f", AUC_mean, AUC_sd),
            "AUC_range" = sprintf("(%.2f,%.2f)", AUC_min, AUC_max)) %>%
  rename(fraction = sampl.frac)

latex.out <- kbl(parentage.prep.xtable, booktabs = T, "latex") %>%
  column_spec(1, bold = T, width = "3.5em", latex_valign= "m", latex_column_spec = "c") %>%
  row_spec(c(1:3, 7:9)-1, extra_latex_after = "\\rowcolor{gray!6}") %>%
  row_spec(0, align="c") %>%
  collapse_rows(1, latex_hline = "none") %>%
  add_header_above(c(" ", " ", "1 - FDR" = 2, "ROC-AUC" = 2),align="c", bold=T)

latex.out.1 <- gsub('\\+/-', "$\\\\pm$", latex.out, perl = T)
latex.out.2 <- gsub('(AUC\\\\_)|(AR\\\\_)', "", latex.out.1, perl = T)
write(latex.out.2,  "notes/table/ch2_poly_parent.tex")

### sibship pairwise analysis

## from this point on forward, we will just stick with the pedfac run with swap and sub moves


standarize_pairwise_tbl <- function(Tib) {
  if ("cluster" %in% colnames(Tib)) {
    Tib <- Tib %>% select(-cluster, -prob.excl) %>%
      rename(prob=prob.incl)
  }
  Tib %>%
    mutate(prob = as.numeric(prob))
}


Plot.SibPair <- function(sib.label = "FS") {

  table.sib.label <- ifelse(sib.label == "FS", "fullsib","halfsib")

compiled.sibship.pairwise.tbl <- Longo %>%
  mutate(
    software = case_when(software == "colony_no_sibprior"~"COLONY",
                         software == "pedFac_w_swap_sub"~"pedFac",
                         TRUE ~ software),
    FS_pairwise_tbl = map(output, paste0(table.sib.label, ".pairwise.tbl")),
    FS_pairwise_tbl = map(FS_pairwise_tbl, standarize_pairwise_tbl)
  ) %>%
  select(-output) %>%
  unnest(cols=c(FS_pairwise_tbl))

## assume two gen ped
total.n.sibs <- unique(c(compiled.sibship.pairwise.tbl$kid.1,
                         compiled.sibship.pairwise.tbl$kid.2)) %>% length

total.n.pairs <- total.n.sibs*(total.n.sibs-1)*0.5

FS.infer.tbl.prep <- compiled.sibship.pairwise.tbl %>%
  ungroup() %>%
  filter(!(is.na(prob) & presence.T ==0)) %>%
  mutate(prob = ifelse(is.na(prob),0,prob),
         prob = ifelse(prob>1,1,prob),
         prob = ifelse(prob<0,0,prob),
         presence.T = ifelse(is.na(presence.T),0,presence.T))

n.low.denom <- FS.infer.tbl.prep %>%
  filter(prob>0, software != "sequoia") %>%
  group_by(rep.indx, sampl.frac, software) %>%
  summarise(n=n()) %>%
  group_by(rep.indx, sampl.frac) %>%
  summarise(cutoff = min(n))

FS.infer.tbl <- left_join(FS.infer.tbl.prep, n.low.denom) %>%
  arrange(desc(prob)) %>%
  group_by(rep.indx, sampl.frac, software) %>%
  mutate(prob = ifelse(prob < 0.5 & software == "pedFac",0,prob)) %>%
  ungroup() %>%
  filter(!(prob == 0 & presence.T ==0))




# # need to conform the 3 methods a bit so that
# FS.infer.tbl %>%
#   select(-presence.T) %>%
#   pivot_wider(id_cols= c(rep.indx, sampl.frac, kid.1, kid.2),
#               names_from = software,
#               values_from = c(prob),
#               values_fill = 0) %>%
#   pivot_longer(
#     cols = sequoia:pedFac,
#     names_to = "software",
#     values_to = "prob"
#     )



#FS.infer.tbl %>%
#  group_by(sampl.frac, rep.indx, software) %>%
#  summarise(n.TP = sum(presence.T == 1),
#            n.FP = sum(presence.T == 0))


prob.rank.FS.tbl <- FS.infer.tbl %>%
  mutate(software = factor(software, levels=c("pedFac",  "COLONY", "sequoia"), labels = c("pedFac",  "COLONY", "sequoia")),
         sampl.frac = factor(sampl.frac, levels=c(1,0.75,0.5, 0.25, 0), labels = c("sf:\n1.0","sf:\n0.75","sf:\n0.5","sf:\n0.25", "sf:\n0"))) %>%
  arrange(desc(prob)) %>%
  group_by(rep.indx, sampl.frac, software) %>%
  mutate(rank = row_number())

prob.rank.FS.err.line.down <- prob.rank.FS.tbl %>%
  filter(presence.T == 1, prob == 0) %>%
  mutate(x=rank, xend=rank, y= prob, yend = prob-0.1)
prob.rank.FS.err.line.up <- prob.rank.FS.tbl %>%
  filter(presence.T == 0) %>%
  mutate(x=rank, xend=rank, y= prob, yend = prob+0.1)

prob.rank.FS.cumsum <- prob.rank.FS.tbl %>%
  arrange(rank) %>%
  group_by(software, sampl.frac, rep.indx) %>%
  mutate(y=cumsum(prob>0 & presence.T == 1)/nrow(full.sib.pair.truth.tbl))

prob.rank.FS.cumsum <- prob.rank.FS.tbl %>%
  arrange(rank) %>%
  group_by(software, sampl.frac, rep.indx) %>%
  mutate(y=cumsum(prob>0 & presence.T == 1)/nrow(full.sib.pair.truth.tbl))

FS.prob.ggplot <- list()

for (i in 1:5) {
  FS.prob.plot <-
    ggplot(data=prob.rank.FS.tbl %>% filter(rep.indx==i)) +
    geom_segment(data=prob.rank.FS.err.line.up %>% filter(rep.indx==i), aes(x=x, y=y, xend=xend, yend=yend), color ="#FF5733" ,alpha=0.8) +
    geom_segment(data=prob.rank.FS.err.line.down %>% filter(rep.indx==i), aes(x=x, y=y, xend=xend, yend=yend),color="orange" ,alpha=0.9)+
    geom_line(aes(x=rank, y=prob))+
    #geom_line(data=prob.rank.FS.cumsum %>% filter(rep.indx==i), aes(x= rank, y= y),color="#5499C7", alpha=0.9)+
    #scale_color_manual(values = c("#5499C7", "#FF5733", "orange"), labels=c("n","d","d1"))+
    facet_grid(sampl.frac~software, switch="y",scales = "free_x")+
    theme_light()+
    scale_y_continuous("Prob", position = "right", breaks=seq(0,1,0.2))+
    theme(#strip.placement = "outside",
      strip.background.y = element_rect(fill= "dark grey"),
      strip.background.x = element_rect(fill= "dark grey"),
      strip.text.y.left = element_text(angle = 0,colour = "white"),
      panel.grid.minor.y = element_blank())

  if (sib.label == "HS") {
    FS.prob.plot <- FS.prob.plot +
      scale_x_continuous("Observed Half-Sibling (HS) Pairs Index")
  } else {
    FS.prob.plot <- FS.prob.plot +
      scale_x_continuous("Observed Full-Sibling (FS) Pairs Index")
  }

  legend.FS <- ggplot()+
    geom_segment(aes(x=0.65, y=0.0, xend=0.65, yend=1), color ="orange" ,alpha=0.8, size=0.9)+
    geom_text(aes(x=0.68, y=0.6, label= "Missing True Assignment"),size=3.5, hjust="left", vjust="middle")+
    geom_segment(aes(x=0.05, y=0.0, xend=0.05, yend=1), color ="#FF5733" ,alpha=0.9, size=0.9)+
    scale_y_continuous("", position = "right", breaks=FALSE,limits = c(0,1),minor_breaks = NULL,labels = NULL)+
    scale_x_continuous("", limits=c(0,1))+
    theme(#strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0,colour = "white"),
      panel.grid.minor.y = element_blank(),
      axis.line=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      legend.position="none",
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())

  if (sib.label == "HS") {
    legend.FS <- legend.FS +
      geom_text(aes(x=0.08, y=0.6, label= "False HS Assignment"),size=3.5, hjust="left", vjust="middle")
  } else {
    legend.FS <- legend.FS +
      geom_text(aes(x=0.08, y=0.6, label= "False FS Assignment"),size=3.5, hjust="left", vjust="middle")
  }


  FS.prob.ggplot[[i]] <- grid.arrange(FS.prob.plot, legend.FS, ncol=1,heights=c(10.5,1),
                                      padding = unit(0, "line"))

  ggsave(paste0("notes/fig/ch2_poly_",sib.label,"Pair_indxplot_",i,".pdf"), FS.prob.plot, device="pdf", height=5, width=9, units="in", dpi=500)

  ggsave(paste0("notes/fig/ch2_poly_",sib.label,"Pair_indxplot_",i,".png"), FS.prob.ggplot[[i]],  device="png", height=5, width=9, units="in", dpi=500)


}

filler.FS.tbl <- prob.rank.FS.tbl %>%
  filter(presence.T == 0) %>%
  group_by(samp_frac_idx, rep.indx, sampl.frac, kid.1, kid.2) %>%
  summarise(has.pedFac = ifelse("pedFac" %in% software, TRUE, FALSE ),
            has.COLONY = ifelse("COLONY" %in% software , TRUE, FALSE ),
            has.sequoia = ifelse("sequoia" %in% software, TRUE, FALSE ),
            n.presence = n())


prob.rank.FS.tbl.1 <- bind_rows(
  prob.rank.FS.tbl,
  filler.FS.tbl %>% filter(!has.pedFac) %>%
    select(samp_frac_idx:kid.2) %>% ungroup() %>%
    mutate(software = "pedFac",
           prob = 0,
           presence.T = 0,
           rank = 10000),
  filler.FS.tbl %>% filter(!has.COLONY) %>%
    select(samp_frac_idx:kid.2) %>% ungroup() %>%
    mutate(software = "COLONY",
           prob = 0,
           presence.T = 0,
           rank = 10000),
  filler.FS.tbl %>% filter(!has.sequoia) %>%
    select(samp_frac_idx:kid.2) %>% ungroup() %>%
    mutate(software = "sequoia",
           prob = 0,
           presence.T = 0,
           rank = 10000)
)%>%
  mutate(software = factor(software, levels=c("pedFac",  "COLONY", "sequoia"), labels = c("pedFac",  "COLONY", "sequoia")),
         sampl.frac = factor(sampl.frac, levels= c("sf:\n1.0","sf:\n0.75","sf:\n0.5","sf:\n0.25", "sf:\n0"))) %>%
  arrange(desc(prob)) %>%
  group_by(rep.indx, sampl.frac, software) %>%
  mutate(rank = row_number())


## ROC plot

roc.FS.mean.line.df <- prob.rank.FS.tbl.1 %>%
  arrange(desc(prob)) %>%
  group_by(software, sampl.frac) %>%
  mutate(n.true.case = cumsum(presence.T == 1),
         n.false.case = cumsum(presence.T != 1),
         tot.true.case = sum(presence.T == 1),
         tot.false.case = sum(presence.T != 1)) %>%
  group_by(software, sampl.frac, prob) %>%
  summarise(tpr = ifelse(tot.true.case[1]==0,0,
                         max(n.true.case)/tot.true.case[1]),
            fpr = ifelse(tot.false.case[1]==0,0,
                         max(n.false.case)/tot.false.case[1]),
            n.true.case = tot.true.case[1],
            n.false.case = tot.false.case[1])

roc.FS.mean.pt.fixed.1 <- roc.FS.mean.line.df %>%
  group_by(software, sampl.frac) %>%
  summarise(prob = 0,
            tpr = 1,
            fpr = 1,
            n.true.case=max(n.true.case),
            n.false.case=max(n.false.case))

roc.FS.mean.pt.fixed.2 <- roc.FS.mean.line.df %>%
  group_by(software,  sampl.frac) %>%
  summarise(prob = 1,
            tpr = 0,
            fpr = 0,
            n.true.case=max(n.true.case[prob==1]),
            n.false.case=max(n.false.case[prob==1]))


roc.FS.line.df <- prob.rank.FS.tbl.1 %>%
  arrange(desc(prob)) %>%
  group_by(software, rep.indx, sampl.frac) %>%
  mutate(n.true.case = cumsum(presence.T == 1),
         n.false.case = cumsum(presence.T != 1),
         tot.true.case = sum(presence.T == 1),
         tot.false.case = sum(presence.T != 1)) %>%
  group_by(software, rep.indx, sampl.frac, prob) %>%
  summarise(tpr = ifelse(tot.true.case[1]==0,0,
                         max(n.true.case)/tot.true.case[1]),
            fpr = ifelse(tot.false.case[1]==0,0,
                         max(n.false.case)/tot.false.case[1]),
            n.true.case = tot.true.case[1],
            n.false.case = tot.false.case[1])

roc.FS.pt.fixed.1 <- roc.FS.line.df %>%
  group_by(software, rep.indx, sampl.frac) %>%
  summarise(prob = 0,
            tpr = 1,
            fpr = 1.0001,
            n.true.case=max(n.true.case),
            n.false.case=max(n.false.case))

roc.FS.pt.fixed.2 <- roc.FS.line.df %>%
  group_by(software, rep.indx, sampl.frac) %>%
  summarise(prob = 1,
            tpr = -0,
            fpr = -0.00001,
            n.true.case=max(n.true.case[prob==1]),
            n.false.case=max(n.false.case[prob==1]))

roc.FS.compiled.df <- bind_rows(roc.FS.pt.fixed.2, roc.FS.line.df, roc.FS.pt.fixed.1)

roc.mean.FS.compiled.df <- bind_rows(roc.FS.mean.pt.fixed.2, roc.FS.mean.line.df, roc.FS.mean.pt.fixed.1) %>%
  arrange(tpr)

Calculate.FS.AUC <- function(prob, categ){
  match.score <- prob[categ==1] # num.correct.calls
  mismatch.score <- prob[categ==0]
  if (length(match.score)==0) return(0)
  if (length(mismatch.score)==0) return(1)
  o <- outer(match.score, mismatch.score, "-")
  auc <- mean((o>0) + .5*(o==0))
  auc
}

AUC.FS.label <- prob.rank.FS.tbl.1 %>% #filter(prob >0) %>%
  arrange(desc(prob)) %>%
  group_by(software, rep.indx, sampl.frac) %>%
  summarise(x=0.8,
            y=0.1,
            auc.score = Calculate.FS.AUC(prob, presence.T)) %>%
  group_by(software, sampl.frac, x, y) %>%
  summarise(auc.label = sprintf("mean AUC: %.3f",signif(mean(auc.score),3)))


#instead of label every minimal threshold,
#it is to space out the labeling
# label.grid.size <- 0.1
# xmin.grid <- seq(0,1-label.grid.size,label.grid.size)
# xmax.grid <- seq(label.grid.size,1,label.grid.size)
# len.grid <- length(xmin.grid)
#
# roc.FS.prob.label <- lapply(1:len.grid,
#                        function(i) roc.FS.line.df %>%
#                          group_by(rep.indx, software, sampl.frac) %>%
#                          filter(fpr>=xmin.grid[i] & fpr<=xmax.grid[i]) %>%
#                          summarise(x=fpr[1],
#                                    y=tpr[1],
#                                    prob=signif(prob[1],2),
#                                    row.number = i%%2)
# ) %>%
#   bind_rows() %>%
#   filter(!is.na(x))

#FS.ROC.ggplot <- list()

#for (i in 1:5) {

#manually removed HS for sequoia since it doesn't call any HS
if(sib.label == "HS") {
  roc.FS.compiled.df <- roc.FS.compiled.df %>%
    filter(software != "sequoia" || sampl.frac != "sf:\n0")

  AUC.FS.label <- AUC.FS.label %>%
    filter(software != "sequoia" || sampl.frac != "sf:\n0")
}

FS.ROC.ggplot <- ggplot(roc.FS.compiled.df %>% arrange(fpr, tpr), aes(x = fpr, y = tpr)) +
  geom_line(aes(group=factor(rep.indx)), alpha=0.5, color="grey", size=0.6) +
  geom_line(data= roc.mean.FS.compiled.df, aes(x = fpr, y = tpr), alpha=1, color="black") +

  #geom_smooth(aes(group=factor(rep.indx)), method = "gam", se = TRUE)+
  geom_abline (intercept = 0, slope = 1, linetype=2) +
  #geom_segment(aes(y=tpr-0.05, yend=tpr+0.05, x=fpr, xend=fpr), color="grey")+
  #geom_text(data=roc.FS.prob.label %>% filter(rep.indx == i),
  #          aes(x=x, y=ifelse(row.number%%2==1,y+0.1,y-0.1),
  #              label=prob), size=2)+
  geom_text(data=AUC.FS.label  ,
            aes(x=x,y=y,label=auc.label), size=2.3)+
  scale_x_continuous("FPR (1-Specif.)",breaks=seq(0,1,0.2)) +
  scale_y_continuous("TPR (Sensit.)", breaks=seq(0,1,0.2), position = "right", limits=c(0,1))+
  facet_grid(sampl.frac~software, switch="y")+
  theme_light()+
  theme(#strip.placement = "outside",
    strip.background.y = element_rect(fill= "dark grey"),
    strip.background.x = element_rect(fill= "dark grey"),
    strip.text.y.left = element_text(angle = 0,colour = "white"),
    panel.grid.minor.y = element_blank())


  ggsave(paste0("notes/fig/ch2_poly_",sib.label,"Pair_ROCplot.pdf"), FS.ROC.ggplot, device="pdf", height=5, width=9, units="in", dpi=500)
  ggsave(paste0("notes/fig/ch2_poly_",sib.label,"Pair_ROCplot.png"), FS.ROC.ggplot,  device="png", height=5, width=9, units="in", dpi=500)

FS.runtime.stat.tbl <-  Longo %>%
    mutate(run.time = map(output, "runtime")) %>%
    group_by(software, rep.indx, sampl.frac) %>%
    summarise (
      stat = "runtime",
      score = as.numeric(run.time[1]))



FS.AUC.stat.tbl <- prob.rank.FS.tbl.1 %>%
  mutate(sampl.frac = c(1,0.75,0.5, 0.25,0)[as.numeric(sampl.frac)]) %>%
  #filter(prob >0) %>%
  arrange(desc(prob)) %>%
  group_by(software, rep.indx, sampl.frac) %>%
  summarise (
    stat = "AUC",
    score = Calculate.FS.AUC(prob, presence.T))

FS.FDR.stat.tbl <-
  prob.rank.FS.tbl %>%
  mutate(sampl.frac = c(1,0.75,0.5, 0.25,0)[as.numeric(sampl.frac)]) %>%
  ungroup() %>%
  filter(prob> 0) %>%
  arrange(desc(prob)) %>%
  group_by(software, rep.indx, sampl.frac) %>%
  summarise (
    stat = "FDR",
    score = sum(presence.T==0)/n())

FS.FNR.stat.tbl <- prob.rank.FS.tbl %>% filter(presence.T ==1) %>%
  mutate(sampl.frac = c(1,0.75,0.5, 0.25,0)[as.numeric(sampl.frac)]) %>%
  arrange(desc(prob)) %>%
  group_by(software, rep.indx, sampl.frac) %>%
  summarise (
    stat = "FNR",
    score = sum(prob==0)/n())

decimalnumcount<-function(x){stopifnot(class(x)=="character")
  x<-gsub("(.*)(\\.)|([0]*$)","",x)
  nchar(x)
}



FS.prep.xtable <-bind_rows(FS.AUC.stat.tbl,
                           FS.FDR.stat.tbl,
                           FS.FNR.stat.tbl,
                           FS.runtime.stat.tbl) %>%
  mutate(software = factor(software, levels=c("pedFac",  "COLONY", "sequoia"))) %>%
  group_by(sampl.frac, software, stat) %>%
  summarise(
    min = min(score),
    max = max(score),
    mean = mean(score),
    sd = sd(score)) %>%
  #CI.95 = signif(1.96*sd(score)/n(), 3) %>% as.character()
  pivot_wider(names_from = stat, names_glue = "{stat}_{.value}", values_from=min:sd) %>%
  group_by(sampl.frac, software) %>%
  summarise("FDR_mean w/ sd" = sprintf("%.4f +/- %.4f", FDR_mean, FDR_sd),
            "FDR_range" = sprintf("(%.3f,%.3f)", FDR_min, FDR_max),
            "FNR_mean w/ sd" = sprintf("%.4f +/- %.4f", FNR_mean, FNR_sd),
            "FNR_range" = sprintf("(%.3f,%.3f)", FNR_min, FNR_max),
            "AUC_mean w/ sd" = sprintf("%.3f +/- %.3f", AUC_mean, AUC_sd),
            "AUC_range" = sprintf("(%.2f,%.2f)", AUC_min, AUC_max),
            "mean w/ sd" = sprintf("%.f +/- %.f", runtime_mean, runtime_sd)) %>%
  arrange(desc(sampl.frac), software) %>%
  rename(fraction = sampl.frac)

if (sib.label != "FS") {
  latex.out <- kbl(FS.prep.xtable %>% select(-"mean w/ sd"), booktabs = T, "latex") %>%
    column_spec(1, bold = T, width = "3.5em", latex_valign= "m") %>%
    row_spec(c(1:3, 7:9, 13:15)-1, extra_latex_after = "\\rowcolor{gray!6}") %>%
    collapse_rows(1, latex_hline = "none") %>%
    row_spec(0, align="c") %>%
    add_header_above(c(" ", " ", "FDR" = 2, "FNR" = 2,"ROC-AUC" = 2), align="c", bold=T)
} else {
  latex.out <- kbl(FS.prep.xtable, booktabs = T, "latex") %>%
    column_spec(1, bold = T, width = "3.5em", latex_valign= "m") %>%
    row_spec(c(1:3, 7:9, 13:15)-1, extra_latex_after = "\\rowcolor{gray!6}") %>%
    collapse_rows(1, latex_hline = "none") %>%
    row_spec(0, align="c") %>%
    add_header_above(c(" ", " ", "FDR" = 2, "FNR" = 2,"ROC-AUC" = 2, "Runtime (sec)" = 1), align="c", bold=T)
}


latex.out.1 <- gsub('\\+/-', "$\\\\pm$", latex.out, perl = T)
latex.out.2 <- gsub('(AUC\\\\_)|(FDR\\\\_)|(FNR\\\\_)', "", latex.out.1, perl = T)
write(latex.out.2, paste0("notes/table/ch2_poly_",sib.label,"pair.tex"))

return(FS.prep.xtable)

## may need to do Halfsib pair
}

fullsib.prep.xtable <- Plot.SibPair()
halfsib.prep.xtable <- Plot.SibPair("HS")

allsib.pre.xtable <- left_join(
fullsib.prep.xtable %>% select(-FDR_range, -FNR_range, -AUC_range, -`mean w/ sd`) %>%
  rename(`FS_FDR(mean,sd)`=`FDR_mean w/ sd`,
         `FS_FNR(mean,sd)`= `FNR_mean w/ sd`,
         `FS_AUC(mean,sd)` = `AUC_mean w/ sd`),

halfsib.prep.xtable %>%
  mutate( `FDR_mean w/ sd` = ifelse(software =="sequoia" & fraction == 0,"Nada", `FDR_mean w/ sd`),
          `FNR_mean w/ sd` = ifelse(software =="sequoia" & fraction == 0,"Nada", `FNR_mean w/ sd`),
          `AUC_mean w/ sd` = ifelse(software =="sequoia" & fraction == 0,"Nada", `AUC_mean w/ sd`)) %>%
  select(-FDR_range, -FNR_range, -AUC_range) %>%
  rename(`HS_FDR(mean,sd)`=`FDR_mean w/ sd`,
         `HS_FNR(mean,sd)`= `FNR_mean w/ sd`,
         `HS_AUC(mean,sd)` = `AUC_mean w/ sd`), by=c("fraction", "software"))

latex.out <- kbl(allsib.pre.xtable, booktabs = T, "latex") %>%
  column_spec(1, bold = T, width = "5em", latex_valign= "m") %>%
  row_spec(c(1:3, 7:9, 13:15)-1, extra_latex_after = "\\rowcolor{gray!6}") %>%
  collapse_rows(1, latex_hline = "none") %>%
  row_spec(0, align="c") %>%
  #add_header_above(c(" ", " ", "FDR" = 1, "FNR" = 1,"AUC" = 1, "FDR" = 1, "FNR" = 1,"AUC" = 1)) %>%
  add_header_above(c(" ", " ", "Assigned Full-Sibling Pair" = 3, "Assigned Half-Sibling Pair" = 3, "Runtime (sec)" = 1), bold= T)

#latex.out.2 <- gsub('(AUC\\\\_)|(FDR\\\\_)|(FNR\\\\_)|(HS\\\\_)|(FS\\\\_)', "", latex.out, perl = T)
latex.out.0 <- gsub('\\+/-', "$\\\\pm$", latex.out, perl = T)
latex.out.1 <- gsub('(HS\\\\_)|(FS\\\\_)', "", latex.out.0, perl = T)
latex.out.2 <- gsub('Nada', "\\\\multicolumn{1}{c}{-}", latex.out.1, perl = T)

write(latex.out.2, paste0("notes/table/ch2_poly_allpair.tex"))





### sibship cluster analysis




standarize_pairwise_tbl <- function(Tib) {
  if ("cluster" %in% colnames(Tib)) {
    Tib <- Tib %>% select(-cluster, -prob.excl) %>%
      rename(prob=prob.incl)
  }
  Tib %>%
    mutate(prob = as.numeric(prob))
}


compiled.truth.cluster.tbl <- Longo %>%
  mutate(
    FS_cluster_tbl = case_when(
      software %in% c("pedFac") ~ map(output, "fullsib.cluster.ls"),
      TRUE ~ map(output, "cluster.ls")),
    ref_truth = map(FS_cluster_tbl, "truth.cluster.assign.ls")) %>%
  select(-output, -FS_cluster_tbl) %>%
  unnest(cols=c(ref_truth))

truth.cluster.compiled.1 <- compiled.truth.cluster.tbl %>%
  mutate(software = factor(software, levels=c("pedFac",  "COLONY", "sequoia")),
         sampl.frac = factor(sampl.frac, levels=c(1,0.75,0.5, 0.25,0),
                             labels = c("sf:\n1.0","sf:\n0.75","sf:\n0.5","sf:\n0.25", "sf:\n0"))) %>%
  arrange(desc(n.sib.truth), truth.grp)%>%
  group_by(rep.indx, sampl.frac, software) %>%
  mutate(index.n = row_number())

true.sib.grp.ggplot = list()

for (i in 1:5) {

  true.sib.group.plot <- ggplot(data= truth.cluster.compiled.1 %>%  filter(rep.indx == i))+
    geom_segment(aes(x=index.n, xend=index.n, y=n.sib.truth, yend=n.sib.truth+n.excess), color ="#FF5733" ,alpha=0.9)+
    geom_segment(aes(x=index.n, xend=index.n, y=n.sib.truth, yend=n.sib.truth-n.missing), color="orange" ,alpha=0.9)+
    geom_point(aes(x=index.n, y=n.sib.truth), size=0.2)+
    facet_grid(sampl.frac~software, switch="y")+
    scale_y_continuous(
      "true sibship size",
      position = "right",
      minor_breaks = seq(0 , 12, 1),
      breaks = seq(0, 12, 4), limits = c(0,12)
    ) +
    theme_light()+
    scale_x_continuous("Index of the true full sibling groups")+
    theme(#strip.placement = "outside",
      strip.background.y = element_rect(fill= "dark grey"),
      strip.background.x = element_rect(fill= "dark grey"),
      strip.text.y.left = element_text(angle = 0,colour = "white"))


  legend.sibgrp <- ggplot()+
    #geom_text(aes(x=0.08, y=0.6, label= "Individual(s) that is/are:"),size=3.5, hjust="left", vjust="middle")+
    geom_segment(aes(x=0.05, y=0.0, xend=0.05, yend=1), color ="#FF5733" ,alpha=0.9, size=0.9)+
    geom_text(aes(x=0.08, y=0.6, label= "excess FS member(s)"),size=3.5, hjust="left", vjust="middle")+
    geom_segment(aes(x=0.55, y=0.0, xend=0.55, yend=1), color ="orange" ,alpha=0.9, size=0.9)+
    geom_text(aes(x=0.58, y=0.6, label= "missing individual(s)"),size=3.5, hjust="left", vjust="middle")+
    scale_y_continuous("", position = "right", breaks=FALSE,limits = c(0,1),minor_breaks = NULL,labels = NULL)+
    scale_x_continuous("", limits=c(0,1))+
    theme(#strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0,colour = "white"),
      panel.grid.minor.y = element_blank(),
      axis.line=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      legend.position="none",
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())

  true.sib.grp.ggplot[[i]] <- grid.arrange(true.sib.group.plot, legend.sibgrp, ncol=1,heights=c(10.5,1),
                                           padding = unit(0, "line"))

  ggsave(paste0("notes/fig/ch2_poly_FSgrp_byTruth_",i,".pdf"), true.sib.group.plot, device="pdf", height=5, width=9, units="in", dpi=500)

  ggsave(paste0("notes/fig/ch2_poly_FSgrp_byTruth_",i,".png"), true.sib.grp.ggplot[[i]],  device="png", height=5, width=9, units="in", dpi=500)
}



compiled.infer.cluster.tbl <- Longo %>%
  mutate(
    FS_cluster_tbl = case_when(
      software %in% c("pedFac") ~ map(output, "fullsib.cluster.ls"),
      TRUE ~ map(output, "cluster.ls")),
    ref_infer = map(FS_cluster_tbl, "posterior.tbl")) %>%
  select(-output, -FS_cluster_tbl) %>%
  unnest(cols=c(ref_infer))


colorBrewer2_scale <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ffff33')

infer.cluster.compiled.1 <- compiled.infer.cluster.tbl %>%
  mutate(software = factor(software, levels=c("pedFac",  "COLONY", "sequoia")),
         sampl.frac = factor(sampl.frac, levels=c(1,0.75,0.5, 0.25,0),
                             labels = c("sf:\n1.0","sf:\n0.75","sf:\n0.5","sf:\n0.25", "sf:\n0"))) %>%
  arrange(desc(prob), desc(n.sib.infer))%>%
  group_by(rep.indx, sampl.frac, software) %>%
  mutate(index.n = row_number()) %>%
  ungroup() %>%
  mutate(n.sib.truth.category = case_when(
    n.sib.truth>5~"8 - 5",
    n.sib.truth==5 ~"5",
    n.sib.truth==4 ~ "4",
    n.sib.truth==3 ~"3",
    n.sib.truth==2 ~ "2",
    n.sib.truth==1 ~ "1",
    TRUE~NA_character_))



infer.sib.grp.ggplot = list()

for (i in 1:5) {


  y.limit <- infer.cluster.compiled.1 %>%  filter(rep.indx == i) %>%
    summarise(max.y = max(prob+n.excess*0.1),
              min.y = min(prob-n.missing*0.1))

  infer.sib.group.plot <-
    ggplot(data= infer.cluster.compiled.1 %>%  filter(rep.indx == i))+
    geom_segment(aes(x=index.n, xend=index.n, y=prob, yend=prob+n.excess*0.1), color ="#FF5733" ,alpha=0.9)+
    geom_segment(aes(x=index.n, xend=index.n, y=prob, yend=prob-n.missing*0.1), color="orange" ,alpha=0.9)+
    geom_rug(aes(x=index.n, color = n.sib.truth.category))+
    scale_color_manual("True Sibship Size",values=colorBrewer2_scale,limits= c("8 - 5","4","3","2","1"))+
    guides(colour = guide_legend(nrow = 1))+
    geom_line(aes(x=index.n, y=prob))+
    facet_grid(sampl.frac~software, switch="y", scale="free_x")+
    scale_y_continuous(
      "prob",
      position = "right",
      minor_breaks = seq(floor(y.limit$min.y/0.1)*0.1 ,
                         ceiling(y.limit$max.y/0.1)*0.1 ,
                         0.1), breaks=seq(0,1,0.5)#, limits = c(-0.2,1.1)
    ) +
    theme_light()+
    scale_x_continuous("Index of inferred full sibling groups")+
    theme(#strip.placement = "outside",
      strip.background.y = element_rect(fill= "dark grey"),
      strip.background.x = element_rect(fill= "dark grey"),
      strip.text.y.left = element_text(angle = 0,colour = "white"),
      legend.position="bottom",
      #panel.grid.minor.y = element_blank(),
    )


  legend.sibgrp <- ggplot()+
    #geom_text(aes(x=0.08, y=0.6, label= "Individual(s) that is/are:"),size=3.5, hjust="left", vjust="middle")+
    geom_segment(aes(x=0.05, y=0.0, xend=0.05, yend=1), color ="#FF5733" ,alpha=0.9, size=0.9)+
    geom_text(aes(x=0.08, y=0.6, label= "excess FS member(s)"),size=3.5, hjust="left", vjust="middle")+
    geom_segment(aes(x=0.55, y=0.0, xend=0.55, yend=1), color ="orange" ,alpha=0.9, size=0.9)+
    geom_text(aes(x=0.58, y=0.6, label= "missing individual(s)"),size=3.5, hjust="left", vjust="middle")+
    scale_y_continuous("", position = "right", breaks=FALSE,limits = c(0,1),minor_breaks = NULL,labels = NULL)+
    scale_x_continuous("", limits=c(0,1))+
    theme(#strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0,colour = "white"),
      panel.grid.minor.y = element_blank(),
      axis.line=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      legend.position="none",
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())
  grid.arrange(infer.sib.group.plot, legend.sibgrp, ncol=1,heights=c(10.5,1),
               padding = unit(0, "line"))

  infer.sib.grp.ggplot[[i]] <- grid.arrange(infer.sib.group.plot, legend.sibgrp, ncol=1,heights=c(10.5,1),
                                            padding = unit(0, "line"))

  ggsave(paste0("notes/fig/ch2_poly_FSgrp_infer_",i,".pdf"), infer.sib.group.plot, device="pdf", height=5, width=9, units="in", dpi=500)

  ggsave(paste0("notes/fig/ch2_poly_FSgrp_infer_",i,".png"), infer.sib.grp.ggplot[[i]],  device="png", height=5, width=9, units="in", dpi=500)
}



# summarise cluster info in table

Calculate.SibCluster.AUC <- function(prob, categ){
  match.score <- prob[categ==1] # num.correct.calls
  mismatch.score <- prob[categ==0]
  if (length(match.score)==0) return(0)
  if (length(mismatch.score)==0) return(1)
  o <- outer(match.score, mismatch.score, "-")
  auc <- mean((o>0) + .5*(o==0))
  auc
}

FS.cluster.AUC.stat.tbl <- compiled.infer.cluster.tbl %>%
  ungroup() %>%
  arrange(desc(prob), desc(n.sib.infer)) %>%
  group_by(rep.indx, sampl.frac, software) %>%
  summarise(stat = "AUC",
            score= Calculate.SibCluster.AUC(prob, is.exact.match))

FS.cluster.FDR.stat.tbl <- compiled.infer.cluster.tbl %>%
  ungroup() %>%
  arrange(desc(prob)) %>%
  group_by(rep.indx, sampl.frac, software) %>%
  summarise(stat = "FDR",
            score=1-sum(is.exact.match == TRUE) /n())

FS.cluster.FNR.stat.tbl <- compiled.truth.cluster.tbl%>%
  ungroup() %>%
  arrange(desc(prob)) %>%
  group_by(rep.indx, sampl.frac, software) %>%
  summarise(stat = "FNR",
            score=1-sum(is.exact.match == TRUE) /n())


FS.cluster.prep.xtable <- bind_rows(FS.cluster.AUC.stat.tbl,
                                    FS.cluster.FDR.stat.tbl,
                                    FS.cluster.FNR.stat.tbl) %>%
  mutate(software = factor(software, levels=c("pedFac",  "COLONY", "sequoia"))) %>%
  group_by(sampl.frac, software, stat) %>%
  summarise(
    min = min(score),
    max = max(score),
    mean = mean(score),
    sd = sd(score)) %>%
  pivot_wider(names_from = stat, names_glue = "{stat}_{.value}", values_from=min:sd) %>%
  group_by(sampl.frac, software) %>%
  summarise("FDR_mean w/ sd" = sprintf("%.4f +/- %.4f", FDR_mean, FDR_sd),
            "FDR_range" = sprintf("(%.3f,%.3f)", FDR_min, FDR_max),
            "FNR_mean w/ sd" = sprintf("%.4f +/- %.4f", FNR_mean, FNR_sd),
            "FNR_range" = sprintf("(%.3f,%.3f)", FNR_min, FNR_max),
            "AUC_mean w/ sd" = sprintf("%.3f +/- %.3f", AUC_mean, AUC_sd),
            "AUC_range" = sprintf("(%.2f,%.2f)", AUC_min, AUC_max)) %>%
  arrange(desc(sampl.frac), software) %>%
  rename(fraction = sampl.frac)

cluster.latex.out <- kbl(FS.cluster.prep.xtable, booktabs = T, "latex") %>%
  column_spec(1, bold = T, width = "3.5em", latex_valign= "m") %>%
  row_spec(c(1:3, 7:9, 13:15)-1, extra_latex_after = "\\rowcolor{gray!6}") %>%
  collapse_rows(1, latex_hline = "none") %>%
  row_spec(0, align="c") %>%
  add_header_above(c(" ", " ", "FDR" = 2, "FNR" = 2,"ROC-AUC" = 2), align="c", bold=T)

cluster.latex.out.1 <- gsub('\\+/-', "$\\\\pm$", cluster.latex.out, perl = T)
cluster.latex.out.2 <- gsub('(AUC\\\\_)|(FDR\\\\_)|(FNR\\\\_)', "", cluster.latex.out.1, perl = T)
write(cluster.latex.out.2, "notes/table/ch2_poly_FSgrp.tex")




##### LL chart


poly.ll.tbl <- expand_grid(
  samp_frac_idx = 1:5,
  rep.indx = 1:5
) %>%
  mutate(
    sampl.frac = c(1, 0.75, 0.5, 0.25, 0)[samp_frac_idx],
    ll = map2(.x = sampl.frac, .y = rep.indx, .f = function(x, y) {
      read.table(paste0("/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/halfSib/sf_",x,"/",y,"/pedfac/ll.txt")) %>%
        tibble() %>%
        mutate(rank = row_number())
      })) %>%
  unnest(cols=c(ll)) %>%
  group_by(samp_frac_idx, rep.indx, V1) %>%
  mutate(first.tick = ifelse(rank == min(rank), 1, 0)) %>%
  ungroup()


ggplot(poly.ll.tbl) +
  geom_line(aes(x=rank, y=V3, color = factor(sampl.frac), linetype = factor(rep.indx))) +
  geom_segment(data = poly.ll.tbl %>% filter(first.tick ==1), aes(y=V3-(400-(V1+1)*4), yend=V3+(400-(V1+1)*4), x=rank, xend=rank, alpha = 1/(V1+1)), color ="black") +
  geom_text(data = poly.ll.tbl %>%
              filter(first.tick ==1, V1==0) %>%
              group_by(sampl.frac) %>%
              summarise(min.post = min(V3)),
            aes(x=1.3, y=min.post-1000, label=paste0("sf: ",sampl.frac), color = factor(sampl.frac)), hjust="left")+
  scale_x_log10()+
  xlab("iterations")+
  ylab("log posterior")+
  scale_linetype_discrete(guide=FALSE)+
  scale_color_discrete(guide=FALSE)+
  scale_alpha_continuous(guide=FALSE)+
  theme_light()+
  annotation_logticks()+
  theme(
    legend.position="bottom",
  )


###

poly.ll.plot <- ggplot(poly.ll.tbl %>% filter(V1>0)) +
  geom_line(aes(x=rank, y=V3, color = factor(sampl.frac), linetype = factor(rep.indx))) +
  geom_segment(data = poly.ll.tbl %>% filter(first.tick ==1, V1>0), aes(y=V3-(400-(V1+1)*4), yend=V3+(400-(V1+1)*4), x=rank, xend=rank, alpha = 1/(V1+1)), color ="black") +
  geom_text(data = poly.ll.tbl %>%
              filter(first.tick ==1, V1==1) %>%
              group_by(sampl.frac) %>%
              summarise(min.post = min(V3),
                        min.iter = min(rank)),
            aes(x=min.iter, y=min.post-800, label=paste0("sf: ",sampl.frac), color = factor(sampl.frac)), hjust="left")+
  geom_text(data = poly.ll.tbl %>%
              filter(first.tick ==1, V1 %in% c(1,2,3,5,10,25,50), sampl.frac == 0) %>%
              group_by(sampl.frac, V1) %>%
              summarise(max.post = max(V3),
                        min.iter = min(rank)),
            aes(x=min.iter, y=max.post+1000, label=V1))+
  scale_x_log10()+
  xlab("iterations")+
  ylab("log posterior")+
  scale_linetype_discrete(guide=FALSE)+
  scale_color_discrete(guide=FALSE)+
  scale_alpha_continuous(guide=FALSE)+
  theme_light()+
  annotation_logticks()+
  theme(
    legend.position="bottom",
  )

ggsave(paste0("notes/fig/ch2_poly_ll.pdf"), poly.ll.plot, device="pdf", height=5, width=9, units="in", dpi=500)

ggsave(paste0("notes/fig/ch2_poly_ll.png"), poly.ll.plot,  device="png", height=5, width=9, units="in", dpi=500)

