subbingID <- function(numeric.id, max.id, sub.ls) ifelse(numeric.id>max.id, -1, sub.ls[numeric.id])

retrieveGrandparents <- function(df, max.id, id.ls) {
  ped.grandpa <- df %>% dplyr::rename(pa=kid,
                               grandpa.pa = pa,
                               grandma.pa = ma)
  ped.grandma <- df %>% dplyr::rename(ma=kid,
                                grandpa.ma = pa,
                                grandma.ma = ma)

  n.iter <- max(df$iter)+1
  ped.join.pa<- dplyr::inner_join(df, ped.grandpa, by=c("iter","pa"))
  if(nrow(ped.join.pa)==0) return()

  dplyr::inner_join(ped.join.pa, ped.grandma, by=c("iter","ma")) %>%
    dplyr::select(-pa, -ma) %>%
    dplyr::filter(kid<=max.id) %>%
    dplyr::group_by(iter, kid) %>%
    dplyr::mutate(grandpa.pa = as.character(subbingID(grandpa.pa, max.id, id.ls)),
                  grandma.pa = as.character(subbingID(grandma.pa, max.id, id.ls)),
                  grandpa.ma = as.character(subbingID(grandpa.ma, max.id, id.ls)),
                  grandma.ma = as.character(subbingID(grandma.ma, max.id, id.ls))) %>%
    dplyr::group_by(kid, grandpa.pa, grandma.pa, grandpa.ma, grandma.ma) %>%
    dplyr::summarise(prob = n()/n.iter) %>%
    dplyr::arrange(desc(prob))  %>%
    dplyr::ungroup() %>%
    dplyr::mutate(kid = subbingID(kid, max.id, id.ls))
}

retrieveParent <- function(df, max.id, id.ls) {

  n.iter <- max(df$iter)+1

  df %>%
    dplyr::filter(kid<=max.id) %>%
    dplyr::group_by(iter, kid) %>%
    dplyr::mutate(pa.id = as.character(subbingID(pa, max.id, id.ls)),
                  ma.id = as.character(subbingID(ma, max.id, id.ls))) %>%
    dplyr::group_by(kid, pa.id, ma.id) %>%
    dplyr::summarise(prob = n()/n.iter) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(prob)) %>%
     dplyr::mutate(kid.id = subbingID(kid, max.id, id.ls)) %>%
    ungroup() %>%
    dplyr::select(kid.id, pa.id, ma.id, prob)
}


RetrieveFullSib <- function(out.ped.tbl, max.id, id.ls) {
  sib.pairs <- out.ped.tbl %>%
    dplyr::filter(kid<=max.id) %>%
    dplyr::group_by(iter, pa, ma) %>%
    dplyr::summarise(full.sib = paste0(kid, collapse = ","),
                     n.sibs = n()) %>%
    dplyr::filter(n.sibs > 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(full.sib) %>%
    unlist()

  n.iter <- max(out.ped.tbl$iter)+1

  if(length(sib.pairs) == 0) return()

  sib.pair.prob.tbl <- lapply(sib.pairs,
         function(x) combn(as.numeric(unlist(strsplit(x,",",perl=T))),
                           m=2,
                           simplify = T)) %>%
    dplyr::bind_cols() %>%
    t %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(kid.1 = ifelse(V1<V2, V1,V2),
           kid.2 = ifelse(V1<V2, V2,V1)) %>%
    dplyr::group_by(kid.1, kid.2) %>%
    dplyr::summarise(prob=n()/n.iter) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(prob)) %>%
    dplyr::mutate(kid.1 = as.character(subbingID(kid.1, max.id, id.ls)),
                  kid.2 = as.character(subbingID(kid.2, max.id, id.ls)))

  return(sib.pair.prob.tbl)
}

RetrieveFullSibTruth <- function(mating.factor.tbl) {
  full.sib.grp.truth.tbl <- mating.factor.tbl %>%
    dplyr::group_by(pa, ma) %>%
    dplyr::summarise(full.sib = paste0(kid, collapse = ","),
                     n.sibs = n()) %>%
    dplyr::filter(n.sibs > 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(full.sib) %>%
    unlist()

  full.sib.pair.truth.tbl <- lapply(full.sib.grp.truth.tbl,
                                    function(x) combn(as.numeric(unlist(strsplit(x,",",perl=T))),
                                                      m=2,
                                                      simplify = T)) %>%
    dplyr::bind_cols() %>%
    t %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(kid.1 = as.character(ifelse(V1<V2, V1,V2)),
                  kid.2 = as.character(ifelse(V1<V2, V2,V1)),
                  presence.T =1) %>%
    select(kid.1, kid.2, presence.T)
  return(full.sib.pair.truth.tbl)
}

RetrieveFullSibGrp <- function(out.ped.tbl, truth.ped.tbl, max.id, id.ls) {
  sib.grp <- out.ped.tbl %>%
    dplyr::filter(kid<=max.id) %>%
    dplyr::group_by(iter, pa, ma) %>%
    dplyr::summarise(full.sib = paste0(kid, collapse = ","),
                     n.sibs = n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(full.sib, n.sibs) %>%
    dplyr::summarise(n.obs = n())

  n.iter <- max(out.ped.tbl$iter)+1

  sib.cluster.prob.tbl <- sib.grp %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(n.obs), desc(n.sibs)) %>%
    mutate(posterior = n.obs/n.iter,
           cluster = row_number()) %>%
    separate_rows(full.sib)

  GreedySearch <- function(sib.process.tbl) {
    min.cluster <- min(sib.process.tbl$cluster)

    sel.cluster <- sib.process.tbl %>% filter(cluster == min.cluster)

    cluster.rm.ls <- sib.process.tbl %>% filter(full.sib %in% sel.cluster$full.sib) %>%
      select(cluster) %>% unique %>% unlist()

    pass.cluster <- sib.process.tbl %>% filter(!cluster %in% cluster.rm.ls)
    if (dim(pass.cluster)[1]==0) {
      return(sel.cluster)
    } else {
    return(bind_rows(sel.cluster, GreedySearch(pass.cluster)))
    }
  }

  greedy.cluster.tbl <- GreedySearch(sib.cluster.prob.tbl)

  greedy.cluster.ls <- greedy.cluster.tbl %>%
    group_by(cluster, posterior, n.sibs) %>%
    summarise(full.sib = paste0(full.sib, sep="", collapse = ","))

  truth.sib.grp <- truth.ped.tbl %>%
    dplyr::filter(kid<=max.id) %>%
    dplyr::group_by(pa, ma) %>%
    dplyr::summarise(full.sib.truth = paste0(kid, collapse = ",")) %>%
    ungroup() %>%
    select(full.sib.truth)

  CountIntersect <- function(infer.cluster, truth.cluster) {
    intersect(unlist(strsplit(infer.cluster, ",")), unlist(strsplit(truth.cluster, ","))) %>%
      length()
  }

  CountNSib <- function(cluster) {
    unlist(strsplit(cluster, ",")) %>% length
  }

  n.matches.infer.truth <- expand_grid(infer.grp = greedy.cluster.ls$full.sib, truth.sib.grp) %>%
    group_by(infer.grp, full.sib.truth) %>%
    mutate(n.match = CountIntersect(infer.grp, full.sib.truth),
           n.sib.infer = CountNSib(infer.grp),
           n.sib.truth = CountNSib(full.sib.truth),
           n.diff = abs(n.sib.infer - n.sib.truth))

  out.ls <- list()

  out.ls$posterior.tbl <- left_join(greedy.cluster.ls %>% rename(infer.grp="full.sib"),
                                    n.matches.infer.truth %>%
    arrange(desc(n.match), n.diff) %>%
    group_by(infer.grp) %>%
    summarise(max.match = n.match[1],
              n.diff = n.diff[1],
              n.sib.infer = n.sib.infer[1],
              n.sib.truth = n.sib.truth[1])) %>%
    group_by(cluster) %>%
    mutate(is.exact.match = (max.match == n.sibs && n.diff == 0),
           # situation in which inferr cluster missing the number of elements
           # situation in which inferr cluster has additional elements
           n.excess = ifelse(n.sib.infer > n.sib.truth, n.sib.infer - max.match, 0),
           # number of elements missing in the inferred cluster
           n.missing = ifelse(n.sib.truth > n.sib.infer,n.sib.truth - max.match, 0)
           )

  max.cluster.indx <- greedy.cluster.ls$cluster %>% max
  out.ls$indiv.cluster.ls <- full_join(greedy.cluster.ls %>% ungroup() %>%
    separate_rows(full.sib) %>%
    rename(indiv.indx = "full.sib", cluster.infer = "cluster") %>%
    select(cluster.infer, indiv.indx),
    truth.sib.grp %>% ungroup() %>%
    mutate(cluster = row_number()) %>%
    separate_rows(full.sib.truth) %>%
    rename(indiv.indx = "full.sib.truth", cluster.truth = "cluster") %>%
    select(cluster.truth, indiv.indx),
  by = "indiv.indx") %>%
    mutate(cluster.infer = ifelse(is.na(cluster.infer), row_number()+max.cluster.indx, cluster.infer))

  return(out.ls)

}

# sequoia verion
RetrieveFullSibGrpSeq <- function(out.ped.tbl, truth.ped.tbl, max.id, id.ls){

  cluster.infer.ls <- out.ped.tbl %>%
    dplyr::filter(kid<=max.id) %>%
    dplyr::group_by(pa, ma) %>%
    dplyr::summarise(full.sib = paste0(kid, collapse = ","),
                     n.sibs = n(),
                     LLR = LLRpair[1]) %>%
    dplyr::ungroup() %>%
    dplyr::select(-pa, -ma) %>%
    dplyr::arrange(desc(LLR), desc(n.sibs)) %>%
    mutate(cluster = row_number())

  truth.sib.grp <- truth.ped.tbl %>%
    dplyr::filter(kid<=max.id) %>%
    dplyr::group_by(pa, ma) %>%
    dplyr::summarise(full.sib.truth = paste0(kid, collapse = ",")) %>%
    ungroup() %>%
    select(full.sib.truth)

  CountIntersect <- function(infer.cluster, truth.cluster) {
    intersect(unlist(strsplit(infer.cluster, ",")), unlist(strsplit(truth.cluster, ","))) %>%
      length()
  }

  CountNSib <- function(cluster) {
    unlist(strsplit(cluster, ",")) %>% length
  }

  n.matches.infer.truth <- expand_grid(infer.grp = cluster.infer.ls$full.sib, truth.sib.grp) %>%
    group_by(infer.grp, full.sib.truth) %>%
    mutate(n.match = CountIntersect(infer.grp, full.sib.truth),
           n.sib.infer = CountNSib(infer.grp),
           n.sib.truth = CountNSib(full.sib.truth),
           n.diff = abs(n.sib.infer - n.sib.truth))

  out.ls <- list()

  out.ls$posterior.tbl <- left_join(cluster.infer.ls %>% rename(infer.grp="full.sib"),
                                    n.matches.infer.truth %>%
                                      arrange(desc(n.match), n.diff) %>%
                                      group_by(infer.grp) %>%
                                      summarise(max.match = n.match[1],
                                                n.diff = n.diff[1],
                                                n.sib.infer = n.sib.infer[1],
                                                n.sib.truth = n.sib.truth[1])) %>%
    group_by(cluster) %>%
    mutate(is.exact.match = (max.match == n.sibs && n.diff == 0),
           # situation in which inferr cluster missing the number of elements
           # situation in which inferr cluster has additional elements
           n.excess = ifelse(n.sib.infer > n.sib.truth, n.sib.infer - max.match, 0),
           # number of elements missing in the inferred cluster
           n.missing = ifelse(n.sib.truth > n.sib.infer,n.sib.truth - max.match, 0)
    )

  max.cluster.indx <- cluster.infer.ls$cluster %>% max
  out.ls$indiv.cluster.ls <- full_join(cluster.infer.ls %>% ungroup() %>%
                                         separate_rows(full.sib) %>%
                                         rename(indiv.indx = "full.sib", cluster.infer = "cluster") %>%
                                         select(cluster.infer, indiv.indx),
                                       truth.sib.grp %>% ungroup() %>%
                                         mutate(cluster = row_number()) %>%
                                         separate_rows(full.sib.truth) %>%
                                         rename(indiv.indx = "full.sib.truth", cluster.truth = "cluster") %>%
                                         select(cluster.truth, indiv.indx),
                                       by = "indiv.indx") %>%
    mutate(cluster.infer = ifelse(is.na(cluster.infer), row_number()+max.cluster.indx, cluster.infer))

  return(out.ls)

}


retrieveHalfSib <- function(df, max.id, id.ls) {
  sib.pa.pairs <- df %>%
    dplyr::filter(kid<=max.id) %>%
    dplyr::group_by(iter, pa) %>%
    dplyr::summarise(sib = paste0(kid, collapse = ","),
                     n.sibs = n()) %>%
    dplyr::filter(n.sibs > 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(sib) %>%
    unlist()

  sib.ma.pairs <- df %>%
    dplyr::filter(kid<=max.id) %>%
    dplyr::group_by(iter, ma) %>%
    dplyr::summarise(sib = paste0(kid, collapse = ","),
                     n.sibs = n()) %>%
    dplyr::filter(n.sibs > 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(sib) %>%
    unlist()

  if(is.null(nrow(sib.ma.pairs)) && is.null(ncol(sib.pa.pairs))) return()

  n.iter <- max(df$iter)+1

  all.sibs <- lapply(c(sib.pa.pairs, sib.ma.pairs),
         function(x) combn(as.numeric(unlist(strsplit(x,",",perl=T))),
                           m=2,
                           simplify = F)) %>%
    dplyr::bind_cols() %>%
    t %>%
    dplyr::tbl_df() %>%
    dplyr::mutate(kid.1 = ifelse(V1<V2, V1,V2),
           kid.2 = ifelse(V1<V2, V2,V1)) %>%
    dplyr::group_by(kid.1, kid.2)%>%
    dplyr::summarise(prob=n()/n.iter)

  dplyr::full_join(all.sibs, retrieveFullSib(df, max.id), by=c("kid.1", "kid.2")) %>%
    dplyr::group_by(kid.1, kid.2) %>%
    dplyr::summarise(prob = ifelse(is.na(prob.y),prob.x,prob.x-(2*prob.y))) %>%
    dplyr::filter(prob>0) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(prob)) %>%
    dplyr::mutate(kid.1 = as.character(subbingID(kid.1, max.id, id.ls)),
                  kid.2 = as.character(subbingID(kid.2, max.id, id.ls)))
}

Plot.SibCluster.plot <- function(cluster.tbl) {

  max.sib.grp <- max(cluster.tbl$n.sibs)
  cluster.tbl <- cluster.tbl %>%
    ungroup() %>%
    arrange(cluster) %>%
    mutate(index = row_number())
  ggplot(data = cluster.tbl) +
    geom_line( aes(x=index, y = posterior)) +
    geom_segment(aes(y=posterior-(n.sibs/max.sib.grp), yend=posterior, x=index, xend=index), color="grey")+
    geom_segment(aes(y=posterior+(n.missing/max.sib.grp), yend=posterior, x=index, xend=index), color="red")+
    geom_segment(aes(y=posterior-(n.excess/max.sib.grp), yend=posterior, x=index, xend=index), color="orange")+
    theme_bw()+
    scale_y_continuous("Posterior",breaks=seq(0,1,0.2))+
    scale_x_continuous("Sibling cluster index")
}

Plot.SibCluster.plot.LLR <- function(cluster.tbl) {

  max.sib.grp <- max(cluster.tbl$n.sibs)
  cluster.tbl <- cluster.tbl %>%
    ungroup() %>%
    arrange(cluster) %>%
    mutate(index = row_number())
  ggplot(data = cluster.tbl) +
    geom_line( aes(x=index, y = LLR)) +
    geom_segment(aes(y=LLR-(n.sibs/max.sib.grp), yend=LLR, x=index, xend=index), color="grey")+
    geom_segment(aes(y=LLR+(n.missing/max.sib.grp), yend=LLR, x=index, xend=index), color="red")+
    geom_segment(aes(y=LLR-(n.excess/max.sib.grp), yend=LLR, x=index, xend=index), color="orange")+
    theme_bw()+
    scale_y_continuous("LLR")+
    scale_x_continuous("Sibling cluster index")
}


Plot.ROC.sibPair.plot <- function(sib.summary, label.grid.size = 0.1) {

  roc.line.df <- sib.summary %>%
    arrange(desc(posterior)) %>%
    group_by(rep) %>%
    mutate(n.true.case = cumsum(is.true == 1),
           n.false.case = cumsum(is.true == 0),
           tot.true.case = max(n.true.case),
           tot.false.case = max(n.false.case)) %>%
    group_by(rep, posterior) %>%
    summarise(tpr = ifelse(tot.true.case[1]==0,0,
                           max(n.true.case)/tot.true.case[1]),
              fpr = ifelse(tot.false.case[1]==0,0,
                           max(n.false.case)/tot.false.case[1]),
              n.true.case = tot.true.case[1],
              n.false.case = tot.false.case[1])


  Calculate.AUC <- function(posterior, labels){
    match.score <- posterior[labels==1] # num.correct.calls
    mismatch.score <- posterior[labels!=1]
    if (length(match.score)==0) return(0)
    if (length(mismatch.score)==0) return(1)
    o <- outer(match.score, mismatch.score, "-")
    auc <- mean((o>0) + .5*(o==0))
  }

  AUC.df <- sib.summary %>%
    group_by(rep) %>%
    summarise(x=0.8,
              y=0.3,
              auc.score = Calculate.AUC(posterior, is.true),
              auc.label = paste0("AUC: ",round(auc.score,2)))

  #instead of label every minimal threshold,
  #it is to space out the labeling
  xmin.grid <- seq(0,1-label.grid.size,label.grid.size)
  xmax.grid <- seq(label.grid.size,1,label.grid.size)
  len.grid <- length(xmin.grid)

  roc.label.df <- lapply(1:len.grid,
                         function(i) roc.line.df %>%
                           group_by(rep) %>%
                           filter(fpr>=xmin.grid[i] & fpr<=xmax.grid[i]) %>%
                           summarise(x=fpr[1],
                                     y=tpr[1],
                                     posterior=round(posterior[1],2),
                                     row.number = i%%2)
  ) %>%
    bind_rows() %>%
    filter(!is.na(x))

  roc.pt.fixed.1 <- roc.line.df %>%
    group_by(rep) %>%
    summarise(posterior = 0,
              tpr = 1,
              fpr = 1.0001,
              n.true.case=max(n.true.case),
              n.false.case=max(n.false.case))

  roc.pt.fixed.2 <- roc.line.df %>%
    group_by(rep) %>%
    summarise(posterior = 1,
              tpr = -0,
              fpr = -0.00001,
              n.true.case=max(n.true.case[posterior==1]),
              n.false.case=max(n.false.case[posterior==1]))


  roc.df <- bind_rows(roc.pt.fixed.2,
                      roc.line.df,
                      roc.pt.fixed.1) %>%
    group_by(rep) %>%
    arrange(posterior)


  ggplot(roc.df, aes(x = fpr, y = tpr)) +
    geom_line(color="red") +
    geom_abline (intercept = 0, slope = 1, linetype=2) +
    geom_segment(aes(y=tpr-0.05, yend=tpr+0.05, x=fpr, xend=fpr), color="grey")+
    geom_text(data=roc.label.df,
              aes(x=x, y=ifelse(row.number%%2==1,y+0.1,y-0.1),
                  label=round(posterior,2)), size=2)+
    geom_text(data=AUC.df,
              aes(x=x,y=y,label=auc.label))+
    theme_bw() +
    scale_x_continuous("FPR (1-Specif.)",breaks=seq(0,1,0.2)) +
    scale_y_continuous("TPR (Sensit.)", breaks=seq(0,1,0.2))+
    facet_grid(rep~.,labeller = label_both)

}


