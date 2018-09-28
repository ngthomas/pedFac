retrieveGrandparents <- function(df, max.id, id.ls) {
  ped.grandpa <- df %>% dplyr::rename(pa=kid,
                               grandpa.pa = pa,
                               grandma.pa = ma)
  ped.grandma <- df %>% dplyr::rename(ma=kid,
                                grandpa.ma = pa,
                                grandma.ma = ma)

  n.iter <- max(df$iter)+1
  ped.join.pa<- dplyr::inner_join(df, ped.grandpa, by=c("iter","pa"))
  dplyr::inner_join(ped.join.pa, ped.grandma, by=c("iter","ma")) %>%
    dplyr::select(-pa, -ma) %>%
    dplyr::filter(kid<=max.id) %>%
    dplyr::group_by(iter, kid) %>%
    dplyr::mutate(grandpa.pa = subbingID(grandpa.pa, max.id, id.ls),
                  grandma.pa = subbingID(grandma.pa, max.id, id.ls),
                  grandpa.ma = subbingID(grandpa.ma, max.id, id.ls),
                  grandma.ma = subbingID(grandma.ma, max.id, id.ls)) %>%
    dplyr::group_by(kid, grandpa.pa, grandma.pa, grandpa.ma, grandma.ma) %>%
    dplyr::summarise(prob = n()/n.iter) %>%
    dplyr::arrange(desc(prob))
}

subbingID <- function(numeric.id, max.id, sub.ls) ifelse(numeric.id>max.id, -1, sub.ls[numeric.id])

retrieveFullSib <- function(df, max.id) {
  sib.pairs <- df %>%
    dplyr::filter(kid<=max.id) %>%
    dplyr::group_by(iter, pa, ma) %>%
    dplyr::summarise(full.sib = paste0(kid, collapse = ","),
                     n.sibs = n()) %>%
    dplyr::filter(n.sibs > 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(full.sib) %>%
    unlist()

  n.iter <- max(df$iter)+1

  lapply(sib.pairs,
         function(x) combn(as.numeric(unlist(strsplit(x,",",perl=T))),
                           m=2,
                           simplify = F)) %>%
    dplyr::bind_cols() %>%
    t %>%
    dplyr::tbl_df() %>%
    mutate(kid.1 = ifelse(V1<V2, V1,V2),
           kid.2 = ifelse(V1<V2, V2,V1)) %>%
    dplyr::group_by(kid.1, kid.2)%>%
    dplyr::summarise(prob=n()/n.iter)
}

retrieveHalfSib <- function(df, max.id) {
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
    dplyr::filter(prob>0)
}
