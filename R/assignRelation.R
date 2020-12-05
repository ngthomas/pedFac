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


retrieveFullSib <- function(df, max.id, id.ls) {
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

  if(is.null(nrow(sib.pairs))) return()

  lapply(sib.pairs,
         function(x) combn(as.numeric(unlist(strsplit(x,",",perl=T))),
                           m=2,
                           simplify = F)) %>%
    dplyr::bind_cols() %>%
    t %>%
    dplyr::tbl_df() %>%
    dplyr::mutate(kid.1 = ifelse(V1<V2, V1,V2),
           kid.2 = ifelse(V1<V2, V2,V1)) %>%
    dplyr::group_by(kid.1, kid.2) %>%
    dplyr::summarise(prob=n()/n.iter) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(prob)) %>%
    dplyr::mutate(kid.1 = as.character(subbingID(kid.1, max.id, id.ls)),
                  kid.2 = as.character(subbingID(kid.2, max.id, id.ls)))
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

