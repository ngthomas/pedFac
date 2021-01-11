
#python src/simPed.py -s 200 -ni 200 -ns 100 -sf 1 --gen 3 -o data/3gen/cyclic/ni_200_ns100_sf1.0_s200

RetrieveGrandparents <- function(df, ref=FALSE) {
  ped.grandpa <- df %>% dplyr::rename(parent.0=kid,
                               grandpa.pa = parent.0,
                               grandmom.pa = parent.1)
  ped.grandmom <- df %>% dplyr::rename(parent.1=kid,
                                grandpa.mom = parent.0,
                                grandmom.mom = parent.1)

  join.label.1 <- "parent.0"
  join.label.2 <- "parent.1"
  if (!ref) {
    join.label.1[2] <- "iter"
    join.label.2[2] <- "iter"
  }

  ped.join.pa<- dplyr::inner_join(df, ped.grandpa, by=join.label.1)
  dplyr::inner_join(ped.join.pa, ped.grandmom, by=join.label.2)
}

#run.path <- "/Users/thomasn/repo/pedigree-factor-graphs/data/3gen/cyclic/ni_200_ns100_sf1.0_s200"
#max.iter <- 20

#impose.acyclic refers to whether simulated ``truth'', not the sampled, pedigree is restricted to be acyclic

ExtractGrandParents <- function(run.path = "",
                                out.folder.name = "out",
                                max.iter = 20,
                                impose.acyclic=FALSE) {

  ped.orig.txt <- read.table(paste0(run.path,"/",out.folder.name,"/ped.txt")) %>% dplyr::tibble()
  colnames(ped.orig.txt) <- c("iter","kid","parent.0","parent.1")
  prior.txt <- read.delim(paste0(run.path,"/prior.txt"), sep="\n")
  maxID.indx <- sapply(1:nrow(prior.txt), function(l) grepl("maxID", prior.txt[l,])) %>% which

  max.id <- prior.txt[maxID.indx,] %>%
    strsplit(., " ") %>%
    unlist() %>% .[2] %>%
    as.integer()

  ## the extra slice action is to remove any extra entry under swapping act
  ped.txt <- ped.orig.txt %>%
    dplyr::mutate(rown = row_number()) %>%
    dplyr::group_by(iter, kid)%>% dplyr::arrange(desc(rown)) %>% dplyr::slice(1) %>%
    dplyr::select(-rown) %>%
    dplyr::ungroup()

  ped.truth <- read.table(paste0(run.path,"/marriage.all")) %>% dplyr::tibble()
  colnames(ped.truth) <- c("kid","parent.0","parent.1",
                           "is.obs.kid",
                           "is.obs.parent.0",
                           "is.obs.parent.1",
                           "is.censor")

  sample.grandparent <- RetrieveGrandparents(ped.txt)
  truth.grandparent <- RetrieveGrandparents(ped.truth, TRUE)

  sample.freq.summary <- sample.grandparent %>%
    dplyr::transmute(kid = kid,
              grandpa.pa = ifelse(grandpa.pa>=max.id,-1,grandpa.pa),
              grandmom.pa = ifelse(grandmom.pa>=max.id,-1,grandmom.pa),
              grandpa.mom = ifelse(grandpa.mom>=max.id,-1,grandpa.mom),
              grandmom.mom = ifelse(grandmom.mom>=max.id,-1,grandmom.mom),
              parent.0 = ifelse(parent.0>=max.id,-1,parent.0),
              parent.1 = ifelse(parent.1>=max.id,-1,parent.1),
              ) %>%
    dplyr::group_by(kid, grandpa.pa, grandmom.pa, grandpa.mom, grandmom.mom, parent.0, parent.1) %>% ##keep parent0 part of the eqn
    dplyr::summarise(n=n(), posterior = n()/max.iter*100) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(n)) %>%
    dplyr::group_by(kid) %>%
    dplyr::top_n(1) %>%
    dplyr::sample_n(1)
  ## taking the best posterior at this point

  ## grabbing the paternal side
  s.1 <- sample.grandparent %>%
    dplyr::transmute(kid = kid,
              grandpa.pa = ifelse(grandpa.pa>=max.id,-1,grandpa.pa),
              grandmom.pa = ifelse(grandmom.pa>=max.id,-1,grandmom.pa),
              grandpa.mom = ifelse(grandpa.mom>=max.id,-1,grandpa.mom),
              grandmom.mom = ifelse(grandmom.mom>=max.id,-1,grandmom.mom)) %>%
    dplyr::select(kid, grandpa.pa, grandmom.pa) %>%
    dplyr::rename(grandpa=grandpa.pa, grandmom = grandmom.pa)

  ## grabbing the maternal side
  s.2 <- sample.grandparent %>%
    dplyr::transmute(kid = kid,
              grandpa.pa = ifelse(grandpa.pa>=max.id,-1,grandpa.pa),
              grandmom.pa = ifelse(grandmom.pa>=max.id,-1,grandmom.pa),
              grandpa.mom = ifelse(grandpa.mom>=max.id,-1,grandpa.mom),
              grandmom.mom = ifelse(grandmom.mom>=max.id,-1,grandmom.mom)) %>%
    dplyr::select(kid, grandpa.mom, grandmom.mom) %>%
    dplyr::rename(grandpa=grandpa.mom, grandmom = grandmom.mom)

  #erase the paternal and maternal labels
  s.3<- rbind(s.1, s.2) %>%
    dplyr::group_by(kid, grandpa, grandmom) %>%
    dplyr::summarise(n=n(), posterior = n()/max.iter*100) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(n)) %>%
    dplyr::group_by(kid) #%>%
  #top_n(2)

  if (impose.acyclic) {
    truth.grandparent <- truth.grandparent %>%
      dplyr::filter(is.obs.kid.x ==1, is.censor.x == 0)
  }
  truth.summary <- truth.grandparent %>%
    dplyr::filter(is.obs.kid.x ==1) %>%
    dplyr::transmute(kid = kid,
              grandpa.pa.T = ifelse(is.obs.parent.0.y == 1, grandpa.pa, -1),
              grandmom.pa.T = ifelse(is.obs.parent.1.y == 1, grandmom.pa, -1),
              grandpa.mom.T = ifelse(is.obs.parent.0 == 1, grandpa.mom, -1),
              grandmom.mom.T = ifelse(is.obs.parent.1 == 1, grandmom.mom, -1),
              pa = ifelse(is.obs.parent.0.x == 1, parent.0, -1*parent.0),
              ma = ifelse(is.obs.parent.1.x == 1, parent.1, -1*parent.1)
              )

  ped.join.truth <- dplyr::left_join(sample.freq.summary, truth.summary, by="kid")


  ped.stat <- ped.join.truth %>%
    dplyr::mutate(num.common = grandpa.pa.T %in% c(grandpa.pa, grandpa.mom,grandmom.pa,grandmom.mom)+
             grandpa.mom.T %in% c(grandpa.pa, grandpa.mom,grandmom.pa,grandmom.mom)+
             grandmom.pa.T %in% c(grandpa.pa, grandpa.mom,grandmom.pa,grandmom.mom)+
             grandmom.mom.T %in% c(grandpa.pa, grandpa.mom,grandmom.pa,grandmom.mom))

  index <- nrow(ped.stat)
  ped.stat %>% cbind(index=1:index)
  # ggplot()+
  #   geom_point(data=ped.stat.plot, aes(x=index, y=posterior))

  # ggplot()+
  #   geom_point(data=ped.stat.plot, aes(x=index, y=posterior))+
  #   geom_segment(aes(x = index,
  #                    y = posterior*1.2,
  #                    xend = index,
  #                    yend = posterior*0.8,
  #                    color=factor(num.common)),
  #                data=ped.stat.plot%>%filter(num.common>0))
}







