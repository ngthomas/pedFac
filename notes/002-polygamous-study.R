library(sequoia)
library(tictoc)
library(tidyverse)

#simulating pedigree for this case study

set.seed(15323)
ped.tbl <- matrix(1, nrow=1000, ncol =3)
ped.tbl.misc <- matrix(1, nrow=1000, ncol =5)
dad.counter <- 1 #index 1 - 50
mom.counter <- 1 #index 51 - 100
m.counter <- 1

cc.counter <- 1
while(dad.counter < 50 || mom.counter < 50) {
  num.partner <- max(rpois(1,3),2)
  for (i in 1:num.partner) {
    ped.tbl[m.counter,] <- c(m.counter, dad.counter, mom.counter)
    ped.tbl.misc[m.counter,] <- c(m.counter, dad.counter, mom.counter, cc.counter, "A")
    m.counter <- m.counter + 1
    mom.counter <- mom.counter + 1
  }
  dad.counter <- dad.counter + 1

  mom.counter <- mom.counter - 1
  num.partner <- max(rpois(1,3),2)
  for (i in 1:num.partner) {
    ped.tbl[m.counter,] <- c(m.counter, dad.counter, mom.counter)
    ped.tbl.misc[m.counter,] <- c(m.counter, dad.counter, mom.counter, cc.counter, "A")
    m.counter <- m.counter + 1
    dad.counter <- dad.counter + 1
  }
  mom.counter <- mom.counter + 1

  cc.counter <- cc.counter +1
}
message("CC counter", cc.counter-1)

while(dad.counter < 100 || mom.counter < 100) {
  num.partner <- max(rpois(1,3),2)
  for (i in 1:num.partner) {
    num.off <- max(rpois(1,3),2)
    for (j in 1:num.off) {
      ped.tbl[m.counter,] <- c(m.counter, dad.counter, mom.counter)
      ped.tbl.misc[m.counter,] <- c(m.counter, dad.counter, mom.counter, cc.counter, "B")
      m.counter <- m.counter + 1
    }
    mom.counter <- mom.counter + 1
  }
  dad.counter <- dad.counter + 1

  mom.counter <- mom.counter - 1
  num.partner <- max(rpois(1,3),2)
  for (i in 1:num.partner) {
    num.off <- max(rpois(1,3),2)
    for (j in 1:num.off) {
      ped.tbl[m.counter,] <- c(m.counter, dad.counter, mom.counter)
      ped.tbl.misc[m.counter,] <- c(m.counter, dad.counter, mom.counter, cc.counter, "B")
      m.counter <- m.counter + 1
    }
    dad.counter <- dad.counter + 1
  }
  mom.counter <- mom.counter + 1

  cc.counter <- cc.counter +1
}
message("CC counter", cc.counter-1)

final.tbl.cut <- ped.tbl[1:(m.counter-1),] %>% as.tibble()
colnames(final.tbl.cut) <- c("kid", "pa", "ma")

final.tbl.more.cut <- ped.tbl.misc[1:(m.counter-1),] %>% as.tibble()
colnames(final.tbl.more.cut) <- c("kid", "pa", "ma", "CC.indx", "group")


final.tbl.randomize <- final.tbl.cut[sample(1:nrow(final.tbl.cut)),]

final.sim.halfSib.ped <- final.tbl.randomize %>%
  mutate(ma = ma+max(pa)) %>%
  mutate(kid = kid+max(ma))

halfSib.ped.info <- final.tbl.more.cut %>%
  mutate(ma = as.numeric(ma)+max(as.numeric(pa))) %>%
  mutate(kid = as.numeric(kid)+max(as.numeric(ma)))

#############

Sim.and.RunPedFac <- function(seed.num = 234,
                        pedigraph.version = 0.185,
                        pedigraph.path = "/Users/thomasn/repo/pedigree-factor-graphs/src",
                        pedfac.opt = "-f 10",

                        outdir.path="/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/halfsib",

                        n.iter = 10,
                        replicate.indx  = 1,
                        sf.rate = 0.75,
                        ped.df = final.sim.halfSib.ped,
                        n.snp =100, geno.err = 0.02, missing.unk = 0.005,

                        alpha.ad = 10, beta.ad = 10,
                        alpha.missing.indiv = NA, beta.missing.indiv = NA,
                        alpha.weight.missing.loci = NA, beta.weight.missing.loci = NA,


                        rerun= TRUE,
                        pastLS = NA,
                        run.label = "pedfac",
                        case.label = "sf_0.75") {

  if (!rerun && is.na(pastLS) ) {stop("need past LS")}

  run.summary <- list()

  run.summary$param <- tibble::lst(seed.num, pedigraph.version, pedigraph.path, outdir.path, replicate.indx, sf.rate, n.snp, geno.err, missing.unk, pedfac.opt, alpha.ad, beta.ad, alpha.missing.indiv,beta.missing.indiv, alpha.weight.missing.loci,beta.weight.missing.loci, run.label, case.label)

  ## factorize ped.df
  mating.factor <- factor(unlist(ped.df))
  mating.lvl <- levels(mating.factor)
  mating.factor.tbl <- matrix(as.numeric(mating.factor), ncol=3) %>% as.tibble()
  colnames(mating.factor.tbl) <- c("kid","pa","ma")

  outfolder.path <- paste0(outdir.path,"/",case.label,"/",replicate.indx)

  pedfac.out.path <- paste0(outdir.path,"/",case.label,"/",replicate.indx,"/",run.label)

  if (rerun) {
  dir.create(file.path(outfolder.path), recursive = TRUE)
  dir.create(file.path(pedfac.out.path), recursive = TRUE)
  }

  set.seed(seed.num)
  random.arr <- ceiling(c(runif(2,0,1e6)))

  if (rerun) {write.table(random.arr, paste0(outfolder.path,"/","rand"),
                         sep = " ",eol = " ",quote = FALSE, col.names = FALSE, row.names = FALSE)

  geno.input.tbl <- simGeno(mating.tbl = mating.factor.tbl, out.path = outfolder.path,
                            n.snp =n.snp, geno.err = geno.err,
                            random.seed = seed.num,
                            sampling.frac.rate = sf.rate,
                            sample.frac.gen = 1,
                            missing.unk = missing.unk,
                            alpha.ad = alpha.ad,
                            beta.ad = beta.ad,
                            alpha.missing.indiv = alpha.missing.indiv,
                            beta.missing.indiv = beta.missing.indiv,
                            alpha.weight.missing.loci = alpha.weight.missing.loci,
                            beta.weight.missing.loci = beta.weight.missing.loci)
  } else {
    geno.input.tbl <- pastLS$geno.input.tbl
  }

  run.summary$geno.input.tbl <- geno.input.tbl
  ## assume indiv ID is preserves


  pedfac.call <- paste0(pedigraph.path,"/pedigraph_",pedigraph.version, " -d ", outfolder.path,
                        " -n ", n.iter,
                        " -r ",outfolder.path,"/rand ", pedfac.opt, collapse = "")

  if (rerun) {
    tic("run pedFac:")
    system(pedfac.call,ignore.stderr = F, ignore.stdout = F)
    runtime <- toc()
    write.table(runtime$toc - runtime$tic, file =paste0(outfolder.path,"/",run.label,"/pedFac.time"), quote = FALSE, row.names = FALSE, col.names = FALSE)
    run.summary$runtime <- runtime$toc - runtime$tic
  } else {
    run.summary$runtime <- pastLS$runtime
  }

  if (rerun) {
    pedfac.outfile.ls <- list.files(paste0(outfolder.path,"/out"), include.dirs=TRUE, full.names = T)
    file.copy(pedfac.outfile.ls,
              paste0(pedfac.out.path), overwrite = T)
  }

  run.summary$pedfac.cmd <- pedfac.call

  out.ped.raw.tbl <- read.table(paste0(outfolder.path, "/out/ped.txt",collapse = ""), stringsAsFactors = FALSE) %>%
    dplyr::tbl_df() %>%
    dplyr::rename("iter"="V1", "kid"="V2", "pa"="V3", "ma"="V4")

  id.ls <-readRDS(paste0(outfolder.path, "/id.rds",collapse = ""))$id
  max.id <- max(id.ls)

  out.ll.tbl <- read.table(paste0(outfolder.path, "/out/ll.txt",collapse = ""), stringsAsFactors = FALSE) %>%
    dplyr::tbl_df() %>%
    dplyr::rename("iter"="V1", "kid"="V2", "ll"="V3")

  run.summary$out.ll.tbl <- out.ll.tbl

  out.ped.tbl <- out.ped.raw.tbl %>%
    mutate(order = row_number())%>%
    arrange(iter, desc(order))%>%
    group_by(iter, kid)%>%
    summarise(pa = pa[1], ma =ma[1])

  run.summary$out.ped.tbl <- out.ped.tbl


  if (sf.rate != 0) {
    message("working out parentage inference:")
    parentage.grp.infer <- retrieveParent(out.ped.tbl, max.id, id.ls)

    observed.id.ls <- geno.input.tbl$sorted.id
    parentage.grp.truth <- mating.factor.tbl %>% ungroup() %>%
      summarise(kid.id = ifelse(kid %in% observed.id.ls, kid, -1),
                pa.T = ifelse(pa %in% observed.id.ls, pa, -1),
                ma.T = ifelse(ma %in% observed.id.ls, ma, -1))

    run.summary$parent.prob.tbl <- left_join(parentage.grp.infer, parentage.grp.truth, by=c("kid.id")) %>%
      ungroup() %>%
      mutate(n.match = (pa.id == pa.T)+(ma.id == ma.T))
  }

  message("writing fullsib_assignment.txt ")
  full.sib.grp.infer <- RetrieveFullSib(out.ped.tbl, max.id, id.ls)
  half.sib.grp.infer <- RetrieveHalfSib(out.ped.tbl, full.sib.grp.infer, max.id, id.ls)

  write.table(full.sib.grp.infer,
              paste0(pedfac.out.path,"/fullsib_assignment.txt"),
              sep = " ",eol = "\n",quote = FALSE, col.names = FALSE, row.names = FALSE)

  write.table(half.sib.grp.infer,
              paste0(pedfac.out.path,"/halfsib_assignment.txt"),
              sep = " ",eol = "\n",quote = FALSE, col.names = FALSE, row.names = FALSE)

  ## perform pairwise sib-pair from the truth ped
  full.sib.pair.truth.tbl <- RetrieveFullSibTruth(mating.factor.tbl)
  half.sib.pair.truth.tbl <- RetrieveHalfSibTruth(mating.factor.tbl, full.sib.pair.truth.tbl)

  ## run cluster analysis

  run.summary$fullsib.pairwise.tbl <- full_join(full.sib.grp.infer, full.sib.pair.truth.tbl, by= c("kid.1", "kid.2"))
  run.summary$halfsib.pairwise.tbl <- full_join(half.sib.grp.infer, half.sib.pair.truth.tbl, by= c("kid.1", "kid.2"))

  run.summary$fullsib.cluster.ls <- RetrieveFullSibGrp(out.ped.tbl, mating.factor.tbl, max.id, id.ls)

  return(run.summary)
}

### polygamous style

RunSequoia <- function(seed.num = 234,
                       outdir.path="/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/halfsib",
                       max.iter = 20,

                       replicate.indx  = 1,
                       sf.rate = 0.75,
                       ped.df = final.sim.halfSib.ped,
                       n.snp =100, geno.err = 0.02, missing.unk = 0.005,

                       alpha.ad = 10, beta.ad = 10,
                       alpha.missing.indiv = NA, beta.missing.indiv = NA,
                       alpha.weight.missing.loci = NA, beta.weight.missing.loci = NA,

                       rerun= TRUE,
                       pastLS = NA,
                       run.label = "pedfac",
                       case.label = "sf_0.75"
                       ) {

  if (!rerun && is.na(pastLS) ) {stop("need past LS")}

  run.summary <- list()

  run.summary$param <- tibble::lst(seed.num, outdir.path, replicate.indx, max.iter, sf.rate, n.snp, geno.err, missing.unk)

  ## require that pedFac has already run Sim
  ## factorize ped.df
  mating.factor <- factor(unlist(ped.df))
  mating.lvl <- levels(mating.factor)
  mating.factor.tbl <- matrix(as.numeric(mating.factor), ncol=3) %>% as.tibble()
  colnames(mating.factor.tbl) <- c("kid","pa","ma")

  outfolder.path <- paste0(outdir.path,"/",case.label,"/",replicate.indx)
  pedFac.geno.tbl <- read.table(paste0(outfolder.path,"/","geno.txt"), stringsAsFactors = F) %>%
    as_tibble()

  lifeHistory.tbl <- pedFac.geno.tbl %>%
    group_by(V1) %>%
    summarise(ID = V1,
              Sex = as.integer(3-V3),
              BirthYear = V4) %>%
    ungroup() %>%
    select(ID, Sex, BirthYear) %>% data.frame()

  geno.matrix <- pedFac.geno.tbl %>% ungroup() %>% select(-V1, -V2, -V3, -V4, -V5) %>%
    as.matrix

  geno.matrix[geno.matrix == 3] <- -9
  rownames(geno.matrix) <- lifeHistory.tbl$ID

  set.seed(seed.num)

  if (rerun) {
    tic("run sequoia:")

    err.rate <- geno.err
    error.rate.0 <- 1-(err.rate)-(err.rate/2)^2
    error.rate.1 <- (err.rate/2)^2
    error.rate.2 <- err.rate /2
    error.rate.3 <- 1-err.rate

    #ErrM <<- matrix(c(error.rate.0,err.rate,error.rate.1,
    #       error.rate.2, 1-err.rate, error.rate.2,
    #      error.rate.1,err.rate,error.rate.0), ncol = 3)

    ErrM <<- matrix(c(error.rate.3,error.rate.2,error.rate.2,
                      error.rate.2, error.rate.3,error.rate.2,
                      error.rate.2,error.rate.2,error.rate.3), ncol = 3)


    ParOUT_Simp <- sequoia(GenoM = as.matrix(geno.matrix), LifeHistData = lifeHistory.tbl, MaxSibIter = 0,Complex = "simp", Err = geno.err)
    SeqOUT <- sequoia(GenoM = geno.matrix,
                      SeqList = ParOUT_Simp,
                      MaxSibIter = max.iter, Err = geno.err, Complex = "simp")

    runtime <- toc()

    save(ParOUT_Simp, SeqOUT, file = paste0(outfolder.path, "/sequoia_out.rds",collapse = ""))

    run.summary$runtime <- runtime$toc - runtime$tic
  } else{
  run.summary$runtime <- pastLS$runtime
  load(paste0(outfolder.path, "/sequoia_out.rds",collapse = ""))
  }


  id.ls <-readRDS(paste0(outfolder.path, "/id.rds",collapse = ""))$id
  max.id <- max(id.ls)

  geno.sim.tbl <-  read.table(paste0(outfolder.path, "/ingeno.txt",collapse = ""), stringsAsFactors = FALSE)

  observed.id.ls <- geno.sim.tbl$V1

  parentage.grp.truth <- mating.factor.tbl %>% ungroup() %>%
    summarise(kid.id = ifelse(kid %in% observed.id.ls, kid, -1),
              pa.T = ifelse(pa %in% observed.id.ls, pa, -1),
              ma.T = ifelse(ma %in% observed.id.ls, ma, -1),
              id = kid.id) %>%
    filter(kid.id != -1)

  if (sf.rate != 0) {

    message("preparing parentage assignment tbl ....")



    run.summary$parent.tbl <- left_join(parentage.grp.truth, SeqOUT$PedigreePar %>% as_tibble() %>% select(id, sire, dam, LLRpair) %>% ungroup() %>% mutate(id =as.numeric(id)) , by="id") %>%
      ungroup() %>%
      rename(pa.id = "sire", ma.id = "dam") %>%
      mutate(pa.id = ifelse(is.na(pa.id), -1, as.numeric(pa.id)),
             ma.id = ifelse(is.na(ma.id), -1, as.numeric(ma.id)),
             n.match = (pa.id == pa.T)+(ma.id == ma.T))
  }


  sequoia.ped.1 <- inner_join(SeqOUT$Pedigree %>% filter(id %in% parentage.grp.truth$kid.id),
                              data.frame(id=rownames(geno.matrix))) %>% as_tibble() %>%
    select(id, sire, dam, LLRpair) %>% ## using LLRpair as the metric of ordering
    ungroup() %>%
    mutate(sire = ifelse(is.na(sire), row_number(), sire),
           dam = ifelse(is.na(dam), row_number(), dam)) %>%
    mutate(norm.LLRpair = ifelse(is.na(LLRpair),
                                 0,
                                 exp(LLRpair)),
           prob=norm.LLRpair/max(norm.LLRpair)) %>%
    select(-LLRpair, -norm.LLRpair)

  sequoia.ped.2 <-
    sequoia.ped.1 %>%
    mutate(sire = as.numeric(factor(sequoia.ped.1$sire))+max.id,
           dam = as.numeric(factor(sequoia.ped.1$dam))+2*max.id) %>%
    rename(pa = "sire", ma = "dam", kid = "id")

  sib.pairs <- sequoia.ped.2 %>%
    dplyr::group_by(pa, ma) %>%
    dplyr::summarise(full.sib = paste0(kid, collapse = ","),
                     n.sibs = n(),
                     prob = prob[1]) %>%
    dplyr::filter(n.sibs > 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(full.sib, prob)


  full.sib.grp.infer <- apply(sib.pairs,1,
                              function(x) {rbind(combn(as.numeric(unlist(strsplit(x[1],",",fixed=T))),
                                                       m=2,
                                                       simplify = T),x[2])}) %>%
    dplyr::bind_cols() %>%
    t %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(V1 = as.numeric(V1),
                  V2 = as.numeric(V2),
                  V3 = as.numeric(V3)) %>%
    dplyr::mutate(kid.1 = ifelse(V1<V2, V1,V2),
                  kid.2 = ifelse(V1<V2, V2,V1)) %>%
    dplyr::group_by(kid.1, kid.2) %>%
    dplyr::summarise(prob=V3[1]) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(prob)) %>%
    dplyr::mutate(kid.1 = as.character(subbingID(as.numeric(kid.1), max.id, id.ls)),
                  kid.2 = as.character(subbingID(as.numeric(kid.2), max.id, id.ls)))


  write.table(full.sib.grp.infer,
              paste0(outfolder.path,"/fullsib_seq_assignment.txt"),
              sep = " ",eol = "\n",quote = FALSE, col.names = FALSE, row.names = FALSE)

  full.sib.pair.truth.tbl <- RetrieveFullSibTruth(mating.factor.tbl)

  run.summary$fullsib.pairwise.tbl <- full_join(full.sib.grp.infer, full.sib.pair.truth.tbl, by= c("kid.1", "kid.2"))


  sib.pa.pairs <- sequoia.ped.2 %>%
    dplyr::filter(as.numeric(kid)<=max.id) %>%
    dplyr::group_by(pa) %>%
    dplyr::summarise(sib = paste0(kid, collapse = ","),
                     n.sibs = n(),
                     prob = max(prob)) %>%
    dplyr::filter(n.sibs > 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(sib, prob)

  sib.ma.pairs <- sequoia.ped.2 %>%
    dplyr::filter(as.numeric(kid)<=max.id) %>%
    dplyr::group_by(ma) %>%
    dplyr::summarise(sib = paste0(kid, collapse = ","),
                     n.sibs = n(),
                     prob = max(prob[1])) %>%
    dplyr::filter(n.sibs > 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(sib, prob)

  all.sibs <- apply(rbind(sib.pa.pairs, sib.ma.pairs),1,
          function(x) {rbind(combn(as.numeric(unlist(strsplit(x[1],",",fixed=T))),
                                   m=2,
                                   simplify = T),x[2])}) %>%
    dplyr::bind_cols() %>%
    t %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(V1 = as.numeric(V1),
                  V2 = as.numeric(V2),
                  V3 = as.numeric(V3)) %>%
    dplyr::mutate(kid.1 = ifelse(V1<V2, V1,V2),
                  kid.2 = ifelse(V1<V2, V2,V1)) %>%
    dplyr::group_by(kid.1, kid.2)%>%
    dplyr::summarise(prob=max(V3)) %>%
    dplyr::mutate(kid.1 = as.character(subbingID(kid.1, max.id, id.ls)),
                  kid.2 = as.character(subbingID(kid.2, max.id, id.ls)))

  has.half.sib <- anti_join(all.sibs, full.sib.grp.infer , by=c("kid.1", "kid.2"))

  half.sib.grp.infer <- has.half.sib %>%
    dplyr::arrange(desc(prob))

  ##it is important the the prob for half-sib relationship wouldn't be making any more sense


  write.table(half.sib.grp.infer,
              paste0(outfolder.path,"/halfsib_assignment.txt"),
              sep = " ",eol = "\n",quote = FALSE, col.names = FALSE, row.names = FALSE)



  half.sib.pair.truth.tbl <- RetrieveHalfSibTruth(mating.factor.tbl, full.sib.pair.truth.tbl)


  run.summary$halfsib.pairwise.tbl <- full_join(half.sib.grp.infer, half.sib.pair.truth.tbl, by= c("kid.1", "kid.2"))

  run.summary$cluster.ls <- RetrieveFullSibGrpSeq(sequoia.ped.2, mating.factor.tbl, max.id, id.ls)
  return(run.summary)

}

RunColony <- function(seed.num = 234,
                      outdir.path="/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/halfsib",

                      replicate.indx  = 1,
                      sf.rate = 0.75,
                      ped.df = final.sim.halfSib.ped,
                      n.snp =100, geno.err = 0.02, missing.unk = 0.005,

                      alpha.ad = 10, beta.ad = 10,
                      alpha.missing.indiv = NA, beta.missing.indiv = NA,
                      alpha.weight.missing.loci = NA, beta.weight.missing.loci = NA,

                      rerun= TRUE,
                      pastLS = NA,
                      case.label = "sf_0.75",
                      colony.label = "colonialTrial",

                      perl.script.path = "/Users/thomasn/repo/pedigree-factor-graphs/src/SimToColony_Poly.pl",
                      colony.path = "/Users/thomasn/Downloads/colony2.mac_.20180730_1/colony2s.out") {

  if (!rerun && is.na(pastLS) ) {stop("need past LS")}

  run.summary <- list()

  run.summary$param <- tibble::lst(seed.num, outdir.path, replicate.indx, case.label, sf.rate, perl.script.path)

  ## require that pedFac has already run Sim
  ## factorize ped.df
  mating.factor <- factor(unlist(ped.df))
  mating.lvl <- levels(mating.factor)
  mating.factor.tbl <- matrix(as.numeric(mating.factor), ncol=3) %>% as.tibble()
  colnames(mating.factor.tbl) <- c("kid","pa","ma")

  outfolder.path <- paste0(outdir.path,"/",case.label,"/",replicate.indx)
  geno.txt.path <- paste0(outfolder.path,"/","geno.txt")
  pedFac.geno.tbl <- read.table(geno.txt.path, stringsAsFactors = F) %>%
    as_tibble()

  set.seed(seed.num)

  colony.prep.cmd <- paste0("perl ", perl.script.path, " -f ", outfolder.path, " -g 0 -n ",colony.label, collapse = "")

  if (rerun) system(colony.prep.cmd,ignore.stderr = F, ignore.stdout = F)

  if (rerun) dir.create(paste0(outfolder.path,"/colony"), recursive = TRUE)
  if (rerun) file.copy(paste0(outfolder.path,"/colony_0.in"),
                       paste0(outfolder.path,"/colony/colony_0.in"), overwrite = T)

  colony.run.cmd <- paste0('faketime "2008-01-01 00:00:00" ', colony.path,  " IFN:",outfolder.path,"/colony_0.in >  ",outfolder.path,"/colony/colony_run.out")

  if (rerun) {
    tic("run colony:")
    system(colony.run.cmd, ignore.stderr = F, ignore.stdout = F)
    runtime <- toc()
    run.summary$runtime <- runtime$toc - runtime$tic
  } else {
    run.summary$runtime <- pastLS$runtime
  }

  run.summary$colony.cmd <- colony.run.cmd

  if (rerun)  colony.run.files <- list.files(".", pattern = paste0(colony.label))
  if (rerun)  file.copy(colony.run.files, paste0(outfolder.path, "/colony/."), recursive = T, copy.date=T, overwrite = T)

  id.ls <-readRDS(paste0(outfolder.path, "/id.rds",collapse = ""))$id
  max.id <- max(id.ls)


  colony.cluster.tbl <-  read.table(paste0(outfolder.path, "/colony/",colony.label,".BestFSFamily",collapse = ""), stringsAsFactors = FALSE, skip = 1) %>% as_tibble()
  colnames(colony.cluster.tbl) <- c("cluster", "prob.incl", "prob.excl", "fullsib")



  if (sf.rate != 0) {
    message("preparing parentage assignment tbl ....")

    colony.infer.tbl <-  read.csv(paste0(outfolder.path, "/colony/",colony.label,".ParentPair",collapse = ""), stringsAsFactors = FALSE) %>% as_tibble()
    colnames(colony.infer.tbl) <- c("kid", "pa", "ma", "prob")
    colony.infer.tbl <- colony.infer.tbl %>% mutate(
      pa = ifelse(pa=="*",max.id+1, pa),
      ma = ifelse(ma=="#",max.id+1, ma)
    )

    colony.id.infer.tbl <- colony.infer.tbl %>% ungroup() %>%
      mutate(id = subbingID(as.numeric(kid), max.id, id.ls),
             pa.id = subbingID(as.numeric(pa), max.id, id.ls),
             ma.id = subbingID(as.numeric(ma), max.id, id.ls))



    geno.sim.tbl <-  read.table(paste0(outfolder.path, "/ingeno.txt",collapse = ""), stringsAsFactors = FALSE)

    observed.id.ls <- geno.sim.tbl$V1

    parentage.grp.truth <- mating.factor.tbl %>% ungroup() %>%
      summarise(kid.id = ifelse(kid %in% observed.id.ls, kid, -1),
                pa.T = ifelse(pa %in% observed.id.ls, pa, -1),
                ma.T = ifelse(ma %in% observed.id.ls, ma, -1),
                id = kid.id) %>%
      filter(kid.id != -1)

    run.summary$parent.tbl <- left_join(parentage.grp.truth, colony.id.infer.tbl %>% select(id, pa.id, ma.id, prob) %>% ungroup() %>% mutate(id =as.numeric(id)) , by="id") %>%
      ungroup() %>%
      mutate(pa.id = ifelse(is.na(pa.id), -1, as.numeric(pa.id)),
             ma.id = ifelse(is.na(ma.id), -1, as.numeric(ma.id)),
             n.match = (pa.id == pa.T)+(ma.id == ma.T))
  }

  colony.FS.tbl <-  read.table(paste0(outfolder.path, "/colony/",colony.label,".FullSibDyad",collapse = ""), stringsAsFactors = FALSE, skip = 1,sep = ",") %>% as_tibble()
  colnames(colony.FS.tbl) <- c("kid.1", "kid.2", "prob.FS")
  colony.HS.tbl <-  read.table(paste0(outfolder.path, "/colony/",colony.label,".HalfSibDyad",collapse = ""), stringsAsFactors = FALSE, skip = 1,sep = ",") %>% as_tibble()
  colnames(colony.HS.tbl) <- c("kid.1", "kid.2", "prob.HS")

  allsib.pairs.ls <- full_join(colony.FS.tbl, colony.HS.tbl, by=c("kid.1", "kid.2")) %>%
    mutate(prob.FS = ifelse(is.na(prob.FS), 0, prob.FS),
           prob.HS = ifelse(is.na(prob.HS), 0, prob.HS),
           stat = ifelse(prob.FS > prob.HS, "FS", "HS"))

  full.sib.pair.infer <- allsib.pairs.ls %>%
    filter(stat == "FS") %>%
    rename(prob = prob.FS) %>%
    select(kid.1, kid.2, prob) %>%
    dplyr::mutate(kid.1 = as.character(subbingID(kid.1, max.id, id.ls)),
                  kid.2 = as.character(subbingID(kid.2, max.id, id.ls)))

  half.sib.pair.infer <- allsib.pairs.ls %>%
    filter(stat == "HS") %>%
    rename(prob = prob.HS) %>%
    select(kid.1, kid.2, prob) %>%
    dplyr::mutate(kid.1 = as.character(subbingID(kid.1, max.id, id.ls)),
                  kid.2 = as.character(subbingID(kid.2, max.id, id.ls)))

  write.table(full.sib.pair.infer,
              paste0(outfolder.path,"/fullsib_colony_assignment.txt"),
              sep = " ",eol = "\n",quote = FALSE, col.names = FALSE, row.names = FALSE)

  write.table(half.sib.pair.infer,
              paste0(outfolder.path,"/halfsib_colony_assignment.txt"),
              sep = " ",eol = "\n",quote = FALSE, col.names = FALSE, row.names = FALSE)


  ## perform pairwise sib-pair from the truth ped
  full.sib.pair.truth.tbl <- RetrieveFullSibTruth(mating.factor.tbl)
  half.sib.pair.truth.tbl <- RetrieveHalfSibTruth(mating.factor.tbl, full.sib.pair.truth.tbl)
  ## run cluster analysis

  run.summary$fullsib.pairwise.tbl <- full_join(full.sib.pair.infer, full.sib.pair.truth.tbl, by= c("kid.1", "kid.2"))
  run.summary$halfsib.pairwise.tbl <- full_join(half.sib.pair.infer, half.sib.pair.truth.tbl, by= c("kid.1", "kid.2"))

  run.summary$cluster.ls <- RetrieveFullSibGrpCOLONY(colony.cluster.tbl, mating.factor.tbl, max.id, id.ls)
  return(run.summary)

}

### sampling frac
seed.array <- c(234, 334, 434, 534, 634)
sf.array <- c(1, 0.75, 0.5, 0.25, 0) ## must do this one
n.snp.array <- c(200, 50)

#gerr.array <- c(0.01, 0.05, 0.1, 0.25)
#create one that centers on common (default) -> flat, rare
alpha.af.array <- c(1,1)
beta.af.array <- c(1,10)

alpha.missing.indiv.array <- c(1,5); beta.missing.indiv.array <- c(100,100);
alpha.weight.missing.loci.array <- c(0.5,0.1); beta.weight.missing.loci.array <- c(10,10);

# case where an expected mean of 10 loci are missing
#alpha.missing.indiv = 5; beta.missing.indiv = 100;
#alpha.weight.missing.loci = 0.1; beta.weight.missing.loci = 10;


halfsib.sf.pedfac <- lapply(1:5, function(sf) {
  lapply(1:5, function (rep) {
    Sim.and.RunPedFac(replicate.indx = rep,
                      n.iter = 100,
                      seed.num = seed.array[rep],
                      sf.rate = sf.array[sf],
                      case.label = paste0("sf_",sf.array[sf]))
  })
})

#save(halfsib.sf.pedfac, file = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/fullsib/halfsib.sf.pedfac")

halfsib.nSNP.pedfac <- lapply(1:2, function(nSNP) {
  lapply(1:5, function (rep) {
    Sim.and.RunPedFac(replicate.indx = rep,
                      n.iter = 100,
                      seed.num = seed.array[rep],
                      n.snp = n.snp.array[nSNP],
                      case.label = paste0("nSNP_",n.snp.array[nSNP]))
  })
})
#save(halfsib.nSNP.pedfac, file = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/fullsib/halfsib.nSNP.pedfac")

halfsib.af.pedfac <- lapply(1:2, function(af) {
  lapply(1:5, function (rep) {
    Sim.and.RunPedFac(replicate.indx = rep,
                      n.iter = 100,
                      seed.num = seed.array[rep],
                      alpha.ad = alpha.af.array[af],
                      beta.ad = beta.af.array[af],
                      case.label = paste0("af_case",af))
  })
})
save(halfsib.af.pedfac, file = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/fullsib/halfsib.af.pedfac")


halfsib.missingR.pedfac <- lapply(1:2, function(missingR) {
  lapply(1:5, function (rep) {
    Sim.and.RunPedFac(replicate.indx = rep,
                      n.iter = 100,
                      seed.num = seed.array[rep],
                      alpha.missing.indiv = alpha.missing.indiv.array[missingR],
                      beta.missing.indiv = beta.missing.indiv.array[missingR],
                      alpha.weight.missing.loci = alpha.weight.missing.loci.array[missingR],
                      beta.weight.missing.loci = beta.weight.missing.loci.array[missingR],
                      case.label = paste0("gerr_case",missingR))
  })
})
save(halfsib.missingR.pedfac, file = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/fullsib/halfsib.missingR.pedfac")

#### Sequoia


# halfsib.sf.sequoia <- lapply(1:5, function(sf) {
#   lapply(1:5, function (rep) {
#     RunSequoia(replicate.indx = rep,
#                       seed.num = seed.array[rep],
#                       sf.rate = sf.array[sf],
#                       case.label = paste0("sf_",sf.array[sf]))
#   })
# })

#save(halfsib.sf.sequoia, file = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/fullsib/halfsib.sf.sequoia")

halfsib.sf.sequoia.rerun <- lapply(1:5, function(sf) {
  lapply(1:5, function (rep) {
    RunSequoia(replicate.indx = rep,
               seed.num = seed.array[rep],
               sf.rate = sf.array[sf],
               case.label = paste0("sf_",sf.array[sf]),
               rerun = FALSE,
               pastLS = halfsib.sf.sequoia[[sf]][[rep]])
  })
})

save(halfsib.sf.sequoia.rerun, file = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/fullsib/halfsib.sf.sequoia.rerun")


halfsib.nSNP.sequoia <- lapply(1:2, function(nSNP) {
  lapply(1:5, function (rep) {
    RunSequoia(replicate.indx = rep,
                      seed.num = seed.array[rep],
                      n.snp = n.snp.array[nSNP],
                      case.label = paste0("nSNP_",n.snp.array[nSNP]))
  })
})
save(halfsib.nSNP.sequoia, file = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/fullsib/halfsib.nSNP.sequoia")

halfsib.af.sequoia <- lapply(1:2, function(af) {
  lapply(1:5, function (rep) {
    RunSequoia(replicate.indx = rep,
                      seed.num = seed.array[rep],
                      alpha.ad = alpha.af.array[af],
                      beta.ad = beta.af.array[af],
                      case.label = paste0("af_case",af))
  })
})

save(halfsib.af.sequoia, file = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/fullsib/halfsib.af.sequoia")


halfsib.missingR.sequoia <- lapply(1:2, function(missingR) {
  lapply(1:5, function (rep) {
    RunSequoia(replicate.indx = rep,
                      seed.num = seed.array[rep],
                      alpha.missing.indiv = alpha.missing.indiv.array[missingR],
                      beta.missing.indiv = beta.missing.indiv.array[missingR],
                      alpha.weight.missing.loci = alpha.weight.missing.loci.array[missingR],
                      beta.weight.missing.loci = beta.weight.missing.loci.array[missingR],
                      case.label = paste0("gerr_case",missingR))
  })
})
save(halfsib.missingR.sequoia, file = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/fullsib/halfsib.missingR.sequoia")


###COLONY

# halfsib.sf.COLONY <- lapply(1:5, function(sf) {
#   lapply(1:5, function (rep) {
#     RunColony(replicate.indx = rep,
#                seed.num = seed.array[rep],
#                sf.rate = sf.array[sf],
#               colony.label = paste0("sf_",sf.array[sf]),
#                case.label = paste0("sf_",sf.array[sf]))
#   })
# })

#save(halfsib.sf.COLONY, file = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/fullsib/halfsib.sf.COLONY")

halfsib.sf.COLONY.rerun <- lapply(1:5, function(sf) {
  lapply(1:5, function (rep) {
    RunColony(replicate.indx = rep,
              seed.num = seed.array[rep],
              sf.rate = sf.array[sf],
              colony.label = paste0("sf_",sf.array[sf]),
              case.label = paste0("sf_",sf.array[sf]),
              rerun = FALSE,
              pastLS = halfsib.sf.COLONY
              )
  })
})

save(halfsib.sf.COLONY.rerun, file = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/fullsib/halfsib.sf.COLONY.rerun")



halfsib.nSNP.COLONY <- lapply(1:2, function(nSNP) {
  lapply(1:5, function (rep) {
    RunColony(replicate.indx = rep,
               seed.num = seed.array[rep],
               n.snp = n.snp.array[nSNP],
               colony.label = paste0("nSNP_",n.snp.array[nSNP],"r",rep),
               case.label = paste0("nSNP_",n.snp.array[nSNP]))
  })
})
#save(halfsib.nSNP.COLONY, file = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/fullsib/halfsib.nSNP.COLONY")

halfsib.af.COLONY <- lapply(1:2, function(af) {
  lapply(1:5, function (rep) {
    RunColony(replicate.indx = rep,
               seed.num = seed.array[rep],
               alpha.ad = alpha.af.array[af],
               beta.ad = beta.af.array[af],
              colony.label = paste0("af_case",af,"r",rep),
               case.label = paste0("af_case",af))
  })
})

save(halfsib.af.COLONY, file = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/fullsib/halfsib.af.COLONY")


halfsib.missingR.COLONY <- lapply(1:2, function(missingR) {
  lapply(1:5, function (rep) {
    RunColony(replicate.indx = rep,
               seed.num = seed.array[rep],
               alpha.missing.indiv = alpha.missing.indiv.array[missingR],
               beta.missing.indiv = beta.missing.indiv.array[missingR],
               alpha.weight.missing.loci = alpha.weight.missing.loci.array[missingR],
               beta.weight.missing.loci = beta.weight.missing.loci.array[missingR],
              colony.label = paste0("gerr_case",missingR,"r",rep),
               case.label = paste0("gerr_case",missingR))
  })
})
save(halfsib.missingR.COLONY, file = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/fullsib/halfsib.missingR.COLONY")
