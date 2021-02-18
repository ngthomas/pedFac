## mock-up 5 generation pedigree, monogamous, mating with unrelated indiv,

# 30 founders, 15 male , 15 female

# mean sib size - pois mean of 2.5
# continue mating with rbin~0.5
library(tidyverse)
library(xtable)
library(kableExtra)
library(gridExtra)
library(tictoc)

ped.tbl <- matrix(1, nrow=10000, ncol =4)

#ped.tbl.misc <- matrix(1, nrow=1000, ncol =5)
max.id <<- 2
counter <<- 1

set.seed(234)
spawn.kid <- function(pa.indx = 1, ma.indx = 2, gen = 0) {

  if (gen > 6) return()

  n.kids <- rpois(1, 2)

  for (i in 1:n.kids) {
    max.id <<- max.id + 1

    #print(c(max.id ,pa.indx,ma.indx))
    ped.tbl[counter,] <<- c(max.id,pa.indx,ma.indx, gen)

    counter <<- counter + 1
    if(rbernoulli(1, 0.5)) {
      max.id <<- max.id + 1
      if (rbernoulli(1)){ spawn.kid(max.id-1, max.id-2, gen+1)}
      else {spawn.kid(max.id-2, max.id-1, gen+1)}
    }
  }
}
spawn.kid()

final.tbl.cut <- ped.tbl[1:(counter-1),] %>% as.tibble()
colnames(final.tbl.cut) <- c("Kid", "Pa", "Ma" , "gen")

# making some sub to connect some broken parts (important not to change any param)

final.tbl.cut.1 <- final.tbl.cut %>%
  mutate(Kid = case_when (
    Kid == 13~28,
    Kid == 49~80,
    TRUE~Kid))

nonfounder.ls <- final.tbl.cut.1$Kid

add.in.layer.indiv <- c(final.tbl.cut.1 %>%
  filter(gen == 4, !Pa %in% nonfounder.ls) %>%
  pull(Pa)  %>% unique,
  final.tbl.cut.1 %>%
    filter(gen == 4, !Ma %in% nonfounder.ls) %>%
    pull(Ma)  %>% unique)


add.in.more.founders <- tibble(
  Kid = add.in.layer.indiv,
  Pa = (1:length(add.in.layer.indiv)) + max.id,
  Ma = (1:length(add.in.layer.indiv)) + (max.id+length(add.in.layer.indiv)),
  gen = 3)

max.id <<- max.id + (length(add.in.layer.indiv)*2)


add.in.more.kids.rep <- rpois(length(add.in.layer.indiv),0.5)
add.in.more.kids <- add.in.more.founders[rep(row.names(add.in.more.founders), add.in.more.kids.rep),] %>%
  mutate(Kid = max.id +row_number())



final.tbl.cut.add <-
  bind_rows(add.in.more.founders,
            add.in.more.kids,
      final.tbl.cut.1) %>%
  filter(gen >2) %>%
  arrange(gen)

refactor.level <- unique(c(final.tbl.cut.add$Pa, final.tbl.cut.add$Ma, final.tbl.cut.add$Kid)[sort.int(c(final.tbl.cut.add$gen,final.tbl.cut.add$gen,final.tbl.cut.add$gen+1), index.return = T)$ix])



crop.tbl <- final.tbl.cut.add  %>%
  mutate(Kid = as.numeric(factor(Kid, levels = refactor.level)),
         Pa = as.numeric(factor(Pa, levels = refactor.level)),
         Ma = as.numeric(factor(Ma, levels = refactor.level))) %>%
  select(-gen) %>%
  mutate(across(.fns=as.character))

crop.indiv.ls <- unique(c(crop.tbl$Kid, crop.tbl$Pa, crop.tbl$Ma))

n.snp <- 200

multigen.folder.path <- "/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/multi_case/ped1"
snp.multigen.path <- paste0(multigen.folder.path,"/snp_",n.snp)

indiv.censor.ls = c("21", "22", "28", "24", "13", "20", "49", "68", "57", "42", "44","34","58","27","33","32","69", "50", "3", "17", "51", "56")

geno.pick.tbl <- simGeno(mating.tbl=crop.tbl %>%
                           mutate(across(.fns=as.numeric)),
        n.snp = n.snp,
        out.path = snp.multigen.path,
        num.unrelated.added.per.gen = 10,
        indiv.censor.ls = indiv.censor.ls,
        #sampling.frac.rate = c(0.5,0.5, 0.5), # treat as r.v.
        #sample.frac.gen = c(1,2,3),
        random.seed = 432,
        max.unobs.layer = 2)

ped2dot(crop.tbl,
        ShowLabelNodes = crop.indiv.ls,
        pfactorNodeStyle = "invis",
        pfactorEdgeStyle = "invis",
        ObsNodes = intersect(as.numeric(geno.pick.tbl$id), crop.indiv.ls) ,
        outf = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/multi_case/ped1/recur")


############### run pedfac

run.multiped1.summary <- list()

## factorize ped.df
set.seed(240)

random.arr <- ceiling(c(runif(2,0,1e6)))

write.table(random.arr, paste0(snp.multigen.path,"/","rand"),
              sep = " ",eol = " ",quote = FALSE, col.names = FALSE, row.names = FALSE)

run.multiped1.summary$geno.input.tbl <- geno.pick.tbl

## set strictly mono
pedfac.call <- paste0("/Users/thomasn/repo/pedigree-factor-graphs/src/pedigraph_0.185 -d ", snp.multigen.path,
                      " -n 100",
                      " -r ",snp.multigen.path,"/rand -j 1 -f 10", collapse = "")

  tic("run pedFac:")
  system(pedfac.call,ignore.stderr = F, ignore.stdout = F)
  runtime <- toc()
  write.table(runtime$toc - runtime$tic, file =paste0(snp.multigen.path,"/pedFac.time"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  run.multiped1.summary$runtime <- runtime$toc - runtime$tic

  run.multiped1.summary$pedfac.cmd <- pedfac.call

out.ped.raw.tbl <- read.table(paste0(snp.multigen.path, "/pedFac/ped.txt",collapse = ""), stringsAsFactors = FALSE) %>%
  dplyr::tbl_df() %>%
  dplyr::rename("iter"="V1", "kid"="V2", "pa"="V3", "ma"="V4")

id.ls <-readRDS(paste0(snp.multigen.path, "/id.rds",collapse = ""))$id
max.id <- max(id.ls)

out.ll.tbl <- read.table(paste0(snp.multigen.path, "/out/ll.txt",collapse = ""), stringsAsFactors = FALSE) %>%
  dplyr::tbl_df() %>%
  dplyr::rename("iter"="V1", "kid"="V2", "ll"="V3")

run.multiped1.summary$out.ll.tbl <- out.ll.tbl


## check out how it looks the last sweep ped
out.ped.tbl.snap <- out.ped.raw.tbl %>%
  filter(iter ==99) %>%
  mutate(order = row_number())%>%
  arrange(desc(order))%>%
  group_by(kid)%>%
  summarise(pa = pa[1], ma =ma[1]) %>%
  mutate(across(.fns=as.character)) %>%
  rename(Kid = kid, Pa = pa, Ma = ma)

crop.snap.indiv.ls <- unique(c(out.ped.tbl.snap$Kid, out.ped.tbl.snap$Pa, out.ped.tbl.snap$Ma)) %>%
  as.numeric

ped2dot(out.ped.tbl.snap,
        ShowLabelNodes = as.character(crop.snap.indiv.ls) ,
        pfactorNodeStyle = "invis",
        pfactorEdgeStyle = "invis",
        ObsNodes = crop.snap.indiv.ls[crop.snap.indiv.ls<=max.id] %>% as.character,
        outf = paste0(snp.multigen.path, "/lastsweep"))



out.ped.tbl <- out.ped.raw.tbl %>%
  mutate(order = row_number())%>%
  arrange(iter, desc(order))%>%
  group_by(iter, kid)%>%
  summarise(pa = pa[1], ma =ma[1])

run.multiped1.summary$out.ped.tbl <- out.ped.tbl

## retrieve all grandparents of any observed individuals of the three recent generation

prior.txt <- read.delim(paste0(snp.multigen.path,"/prior.txt"), sep="\n")
maxID.indx <- sapply(1:nrow(prior.txt), function(l) grepl("maxID", prior.txt[l,])) %>% which

max.id <- prior.txt[maxID.indx,] %>%
  strsplit(., " ") %>%
  unlist() %>% .[2] %>%
  as.integer()

RetrieveGrandparents.simple <- function(df, ref=FALSE, full.join = FALSE) {
  ped.grandpa <- df %>% dplyr::rename(pa=kid,
                                      grandpa.pa = pa,
                                      grandma.pa = ma)
  ped.grandma <- df %>% dplyr::rename(ma=kid,
                                       grandpa.ma = pa,
                                       grandma.ma = ma)

  join.label.1 <- "pa"
  join.label.2 <- "ma"
  if (!ref) {
    join.label.1[2] <- "iter"
    join.label.2[2] <- "iter"
  }

  if (full.join) {
    ped.join.pa<- dplyr::full_join(df, ped.grandpa, by=join.label.1)
    dplyr::full_join(ped.join.pa, ped.grandma, by=join.label.2)

  } else {
  ped.join.pa<- dplyr::inner_join(df, ped.grandpa, by=join.label.1)
  dplyr::inner_join(ped.join.pa, ped.grandma, by=join.label.2)
  }
}


## identify any founder that held in between ped and add into the -1 category


intermed.founder.ls <- final.tbl.cut.add  %>%
  mutate(Kid = as.numeric(factor(Kid, levels = refactor.level)),
         Pa = as.numeric(factor(Pa, levels = refactor.level)),
         Ma = as.numeric(factor(Ma, levels = refactor.level))) %>%
  filter(gen >=5) %>%
  group_by(Kid) %>%
  mutate(
    founder.pa = (!Pa %in% as.numeric(crop.tbl$Kid)) &  (!Pa %in% as.numeric(indiv.censor.ls)),
    founder.ma = (!Ma %in% as.numeric(crop.tbl$Kid)) &  (!Ma %in% as.numeric(indiv.censor.ls))
    ) %>%
  filter(founder.pa || founder.ma)


prep.truth.tbl <- bind_rows(
intermed.founder.ls %>%
  filter(founder.pa) %>%
  mutate(Kid = Pa,
         Pa = -1,
         Ma = -1) %>%
  group_by(Kid) %>%
  slice_head() %>%
  select(Kid, Pa , Ma),
intermed.founder.ls %>%
  filter(founder.ma) %>%
  mutate(Kid = Ma,
         Pa = -1,
         Ma = -1)%>%
  group_by(Kid) %>%
  slice_head() %>%
  select(Kid, Pa , Ma),
tibble(Kid= -1,
       Pa = -1,
       Ma = -1),
crop.tbl %>%
  mutate(across(.fns=as.numeric))
)%>%
  rename(pa=Pa, kid=Kid, ma=Ma)




infer.grandparent <- RetrieveGrandparents.simple(out.ped.tbl)
truth.grandparent <- RetrieveGrandparents.simple(prep.truth.tbl, TRUE)

## let's focus on grandkids that are observed first, then we can look at special cases

infer.grandparent.summary <- infer.grandparent %>%
  filter(kid <max.id ) %>%
  dplyr::transmute(kid = kid,
                   grandpa.pa = ifelse(grandpa.pa>=max.id,-1,grandpa.pa),
                   grandma.pa = ifelse(grandma.pa>=max.id,-1,grandma.pa),
                   grandpa.ma = ifelse(grandpa.ma>=max.id,-1,grandpa.ma),
                   grandma.ma = ifelse(grandma.ma>=max.id,-1,grandma.ma),
                   pa = ifelse(pa>=max.id,-1,pa),
                   ma = ifelse(ma>=max.id,-1,ma),
                   flip = ifelse(pa == -1 & ma == -1 & paste0(grandpa.pa,grandma.pa) > paste0(grandpa.ma,grandma.ma),
                                 T, F))

infer.grandparent.summary.1 <- bind_rows(
  infer.grandparent.summary %>%
    filter(flip) %>%
    mutate(tem.1 = grandpa.pa,
           tem.2 = grandma.pa,
           grandpa.pa = grandpa.ma,
           grandma.pa = grandma.ma,
           grandpa.ma = tem.1,
           grandma.ma = tem.2) %>%
    select(-tem.1, -tem.2),
  infer.grandparent.summary %>%
    filter(!flip)) %>%
  dplyr::group_by(kid, grandpa.pa, grandma.pa, grandpa.ma, grandma.ma) %>% ##keep parent0 part of the eqn
  dplyr::summarise(n=n(), posterior = n()/100) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::group_by(kid) %>%
  dplyr::top_n(1) %>%
  dplyr::sample_n(1)



truth.grandparent.summary <- truth.grandparent %>%
  dplyr::filter(!kid %in% as.numeric(indiv.censor.ls)) %>%
  dplyr::transmute(kid = kid,
                   grandpa.pa = ifelse(grandpa.pa %in% as.numeric(indiv.censor.ls),-1,grandpa.pa),
                   grandma.pa = ifelse(grandma.pa %in% as.numeric(indiv.censor.ls),-1,grandma.pa),
                   grandpa.ma = ifelse(grandpa.ma %in% as.numeric(indiv.censor.ls),-1,grandpa.ma),
                   grandma.ma = ifelse(grandma.ma %in% as.numeric(indiv.censor.ls),-1,grandma.ma),
                   pa = ifelse(pa %in% as.numeric(indiv.censor.ls),-1,pa),
                   ma = ifelse(ma %in% as.numeric(indiv.censor.ls),-1,ma),
                   flip = ifelse(pa == -1 & ma == -1 & paste0(grandpa.pa,grandma.pa) > paste0(grandpa.ma,grandma.ma),
                                 T, F))

truth.grandparent.summary.1 <- bind_rows(
  truth.grandparent.summary %>%
    filter(flip) %>%
    mutate(tem.1 = grandpa.pa,
           tem.2 = grandma.pa,
           grandpa.pa = grandpa.ma,
           grandma.pa = grandma.ma,
           grandpa.ma = tem.1,
           grandma.ma = tem.2) %>%
    select(-tem.1, -tem.2),
  truth.grandparent.summary %>%
    filter(!flip)) %>%
  select(kid, grandpa.pa, grandma.pa, grandpa.ma, grandma.ma)


ped.join.truth <- dplyr::left_join(infer.grandparent.summary.1, truth.grandparent.summary.1, by="kid") %>%
  mutate(across(.fn=~ifelse(is.na(.),-1,.))) %>% group_by(kid) %>%
  mutate(n.matches =
           (grandpa.pa.x == grandpa.pa.y) +
           (grandpa.ma.x == grandpa.ma.y) +
           (grandma.pa.x == grandma.pa.y) +
           (grandma.ma.x == grandma.ma.y) )


#View(ped.join.truth)


highlight.unobs.ls <- c("21", "22", "28", "24", "49", "68", "57", "42", "44","34","58","27","33","32", "50", "51")

highlight.kids.ls <- RetrieveGrandparents.simple(crop.tbl %>%
                                                       rename(pa=Pa, kid=Kid, ma=Ma) %>%
                                                       mutate(across(.fns=as.numeric)), TRUE, full.join = T) %>%
  filter(pa %in% highlight.unobs.ls | ma %in% highlight.unobs.ls, !is.na(kid))



highlight.unobs.tbl <- bind_rows(
 left_join(highlight.kids.ls , ped.join.truth, by="kid") %>%
   filter(pa %in% highlight.unobs.ls) %>%
   filter(!is.na(n.matches)) %>%
   group_by(pa) %>%
   summarise(posterior = max(posterior, na.rm=T),
             n.match = max(n.matches, na.rm=T)) %>%
   rename(indiv = pa),
 left_join(highlight.kids.ls , ped.join.truth, by="kid") %>%
   filter(ma %in% highlight.unobs.ls) %>%
   filter(!is.na(n.matches)) %>%
   group_by(ma) %>%
   summarise(posterior = max(posterior, na.rm=T),
             n.match = max(n.matches, na.rm=T)) %>%
   rename(indiv = ma)
)

### specifically focus on indiv  51

case_51.result <- left_join(
out.ped.tbl %>%
  filter(kid == 70) %>%
  group_by(iter) %>%
  summarise(kid = ma),
infer.grandparent ) %>%
  group_by(iter) %>%
  mutate(flip = ifelse(paste0(grandpa.pa,grandma.pa) > paste0(grandpa.ma,grandma.ma),
                T, F))

case_51.result.b <- bind_rows(
  case_51.result %>%
    filter(flip) %>%
    mutate(tem.1 = grandpa.pa,
           tem.2 = grandma.pa,
           grandpa.pa = grandpa.ma,
           grandma.pa = grandma.ma,
           grandpa.ma = tem.1,
           grandma.ma = tem.2) %>%
    select(-tem.1, -tem.2),
  case_51.result %>%
    filter(!flip)) %>%
  dplyr::group_by(grandpa.pa, grandma.pa, grandpa.ma, grandma.ma) %>% ##keep parent0 part of the eqn
  dplyr::summarise(n=n(), posterior = n()/100) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::top_n(1) %>%
  dplyr::sample_n(1)

highlight.unobs.tbl$posterior[highlight.unobs.tbl$indiv %in% 51] <- case_51.result.b$posterior

highlight.unobs.style.tbl <-  full_join(
  highlight.unobs.tbl, tibble(indiv=as.numeric(highlight.unobs.ls))) %>%
  group_by(indiv) %>%
  summarise(label.str = ifelse(is.na(posterior), "", paste('xlabel="',posterior,'"',",label_scheme=1,", sep = "", collapse = ""))) %>%
  mutate(indiv = as.character(indiv))

intermed.parent.ls <- tibble(indiv= c("69","3", "13", "17", "17"),
                              kid = c(83, 36, 23, 31, 37))
## add in immediate parent
posterior.intermed.result <- out.ped.tbl %>%
  filter(kid %in% c(83, 36, 23, 31, 37)) %>%
  mutate(pa = ifelse(pa > max.id, -1, pa),
         ma = ifelse(ma > max.id, -1, ma)) %>%
  group_by(kid, pa, ma ) %>%
  summarise(posterior = n()/100) %>%
  arrange(desc(posterior)) %>%
  left_join(., intermed.parent.ls) %>%
  group_by(indiv) %>%
  summarise(label.str = ifelse(is.na(posterior[1]), "", paste('xlabel="',max(posterior),'"',",label_scheme=1,", sep = "", collapse = ""))) %>%
  mutate(indiv = as.character(indiv))





ped2dot(crop.tbl,
        ShowLabelNodes = crop.indiv.ls,
        pfactorNodeStyle = "invis",
        pfactorEdgeStyle = "invis",
        ObsNodes = intersect(as.numeric(geno.pick.tbl$id), crop.indiv.ls),,
        outf = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/multi_case/ped1/multiGen_pedFac",
        highlight.unobs.edge = TRUE,
        highlightNodes = c("21", "22", "28", "24", "49", "68", "57", "42", "44","34","58","27","33","32", "50", "51", "3", "13", "17", "69", "20"),
        highlight.unobs.style = bind_rows(highlight.unobs.style.tbl, posterior.intermed.result)
)


### add 36 & 23

write_rds(highlight.unobs.style.tbl, "~/Downloads/highlight.unobs.style.tbl")

write_rds(run.multiped1.summary, "/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/multi_case/ped1/snp_200/pedFac.summary")




###### check out result for sequoia

seq.summary <- list()


## require that pedFac has already run Sim
## factorize ped.df


pedFac.geno.tbl <- read.table(paste0("/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/multi_case/ped1/snp_200/geno.txt"), stringsAsFactors = F) %>%
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

set.seed(234)

  tic("run sequoia:")

  err.rate <- 0.02
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
  ParOUT_mono <- sequoia(GenoM = as.matrix(geno.matrix), LifeHistData = lifeHistory.tbl, MaxSibIter = 0,Complex = "mono", Err = geno.err)

  ParOUT_Simp <- sequoia(GenoM = as.matrix(geno.matrix), LifeHistData = lifeHistory.tbl, MaxSibIter = 0,Complex = "simp", Err = geno.err)
  SeqOUT <- sequoia(GenoM = geno.matrix,
                    SeqList = ParOUT_Simp,
                    MaxSibIter = 20, Err = geno.err, Complex = "simp")

  runtime <- toc()

  outfolder.path <- "/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/multi_case/ped1"
  save(ParOUT_Simp, SeqOUT, file = paste0(outfolder.path, "/sequoia_out.rds",collapse = ""))

  seq.summary <- runtime$toc - runtime$tic

# fill in unobs node

  all.indiv.uncensor.ls <- as.character(pedFac.geno.tbl$V1)

seq.ped.conform.1 <- SeqOUT$Pedigree %>%
  as.tibble() %>%
  select(id, dam, sire, LLRpair) %>%
  rename(Kid = id) %>%
  inner_join(., crop.tbl) %>%
  group_by(Kid) %>%
  mutate(
    Ma.rep = ifelse( (!dam %in% all.indiv.uncensor.ls ) && Ma %in% indiv.censor.ls ,Ma,dam),
    Pa.rep = ifelse( (!sire %in% all.indiv.uncensor.ls ) && Pa %in% indiv.censor.ls ,Pa,sire),
    is.Ma.sub = (Ma.rep == Ma && Ma %in% indiv.censor.ls),
    is.Pa.sub = (Pa.rep == Pa && Pa %in% indiv.censor.ls)
  )

seq.id.sub.ls <- bind_rows(
seq.ped.conform.1 %>%
  filter(is.Ma.sub) %>%
  select(dam, Ma.rep) %>%
  filter(!is.na(dam)) %>%
  select(dam, Ma.rep) %>%
  rename(id = dam, rep.id = Ma.rep),
seq.ped.conform.1 %>%
  filter(is.Pa.sub) %>%
  select(sire, Pa.rep) %>%
  filter(!is.na(sire)) %>%
  select(sire, Pa.rep) %>%
  rename(id = sire, rep.id = Pa.rep)
) %>%
  group_by(id) %>%
  summarise(rep.id = rep.id[1])

seq.ped.conform.2 <- SeqOUT$DummyIDs %>%
  as.tibble() %>%
  select(id, dam, sire, LLRpair) %>%
  left_join(., seq.id.sub.ls) %>%
  filter(!is.na(rep.id)) %>%
  mutate(id = rep.id) %>%
  select(-rep.id) %>%
  rename(Kid = id) %>%
  inner_join(., crop.tbl) %>%
  group_by(Kid) %>%
  mutate(
    Ma.rep = ifelse( (!dam %in% all.indiv.uncensor.ls ) && Ma %in% indiv.censor.ls ,Ma,dam),
    Pa.rep = ifelse( (!sire %in% all.indiv.uncensor.ls ) && Pa %in% indiv.censor.ls ,Pa,sire),
    is.Ma.sub = (Ma.rep == Ma && Ma %in% indiv.censor.ls),
    is.Pa.sub = (Pa.rep == Pa && Pa %in% indiv.censor.ls)
  )


seq.ped.combn.result <- bind_rows(seq.ped.conform.1,
          seq.ped.conform.2) %>%
  mutate(n.match = (Pa == Pa.rep) + (Ma == Ma.rep)) %>%
  select(Kid, Ma.rep, Pa.rep, LLRpair) %>%
  left_join(crop.tbl, .)

# tease out what's being omitted
kid.rm.ls <- seq.ped.combn.result %>%
  filter(( is.na(Ma.rep) | Ma.rep != Ma) | (is.na(Pa.rep) | Pa.rep != Pa)) %>%
  pull(Kid)

parent.rm.ls <- seq.ped.combn.result %>%
  filter(( is.na(Ma.rep) | Ma.rep != Ma) | (is.na(Pa.rep) | Pa.rep != Pa)) %>%
  group_by(Pa, Ma) %>%
  summarise(n.rm = n()) %>%
  ungroup() %>%
  left_join(.,
            crop.tbl %>%
              group_by(Pa, Ma) %>%
              summarise(n.kid = n()) ) %>%
  filter(n.rm == n.kid) %>%
  mutate(mnode = paste0(Pa, "x", Ma))

sequoia.unobs.style.tbl <-
  seq.ped.combn.result %>%
  filter(!is.na(Ma.rep), !is.na(Pa.rep)) %>%
  # mutate(norm.LLRpair = ifelse(is.na(LLRpair),
  #                              0,
  #                              exp(LLRpair)),
  #        prob=norm.LLRpair/max(norm.LLRpair)) %>%
  pivot_longer(cols=Pa.rep:Ma.rep) %>%
  filter(value %in% indiv.censor.ls) %>%
  rename(indiv = value) %>%
  group_by(indiv) %>%
  summarise(label.str = paste('xlabel="',max(LLRpair),'"',",label_scheme=1,", sep = "", collapse = "")) %>%
  mutate(indiv = as.character(indiv))




indiv.censor.ls <- c("21", "22", "28", "24", "13", "20", "49", "68", "57", "42", "44","34","58","27","33","32","69", "50", "3", "17", "51", "56")
all.indiv.uncensor.ls <- pedFac.geno.tbl$V1
highlight.unobs.ls <- c("21", "22", "28", "24", "49", "68", "57", "42", "44","34","58","27","33","32", "50", "51", "3", "13","17")
crop.indiv.ls <- unique(c(crop.tbl$Kid, crop.tbl$Pa, crop.tbl$Ma))



highlight.unobs.style.tbl <-  full_join(
  highlight.unobs.tbl,
  tibble(indiv=as.numeric(highlight.unobs.ls))) %>%
  group_by(indiv) %>%
  summarise(label.str = ifelse(is.na(posterior), "", paste('xlabel="',posterior,'"',",label_scheme=1,", sep = "", collapse = ""))) %>%
  mutate(indiv = as.character(indiv))


ped2dot(crop.tbl,
        ShowLabelNodes = crop.indiv.ls,
        pfactorNodeStyle = "invis",
        pfactorEdgeStyle = "invis",
        ObsNodes = intersect(as.numeric(geno.pick.tbl$id), crop.indiv.ls),
        outf = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/multi_case/ped1/multiGen_sequoia",
        highlight.unobs.edge = TRUE,
        highlightNodes = c("21", "22", "28", "24", "49", "68", "57", "42", "44","34","58","27","33","32", "50", "51", "3", "13", "17", "69", "20"),
        highlight.unobs.style =sequoia.unobs.style.tbl,
        opt.turnoff.edge = TRUE,
        m_node.invis.ls =  parent.rm.ls$mnode,
        invis.ls.kid = kid.rm.ls,
        invis.ls.pa =  parent.rm.ls$Pa,
        invis.ls.ma =  parent.rm.ls$Ma
)


### Work on franz
# i have update simtoFranz_v200.pl

# run perl SimToFranz_v200.pl -f /Users/thomasn/repo/pedigree-factor-graphs/data/ch2/multi_case/ped1/snp_200
#Run the following:
 #mkdir -p /Users/thomasn/repo/pedigree-factor-graphs/data/ch2/multi_case/ped1/snp_200/FRANz_run; cd /Users/thomasn/repo/pedigree-factor-graphs/data/ch2/multi_case/ped1/snp_200/FRANz_run

#FRANz --N 50 --femrepro 1:100 --malerepro 1:100 --updatefreqs --typingerror 0.02 --seed 234 --pedigreeoutformat 1,2,3 ../franz_100.in

#FRANz --N 50 --femrepro 1:100 --malerepro 1:100 --updatefreqs --typingerror    109.34s user 0.10s system 99% cpu 1:49.50 total

#I did edit add extra column header for parentage.csv!!!

franz.result.ls <- list()
franz.result.ls$runtime <- 110
all.female.ls <- crop.tbl$Ma %>% unique
all.male.ls <- crop.tbl$Pa %>% unique

FRANz.parentage.result.0 <- read.csv("/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/multi_case/ped1/snp_200/FRANz_run/parentage.csv", row.names = NULL) %>% as.tibble() %>%
  select(Offspring, Parent.1, Parent.2, Posterior) %>%
  rename(Kid = Offspring, dam = Parent.1, sire =Parent.2) %>%
  ungroup() %>%
  mutate(across(.cols=Kid:sire, .fns=as.character)) %>%
  mutate(need.switch = ifelse (sire %in% all.female.ls | dam %in% all.male.ls, 1, 0 ))

FRANz.parentage.result <- bind_rows(
  FRANz.parentage.result.0 %>%
  filter(need.switch == 0),
  FRANz.parentage.result.0 %>%
    filter(need.switch == 1) %>%
    rename(sire = dam, dam= sire))


FRANz.ped.combn.result <- FRANz.parentage.result %>%
  inner_join(., crop.tbl) %>%
  group_by(Kid) %>%
  mutate(
    Ma.rep = ifelse( (!dam %in% all.indiv.uncensor.ls ) && Ma %in% indiv.censor.ls ,Ma,dam),
    Pa.rep = ifelse( (!sire %in% all.indiv.uncensor.ls ) && Pa %in% indiv.censor.ls ,Pa,sire),
    is.Ma.sub = (Ma.rep == Ma && Ma %in% indiv.censor.ls),
    is.Pa.sub = (Pa.rep == Pa && Pa %in% indiv.censor.ls)
  )  %>%
  mutate(n.match = (Pa == Pa.rep) + (Ma == Ma.rep)) %>%
  select(Kid, Ma.rep, Pa.rep, Posterior) %>%
  left_join(crop.tbl, .)

# tease out what's being omitted
FRANz.kid.rm.ls <- FRANz.ped.combn.result %>%
  filter(( is.na(Ma.rep) | Ma.rep != Ma) | (is.na(Pa.rep) | Pa.rep != Pa)) %>%
  pull(Kid)

FRANz.parent.rm.ls <- FRANz.ped.combn.result %>%
  filter(( is.na(Ma.rep) | Ma.rep != Ma) | (is.na(Pa.rep) | Pa.rep != Pa)) %>%
  group_by(Pa, Ma) %>%
  summarise(n.rm = n()) %>%
  ungroup() %>%
  left_join(.,
            crop.tbl %>%
              group_by(Pa, Ma) %>%
              summarise(n.kid = n()) ) %>%
  filter(n.rm == n.kid) %>%
  mutate(mnode = paste0(Pa, "x", Ma))

FRANz.unobs.style.tbl <-
  FRANz.ped.combn.result %>%
  filter(!is.na(Ma.rep), !is.na(Pa.rep)) %>%
  # mutate(norm.LLRpair = ifelse(is.na(LLRpair),
  #                              0,
  #                              exp(LLRpair)),
  #        prob=norm.LLRpair/max(norm.LLRpair)) %>%
  pivot_longer(cols=Pa.rep:Ma.rep) %>%
  filter(value %in% indiv.censor.ls) %>%
  rename(indiv = value) %>%
  group_by(indiv) %>%
  summarise(label.str = paste('xlabel="',max(Posterior),'"',",label_scheme=1,", sep = "", collapse = "")) %>%
  mutate(indiv = as.character(indiv))


ped2dot(crop.tbl,
        ShowLabelNodes = crop.indiv.ls,
        pfactorNodeStyle = "invis",
        pfactorEdgeStyle = "invis",
        ObsNodes = intersect(as.numeric(geno.pick.tbl$id), crop.indiv.ls),
        outf = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/multi_case/ped1/multiGen_FRANz",
        highlight.unobs.edge = TRUE,
        highlightNodes = c("21", "22", "28", "24", "49", "68", "57", "42", "44","34","58","27","33","32", "50", "51", "3", "13", "17", "69", "20"),
        highlight.unobs.style =FRANz.unobs.style.tbl,
        opt.turnoff.edge = TRUE,
        m_node.invis.ls =  FRANz.parent.rm.ls$mnode,
        invis.ls.kid = FRANz.kid.rm.ls,
        invis.ls.pa =  FRANz.parent.rm.ls$Pa,
        invis.ls.ma =  FRANz.parent.rm.ls$Ma
)
