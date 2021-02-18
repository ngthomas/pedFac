

library(tidyverse)
library(tictoc)


simMendelGeno <- function(
                    n.snp = 10,
                    alpha.ad = 10, beta.ad = 10,
                    geno.err = 0.02,
                    missing.unk = 0.005,
                    random.seed = 46,
                    alpha.missing.indiv = NA, beta.missing.indiv = NA,
                    alpha.weight.missing.loci = NA, beta.weight.missing.loci = NA,
                    num.unrelated.added.per.gen = 0,
                    indiv.censor.ls  = NA,
                    sampling.frac.rate = 1, # treat as r.v.
                    sample.frac.gen = -1, # which generation sample rate for - 0 being the most recent generation
                    one.side.censorship = 0, # 1 - remove all male; 2 - remove all female
                    out.path = tempdir(),
                    max.unobs.layer = 1) {
  set.seed(random.seed)

  # sampling background alllelic frequency
  af <- rbeta(n.snp, alpha.ad, beta.ad)
  maf <- pmin(1-af,af) #minor
  major.af <- 1-maf


  m1 <- tibble(
    Chrom = "Unk",
    Locus = paste0("Unk-", 1:n.snp),
    Pos = NA_real_,
    LocIdx = 1:n.snp,
    a1 = major.af,
    a2 = maf
  ) %>%
    pivot_longer(
      cols = c(a1, a2),
      names_to = "Allele",
      values_to = "Freq"
    ) %>%
    group_by(Locus) %>%
    mutate(AlleIdx = 1:n()) %>%
    ungroup() %>%
    select(Chrom, Locus, Pos, Allele, LocIdx, AlleIdx, Freq)

  markers <- gdropR::reindex_markers(m1)

  genome <- gdropR::chinook_chromosomes %>%
    mutate(
      idx = 1:n(),
      chrom = name1,
      num_bases = length,
      scaled_length = num_bases / max(num_bases)
    ) %>%
    select(idx, chrom, scaled_length, num_bases)

  mapped_markers <- gdropR::sprinkle_markers_into_genome(
    markers,
    genome
  ) %>%
    group_by(Chrom) %>%
    mutate(Locus = paste0(Chrom, "-", rep(1:(n()/2), each = 2))) %>%
    ungroup()

  #sex:1 male 2:female
  pedigree <- gdropR::ped5gen_acyclic %>%
    gdropR::add_explicit_founder_parents() %>%
    gdropR::add_pedigree_sex_column(fem_prob = 1.0)

  long_genos <- gdropR::simulate_linked_genotypes(
    mapped_markers,
    pedigree
  ) %>%
    pivot_longer(
      cols = c(-indiv, -sex),
      names_to = "Locus",
      values_to = "Genotype"
    )

  err_genos <- long_genos %>%
    mutate(
      obs_geno = gdropR::biallelic_geno_error(
        genos = Genotype,
        het_miscall = geno.err,
        hom_miscall = 0.005
      )
    )

  ped_with_gens <- pedigree %>%
    mutate(
      generation = case_when(
        !is.na(generation) ~ as.integer(generation),
        (Pa == "0" | Ma == "0") & as.integer(Kid) < 1000 ~ 0L,
        as.integer(Kid) > 1000 ~ as.integer(str_sub(Kid, 1, 1)) - 2L,
        TRUE ~ NA_integer_
      )
    )

  geno.1 <- err_genos %>%
    mutate(obs012 = case_when(
      obs_geno %in% c("1/1","1/2") ~ 0L,
      obs_geno %in% c("2/1", "2/2") ~ 1L,
      TRUE ~ NA_integer_
    )) %>%
    select(-Genotype, -obs_geno) %>%
    pivot_wider(
      names_from = Locus,
      values_from = obs012
    ) %>%
    select(-indiv, -sex)

  geno.2 <- err_genos %>%
    mutate(obs012 = case_when(
      obs_geno %in% c("1/1","2/1") ~ 0L,
      obs_geno %in% c("1/2", "2/2") ~ 1L,
      TRUE ~ NA_integer_
    )) %>%
    select(-Genotype, -obs_geno) %>%
    pivot_wider(
      names_from = Locus,
      values_from = obs012
    )%>%
    select(-indiv, -sex)

  id.ls <- ped_with_gens$Kid
  n.id <- length(id.ls)
  indiv.sex <- ped_with_gens$Sex

  is.missing <- rbinom(n.id*n.snp, 1, missing.unk) %>% matrix(ncol=n.snp, by=T)
  is.avail <- 1-is.missing

  final.geno.1 <- (abs(geno.1-make.err.1)*is.avail) +(-1*is.missing)
  final.geno.2 <- (abs(geno.2-make.err.2)*is.avail) +(-1*is.missing)

  geno.long <- cbind(final.geno.1, final.geno.2)[,rbind(1:n.snp, (1:n.snp)+n.snp) %>% as.vector()]

  geno.input.tbl <- cbind(id = id.ls,
                          obs = rep(1, n.id),
                          sex = ped_with_gens$Sex,
                          gen =  ped_with_gens$generation,
                          geno.long) %>% as.tibble()

  gen.indx <- ped_with_gens$generation
  total.gen <- max(gen.indx)+1

  ## add in unrelated indiv to any generation, except to the most recent one
  sample.frac <- rep(1, total.gen)

  if (one.side.censorship == 0) {
  # process sampling frac
  if(sample.frac.gen[1] != -1) {
    sample.frac[sample.frac.gen+1] <- sampling.frac.rate
    gen.indx <- ped_with_gens$generation
    is.avail <- rbinom(length(gen.indx),1,sample.frac[gen.indx+1])==1
    geno.input.tbl <- geno.input.tbl[is.avail, ]
  }
  } else {
    sample.frac.tbl <- ped_with_gens %>%
      group_by(generation) %>%
      mutate(n.indiv = n()) %>%
      ungroup() %>%
      filter(Sex != one.side.censorship) %>%
      group_by(generation) %>%
      summarise(sample.frac = n()/n.indiv[1])

    sample.frac[sample.frac.tbl$generation+1] <- sample.frac.tbl$sample.frac

    geno.input.tbl <- geno.input.tbl %>%
      filter(sex != one.side.censorship)
  }

  if (!is.na(indiv.censor.ls )) {
    geno.input.tbl <- geno.input.tbl %>%
      filter(!id %in% indiv.censor.ls)
  }

  dir.create(file.path(out.path), recursive = TRUE)
  geno.path <- paste0(out.path, "/ingeno.txt")
  write.table(geno.input.tbl,
              geno.path,sep = " ", eol = "\n", quote = FALSE, col.names = FALSE, row.names = FALSE)

  param<-list(geno.path=geno.path,
              output.path=out.path,
              random.seed=random.seed,
              observe.frac= sample.frac,
              af = maf,
              max.unobs=max.unobs.layer, max.gen=0, #max.gen as the number extend past the current generation grp
              min.age=1, max.age=1, geno.err=geno.err,
              n.marr = 0)

  message("writing final files to ", out.path)

  writeIntermedGeno(param, preservesID = TRUE, preservesGeno = TRUE)

  return(geno.input.tbl)
}




### run pedFac

run.pedFac.multiGen <- function(n.snp = 100) {

  snp.multigen.path <- paste0("/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/multi_case/ped2/snp", n.snp)

  simMendelGeno(n.snp = n.snp, out.path =snp.multigen.path, one.side.censorship = 1)


run.multiped.summary <- list()
set.seed(240)
random.arr <- ceiling(c(runif(2,0,1e6)))


write.table(random.arr, paste0(snp.multigen.path,"/","rand"),
            sep = " ",eol = " ",quote = FALSE, col.names = FALSE, row.names = FALSE)

## set strictly mono
pedfac.call <- paste0("/Users/thomasn/repo/pedigree-factor-graphs/src/pedigraph_0.185 -d ", snp.multigen.path,
                      " -n 100",
                      " -r ",snp.multigen.path,"/rand -j 1 -f 10", collapse = "")

tic("run pedFac:")
system(pedfac.call,ignore.stderr = F, ignore.stdout = F)
runtime <- toc()
write.table(runtime$toc - runtime$tic, file =paste0(snp.multigen.path,"/pedFac.time"), quote = FALSE, row.names = FALSE, col.names = FALSE)
run.multiped.summary$runtime <- runtime$toc - runtime$tic
run.multiped.summary$pedfac.cmd <- pedfac.call

return(run.multiped.summary)
}

snp.array <- c(100, 200, 400, 800, 1200)

multiped.summary <- lapply(snp.array, function(x) run.pedFac.multiGen(x))


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
