


library(tidyverse)
library(tictoc)
library(sequoia)
library(ggh4x)
library(kableExtra)
library(gdropR)


fourgen.ped <- gdropR::FourGen_555_pedigree %>%
  mutate(Kid = as.numeric(factor(Kid,levels = unique(c(FourGen_555_pedigree$Pa, FourGen_555_pedigree$Ma, FourGen_555_pedigree$Kid)))),
         Pa = as.numeric(factor(Pa,levels = unique(c(FourGen_555_pedigree$Pa, FourGen_555_pedigree$Ma, FourGen_555_pedigree$Kid)))),
         Ma = as.numeric(factor(Ma,levels = unique(c(FourGen_555_pedigree$Pa, FourGen_555_pedigree$Ma, FourGen_555_pedigree$Kid)))))



simMendelGeno <- function(
  starting.ped = NA,
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
  max.num.chr = NA, # set the upper limit - the smaller the number, the higher chance of being physially linked
  ref.genome = gdropR::chinook_chromosomes,
                    max.unobs.layer = 1) {

  if (is.na(starting.ped )) return()
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


  if(!is.na(max.num.chr)) {
    ref.genome <- head(ref.genome, max.num.chr)
  }

  genome <- ref.genome %>%
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
  pedigree <- starting.ped %>%
    gdropR::add_explicit_founder_parents() %>%
    gdropR::add_pedigree_sex_column(fem_prob = 1.0)

  generation.filler <-
    starting.ped %>%
    pivot_longer(cols=Kid:Ma) %>%
    mutate(generation.update = ifelse(name == "Kid", generation, generation -1)) %>%
    rename(Kid = value) %>%
    group_by(Kid) %>%
    summarise(generation.update = generation.update[1])

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

  ped_with_gens <- left_join(pedigree, generation.filler, by ="Kid") %>%
    select(-generation) %>%
    rename(generation = generation.update)

  genos012 <- err_genos %>%
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
    mutate(indiv = case_when(
      str_detect(indiv, "_$") ~ str_c(indiv, "f"),
      TRUE ~ indiv
    ))

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


  is.missing <- rbinom(n.id*n.snp, 1, missing.unk) %>% matrix(ncol=n.snp, by=T)
  is.avail <- 1-is.missing

  final.geno.1 <- (abs(geno.1-make.err.1)*is.avail) +(-1*is.missing)
  final.geno.2 <- (abs(geno.2-make.err.2)*is.avail) +(-1*is.missing)

  geno.long <- cbind(final.geno.1, final.geno.2)[,rbind(1:n.snp, (1:n.snp)+n.snp) %>% as.vector()]


  geno.meta.info <- left_join(tibble(Kid =genos012$indiv), ped_with_gens)

  geno.input.tbl <- cbind(id = geno.meta.info$Kid,
                          obs = rep(1, length(geno.meta.info$Kid)),
                          sex = geno.meta.info$Sex,
                          gen =  geno.meta.info$generation,
                          geno.long) %>% as.tibble()

  gen.indx <- geno.meta.info$generation
  total.gen <- max(gen.indx)+1

  ## add in unrelated indiv to any generation, except to the most recent one
  sample.frac <- rep(1, total.gen)

  if (one.side.censorship == 0) {
  # process sampling frac
  if(sample.frac.gen[1] != -1) {
    sample.frac[sample.frac.gen+1] <- sampling.frac.rate
    is.avail <- rbinom(length(gen.indx),1,sample.frac[gen.indx+1])==1
    geno.input.tbl <- geno.input.tbl[is.avail, ]
  }
  } else {
    sample.frac.tbl <- geno.meta.info %>%
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

run.pedFac.multiGen <- function(input.ped = NA, n.snp = 100, n.iter = 10, folder.label = "", max.num.chr = NA, random.seed = 234) {

  snp.multigen.path <- paste0("/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/multi_case/ped2/",folder.label)

  dir.create(file.path(snp.multigen.path), recursive = TRUE, showWarnings = FALSE)

  if(is.na(input.ped)) {
  starting.ped  <- gdropR::FourGen_555_pedigree %>%
    filter(generation != 1) %>%
    mutate(Kid = as.numeric(factor(Kid,levels = unique(c(FourGen_555_pedigree$Pa, FourGen_555_pedigree$Ma, FourGen_555_pedigree$Kid)))),
           Pa = as.numeric(factor(Pa,levels = unique(c(FourGen_555_pedigree$Pa, FourGen_555_pedigree$Ma, FourGen_555_pedigree$Kid)))),
           Ma = as.numeric(factor(Ma,levels = unique(c(FourGen_555_pedigree$Pa, FourGen_555_pedigree$Ma, FourGen_555_pedigree$Kid))))) %>%
    mutate(across(.cols=Kid:Ma, .fns = as.character)) %>%
    mutate(generation = generation -1)
  } else {
    starting.ped <- input.ped
  }

  simMendelGeno(starting.ped = starting.ped,
                  n.snp = n.snp,
                  out.path =snp.multigen.path,
                  one.side.censorship = 1,
                  max.num.chr = max.num.chr,
                random.seed = random.seed)



run.multiped.summary <- list()
set.seed(random.seed)
random.arr <- ceiling(c(runif(2,0,1e6)))


write.table(random.arr, paste0(snp.multigen.path,"/","rand"),
            sep = " ",eol = " ",quote = FALSE, col.names = FALSE, row.names = FALSE)

## set strictly mono
pedfac.call <- paste0("/Users/thomasn/repo/pedigree-factor-graphs/src/pedigraph_0.187 -d ", snp.multigen.path,
                      " -n ",n.iter,
                      " -r ",snp.multigen.path,"/rand -j 1 -f 10", collapse = "")

tic("run pedFac:")
system(pedfac.call,ignore.stderr = F, ignore.stdout = F)
runtime <- toc()
write.table(runtime$toc - runtime$tic, file =paste0(snp.multigen.path,"/pedFac.time"), quote = FALSE, row.names = FALSE, col.names = FALSE)
run.multiped.summary$runtime <- runtime$toc - runtime$tic
run.multiped.summary$pedfac.cmd <- pedfac.call

return(run.multiped.summary)
}

run.sequoia.multiGen <- function(input.ped = NA, n.snp = 100, n.iter = 10, folder.label = "", max.num.chr = NA, random.seed = 234) {
  snp.multigen.path <- paste0("/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/multi_case/ped2/",folder.label)
  if(is.na(input.ped)) {
    starting.ped  <- gdropR::FourGen_555_pedigree %>%
      filter(generation != 1) %>%
      mutate(Kid = as.numeric(factor(Kid,levels = unique(c(FourGen_555_pedigree$Pa, FourGen_555_pedigree$Ma, FourGen_555_pedigree$Kid)))),
             Pa = as.numeric(factor(Pa,levels = unique(c(FourGen_555_pedigree$Pa, FourGen_555_pedigree$Ma, FourGen_555_pedigree$Kid)))),
             Ma = as.numeric(factor(Ma,levels = unique(c(FourGen_555_pedigree$Pa, FourGen_555_pedigree$Ma, FourGen_555_pedigree$Kid))))) %>%
      mutate(across(.cols=Kid:Ma, .fns = as.character)) %>%
      mutate(generation = generation -1)
  } else {
    starting.ped <- input.ped
  }


  pedFac.geno.tbl <- read.table(paste0("/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/multi_case/ped2/", folder.label,"/geno.txt"), stringsAsFactors = F) %>%
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

  set.seed(random.seed)

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

  ParOUT_Simp <- sequoia(GenoM = as.matrix(geno.matrix), LifeHistData = lifeHistory.tbl, MaxSibIter = 0,Complex = "simp", Err = geno.err)
  SeqOUT <- sequoia(GenoM = geno.matrix,
                    SeqList = ParOUT_Simp,
                    MaxSibIter = 20, Err = geno.err, Complex = "simp")

  runtime <- toc()
  write.table(runtime$toc - runtime$tic, file =paste0(snp.multigen.path,"/sequoia.time"), quote = FALSE, row.names = FALSE, col.names = FALSE)

  save(SeqOUT, file = paste0(snp.multigen.path, "/sequoia_out.rds",collapse = ""))

  seq.summary <- runtime$toc - runtime$tic
}

report.sequoia.multiGen <- function(input.ped = NA, n.snp = 100, n.iter = 10, folder.label = "", max.num.chr = NA, random.seed = 234) {
  # fill in unobs node

  snp.multigen.path <- paste0("/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/multi_case/ped2/",folder.label)

  if(is.na(input.ped)) {
    starting.ped  <- gdropR::FourGen_555_pedigree %>%
      filter(generation != 1) %>%
      mutate(Kid = as.numeric(factor(Kid,levels = unique(c(FourGen_555_pedigree$Pa, FourGen_555_pedigree$Ma, FourGen_555_pedigree$Kid)))),
             Pa = as.numeric(factor(Pa,levels = unique(c(FourGen_555_pedigree$Pa, FourGen_555_pedigree$Ma, FourGen_555_pedigree$Kid)))),
             Ma = as.numeric(factor(Ma,levels = unique(c(FourGen_555_pedigree$Pa, FourGen_555_pedigree$Ma, FourGen_555_pedigree$Kid))))) %>%
      mutate(across(.cols=Kid:Ma, .fns = as.character)) %>%
      mutate(generation = generation -1)
  } else {
    starting.ped <- input.ped
  }

  pedFac.geno.tbl <- read.table(paste0("/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/multi_case/ped2/", folder.label,"/geno.txt"), stringsAsFactors = F) %>%
    as_tibble()

  all.indiv.uncensor.ls <- as.character(pedFac.geno.tbl$V1)
  indiv.censor.ls <- setdiff(unique(c(starting.ped$Kid, starting.ped$Pa, starting.ped$Ma)), all.indiv.uncensor.ls)

  load(paste0(snp.multigen.path, "/sequoia_out.rds"))

  seq.ped.1 <- bind_rows(
  SeqOUT$Pedigree %>%
    as.tibble() %>%
    select(id, dam, sire, LLRpair) %>%
    rename(kid = id, pa = sire, ma = dam) %>%
    mutate(ma = ifelse(is.na(ma), -1, ma),
           pa = ifelse(is.na(pa), -1, pa)),
  tibble(kid = "-1",
         ma = "-1",
         pa = "-1",
         LLRpair = NA)
  )

  infer.grandparent <- RetrieveGrandparents.simple(seq.ped.1,ref = TRUE) %>%
    mutate(norm.LLRpair = ifelse(is.na(LLRpair),
                                 0,
                                 exp(LLRpair)),
           norm.LLRpair.x = ifelse(is.na(LLRpair.x),
                                 0,
                                 exp(LLRpair.x)),
           prob.x=norm.LLRpair.x/max(norm.LLRpair.x, norm.LLRpair),
           prob.y=norm.LLRpair/max(norm.LLRpair.x, norm.LLRpair)
           ) %>%
    group_by(kid) %>%
    mutate(prob = max(prob.x, prob.y)) %>%
    select(-norm.LLRpair, -norm.LLRpair.x, -prob.x, -prob.y, -LLRpair.x, -LLRpair) %>%
    ungroup()
  truth.grandparent <- RetrieveGrandparents.simple(prep.truth.tbl, TRUE)

  ## let's focus on grandkids that are observed first, then we can look at special cases

  infer.grandparent.summary <- infer.grandparent %>%
    filter(kid %in% all.indiv.uncensor.ls) %>%
    dplyr::transmute(kid = kid,
                     grandpa.pa = ifelse(!grandpa.pa %in% all.indiv.uncensor.ls,-1,grandpa.pa),
                     grandma.pa = ifelse(!grandma.pa %in% all.indiv.uncensor.ls,-1,grandma.pa),
                     grandpa.ma = ifelse(!grandpa.ma %in% all.indiv.uncensor.ls,-1,grandpa.ma),
                     grandma.ma = ifelse(!grandma.ma %in% all.indiv.uncensor.ls,-1,grandma.ma),
                     pa = ifelse(!pa %in% all.indiv.uncensor.ls,-1,pa),
                     ma = ifelse(!ma %in% all.indiv.uncensor.ls,-1,ma),
                     prob = prob) %>%
    mutate(across(.fn=as.numeric),
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
    select(-flip)


  grandkids.ls <- ThreeGen_ped %>% filter(generation != 1) %>% pull(Kid)

  truth.grandparent.summary <- truth.grandparent %>%
    dplyr::filter(kid %in% as.numeric(grandkids.ls)) %>%
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



  ped.join.truth <- dplyr::left_join(infer.grandparent.summary.1 %>%
                                     filter(kid %in% truth.grandparent.summary.1$kid) %>%
                                       mutate(across(.fn=as.numeric)),
                                     truth.grandparent.summary.1, by="kid") %>%
    mutate(across(.fn=~ifelse(is.na(.),-1,.))) %>% group_by(kid) %>%
    mutate(n.matches =
             #(grandpa.pa.x == grandpa.pa.y) +
             #(grandpa.ma.x == grandpa.ma.y) +
             (grandma.pa.x == grandma.pa.y) +
             (grandma.ma.x == grandma.ma.y) ) %>%
    select(-pa, -ma)


}




###

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

prepping.complete.truth.tbl <- function () {
all.members.ls <- unique(c(gdropR::FourGen_555_pedigree$Pa, gdropR::FourGen_555_pedigree$Ma, gdropR::FourGen_555_pedigree$Kid))
ThreeGen_ped<- gdropR::FourGen_555_pedigree %>%
  filter(generation != 1) %>%
  mutate(Kid = as.numeric(factor(Kid,levels = all.members.ls)),
         Pa = as.numeric(factor(Pa,levels = all.members.ls)),
         Ma = as.numeric(factor(Ma,levels = all.members.ls))) %>%
  mutate(across(.cols=Kid:Ma, .fns = as.character)) %>%
  mutate(generation = generation -1)

## ID' any founder in the second gen
intermed.founder.ls <- ThreeGen_ped %>%
  filter(generation != 1) %>%
  mutate(founder.pa = (!Pa %in% as.numeric(ThreeGen_ped$Kid)),
         founder.ma = (!Ma %in% as.numeric(ThreeGen_ped$Kid)) ) %>%
  filter(founder.pa || founder.ma)

# intermed.founder.ls <- final.tbl.cut.add  %>%
#   mutate(Kid = as.numeric(factor(Kid, levels = refactor.level)),
#          Pa = as.numeric(factor(Pa, levels = refactor.level)),
#          Ma = as.numeric(factor(Ma, levels = refactor.level))) %>%
#   filter(gen >=5) %>%
#   group_by(Kid) %>%
#   mutate(
#     founder.pa = (!Pa %in% as.numeric(crop.tbl$Kid)) &  (!Pa %in% as.numeric(indiv.censor.ls)),
#     founder.ma = (!Ma %in% as.numeric(crop.tbl$Kid)) &  (!Ma %in% as.numeric(indiv.censor.ls))
#   ) %>%
#   filter(founder.pa || founder.ma)


prep.truth.tbl <- bind_rows(
  intermed.founder.ls %>%
    filter(founder.pa) %>%
    mutate(Kid = as.numeric(Pa),
           Pa = -1,
           Ma = -1) %>%
    group_by(Kid) %>%
    slice_head() %>%
    select(Kid, Pa , Ma),
  intermed.founder.ls %>%
    filter(founder.ma) %>%
    mutate(Kid = as.numeric(Ma),
           Pa = -1,
           Ma = -1)%>%
    group_by(Kid) %>%
    slice_head() %>%
    select(Kid, Pa , Ma),
  tibble(Kid= -1,
         Pa = -1,
         Ma = -1),
  ThreeGen_ped %>% select(-generation) %>%
    mutate(across(.fns=as.numeric))
)%>%
  rename(pa=Pa, kid=Kid, ma=Ma)

return(prep.truth.tbl)
}

####

### exclusive for 3-gen pedigree (if more gen, need to rexamine some assumption)

reportGrandparentage.multi <- function(tag, base.path = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/multi_case/ped2/" , pedFac.folder.label = "out", truth.tbl = NA, indiv.censor.ls = NA) {

  if (is.na(truth.tbl)) {
    message("missing truth.tbl ")
    return()
  }
snp.multigen.path <- paste0(base.path, tag)

out.ped.raw.tbl <- read.table(paste0(snp.multigen.path, "/",pedFac.folder.label,"/ped.txt",collapse = ""), stringsAsFactors = FALSE) %>%
  dplyr::tbl_df() %>%
  dplyr::rename("iter"="V1", "kid"="V2", "pa"="V3", "ma"="V4")

id.ls <-readRDS(paste0(snp.multigen.path, "/id.rds",collapse = ""))$id
max.id <- max(id.ls)

out.ll.tbl <- read.table(paste0(snp.multigen.path, "/",pedFac.folder.label,"/ll.txt",collapse = ""), stringsAsFactors = FALSE) %>%
  dplyr::tbl_df() %>%
  dplyr::rename("iter"="V1", "kid"="V2", "ll"="V3")

out.ped.tbl <- out.ped.raw.tbl %>%
  mutate(order = row_number())%>%
  arrange(iter, desc(order))%>%
  group_by(iter, kid)%>%
  summarise(pa = pa[1], ma =ma[1])

## retrieve all grandparents of any observed individuals of the three recent generation
prior.txt <- read.delim(paste0(snp.multigen.path,"/prior.txt"), sep="\n")
maxID.indx <- sapply(1:nrow(prior.txt), function(l) grepl("maxID", prior.txt[l,])) %>% which

max.id <- prior.txt[maxID.indx,] %>%
  strsplit(., " ") %>%
  unlist() %>% .[2] %>%
  as.integer()

## identify any founder that held in between ped and add into the -1 category





infer.grandparent <- RetrieveGrandparents.simple(out.ped.tbl)
truth.grandparent <- RetrieveGrandparents.simple(truth.tbl, TRUE)

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

n.iter <- length(unique(infer.grandparent$iter))
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
  dplyr::summarise(n=n(), posterior = n()/n.iter) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::group_by(kid) %>%
  dplyr::top_n(1) %>%
  dplyr::sample_n(1)


min.gen <- min(truth.tbl$generation)
grandkids.ls <- truth.tbl %>% filter(min.gen < generation) %>% pull(kid)

truth.grandparent.summary <- truth.grandparent %>% ungroup() %>% mutate(across(.fns=as.numeric)) %>%
  dplyr::filter(kid %in% as.numeric(grandkids.ls)) %>%
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


ped.join.truth <- dplyr::left_join(infer.grandparent.summary.1 %>%
                                     filter(kid %in% truth.grandparent.summary.1$kid)
                                   , truth.grandparent.summary.1, by="kid") %>%
  mutate(across(.fn=~ifelse(is.na(.),-1,.))) %>% group_by(kid) %>%
  mutate(n.matches.a =
           (grandpa.pa.x == grandpa.pa.y) +
           (grandpa.ma.x == grandpa.ma.y) +
           (grandma.pa.x == grandma.pa.y) +
           (grandma.ma.x == grandma.ma.y),
         n.matches.b =
           (grandpa.ma.x == grandpa.pa.y) +
           (grandpa.pa.x == grandpa.ma.y) +
           (grandma.ma.x == grandma.pa.y) +
           (grandma.pa.x == grandma.ma.y),
         n.matches = max(n.matches.a, n.matches.b)
         ) %>%
  select(-n, -n.matches.a, -n.matches.b)

}


###
snp.array <- c(100, 200, 400, 800)
random.seed.array <-c(234, 334, 432, 543, 623)

# multiped.summary.34chr <-lapply(1:4, function(snp.indx) {
#   lapply(1:5, function (rep) {
#     run.pedFac.multiGen(n.snp = snp.array[snp.indx],
#                         random.seed = random.seed.array[rep],
#                         n.iter = 100,
#                         folder.label = paste0("34chr/snp",snp.array[snp.indx],"/",rep))
#   })
# })
#
# multiped.summary.3chr <-lapply(1:4, function(snp.indx) {
#   lapply(1:5, function (rep) {
#     run.pedFac.multiGen(n.snp = snp.array[snp.indx],
#                         random.seed = random.seed.array[rep],
#                         n.iter = 100,
#                         folder.label = paste0("3chr/snp",snp.array[snp.indx],"/",rep),
#                         max.num.chr = 3)
#   })
# })



spreadLD.gp.pedfac.result <- lapply(1:4, function(snp.indx) {
  lapply(1:5, function (rep) {
    reportGrandparentage.multi(tag = paste0("34chr/snp",snp.array[snp.indx],"/",rep))
  })
})

save(spreadLD.gp.pedfac.result, file = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/multi_case/ped2/spreadLD.gp.pedfac.result")

highLD.gp.pedfac.result <- lapply(1:4, function(snp.indx) {
  lapply(1:5, function (rep) {
    reportGrandparentage.multi(tag = paste0("3chr/snp",snp.array[snp.indx],"/",rep))
  })
})

save(highLD.gp.pedfac.result, file = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/multi_case/ped2/highLD.gp.pedfac.result")

# spreadLD.gp.sequoia.result <- lapply(1:4, function(snp.indx) {
#   lapply(1:5, function (rep) {
#     run.sequoia.multiGen(n.snp = snp.array[snp.indx],
#                          random.seed = random.seed.array[rep],
#                          n.iter = 100,
#                          folder.label = paste0("34chr/snp",snp.array[snp.indx],"/",rep))
#   })
# })
#
# highLD.gp.sequoia.result <- lapply(1:4, function(snp.indx) {
#   lapply(1:5, function (rep) {
#     run.sequoia.multiGen(n.snp = snp.array[snp.indx],
#                          random.seed = random.seed.array[rep],
#                          n.iter = 100,
#                          folder.label = paste0("3chr/snp",snp.array[snp.indx],"/",rep))
#   })
# })


spreadLD.gp.sequoia.result <- lapply(1:4, function(snp.indx) {
  lapply(1:5, function (rep) {
    report.sequoia.multiGen(n.snp = snp.array[snp.indx],
                         random.seed = random.seed.array[rep],
                         n.iter = 100,
                         folder.label = paste0("34chr/snp",snp.array[snp.indx],"/",rep))
  })
})

save(spreadLD.gp.sequoia.result, file = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/multi_case/ped2/spreadLD.gp.sequoia.result")

highLD.gp.sequoia.result <- lapply(1:4, function(snp.indx) {
  lapply(1:5, function (rep) {
    report.sequoia.multiGen(n.snp = snp.array[snp.indx],
                         random.seed = random.seed.array[rep],
                         n.iter = 100,
                         folder.label = paste0("3chr/snp",snp.array[snp.indx],"/",rep))
  })
})

save(highLD.gp.sequoia.result, file = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/multi_case/ped2/highLD.gp.sequoia.result")


###############

load("/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/multi_case/ped2/spreadLD.gp.pedFac.result")
load("/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/multi_case/ped2/highLD.gp.pedFac.result")


load("/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/multi_case/ped2/spreadLD.gp.sequoia.result")
load("/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/multi_case/ped2/highLD.gp.sequoia.result")

Longo <-
  bind_rows(
expand_grid(
  snp.indx = 1:4,
  rep.indx = 1:5
) %>%
  mutate(
    n.snp = c(100, 200, 400, 800)[snp.indx],
    chr.category = "N:34",
    pedFac = map2(.x = snp.indx, .y = rep.indx, .f = function(x, y) spreadLD.gp.pedfac.result[[x]][[y]]),
    sequoia = map2(.x = snp.indx, .y = rep.indx, .f = function(x, y) spreadLD.gp.sequoia.result[[x]][[y]]),
  ) %>%
  pivot_longer(
    cols = pedFac:sequoia,
    names_to = "software",
    values_to = "output"
  ),
expand_grid(
  snp.indx = 1:4,
  rep.indx = 1:5
) %>%
  mutate(
    n.snp = c(100, 200, 400, 800)[snp.indx],
    chr.category = "N:3",
    pedFac = map2(.x = snp.indx, .y = rep.indx, .f = function(x, y) highLD.gp.pedfac.result[[x]][[y]]),
    sequoia = map2(.x = snp.indx, .y = rep.indx, .f = function(x, y) highLD.gp.sequoia.result[[x]][[y]])
  ) %>%
  pivot_longer(
    cols = pedFac:sequoia,
    names_to = "software",
    values_to = "output"
  ))

cleanOutGrandParentTbl <- function(Tib) {
  if ("posterior" %in% colnames(Tib)) {
    Tib <- Tib %>% rename(prob = posterior)}

  Tib %>%
    arrange(desc(prob)) %>%
    group_by(kid) %>% slice_head(n=1) %>%
    ungroup() %>%
    arrange(desc(prob)) %>%
    mutate(rank = row_number()) %>%
    rename(num.correct.calls = n.matches)
}



grandparentage.stat <- Longo %>%
  group_by(snp.indx, rep.indx, chr.category, software) %>%
  mutate(output.tbl = map(output, cleanOutGrandParentTbl)) %>%
  select(-output) %>%
  unnest(cols=c(output.tbl))


prob.rank.plot.tbl <- grandparentage.stat %>%
  mutate(software = factor(software,
                           levels=c("pedFac",  "sequoia")),
         chr.category = factor(chr.category, levels = c("N:3", "N:34"), labels = c("N=3","N=34")),
         n.snp = factor(n.snp, levels=rev(c(100, 200, 400, 800)),
                        labels = rev(c("snp:\n100","snp:\n200","snp:\n400","snp:\n800"))),
         n.obs.grandma = as.character((grandma.ma.y != -1)+(grandma.pa.y != -1))
         )


colorBrewer2_scale <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ffff33')


prob.rank.err.line2.tbl <- prob.rank.plot.tbl %>%
  filter( num.correct.calls != 2) %>%
  mutate(
    n.miss.obs.grandparent = (grandma.ma.y >= 0 & grandma.ma.y != grandma.ma.x)+
      (grandma.pa.y >= 0 & grandma.pa.y != grandma.pa.x),
    x=rank, xend=rank, y= prob, yend = prob+0.1*(n.miss.obs.grandparent))
prob.rank.err.line3.tbl <- prob.rank.plot.tbl %>%
  filter( num.correct.calls != 2) %>%
  mutate(
    n.miss.unobs.grandparent = (grandma.ma.y == -1 & grandma.ma.y != grandma.ma.x)+
      (grandma.pa.y == -1 & grandma.pa.y != grandma.pa.x),
    x=rank, xend=rank, y= prob, yend = prob-0.1*(n.miss.unobs.grandparent))


grandparentage.prob.ggplot <- list()

for (i in 1:5) {

  grandparentage.prob.plot <-
    ggplot(data=prob.rank.plot.tbl %>% filter(rep.indx==i)) +
    #geom_segment(data=prob.rank.nobs.line.tbl, aes(x=x, y=y, xend=xend, yend=yend), color =
    #              "grey", alpha=0.4) +
    geom_segment(data=prob.rank.err.line2.tbl %>% filter(rep.indx==i), aes(x=x, y=y, xend=xend, yend=yend), color ="#a45eff" ,alpha=0.8) +
    geom_segment(data=prob.rank.err.line3.tbl %>% filter(rep.indx==i), aes(x=x, y=y, xend=xend, yend=yend),color="#a6d7ff" ,alpha=0.9)+
    geom_line(aes(x=rank, y=prob))+
    geom_rug(aes(x=rank, color = n.obs.grandma))+
    scale_color_manual("Number of grandmothers in the sample",values=c("red", "orange"), limits=c("2", "1"))+
    guides(colour = guide_legend(nrow = 1))+
    #geom_line(data=prob.err.cum.sum %>% filter(rep.indx==i), aes(x= rank, y= y),color="#5499C7", alpha=0.9)+
    #scale_color_manual(values = c("#5499C7", "#FF5733", "orange"), labels=c("n","d","d1"))+
    facet_nested(n.snp~chr.category + software, switch="y")+
    theme_light()+
    scale_y_continuous("Prob", position = "right", breaks=seq(0,1,0.2), minor_breaks = seq(-0.1,1.1,0.1))+
    scale_x_continuous("Individual index")+
    theme(#strip.placement = "outside",
      panel.spacing=unit(0,"lines"),
      strip.background.y = element_rect(fill= "dark grey"),
      strip.background.x = element_rect(fill= "dark grey"),
      strip.text.y.left = element_text(angle = 0,colour = "white"),
      legend.position="bottom")

  #panel.grid.minor.y = element_blank())

  legend.parentage <- ggplot()+
    geom_text(aes(x=0.08, y=0.6, label= "Grandmother that is misassigned is:"),size=3.5, hjust="left", vjust="middle")+
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

  ggsave(paste0("notes/fig/ch2_multi_grandparent_indxplot_",i,".pdf"), grandparentage.prob.plot,  device="pdf", height=5, width=9, units="in", dpi=500)


  grandparentage.prob.ggplot[[i]] <- grid.arrange(grandparentage.prob.plot, legend.parentage, ncol=1,heights=c(10.5,1),
                                             padding = unit(0, "line"))
  ggsave(paste0("notes/fig/ch2_multi_grandparent_indxplot_",i,".png"), grandparentage.prob.ggplot[[i]],  device="png", height=5, width=9, units="in", dpi=500)

}



## prepare for the ROC curve

roc.line.df <- grandparentage.stat %>%
  arrange(rank, desc(prob)) %>%
  group_by(software, chr.category, rep.indx, n.snp) %>%
  mutate(n.true.case = cumsum(num.correct.calls == 2),
         n.false.case = cumsum(num.correct.calls != 2),
         tot.true.case = sum(n.true.case),
         tot.false.case = sum(n.false.case)) %>%
  group_by(software, chr.category, rep.indx, n.snp, prob) %>%
  summarise(tpr = ifelse(tot.true.case[1]==0,0,
                         max(n.true.case)/tot.true.case[1]),
            fpr = ifelse(tot.false.case[1]==0,0,
                         max(n.false.case)/tot.false.case[1]),
            n.true.case = tot.true.case[1],
            n.false.case = tot.false.case[1])

roc.pt.fixed.1 <- roc.line.df %>%
  group_by(software, chr.category, rep.indx, n.snp) %>%
  summarise(posterior = 0,
            tpr = 1,
            fpr = 1.0001,
            n.true.case=max(n.true.case),
            n.false.case=max(n.false.case))

roc.pt.fixed.2 <- roc.line.df %>%
  group_by(software, chr.category, rep.indx, n.snp) %>%
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

grandparentage.stat.tbl <-
  grandparentage.stat %>%
  mutate(software = factor(software,
                           levels=c("pedFac",  "sequoia")),
         chr.category = factor(chr.category, levels = c("N:3", "N:34")),
         n.snp = factor(n.snp, levels=c(100, 200, 400, 800)),
         n.obs.grandma = as.character((grandma.ma.y != -1)+(grandma.pa.y != -1))
  ) %>%
  arrange(rank, desc(prob)) %>%
  group_by(software, chr.category, rep.indx, n.snp) %>%
  summarise (AUC = Calculate.AUC(prob, num.correct.calls),
             AR = sum(num.correct.calls==2)/n(),
             maOnpat = sum(grandma.pa.x == grandma.pa.y)/n(),
             maOnmat = sum(grandma.ma.x == grandma.ma.y)/n(),
             ) %>%
  pivot_longer(
    cols = AUC:maOnmat,
    names_to = "stat",
    values_to = "score"
  )


grandparentage.boxstat.tbl <- grandparentage.stat.tbl %>%
  group_by(software,chr.category,  n.snp, stat) %>%
  summarise(min = min(score),
            max = max(score),
            mean = mean(score),
            sd = sd(score))

###create xtable in respect to graph
grandparentage.prep.xtable <-  grandparentage.boxstat.tbl %>% ungroup() %>%
  filter(stat != "AUC") %>%
  mutate(
    chr.category = factor(chr.category, levels = c("N:34", "N:3"), labels=c("low", "high")),
    software = factor(software,  levels=c("pedFac", "sequoia")),
    n.snp = factor(n.snp, c(100,200,400,800)),
    min = round(min, 2),
    max = round(max, 2),
    mean = round(mean, 2),
    sd = round(sd, 2)) %>%
  pivot_wider(names_from = stat, names_glue = "{stat}_{.value}", values_from=min:sd) %>%
  arrange(chr.category, n.snp,  desc(software)) %>%
  group_by(chr.category, n.snp, software) %>%
  summarise("both" = sprintf("%.2f +/- %.2f", AR_mean, AR_sd),
            #"AR_range" = sprintf("(%.2f,%.2f)", AR_min, AR_max),
            "paternal" = sprintf("%.2f +/- %.2f", maOnpat_mean, maOnpat_sd),
            "maternal" = sprintf("%.2f +/- %.2f", maOnmat_mean, maOnmat_sd)) %>%
  pivot_wider(names_from=chr.category, names_glue = "{chr.category}_{.value}", values_from= 4:6) %>%
  relocate(c(2,1,8,6,4,7,5,3)) %>%
  rename("nSNPs"="n.snp") %>%
  arrange(software, nSNPs)



latex.out <-
  kbl( grandparentage.prep.xtable, booktabs = T,format = "latex") %>%
  column_spec(1, bold = T) %>%
  row_spec(c(1:4)-1, extra_latex_after = "\\rowcolor{gray!6}") %>%
  row_spec(0, align="c") %>%
  collapse_rows(1, latex_hline = "none", valign = "middle") %>%
  #add_header_above(c(" ", " ", "paternal" = 1, "maternal"=1, "both " = 1, "paternal" = 1, "maternal"=1,"both " = 1),align="c", bold=F) %>%
  add_header_above(c(" ", " ", "N=3" = 3, "N=34" = 3),align="c", bold=F) %>%
add_header_above(c(" ", " ", "1 - FDR (mean +/- sd)" = 6),align="c", bold=F)

latex.out.0 <- gsub('\\+/-', "$\\\\pm$", latex.out, perl = T)
latex.out.1 <- gsub('nSNPs', "\\\\textnumero \\$ \\$ SNPs", latex.out.0, perl = T)
latex.out.2 <- gsub('(AUC\\\\_)|(AR\\\\_)|(high\\\\_)|(low\\\\_)|(maOnpat\\\\_)|(maOnmat\\\\_)', "", latex.out.1, perl = T)
write(latex.out.2,  "notes/table/ch2_multi_3gen_grandparents.tex")


###

## check out how it looks the last sweep ped
out.ped.tbl.snap <- out.ped.raw.tbl %>%
  filter(iter ==max(iter)) %>%
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


starting.ped  <- gdropR::FourGen_555_pedigree %>%
  filter(generation != 1) %>%
  mutate(Kid = as.numeric(factor(Kid,levels = unique(c(FourGen_555_pedigree$Pa, FourGen_555_pedigree$Ma, FourGen_555_pedigree$Kid)))),
         Pa = as.numeric(factor(Pa,levels = unique(c(FourGen_555_pedigree$Pa, FourGen_555_pedigree$Ma, FourGen_555_pedigree$Kid)))),
         Ma = as.numeric(factor(Ma,levels = unique(c(FourGen_555_pedigree$Pa, FourGen_555_pedigree$Ma, FourGen_555_pedigree$Kid))))) %>%
  mutate(across(.cols=Kid:Ma, .fns = as.character)) %>%
  mutate(generation = generation -1)



observ.node.ls <- unique(c(starting.ped$Ma,
                           starting.ped %>% filter(generation ==2) %>% pull(Kid)))

ped2dot(starting.ped,
        ShowLabelNodes = unique(c(starting.ped$Kid, starting.ped$Pa, starting.ped$Ma)),
        pfactorNodeStyle = "invis",
        pfactorEdgeStyle = "invis",
        ObsNodes = observ.node.ls,
        outf = paste0(snp.multigen.path, "/truth"))


### making a smaller ped version


indiv.ls.incl <- ls

branch.backward <- function(gen = 5, indiv.list = c("712", "713", "714")) {

  select.ped <- gdrop.ped %>%
    filter(generation == gen, Kid %in% indiv.list)

  new.ls <- c(select.ped$Kid, select.ped$Pa, select.ped$Ma, indiv.list) %>%
    unique


  if (as.numeric(gen) != 0) {
    final.ls <- branch.backward(gen - 1, new.ls)
  } else  {
    final.ls <- new.ls
  }

  return(final.ls)
}

branch.forward <- function(gen = 1, indiv.list = c("78")) {

  select.ped <- gdrop.ped %>%
    filter(generation == gen, Pa %in% indiv.list | Ma %in% indiv.list)

  new.ls <- c(select.ped$Kid, select.ped$Pa, select.ped$Ma, indiv.list) %>%
    unique


  if (as.numeric(gen) != 5) {
    final.ls <- branch.forward(gen + 1, new.ls)
  } else  {
    final.ls <- new.ls
  }

  return(final.ls)
}


sub.indiv.ls <- branch.backward() %>% unique()
sub.gdrop.ped <- gdrop.ped %>%
  filter(Kid %in% sub.indiv.ls | Pa %in% sub.indiv.ls | Ma %in% sub.indiv.ls)


sub.indiv.ls.1 <- branch.forward(gen = 1, sub.indiv.ls) %>% unique()

sub.indiv.ls.1 <- branch.forward() %>% unique()

sub.gdrop.ped <- gdrop.ped %>%
  filter(Kid %in% sub.indiv.ls.1 | Pa %in% sub.indiv.ls.1 | Ma %in% sub.indiv.ls.1)

ped2dot(sub.gdrop.ped,
        ShowLabelNodes = unique(c(sub.gdrop.ped$Kid, sub.gdrop.ped$Pa, sub.gdrop.ped$Ma)),
        pfactorNodeStyle = "invis",
        pfactorEdgeStyle = "invis",
        ObsNodes = unique(c(sub.gdrop.ped$Ma, sub.gdrop.ped %>% filter(Sex == 2) %>% pull(Kid))),
        outf = paste0(snp.multigen.path, "/subped"))



run.pedFac.multiGen(input.ped =  sub.gdrop.ped %>% select(-Sex), n.snp = 200, n.iter = 50)


### what i learn is that the sample ped Eric provide has loops, need to revise the founder fn since
# it spits out gen not properly if i try to create a samller pedigree





starter <- gdrop.ped %>%
  filter(generation == 1, Ma %in% c("61","62", "63"))

gdrop.ped %>%
  filter(generation == 1, Ma %in% c("61","62", "63"))



