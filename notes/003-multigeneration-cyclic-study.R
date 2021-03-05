#
#
# simulation of three generation pedigree
#
# i think it will be best to use simPed.python with s.f. set as 1
# 200 founders (w/3 generations)
#
# mean of fecundity (2,2) - multiple partner and potential sibling
#
# we will be playing around with sampling fraction and nSNPs
# sampling frac (on middle generation first )
#
# and also sampling frac (on oldest gen with some missing indiv)
#
# num of SNP: 200 markers, 400 markers (MENDEL = N:36)
# pedfac mode: polygamous acyclic and cyclic
# sequoia??
# grandparentage inference

#python simPed.py -o /Users/thomasn/repo/pedigree-factor-graphs/data/ch3/poly1/truth -g 3 -ni 200 -sf 1 -s 400

true.ped.path <- "/Users/thomasn/repo/pedigree-factor-graphs/data/ch3/poly1/truth"
true.ped <- read.table("/Users/thomasn/repo/pedigree-factor-graphs/data/ch3/poly1/truth/marriage.all") %>%
  as.tibble %>%
  select(V1, V2, V3, V8) %>%
  rename(Kid = V1, Pa = V2, Ma= V3, generation = V8) %>%
  mutate(across(Kid:Ma, .fns=as.character))

all.indiv.ls <- c(true.ped$Kid, true.ped$Pa, true.ped$Ma) %>% unique
not.middle.gen.indiv.ls <- c(true.ped$Kid[true.ped$generation ==2 ],
                         true.ped$Pa[true.ped$generation ==1 ],
                         true.ped$Ma[true.ped$generation ==1 ]) %>% unique

loop.indiv.ls <- unique(c(5,9,20,24,34,39,45,47,63,64,73,74,81,84,101,115,117,118,138,142,144,146,164,168,169,174,175,176,177,194,196,270,271,272,273,203,274,205,275,209,210,276,277,214,281,216,283,217,220,221,226,228,288,230,290,232,291,233,234,294,296,239,240,241,242,243,298,245,299,300,301,248,250,303,252,305,253,306,307,308,255,256,258,310,311,313,260,261,262,263,316,265,317,318,269)) %>% as.character()

left_join(true.ped, tibble(Kid = loop.indiv.ls, in.loop = 1))

ped2dot(true.ped,
        ShowLabelNodes = NA,
        pfactorNodeStyle = "invis",
        pfactorEdgeStyle = "invis",
        RankSep = 2,
        NodeSep = 0.25,
        #ObsStyle = list(style="filled", fillcolor="#006bce"),
        ObsNodes = all.indiv.ls,
        #outf = paste0(true.ped.path, "/truth"),
        outf = paste0("notes/fig/ch3_3gen_truth"),
        highlight.unobs.edge.2 = TRUE,
        highlightNodes = unique(c(5,9,20,24,34,39,45,47,63,64,73,74,81,84,101,115,117,118,138,142,144,146,164,168,169,174,175,176,177,194,196,270,271,272,273,203,274,205,275,209,210,276,277,214,281,216,283,217,220,221,226,228,288,230,290,232,291,233,234,294,296,239,240,241,242,243,298,245,299,300,301,248,250,303,252,305,253,306,307,308,255,256,258,310,311,313,260,261,262,263,316,265,317,318,269)) %>% as.character()
        )

# caption: Three generation pedigrees of 335 individuals with varying mating patterns and 33 simple cycles (all complex loop type). Individuals (31 & 64 in the 1st and 2nd generation) and edges that are highlighted in blue are part of the loops. Around 1/3 of them are involved in the loop.


### need to run Mendel assume 34 chr, 200 (400?) markers with middle gen of sf of 0, 0.5, 1
# uses simMendelGeno from 002d doc


run.pedFac.3Gen.loop <- function(input.ped = NA,
                                 n.snp = 200,
                                 n.iter = 10,
                                 folder.label = "demo",
                                 max.num.chr = NA,
                                 random.seed = 234,
                                 pedfac.opts= "-f 10",
                                 pedfac.folder.label = "pedFac",
                                 sampling.frac.rate = 1, # treat as r.v.
                                 sample.frac.gen = -1) {

  snp.multigen.path <- paste0("/Users/thomasn/repo/pedigree-factor-graphs/data/ch3/poly1/",folder.label)

  dir.create(file.path(snp.multigen.path), recursive = TRUE, showWarnings = FALSE)

  if(is.na(input.ped)) {
    starting.ped <- read.table("/Users/thomasn/repo/pedigree-factor-graphs/data/ch3/poly1/truth/marriage.all") %>%
      as.tibble %>%
      select(V1, V2, V3, V8) %>%
      rename(Kid = V1, Pa = V2, Ma= V3, generation = V8) %>%
      mutate(across(Kid:Ma, .fns=as.character))
  } else {
    starting.ped <- input.ped
  }

  simMendelGeno(starting.ped = starting.ped,
                n.snp = n.snp,
                out.path =snp.multigen.path,
                max.num.chr = max.num.chr,
                random.seed = random.seed,
                sampling.frac.rate = sampling.frac.rate, # treat as r.v.
                sample.frac.gen = sample.frac.gen
                )

  run.multiped.summary <- list()
  set.seed(random.seed)
  random.arr <- ceiling(c(runif(2,0,1e6)))


  write.table(random.arr, paste0(snp.multigen.path,"/","rand"),
              sep = " ",eol = " ",quote = FALSE, col.names = FALSE, row.names = FALSE)

  ## set strictly mono
  pedfac.call <- paste0("/Users/thomasn/repo/pedigree-factor-graphs/src/pedigraph_0.187.1 -d ", snp.multigen.path,
                        " -n ",n.iter,
                        " -r ",snp.multigen.path,"/rand ",pedfac.opts, collapse = "")

  pedfac.out.path <- paste0(snp.multigen.path, "/",pedfac.folder.label)
  dir.create(file.path(pedfac.out.path), recursive = TRUE, showWarnings = FALSE)

  print(pedfac.call)
  tic("run pedFac:")
  system(pedfac.call,ignore.stderr = F, ignore.stdout = F)
  runtime <- toc()

  pedfac.outfile.ls <- list.files(paste0(snp.multigen.path,"/out"), include.dirs=TRUE, full.names = T)
  file.copy(pedfac.outfile.ls,
            paste0(pedfac.out.path), overwrite = T)


  write.table(runtime$toc - runtime$tic, file =paste0(pedfac.out.path,"/pedFac.time"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  run.multiped.summary$runtime <- runtime$toc - runtime$tic
  run.multiped.summary$pedfac.cmd <- pedfac.call

  return(run.multiped.summary)
}


##run sequoia


run.sequoia.3Gen.loop <- function(input.ped = NA,
                                                                   n.snp = 200,
                                                                   n.iter = 10,
                                                                   folder.label = "demo",
                                                                   max.num.chr = NA,
                                                                   random.seed = 234,
                                                                   pedfac.opts= "-f 10",
                                                                   pedfac.folder.label = "pedFac",
                                                                   sampling.frac.rate = 1, # treat as r.v.
                                                                   sample.frac.gen = -1) {

  snp.multigen.path <- paste0("/Users/thomasn/repo/pedigree-factor-graphs/data/ch3/poly1/",folder.label)

  if(is.na(input.ped)) {
    starting.ped <- read.table("/Users/thomasn/repo/pedigree-factor-graphs/data/ch3/poly1/truth/marriage.all") %>%
      as.tibble %>%
      select(V1, V2, V3, V8) %>%
      rename(Kid = V1, Pa = V2, Ma= V3, generation = V8) %>%
      mutate(across(Kid:Ma, .fns=as.character))
  } else {
    starting.ped <- input.ped
  }


  pedFac.geno.tbl <- read.table(paste0("/Users/thomasn/repo/pedigree-factor-graphs/data/ch3/poly1/", folder.label,"/geno.txt"), stringsAsFactors = F) %>%
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

  ParOUT_Simp <- sequoia(GenoM = as.matrix(geno.matrix), LifeHistData = lifeHistory.tbl, Err = geno.err, MaxSibIter = 0) #Complex = "simp"
  SeqOUT <- sequoia(GenoM = geno.matrix,
                    SeqList = ParOUT_Simp,Err = geno.err)
                    #MaxSibIter = 20, Err = geno.err, Complex = "simp")#MaxSibIter = 20,

  runtime <- toc()
  write.table(runtime$toc - runtime$tic, file =paste0(snp.multigen.path,"/sequoia.time"), quote = FALSE, row.names = FALSE, col.names = FALSE)

  save(SeqOUT, file = paste0(snp.multigen.path, "/sequoia_out.rds",collapse = ""))

  seq.summary <- runtime$toc - runtime$tic
}

report.sequoia.3Gen.loop <- function(tag, base.path = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/multi_case/ped2/" , folder.label = "out", truth.tbl = NA, indiv.censor.ls = NA) {
  # fill in unobs node

  snp.multigen.path <- paste0(base.path, tag)

  if(is.na(truth.tbl)) {
    truth.tbl <- read.table("/Users/thomasn/repo/pedigree-factor-graphs/data/ch3/poly1/truth/marriage.all") %>%
      as.tibble %>%
      select(V1, V2, V3, V8) %>%
      rename(kid = V1, pa = V2, ma= V3, generation = V8) %>%
      mutate(across(kid:ma, .fns=as.character))
  }

  pedFac.geno.tbl <- read.table(paste0(snp.multigen.path,"/geno.txt"), stringsAsFactors = F) %>%
    as_tibble()

  all.indiv.uncensor.ls <- as.character(pedFac.geno.tbl$V1)
  indiv.censor.ls <- setdiff(unique(c(truth.tbl$kid, truth.tbl$pa, truth.tbl$ma)), all.indiv.uncensor.ls)

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

  truth.grandparent <- RetrieveGrandparents.simple(truth.tbl, TRUE)

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


  min.gen <- min(truth.tbl$generation)
  grandkids.ls <- truth.tbl %>% filter(min.gen < generation) %>% pull(kid)

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
                                     truth.grandparent.summary.1 %>%
                                       mutate(across(.fn=as.numeric)), by="kid") %>%
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
    select( -n.matches.a, -n.matches.b, -pa, -ma)


}

# run.pedFac.3Gen.loop(pedfac.opts = "-f 10", pedfac.folder.label = "pedFac.v", n.iter = 25)
# demo.ac.pedFac.result <- reportGrandparentage.multi(base.path = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch3/poly1/",
#                            tag= "demo",
#                             pedFac.folder.label = "pedFac.v", truth.tbl = true.ped %>% rename(kid =Kid, pa =Pa, ma=Ma))
# demo.ac.pedFac.result %>% ungroup() %>% count(n.matches)
#
#
# run.pedFac.3Gen.loop(pedfac.opts = "-f 10 -c 1", pedfac.folder.label = "cpedFac", n.iter = 25)
# demo.c.pedFac.result <- reportGrandparentage.multi(base.path = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch3/poly1/",
#                                                     tag= "demo",
#                                                     pedFac.folder.label = "cpedFac", truth.tbl = true.ped %>% rename(kid =Kid, pa =Pa, ma=Ma))
# demo.c.pedFac.result %>% ungroup() %>% count(n.matches)


sf.array <- c(1, 0.75, 0.5, 0.25, 0)
random.seed.array <-c(234, 334, 432, 543, 623)

# multigen.ac.summary <-lapply(1:5, function(sf.indx) {
#   lapply(1:5, function (rep) {
#     run.pedFac.3Gen.loop(pedfac.opts = "-f 10",
#                          folder.label = paste0("snp100/sf",sf.array[sf.indx],"/",rep),
#                          pedfac.folder.label = "pedFac_ac",
#                          n.iter = 25,
#                          random.seed = random.seed.array[rep],
#                          sampling.frac.rate = sf.array[sf.indx], # treat as r.v.
#                          sample.frac.gen = 1,
#                          n.snp = 100)
#   })
# })

multigen.c.summary <-lapply(1:5, function(sf.indx) {
  lapply(1:5, function (rep) {
    run.pedFac.3Gen.loop(pedfac.opts = "-f 10 -c 1",
                         folder.label = paste0("snp100/sf",sf.array[sf.indx],"/",rep),
                         pedfac.folder.label = "pedFac_c",
                         n.iter = 25,
                         random.seed = random.seed.array[rep],
                         sampling.frac.rate = sf.array[sf.indx], # treat as r.v.
                         sample.frac.gen = 1,
                         n.snp = 100)
  })
})

multigen.seq.summary <-lapply(1:5, function(sf.indx) {
  lapply(1:5, function (rep) {
    run.sequoia.3Gen.loop(pedfac.opts = "-f 10 -c 1",
                         folder.label = paste0("snp100/sf",sf.array[sf.indx],"/",rep),
                         pedfac.folder.label = "pedFac_c",
                         n.iter = 25,
                         random.seed = random.seed.array[rep],
                         sampling.frac.rate = sf.array[sf.indx], # treat as r.v.
                         sample.frac.gen = 1,
                         n.snp = 100)
  })
})


# ## did the one earlier
# multigen.ac.summary <-lapply(1:5, function(sf.indx) {
#   lapply(1:5, function (rep) {
#     run.pedFac.3Gen.loop(pedfac.opts = "-f 10",
#                          folder.label = paste0("snp200/sf",sf.array[sf.indx],"/",rep),
#                          pedfac.folder.label = "pedFac_ac",
#                          n.iter = 25,
#                          random.seed = random.seed.array[rep],
#                          sampling.frac.rate = sf.array[sf.indx], # treat as r.v.
#                          sample.frac.gen = 1,
#                          n.snp = 200)
#   })
# })
#
multigen.c.summary <-lapply(5, function(sf.indx) {
  lapply(5, function (rep) {
    run.pedFac.3Gen.loop(pedfac.opts = "-f 10 -c 1",
                         folder.label = paste0("snp200/sf",sf.array[sf.indx],"/",rep),
                         pedfac.folder.label = "pedFac_c",
                         n.iter = 25,
                         random.seed = random.seed.array[rep],
                         sampling.frac.rate = sf.array[sf.indx], # treat as r.v.
                         sample.frac.gen = 1,
                         n.snp = 200)
  })
})

multigen.seq.summary <-lapply(5, function(sf.indx) {
  lapply(1:5, function (rep) {
    run.sequoia.3Gen.loop(pedfac.opts = "-f 10 -c 1",
                         folder.label = paste0("snp200/sf",sf.array[sf.indx],"/",rep),
                         pedfac.folder.label = "pedFac_c",
                         n.iter = 25,
                         random.seed = random.seed.array[rep],
                         sampling.frac.rate = sf.array[sf.indx], # treat as r.v.
                         sample.frac.gen = 1,
                         n.snp = 200)
  })
})

multigen.ac.snp200.result <- lapply(1:5, function(sf.indx) {
  lapply(1:5, function (rep) {
    reportGrandparentage.multi(base.path = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch3/poly1/",
                               tag= paste0("snp200/sf",sf.array[sf.indx],"/",rep),
                               pedFac.folder.label = "pedFac_ac", truth.tbl = true.ped %>% rename(kid =Kid, pa =Pa, ma=Ma))
  })
})
save(multigen.ac.snp200.result, file = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch3/poly1/multigen.ac.snp200.result")


multigen.c.snp200.result <- lapply(1:5, function(sf.indx) {
  lapply(1:5, function (rep) {
    reportGrandparentage.multi(base.path = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch3/poly1/",
                               tag= paste0("snp200/sf",sf.array[sf.indx],"/",rep),
                               pedFac.folder.label = "pedFac_c", truth.tbl = true.ped %>% rename(kid =Kid, pa =Pa, ma=Ma))
  })
})

save(multigen.c.snp200.result, file = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch3/poly1/multigen.c.snp200.result")

multigen.seq.snp200.result <- lapply(1:5, function(sf.indx) {
  lapply(1:5, function (rep) {
    report.sequoia.3Gen.loop(base.path = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch3/poly1/",
                               tag= paste0("snp200/sf",sf.array[sf.indx],"/",rep))
  })
})

save(multigen.seq.snp200.result, file = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch3/poly1/multigen.seq.snp200.result")


######100 SNPs

multigen.ac.snp100.result <- lapply(1:5, function(sf.indx) {
  lapply(1:5, function (rep) {
    reportGrandparentage.multi(base.path = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch3/poly1/",
                               tag= paste0("snp100/sf",sf.array[sf.indx],"/",rep),
                               pedFac.folder.label = "pedFac_ac", truth.tbl = true.ped %>% rename(kid =Kid, pa =Pa, ma=Ma))
  })
})
save(multigen.ac.snp100.result, file = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch3/poly1/multigen.ac.snp100.result")


multigen.c.snp100.result <- lapply(1:5, function(sf.indx) {
  lapply(1:5, function (rep) {
    reportGrandparentage.multi(base.path = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch3/poly1/",
                               tag= paste0("snp100/sf",sf.array[sf.indx],"/",rep),
                               pedFac.folder.label = "pedFac_c", truth.tbl = true.ped %>% rename(kid =Kid, pa =Pa, ma=Ma))
  })
})

save(multigen.c.snp100.result, file = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch3/poly1/multigen.c.snp100.result")

multigen.seq.snp100.result <- lapply(1:5, function(sf.indx) {
  lapply(1:5, function (rep) {
    report.sequoia.3Gen.loop(base.path = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch3/poly1/",
                             tag= paste0("snp100/sf",sf.array[sf.indx],"/",rep))
  })
})

save(multigen.seq.snp100.result, file = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch3/poly1/multigen.seq.snp100.result")



Longo <-
  expand_grid(
    sf.indx = 1:5,
    rep.indx = 1:5
  ) %>%
  mutate(
    sf = c(1, 0.75, 0.5, 0.25, 0)[sf.indx],
    pedFac.ac = map2(.x = sf.indx, .y = rep.indx, .f = function(x, y) multigen.ac.snp100.result[[x]][[y]]),
    pedFac.c = map2(.x = sf.indx, .y = rep.indx, .f = function(x, y) multigen.c.snp100.result[[x]][[y]]),
    sequoia = map2(.x = sf.indx, .y = rep.indx, .f = function(x, y) multigen.seq.snp100.result[[x]][[y]]),
  ) %>%
  pivot_longer(
    cols = pedFac.ac:sequoia,
    names_to = "software",
    values_to = "output"
  )




Longo <-
    expand_grid(
      sf.indx = 1:5,
      rep.indx = 1:5
    ) %>%
      mutate(
        sf = c(1, 0.75, 0.5, 0.25, 0)[sf.indx],
        pedFac.ac = map2(.x = sf.indx, .y = rep.indx, .f = function(x, y) multigen.ac.result[[x]][[y]]),
        pedFac.c = map2(.x = sf.indx, .y = rep.indx, .f = function(x, y) multigen.c.snp200.result[[x]][[y]]),
        sequoia = map2(.x = sf.indx, .y = rep.indx, .f = function(x, y) multigen.seq.snp200.result[[x]][[y]]),
      ) %>%
      pivot_longer(
        cols = pedFac.ac:sequoia,
        names_to = "software",
        values_to = "output"
      )

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
  group_by(sf.indx, rep.indx, software) %>%
  mutate(output.tbl = map(output, cleanOutGrandParentTbl)) %>%
  select(-output) %>%
  unnest(cols=c(output.tbl))


prob.rank.plot.tbl <- grandparentage.stat %>%
  ungroup() %>%
  mutate(software = factor(software,
                           levels=c("pedFac.c","pedFac.ac",  "sequoia"), labels=c("pedFac (cyclic sampler)", "pedFac  (acyclic sampler)", "sequoia")),
         sf = factor(sf, levels=c(1, 0.75, 0.5, 0.25, 0),
                        labels = c("sf:\n1","sf:\n0.75", "sf:\n0.5","sf:\n0.25", "sf:\n0"))) %>%
  group_by(sf, software, rep.indx, kid) %>%
  mutate(involves.in.loops = sum((grandpa.pa.y %in% as.numeric(loop.indiv.ls))+
           (grandpa.ma.y %in% as.numeric(loop.indiv.ls))+
           (grandma.pa.y %in% as.numeric(loop.indiv.ls))+
           (grandma.ma.y %in% as.numeric(loop.indiv.ls)))>0
         ) %>% ungroup()


#colorBrewer2_scale <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ffff33')

prob.rank.err.line2.tbl <- prob.rank.plot.tbl %>%
  filter( num.correct.calls != 4) %>%
  mutate(
    x=rank, xend=rank, y= prob, yend = prob+0.1*(4-num.correct.calls))
#
# prob.rank.err.line3.tbl <- prob.rank.plot.tbl %>%
#   filter( num.correct.calls != 2) %>%
#   mutate(
#     n.miss.unobs.grandparent = (grandma.ma.y == -1 & grandma.ma.y != grandma.ma.x)+
#       (grandma.pa.y == -1 & grandma.pa.y != grandma.pa.x),
#     x=rank, xend=rank, y= prob, yend = prob-0.1*(n.miss.unobs.grandparent))

grandparentage.prob.ggplot <- list()

for (i in 1:5) {

  grandparentage.prob.plot <-
    ggplot(data=prob.rank.plot.tbl %>% filter(rep.indx==i)) +
    #geom_segment(data=prob.rank.nobs.line.tbl, aes(x=x, y=y, xend=xend, yend=yend), color =
    #              "grey", alpha=0.4) +
    geom_segment(data=prob.rank.err.line2.tbl %>% filter(rep.indx==i, num.correct.calls <4), aes(x=x, y=y, xend=xend, yend=yend, color = as.factor(4-num.correct.calls)) ,alpha=1) +
    #geom_segment(data=prob.rank.err.line3.tbl %>% filter(rep.indx==i), aes(x=x, y=y, xend=xend, yend=yend),color="#a6d7ff" ,alpha=0.9)+
    geom_line(aes(x=rank, y=prob))+
    #geom_rug(aes(x=rank, color = involves.in.loops))+
    #scale_color_manual("grandparents involved in loop",values=c("red", "blue"), limits=c(TRUE, FALSE))+
    guides(colour = guide_legend(nrow = 1))+
    scale_color_manual("No. incorrect assigned grandparents",values = rev(c('#ff7512','#fdb863','#b2abd2','#5e3c99')))+
    facet_grid(sf~software,  switch="y")+
    theme_light()+
    scale_y_continuous("Prob", position = "right", breaks=seq(0,1,0.2), minor_breaks = seq(-0.1,1.1,0.1))+
    scale_x_continuous("Individual index")+
    theme(
      strip.placement = "outside",
      panel.spacing=unit(0,"lines"),
      strip.background.y = element_rect(fill= "dark grey"),
      strip.background.x = element_rect(fill= "dark grey"),
      strip.text.y.left = element_text(angle = 0,colour = "white",size = 10),
      strip.text.x.top = element_text(angle = 0,colour = "white",size = 10),
      legend.position="bottom")

  #panel.grid.minor.y = element_blank())

  # legend.parentage <- ggplot()+
  #   geom_text(aes(x=0.08, y=0.6, label= "Grandmother that is misassigned is:"),size=3.5, hjust="left", vjust="middle")+
  #   geom_segment(aes(x=0.55, y=0.0, xend=0.55, yend=1), color ="#a45eff" ,alpha=0.8, size=0.9)+
  #   geom_text(aes(x=0.58, y=0.6, label= "observed"),size=3.5, hjust="left", vjust="middle")+
  #   geom_segment(aes(x=0.75, y=0.0, xend=0.75, yend=1), color ="#a6d7ff" ,alpha=0.9, size=0.9)+
  #   geom_text(aes(x=0.78, y=0.6, label= "unobserved"),size=3.5, hjust="left", vjust="middle")+
  #   scale_y_continuous("", position = "right", breaks=FALSE,limits = c(0,1),minor_breaks = NULL,labels = NULL)+
  #   scale_x_continuous("", limits=c(0,1))+
  #   theme(#strip.placement = "outside",
  #     strip.text.y.left = element_text(angle = 0,colour = "white"),
  #     panel.grid.minor.y = element_blank(),
  #     axis.line=element_blank(),
  #     axis.text.x=element_blank(),
  #     axis.text.y=element_blank(),
  #     axis.ticks=element_blank(),
  #     axis.title.x=element_blank(),
  #     axis.title.y=element_blank(),
  #     legend.position="none",
  #     panel.background=element_blank(),
  #     panel.border=element_blank(),
  #     panel.grid.major=element_blank(),
  #     panel.grid.minor=element_blank(),
  #     plot.background=element_blank())

  ggsave(paste0("notes/fig/ch3_3gen_100_grandparent_indxplot_",i,".pdf"), grandparentage.prob.plot,  device="pdf", height=5, width=9, units="in", dpi=500)

  ggsave(paste0("notes/fig/ch3_3gen_100_grandparent_indxplot_",i,".png"), grandparentage.prob.plot,  device="png", height=5, width=9, units="in", dpi=500)

}

roc.line.df <- grandparentage.stat %>%
  arrange(rank, desc(prob)) %>%
  group_by(software, rep.indx, sf) %>%
  mutate(n.true.case = cumsum(num.correct.calls == 4),
         n.false.case = cumsum(num.correct.calls != 4),
         tot.true.case = sum(n.true.case),
         tot.false.case = sum(n.false.case)) %>%
  group_by(software, rep.indx, sf, prob) %>%
  summarise(tpr = ifelse(tot.true.case[1]==0,0,
                         max(n.true.case)/tot.true.case[1]),
            fpr = ifelse(tot.false.case[1]==0,0,
                         max(n.false.case)/tot.false.case[1]),
            n.true.case = tot.true.case[1],
            n.false.case = tot.false.case[1])

roc.pt.fixed.1 <- roc.line.df %>%
  group_by(software, rep.indx, sf) %>%
  summarise(posterior = 0,
            tpr = 1,
            fpr = 1.0001,
            n.true.case=max(n.true.case),
            n.false.case=max(n.false.case))

roc.pt.fixed.2 <- roc.line.df %>%
  group_by(software, rep.indx, sf) %>%
  summarise(posterior = 1,
            tpr = -0,
            fpr = -0.00001,
            n.true.case=max(n.true.case[posterior==1]),
            n.false.case=max(n.false.case[posterior==1]))

roc.compiled.df <- bind_rows(roc.pt.fixed.2, roc.line.df, roc.pt.fixed.1)

Calculate.AUC <- function(prob, categ){
  match.score <- prob[categ==4] # num.correct.calls
  mismatch.score <- prob[categ!=4]
  if (length(match.score)==0) return(0)
  if (length(mismatch.score)==0) return(1)
  o <- outer(match.score, mismatch.score, "-")
  auc <- mean((o>0) + .5*(o==0))
}

grandparentage.stat.tbl <-
  grandparentage.stat %>%
  mutate(software = factor(software,
                           levels=c("pedFac.c","pedFac.ac",  "sequoia"), labels=c("pedac", "pedc", "sequoia")),
         sf = factor(sf, levels=c(1, 0.75, 0.5, 0.25, 0))
  ) %>%
  arrange(rank, desc(prob)) %>%
  group_by(software, rep.indx, sf) %>%
  summarise (AUC = Calculate.AUC(prob, num.correct.calls),
             AR = sum(num.correct.calls)/(4*n())) %>% #more of an call rate
  pivot_longer(
    cols = AUC:AR,
    names_to = "stat",
    values_to = "score"
  )


grandparentage.boxstat.tbl <- grandparentage.stat.tbl %>%
  group_by(software, sf, stat) %>%
  summarise(min = min(score),
            max = max(score),
            mean = mean(score),
            sd = sd(score))

###create xtable in respect to graph
gp.prep.xtable <- grandparentage.boxstat.tbl %>% ungroup() %>%
  filter(stat != "AUC") %>%
  mutate(
    min = round(min, 2),
    max = round(max, 2),
    mean = round(mean, 2),
    sd = round(sd, 2)) %>%
  pivot_wider(names_from = stat, names_glue = "{stat}_{.value}", values_from=min:sd) %>%
  arrange(sf, desc(software)) %>%
  group_by(sf, software) %>%
  summarise("AR_mean +/- sd" = sprintf("%.2f +/- %.2f", AR_mean, AR_sd),
            "AR_range" = sprintf("(%.2f,%.2f)", AR_min, AR_max),
            #"AUC_mean w/ sd" = sprintf("%.2f +/- %.2f", AUC_mean, AUC_sd),
            #"AUC_range" = sprintf("(%.2f,%.2f)", AUC_min, AUC_max)
            ) %>%
  pivot_wider(names_from=software, names_glue = "{software}_{.value}", values_from= 3:4) %>%
  relocate(c(1,2,5,3,6,4,7)) %>%
  rename(fraction = sf)

latex.out <- kbl(gp.prep.xtable, booktabs = T, "latex") %>%
  column_spec(1, bold = T, width = "3.5em", latex_valign= "m", latex_column_spec = "c") %>%
  row_spec(c(1,3,5)-1, extra_latex_after = "\\rowcolor{gray!6}") %>%
  row_spec(0, align="c") %>%
  collapse_rows(1, latex_hline = "none") %>%
  add_header_above(c(" ", "pedFac (cyclic sampler)" = 2, "pedFac  (acyclic sampler)" = 2, "sequoia" = 2),align="c", bold=F) %>%
  add_header_above(c(" ", "fraction of correctly assigned grandparents" = 6),align="c", bold=T)

latex.out.1 <- gsub('\\+/-', "$\\\\pm$", latex.out, perl = T)
latex.out.2 <- gsub('(AUC\\\\_)|(AR\\\\_)|(pedc\\\\_)|(pedac\\\\_)|(sequoia\\\\_)', "", latex.out.1, perl = T)
write(latex.out.2,  "notes/table/ch3_3gen_gp_100.tex")




### work with the mini ped sample
