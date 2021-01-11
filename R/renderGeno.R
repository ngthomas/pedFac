#'check to see whether all the entries of the input genotype file are sound
#'
#'The function \code{checkGeno} examines whether the input genotype file is accepteable
#'
#'@param geno.path string. Path to the input genotype file. Required.
#'
#'About the genotype file:
#'The genotype file is a space separate file that contains individual's genotype and meta information.
#'Each row is an individual entry with its associate genotype information, in the order as follows:
#'unique indiv ID | is the indiv observed? | sex of individual | birth year | genotype(s) information.
#'
#'geno.path = "/Users/thomasn/repo/pedfac/example/case_0/genotype.txt"
checkGeno <- function(geno.path) {

  if(!file.exists(geno.path)) stop("the genotype file 'geno.path' provided - ", geno.path, " does not exist")
  geno.tbl <- read.table(geno.path, stringsAsFactors = FALSE) %>%
    dplyr::tbl_df()

  #check whether first of all whether any of the ids have duplicated entries
  n.id <- unique(geno.tbl$V1) %>% length()
  if (n.id != nrow(geno.tbl)) stop("Duplicate ID entries are found in ",geno.path)

  # column 2: only allow 0 or 1
  geno.checker <- geno.tbl %>%
    dplyr::group_by(V1) %>%
    dplyr::summarise(obs.check = !V2 %in% c(0,1),
                     sex.check = !V3 %in% c(0,1,2),
                     birth.year.check = is.character(V4)) %>%
    dplyr::ungroup() %>%
    dplyr::summarise(obs.pass = sum(obs.check)==0,
                     sex.pass = sum(sex.check)==0,
                     birth.pass = sum(birth.year.check)==0)

  if(!geno.checker$obs.pass) stop("The second column of the genotype files - whether indiv is observed - can accept value 0 or 1")
  if(!geno.checker$sex.pass) stop("The third column of the genotype files - sex of indiv - can accept value 0, 1, or 2")
  if(!geno.checker$birth.pass) stop("The fourth column of the genotype files - birth years - does not allow string character")


  colnames(geno.tbl) <- 1:ncol(geno.tbl)
  colnames(geno.tbl)[1] <- "id"

  genotype.chk <- dplyr::left_join(geno.tbl %>%
     tidyr::gather("col","snp", seq(5,ncol(geno.tbl),2)) %>%
       dplyr::select(id, col, snp) %>%
       dplyr::group_by(id, col) %>%
       dplyr::mutate(snp.pos = floor((as.numeric(col)-5)/2)),
    geno.tbl %>%
      tidyr::gather("col","snp", seq(6,ncol(geno.tbl),2)) %>%
      dplyr::select(id, col, snp) %>%
      dplyr::group_by(id,col) %>%
      dplyr::mutate(snp.pos = floor((as.numeric(col)-5)/2)),
    by=c("id", "snp.pos" )) %>%
    dplyr::group_by(id, snp.pos) %>%
    dplyr::mutate(
      snp.x = as.character(snp.x),
      snp.y = as.character(snp.y),
      comma.entries.1 = length(unlist(strsplit(snp.x, ",",fixed = T))) ==2,
                  comma.entries.2 = length(unlist(strsplit(snp.y, ",",fixed = T))) ==2,
                  not.match.up = ((snp.x %in% c("N",-1))+(snp.y %in% c("N",-1))) ==1) %>%
    dplyr::ungroup() %>%
    dplyr::summarise(missing.entries = sum(comma.entries.1)+sum(comma.entries.2),
                     not.match.NA = sum(not.match.up))

  if(genotype.chk$missing.entries>0) stop("Detect missing entries in genotype likelihood format")
  if(genotype.chk$not.match.NA>0) stop("Genotype are reported inconsistently. Under missing genotype, both of allelic columns must either be set as N or -1.")
}

writeIntermedGeno <- function(param, preservesID=F, preservesGeno = F) {
  #param$geno.path = "/Users/thomasn/repo/pedfac/example/case_0/genotype.txt"
  #param$geno.path ="/var/folders/6y/f_3ctjk96fb17pkn2ct2x3w00000gn/T//RtmpBgiZIO/ingeno.txt"

  geno.tbl <- read.table(param$geno.path, stringsAsFactors = FALSE) %>%
    dplyr::tbl_df()

  geno.ls <- list()
  geno.ls$max.id <- max(geno.tbl[,1])
  geno.ls$n.indiv <- nrow(geno.tbl)
  geno.ls$n.obs.indiv <- sum(geno.tbl[,2]==1)
  geno.ls$min.yr <- min(geno.tbl[,4])
  geno.ls$max.yr <- max(geno.tbl[,4])

  geno.ls$yr <- geno.tbl[,4]

  geno.ls$gen <- floor((geno.ls$max.yr-geno.tbl[,4])/param$min.age) %>% unlist()

  geno.ls$n.gen <- max(geno.ls$gen, 2)+param$max.gen
  geno.ls$is.founder <- rep(1, geno.ls$n.indiv)

  ##this opt may change if and only we feed the pedfac with initial pedigree config
  #if(param$max.gen == 0) geno.ls$is.founder <- 1*(geno.ls$gen ==max(geno.ls$gen))

  if(param$observe.frac == -1) {
    geno.ls$observe.frac <- rep(-1, geno.ls$n.gen)
  } else {
    geno.ls$observe.frac <- rep(0, geno.ls$n.gen)
    geno.ls$observe.frac[(geno.ls$gen)+1] <- param$observe.frac
  }

  geno.ls$n.SNP <- (ncol(geno.tbl)-4)/2

  colnames(geno.tbl) <- 1:ncol(geno.tbl)
  colnames(geno.tbl)[1] <- "id"

  geno.left <- geno.tbl %>%
    tidyr::gather("col","snp", seq(5,ncol(geno.tbl),2)) %>%
    dplyr::select(id, col, snp) %>%
    dplyr::group_by(id, col) %>%
    dplyr::mutate(snp.pos = floor((as.numeric(col)-5)/2)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-col)

  # normalize the comma separate geno likelihood
  geno.right <- geno.tbl %>%
    tidyr::gather("col","snp", seq(6,ncol(geno.tbl),2)) %>%
    dplyr::select(id, col, snp) %>%
    dplyr::group_by(id,col) %>%
    dplyr::mutate(snp.pos = floor((as.numeric(col)-5)/2),
                  snp = ifelse(grepl(",", snp, fixed=T),
                               paste0(round(as.numeric(unlist(strsplit(snp,",",fixed = T)))/
                                 sum(as.numeric(unlist(strsplit(snp,",",fixed = T)))),4),
                                 collapse = ","),
                               snp)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-col)

  geno.ct.1 <- dplyr::bind_rows(geno.left, geno.right) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!snp %in% c(-1, "N"), !grepl(",", snp , fixed=T)) %>%
    dplyr::group_by(snp.pos,snp) %>%
    dplyr::summarise(ct = n())


  # in working with genotype ll
  if(sum(grepl(",", geno.right$snp, fixed=T)) > 0) {
  geno.ct.2 <- geno.right %>%
    dplyr::filter(!snp %in% c(-1, "N"), grepl(",", snp , fixed=T)) %>%
    dplyr::group_by(snp.pos, id) %>%
    dplyr::summarise(`0`=sum(as.numeric(unlist(strsplit(snp,",",fixed = T)))*c(1,0.5,0)),
                     `1`=1-`0`) %>%
    tidyr::gather("snp","ct", 3:4) %>%
    dplyr::group_by(snp.pos,snp) %>%
    dplyr::summarise(ct = sum(ct))

    geno.ct.1 <- dplyr::bind_rows(geno.ct.1, geno.ct.2)

  }

  geno.freq.tbl <- geno.ct.1 %>%
    dplyr::group_by(snp.pos,snp) %>%
    dplyr::summarise(ct = sum(ct)) %>%
    dplyr::group_by(snp.pos) %>%
    dplyr::mutate(tot = sum(ct)) %>%
    dplyr::group_by(snp.pos, snp) %>%
    dplyr::summarise(freq = ct/tot) %>%
    dplyr::arrange(snp.pos, desc(freq)) %>%
    dplyr::group_by(snp.pos) %>%
    dplyr::mutate(rank=dplyr::row_number(),
                  base = freq[1])

  # by organizing multiallelic marker - method 1 : base, jam, crust - we want to group a subset of markers as the "ancestral/dominant" allele

  # caution: this removes any snp pos with all in
  geno.cum.tbl <- geno.freq.tbl %>%
    dplyr::filter(rank > 1) %>%
    dplyr::arrange(snp.pos, freq) %>%
    dplyr::group_by(snp.pos) %>%
    dplyr::mutate(cfreq = cumsum(freq)+base[1])

  zero.allele.tbl <- dplyr::bind_rows(
    geno.freq.tbl %>% dplyr::filter(rank == 1),
    geno.cum.tbl %>% dplyr::filter(base < 0.5, cfreq<0.5) %>% dplyr::select(-cfreq),
    geno.cum.tbl %>% dplyr::filter(base < 0.5, cfreq>=0.5) %>%
      dplyr::group_by(snp.pos) %>% dplyr::top_n(1, cfreq) %>% dplyr::select(-cfreq))

  zero.freq <- zero.allele.tbl %>%
    dplyr::group_by(snp.pos) %>%
    dplyr::summarise(freq=sum(freq))

  if (preservesGeno) {
    afreq <- paste0(round(param$af,4), collapse = " ")
    } else {
  afreq <- paste0(round(1-(zero.freq %>%
                             dplyr::arrange(snp.pos) %>%
                             dplyr::select(freq)),3) %>% unlist(),
                  collapse = " ")
    }

  zero.rep.key <- zero.allele.tbl %>% dplyr::select(snp.pos, snp) %>% dplyr::mutate(allele=0)

  geno.left.mod <- dplyr::left_join(geno.left %>% dplyr::filter(!grepl(",", snp , fixed=T)), zero.rep.key, by= c("snp.pos", "snp")) %>%
    dplyr::group_by(id, snp.pos) %>%
    dplyr::mutate(allele = ifelse(is.na(allele),
                                  ifelse(snp %in% c("N", -1, "0,1,2"), -1, 1),
                                  0))

  geno.right.mod <- dplyr::left_join(geno.right %>% dplyr::filter(!grepl(",", snp , fixed=T)), zero.rep.key, by= c("snp.pos", "snp")) %>%
    dplyr::group_by(id, snp.pos) %>%
    dplyr::mutate(allele = ifelse(is.na(allele),
                                  ifelse(snp %in% c("N", -1, "0,1,2"), -1, 1),
                                  0))


  if (nrow(geno.right %>% dplyr::filter(grepl(",", snp , fixed=T))) != 0 ) {
  geno.spread.tbl <- dplyr::bind_rows(
    dplyr::left_join(geno.right %>% dplyr::filter(grepl(",", snp , fixed=T)),
                     zero.rep.key %>% dplyr::rename("anc"= "snp"),
                     by= c("snp.pos")) %>%
      dplyr::group_by(id, snp.pos) %>%
      dplyr::summarise(geno.cl = ifelse(as.numeric(anc)==allele,
                                        snp,
                                        paste0(rev(unlist(strsplit(snp,",",fixed = T))), collapse = ","))),
    dplyr::left_join(geno.left.mod, geno.right.mod, by=c("id", "snp.pos" )) %>%
      dplyr::group_by(id, snp.pos) %>%
      dplyr::mutate(geno = allele.x + allele.y,
                    geno.cl = ifelse(geno==-2,"3",as.character(geno))) %>%
      dplyr::select(id, snp.pos, geno.cl)) %>%
    tidyr::spread(snp.pos, geno.cl)
  } else {
    geno.spread.tbl <- dplyr::left_join(geno.left.mod, geno.right.mod, by=c("id", "snp.pos" )) %>%
        dplyr::group_by(id, snp.pos) %>%
        dplyr::mutate(geno = allele.x + allele.y,
                      geno.cl = ifelse(geno==-2,"3",as.character(geno))) %>%
        dplyr::select(id, snp.pos, geno.cl) %>%
      tidyr::spread(snp.pos, geno.cl)
  }

  if (preservesGeno) {
    geno.spread.tbl <- dplyr::left_join(geno.left, geno.right, by=c("id", "snp.pos" )) %>%
      dplyr::group_by(id, snp.pos) %>%
      dplyr::mutate(geno = snp.x + snp.y,
                    geno.cl = ifelse(geno==-2,"3",as.character(geno))) %>%
      dplyr::select(id, snp.pos, geno.cl) %>%
      tidyr::spread(snp.pos, geno.cl)
  }
  # save RDS for id.factor = id and geno.allele class rep
  # format and write intermed geno file

  geno.prefix <- dplyr::bind_cols(id=geno.tbl$id,
                                  is.obs = geno.tbl$`2`,
                                  sex = geno.tbl$`3`,
                                  yr = geno.ls$yr,
                                  is.founder = geno.ls$is.founder)

  geno.join.tbl <- dplyr::left_join(geno.prefix, geno.spread.tbl, by="id")

  if (preservesID) {

    saveRDS(tibble(id.num=1:geno.ls$max.id, id=1:geno.ls$max.id),
            paste0(param$output.path,"/id.rds"))
    } else {

      saveRDS(geno.join.tbl %>% dplyr::select(id) %>% dplyr::ungroup() %>% dplyr::mutate(id.num = dplyr::row_number()),
              paste0(param$output.path,"/id.rds"))

      geno.join.tbl <- geno.join.tbl %>%
        dplyr::ungroup() %>%
        dplyr::mutate(id = dplyr::row_number())

  }
  saveRDS(zero.rep.key, paste0(param$output.path,"/zero_rep.rds"))

  write.table(geno.join.tbl %>% arrange(id), paste0(param$output.path,"/geno.txt"),
              sep = " ",eol = "\n",quote = FALSE, col.names = FALSE, row.names = FALSE)

  write(paste0(c("nIndiv ", geno.ls$n.indiv, "\n",
                 "nObsIndiv ", geno.ls$n.obs.indiv, "\n",
                 "nGen ", geno.ls$n.gen, "\n",
                 "nSNP ", geno.ls$n.SNP, "\n",
                 "nMar 0\n",# param$n.marr,"\n",
                 "maxID ", geno.ls$max.id + 1, "\n",
                 "aFreq ", afreq, "\n",
                 "epsilon ", paste0(rep(param$geno.err, geno.ls$n.SNP), collapse=" "), "\n",
                 "obsFrac ", paste0(round(geno.ls$observe.frac,3), collapse=" "), "\n",
                 "minAge 1\n", # need to edit this
                 "maxAge 1\n", # need to edit this
                 #"semelGap \n",
                 "recentYr ", geno.ls$max.yr, "\n",
                 #"maxMarrGap ", ceiling(param$max.age/param$min.age), "\n",
                 "maxUnobsLayer ", param$max.unobs ),collapse = ""),
        paste0(param$output.path,"/prior.txt")
        )
}

#' simulate genotype entry based on mating table (kid, pa, ma)
#'
#' The function \code{simGeno} simulates genotype of individuals given a marriage table, with the working assumption that the parents of the offspring are from the same generation
#'@param mating.path string. Path to a 3-column mating file. Required if mating.tbl is missing
#'@param mating.tbl data frame. Data frame of 3-column mating (kid, pa ,ma). Required if mating.path is missing.
#'@param n.snp positive integer. number of SNPs. 10 default. Optional
#'@param alpha.ad positive numeric. parameter of a beta distribution - model for allelic density map. 10 default. Optional
#'@param beta.ad positive numeric. parameter of a beta distribution - model for allelic density map. 20 default. Optional
#'@param geno.err positive value between 0 and 1. the rate of observing error in genotype. 0.02 default. Optional
#'@param missing.unk positive value between 0 and 1. Dictates the rate of genotype of a SNP site is missing. 0.005  default. Optional
#'@param skip.gen.ls numeric value or or a list of generations of the individuals being removed. Note that the generation level we use increments by value of 1, with 0 being the founder generation. NA default. Optional
#'@param random.seed positive numeric. random seed set for genotype sampling model. 46 default. Optional
#'@param out.path string. Path to hold the geno info
#'
#' library(tidyverse)
#' mating.path <-"/Users/thomasn/repo/pedigree-factor-graphs/data/loop_test1/marriage.txt"
#' simGeno(mating.path)
#'
simGeno <- function(mating.path = "",
                    mating.tbl ="",
                    n.snp = 10,
                    alpha.ad = 10, beta.ad = 10,
                    geno.err = 0.02,
                    missing.unk = 0.005,
                    skip.gen.ls = NA,
                    random.seed = 46,
                    out.path = tempdir()) {
  set.seed(random.seed)

  if(file.exists(mating.path)) {
    # reading matable tbl
    mating.tbl <- read.table(mating.path, stringsAsFactors = FALSE) %>%
      dplyr::tbl_df()

    colnames(marr.tbl) <- c("kid", "pa", "ma")
  }

  if (is.null(mating.tbl)) {
    stop("the mating file 'mating.path' provided - ", geno.path, " does not exist")
  }

  if (ncol(mating.tbl) != 3) {
    stop("the mating table is not a 3-column table. ")
  }

  mating.factor <- factor(unlist(mating.tbl))
  mating.lvl <- levels(mating.factor)

  # set up param (number id)
  n.id <- factor(unlist(mating.tbl)) %>% levels %>% length()
  cc.indx <- 1:n.id # associate connecting components
  gen.indx <- rep(0, n.id) # assign generation level

  # sampling background alllelic frequency
  af <- rbeta(n.snp, alpha.ad, beta.ad)
  maf <- pmin(1-af,af)

  # create first copy and ..
  geno.1 <- rbinom(n.id*n.snp, 1,maf) %>% matrix(ncol=n.snp, by=T)
  # second copy
  geno.2 <- rbinom(n.id*n.snp, 1,maf) %>% matrix(ncol=n.snp, by=T)
  # indiv sex (default, update during part 2)
  indiv.sex <- rbinom(n.id, 1, 0.5)+1

  # 1st step: estimate generation of each indiv (assume everyone is gen 0, until
  # proven otherwise)
  factor.tbl <- matrix(as.numeric(mating.factor), ncol=3)
  for(i in 1:nrow(factor.tbl)) {
    x<-factor.tbl[i,]
    if(cc.indx[x[2]] != cc.indx[x[3]]) {
      diff.rank <- gen.indx[x[2]] - gen.indx[x[3]]
      if(diff.rank != 0) {
        if(diff.rank > 0) {
          gen.indx = gen.indx + diff.rank*(cc.indx==cc.indx[x[3]])
        } else {
          gen.indx = gen.indx + -1*diff.rank*(cc.indx==cc.indx[x[2]])
        }
      }
      cc.indx[(cc.indx == cc.indx[x[3]])] <- cc.indx[x[2]]
    }

    if(cc.indx[x[2]] != cc.indx[x[1]]){
      diff.rank <- gen.indx[x[2]] - (gen.indx[x[1]]+1)
      if(diff.rank != 0) {
        if(diff.rank > 0) {
          gen.indx = gen.indx + diff.rank*(cc.indx==cc.indx[x[1]])
        } else {
          gen.indx = gen.indx + -1*diff.rank*(cc.indx==cc.indx[x[2]])
        }
      }
      cc.indx[(cc.indx == cc.indx[x[1]])] <- cc.indx[x[2]]
    }

  }

  # 2nd step : rank the marr table from least recent to fill in offspring genotype
  sorted.order <- gen.indx[factor.tbl[,1]] %>% sort(decreasing = T, index.return=T) %>% .$ix
  sorted.id <- 1:n.id
  for(i in 1:nrow(factor.tbl)) {
    x<-factor.tbl[sorted.order[i],]
    indiv.sex[x[2]] <- 1
    indiv.sex[x[3]] <- 2
    sorted.id[x[1]] <- mating.lvl[x[1]]
    sorted.id[x[2]] <- mating.lvl[x[2]]
    sorted.id[x[3]] <- mating.lvl[x[3]]

    parent.copy <- rbinom(2,1,0.5)

    select.choice <- rbinom(n.snp, 1, 0.5)
    geno.1[x[1],] <- geno.1[x[2],]*select.choice + geno.2[x[2],]*(1-select.choice)
    select.choice <- rbinom(n.snp, 1, 0.5)
    geno.2[x[1],] <- geno.1[x[3],]*select.choice + geno.2[x[3],]*(1-select.choice)

  }
  # impose genotype error at the end
  make.err.1 <- rbinom(n.id*n.snp, 1, geno.err/2) %>% matrix(ncol=n.snp, by=T)
  make.err.2 <- rbinom(n.id*n.snp, 1, geno.err/2) %>% matrix(ncol=n.snp, by=T)

  is.missing <- rbinom(n.id*n.snp, 1, missing.unk) %>% matrix(ncol=n.snp, by=T)
  is.avail <- 1-is.missing

  final.geno.1 <- (abs(geno.1-make.err.1)*is.avail) +(-1*is.missing)
  final.geno.2 <- (abs(geno.2-make.err.2)*is.avail) +(-1*is.missing)

  geno.long <- cbind(final.geno.1, final.geno.2)[,rbind(1:n.snp, (1:n.snp)+n.snp) %>% as.vector()]

  geno.input.tbl <- cbind(sorted.id,
                          rep(1, n.id),
                          indiv.sex,
                          gen = max(gen.indx) - gen.indx,
                          geno.long) %>% as.tibble()

  ##remove any individuals from the generation censor blacklist
  if (!is.na(skip.gen.ls[1])) {
    geno.input.tbl <- geno.input.tbl %>%
      dplyr::filter(!gen %in% skip.gen.ls)
  }

  dir.create(file.path(out.path), recursive = TRUE)
  geno.path <- paste0(out.path, "/ingeno.txt")
  write.table(geno.input.tbl,
              geno.path,sep = " ", eol = "\n", quote = FALSE, col.names = FALSE, row.names = FALSE)

  param<-list(geno.path=geno.path,
              output.path=out.path,
              random.seed=random.seed,
              observe.frac= 1,
              af = maf,
              max.unobs=1, max.gen=0, #max.gen as the number extend past the current generation grp
              min.age=1, max.age=1, geno.err=geno.err,
              n.marr = paste0(mating.tbl$pa, "_",mating.tbl$ma) %>% unique %>% length)

  message("writing final files to ", out.path)

  writeIntermedGeno(param, preservesID = TRUE, preservesGeno = TRUE)
}

