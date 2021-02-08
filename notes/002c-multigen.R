## mock-up 5 generation pedigree, monogamous, mating with unrelated indiv,

# 30 founders, 15 male , 15 female

# mean sib size - pois mean of 2.5
# continue mating with rbin~0.5


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
  Pa = (1:length(add.in.layer.indiv)) + max.ID,
  Ma = (1:length(add.in.layer.indiv)) + (max.ID+length(add.in.layer.indiv)),
  gen = 3)

max.ID <<- max.ID + (length(add.in.layer.indiv)*2)


add.in.more.kids.rep <- rpois(length(add.in.layer.indiv),0.5)
add.in.more.kids <- add.in.more.founders[rep(row.names(add.in.more.founders), add.in.more.kids.rep),] %>%
  mutate(Kid = max.ID +row_number())



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


geno.pick.tbl <- simGeno(mating.tbl=crop.tbl %>%
                           mutate(across(.fns=as.numeric)),
        n.snp = 300,
        out.path = "/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/multi_case/ped1",
        num.unrelated.added.per.gen = 10,
        indiv.censor.ls = c("21", "22", "28", "24", "13", "20", "49", "68", "57", "42", "44","34","58","27","33","32","69", "50", "3", "17", "51", "56"),
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


## run pedfac

run.multiped1.summary <- list()

## factorize ped.df
set.seed(234)
outfolder.path <- "/Users/thomasn/repo/pedigree-factor-graphs/data/ch2/multi_case/ped1"
random.arr <- ceiling(c(runif(2,0,1e6)))

write.table(random.arr, paste0(outfolder.path,"/","rand"),
              sep = " ",eol = " ",quote = FALSE, col.names = FALSE, row.names = FALSE)

run.summary$geno.input.tbl <- geno.pick.tbl

## set strictly mono
pedfac.call <- paste0("/Users/thomasn/repo/pedigree-factor-graphs/src/pedigraph_0.185 -d ", outfolder.path,
                      " -n 100",
                      " -r ",outfolder.path,"/rand -j 1 -f 10", collapse = "")

  tic("run pedFac:")
  system(pedfac.call,ignore.stderr = F, ignore.stdout = F)
  runtime <- toc()
  write.table(runtime$toc - runtime$tic, file =paste0(outfolder.path,"/pedFac.time"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  run.summary$runtime <- runtime$toc - runtime$tic

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
