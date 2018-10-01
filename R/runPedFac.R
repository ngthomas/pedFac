#' Initialize pedigree run
#'
#' The function \code{runPedFac} prepares and runs a multigeneration pedigree sampler given input genetic data
#'
#' @param geno.path string. Path to the input genotype file (see checkGeno help for formatting). Required.
#' @param marker.path string. Directed path to the option marker info file.
#' @param output.path string. Path to store intermed and final output files. By default, outputPath will be set in an upper folder of the current user's working directory.
#'
#' Regarding sampling:
#' @param random.seed postive integer. Random seed to pass on pedigree sampler. Default:26.
#' @param n.iter positive integer. Number of sampling iteration. Default: 1.
#' @param cyclic.choice integer from 0 to 2. Choices of handling loops or cyclic path in pedigree. Default: 0. 0 - not allowing loops; 1 - throttle method; 2 - decimation method
#' @param observe.frac float value from 0 to 1. Assumed sampling fraction. Default: 0.8; use -1 for unknown.
#' @param max.unobs nonnegative integer. Maximum number of unobserved individuals allowed in between any two individuals. Default: 1.
#' @param max.gen nonnegative integer. Number of predecessor generation(s) considered beyond the earliest observed generation. Default: 0. Setting it as 0 means that individuals of the earliest observed generation are treated as founders.
#'
#' Regarding specie life history:
#' @param min.age positive float value. Minimal age of sexual maturation or fecundity (in year). Default: 1.
#' @param max.age positive float value. Maximum age of sexual maturation or fecundity (in year). Default: 1.
#'
#' Regarding genotype marker:
#' @param haplo.method positive integer 0 - 2. Selected method in the case of handling multiallelic markers. Default: 0. 0 - taking the most informative allele whose frequency is closest to 0.5; 1 - (not avail) deconstructing haplotype into a set of nucleotide units; 2 - (not avail) reduce the multiallelic basis into n class of binomial switches
#' @param geno.err float value from 0 to 1. Assumed background genotype error rate in the form of epsilon. Default: 0.02. If the genotype error row - 'gerror' of marker_info.txt is provided, this param will be overridden.
#'
#' @export
#' @examples
#' # runHaplot(run.label, sam.path, label.path, vcf.path)
#' runPedFac(geno.path="/Users/thomasn/repo/pedfac/example/case_0/genotype.txt")
runPedFac <- function(geno.path,
                      marker.path="",
                      output.path="",
                      random.seed=26,
                      n.iter =1,
                      cyclic.choice=0,
                      observe.frac=0.8,
                      max.unobs =1,
                      max.gen =0,
                      min.age =1,
                      max.age =1,
                      haplo.method =0,
                      geno.err =0.02){

  param<-list(geno.path=geno.path, marker.path=marker.path,
              output.path=output.path,
              random.seed=random.seed,
              n.iter=n.iter, cyclic.choice=cyclic.choice,
              observe.frac=observe.frac,
              max.unobs=max.unobs, max.gen=max.gen,
              min.age=min.age, max.age=max.age,
              haplo.method=haplo.method, geno.err=geno.err)


  if(!file.exists(geno.path)) stop("the genotype file 'geno.path' provided - ", geno.path, " does not exist")
  if(!file.exists(marker.path) && marker.path != "") stop("the marker file 'marker.path' provided - ", marker.path, " does not exist")
  if(!file.exists(output.path) && output.path != "") stop("the output path 'output.path' provided - ", output.path, " does not exist")

  #return(param$output.path)
  if(param$output.path=="") param$output.path = getwd()

  param$output.path <- paste0(param$output.path,"/","run")

  if(!file.exists(param$output.path)) {
    message("creating directory:", param$output.path)
    dir.create(param$output.path)
  }

  if(is.character(random.seed)) stop("Random seed cannot be a string")
  set.seed(random.seed)

  if (n.iter <1) stop("Number of iteration, n.iter, must be equal or greater than 1")
  if (!(cyclic.choice %in% c(0,1,2))) stop("Out of acceptable range for cyclic.choice")
  if (observe.frac > 1 | (observe.frac < 0 & observe.frac != -1)) stop("Out of acceptable range for observe.frac")
  if (max.unobs <0) stop("max.unobs must be equal or greater than 0")
  if (max.gen <0) stop("max.gen must be equal or greater than 0")

  if (min.age <=0) stop("min.age must be greater than 0")
  if (max.age <=0) stop("max.age must be greater than 0")

  if (!(haplo.method %in% c(0,1,2))) stop("Out of acceptable range for haplo.method")
  if (geno.err > 1 | geno.err < 0) stop("Out of acceptable range for geno.err")

  # create 2 random seeds
  random.arr <- ceiling(c(runif(2,0,1e6)))
  write.table(random.arr, paste0(param$output.path,"/","rand"),
              sep = " ",eol = " ",quote = FALSE, col.names = FALSE, row.names = FALSE)

  # 1. check whether geno.path file is acceptable
  checkGeno(geno.path)

  # 2. if so, write intermed files
  message("Creating intermediate files ...")
  writeIntermedGeno(param)

  # 3. run pedFactory
  message("Running pedFactory ...")

  pedfac.call <- paste0(pedigraph_binary(), " -d ", param$output.path,
                         " -n ",param$n.iter,
                         " -r ",param$output.path,"/rand", collapse = "")

  if (Sys.info()['sysname'] == "Windows") {
    shell(pedfac.call)
  }
  else{
    system(pedfac.call)
  }

  # 4. clean up pedFactory
  # param$output.path = "/Users/thomasn/repo/pedigree-factor-graphs/data/unit_norm_1"
  # param$output.path = "/Users/thomasn/repo/pedigree-factor-graphs/data/unit_threegen_1"
  ped.tbl <- read.table(paste0(param$output.path, "/out/ped.txt",collapse = ""), stringsAsFactors = FALSE) %>%
    dplyr::tbl_df() %>%
    dplyr::rename("iter"="V1", "kid"="V2", "pa"="V3", "ma"="V4") %>%
    dplyr::group_by(iter, kid) %>%
    dplyr::summarise(pa = pa[dplyr::n()],
                     ma = ma[dplyr::n()])

  id.ls <-readRDS(paste0(param$output.path, "/id.rds",collapse = ""))$id
  max.id <- length(id.ls)
  # report complete pedigree, parentage, half-sib, full-sib, grandparentage posterior assignment

  #converting the table into id
  ped.id.tbl <- ped.tbl %>%
    dplyr::group_by(iter, kid, pa , ma) %>%
    dplyr::mutate(kid.id = ifelse(kid > length(id.ls), paste0("unobs.", kid, collapse = ""),
                                  as.character(id.ls[kid])),
                  pa.id = ifelse(pa > length(id.ls), paste0("unobs.", pa, collapse = ""),
                                 as.character(id.ls[pa])),
                  ma.id = ifelse(ma > length(id.ls), paste0("unobs.", ma, collapse = ""),
                                 as.character(id.ls[ma])))

  message("writing ped_sample.txt")
  write.table(ped.id.tbl %>% dplyr::ungroup() %>% dplyr::select(-kid, -pa, -ma),
              paste0(param$output.path,"/ped_sample.txt"),
              sep = " ",eol = "\n",quote = FALSE, col.names = FALSE, row.names = FALSE)

  message("writing parent_assignment.txt")
  write.table(retrieveParent(ped.tbl, max.id, id.ls),
              paste0(param$output.path,"/parent_assignment.txt"),
              sep = " ",eol = "\n",quote = FALSE, col.names = FALSE, row.names = FALSE)

  message("writing fullsib_assignment.txt halfsib_assignment.txt ")
  write.table(retrieveFullSib(ped.tbl, max.id, id.ls),
              paste0(param$output.path,"/fullsib_assignment.txt"),
              sep = " ",eol = "\n",quote = FALSE, col.names = FALSE, row.names = FALSE)
  write.table(retrieveHalfSib(ped.tbl, max.id, id.ls),
              paste0(param$output.path,"/halfsib_assignment.txt"),
              sep = " ",eol = "\n",quote = FALSE, col.names = FALSE, row.names = FALSE)

  message("writing grandparent_assignment.txt")
  grandparent.tbl <- retrieveGrandparents(ped.tbl, max.id, id.ls)
  write.table(grandparent.tbl,
              paste0(param$output.path,"/grandparent_assignment.txt"),
              sep = " ",eol = "\n",quote = FALSE, col.names = FALSE, row.names = FALSE)

  #report grandparentage assignment disregarding whether it is belong to maternal side or paternal side

  if(nrow(grandparent.tbl)>0) {
  grandparent.v2 <- grandparent.tbl %>%
    dplyr::group_by(kid, grandpa.pa, grandma.pa, grandpa.ma, grandma.ma) %>%
    dplyr::mutate(univ.id = paste0(sort(c(paste0(c(grandpa.pa,grandma.pa),collapse = "-"),
                                          paste0(c(grandpa.ma,grandma.ma),collapse = "-"))),
                                   collapse = "-")) %>%
    dplyr::group_by(kid, univ.id) %>%
    dplyr::mutate(prob = sum(prob),
                  rnum = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(rnum ==1) %>%
    dplyr::select(-univ.id, -rnum)

  write.table(grandparent.v2,
              paste0(param$output.path,"/grandparent_assignment_v2.txt"),
              sep = " ",eol = "\n",quote = FALSE, col.names = FALSE, row.names = FALSE)
  }
}
