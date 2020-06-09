###################################################
######## Automate ABSOLUTE Calling ################
###################################################
# Copyright @Shixiang Wang <w_shixiang@163.com> ###
###################################################
#' Automate ABSOLUTE calling for multiple samples in parallel way
#'
#' An example can be found at [README](https://github.com/ShixiangWang/DoAbsolute#example).
#'
#' [ABSOLUTE](https://www.nature.com/articles/nbt.2203) is a famous software
#' developed by Broad Institute, however the `RunAbsolute` function is designed for
#' computing one sample each time and set no default values. **DoAbsolute** help
#' user set default parameters according to [ABSOLUTE documentation](http://software.broadinstitute.org/cancer/software/genepattern/modules/docs/ABSOLUTE), provide an uniform interface to
#' input data easily and run RunAbsolute parallelly.
#'
#' More detail about how to analyze ABSOLUTE results please see [this link](http://software.broadinstitute.org/cancer/software/genepattern/analyzing-absolute-data).
#' @param Seg a `data.frame` or a file (path) contains columns
#' "Sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean".
#' @param Maf MAF, default is `NULL`, can provided as `data.frame` or file path.
#' @param sigma.p Provisional value of excess sample level variance used for mode search. Default: 0
#' @param max.sigma.h Maximum value of excess sample level variance (Eq. 6). Default: 0.2
#' @param min.ploidy Minimum ploidy value to consider. Solutions implying lower ploidy values will be discarded. Default: 0.5
#' @param max.ploidy Maximum ploidy value to consider. Solutions implying greater ploidy values will be discarded. Default: 10
#' @param primary.disease Primary disease of the sample. Default: `NA`
#' @param platform one of "SNP_6.0", "Illumina_WES", "SNP_250K_STY". Default: "SNP_6.0"
#' @param temp.dir directory path used to store tempory files. Default: Absolute subdirectory under `tempdir()`
#' @param clean.temp if `TRUE`, auto-clean temp dir at the end. Default: `FALSE`
#' @param results.dir directory path used to store result files. Default: work directory
#' @param max.as.seg.count Maximum number of allelic segments. Samples with a higher segment count will be flagged as 'failed'. Default: 1500
#' @param max.non.clonal Maximum genome fraction that may be modeled as non-clonal (subclonal SCNA). Solutions implying greater values will be discarded. Default: 0.05
#' @param max.neg.genome Maximum genome fraction that may be modeled as non-clonal with copy-ratio below that of clonal homozygous deletion. Solutions implying greater values will be discarded. Default: 0.005
#' @param copy.num.type The type of copy number to be handled. Either total or allelic. Total is what this package for. Default: "total"
#' @param min.mut.af Minimum mutation allelic fraction. Mutations with lower allelic fractions will be filtered out before analysis. Default: 0.1
#' @param min.no.mut Minor allele frequency file, or NULL if one is not available. This specifies the data for somatic point mutations to be used by ABSOLUTE. Default: 5
#' @param verbose if `TRUE`, print extra info. Default: `FALSE`
#' @param nThread number of cores used for computation. Default: 1L
#' @param keepAllResult if `TRUE`, clean all results, otherwise clean result directory and keep most important results. Default: `TRUE`
#' @param recover if `TRUE`, recover previous unfinished work.
#' This is helpful when program stop unexpectedly when `clean.temp` is FALSE. Default: `FALSE`
#' @author Shixiang Wang <w_shixiang@163.com>
#' @return NULL
#' @import foreach doParallel data.table utils parallel
#' @export
#' @references Carter, Scott L., et al. "Absolute quantification of somatic DNA alterations in human cancer." Nature biotechnology 30.5 (2012): 413.
#'
DoAbsolute <- function(Seg, Maf = NULL,
                       sigma.p = 0, max.sigma.h = 0.2, min.ploidy = 0.5, max.ploidy = 10,
                       primary.disease = NA, platform = c("SNP_6.0", "Illumina_WES", "SNP_250K_STY"),
                       temp.dir = file.path(tempdir(), "Absolute"), clean.temp = FALSE,
                       results.dir = getwd(), max.as.seg.count = 1500,
                       max.non.clonal = 0.05, max.neg.genome = 0.005, copy.num.type = c("total", "allelic"),
                       min.mut.af = 0.1, min.no.mut = 5, verbose = FALSE, nThread = 1L, keepAllResult = TRUE,
                       recover = FALSE) {
  # if (dir.exists(temp.dir) & length(dir(path = temp.dir)) != 0) {
  #     stop("Your 'temp.dir' is not empty, please reset it.")
  # }

  if (!suppressMessages(requireNamespace("ABSOLUTE"))) {
    stop("Find no package called 'ABSOLUTE', please install it...")
  }

  if (!suppressMessages(requireNamespace("foreach"))) {
    warning("Find no package called 'foreach', try to install it...")
    install.packages("foreach", dependencies = TRUE)
  }

  if (!suppressMessages(requireNamespace("doParallel"))) {
    warning("Find no package called 'doParallel', try to install it...")
    install.packages("doParallel", dependencies = TRUE)
  }

  if (!suppressMessages(requireNamespace("data.table"))) {
    warning("Find no package called 'data.table', try to install it...")
    install.packages("data.table", dependencies = TRUE)
  }

  if (verbose) cat("-> Loading segmentation data...\n")
  if (is.character(Seg)) {
    if (file.exists(Seg)) {
      Seg <- data.table::fread(input = Seg)
    } else {
      stop("file ", Seg, " does not exist")
    }
  } else {
    if (inherits(Seg, "data.frame")) {
      data.table::setDT(Seg)
    } else {
      Stop("Unsupport Segmentation Format!")
    }
  }

  if (!is.null(Maf)) {
    if (verbose) cat("-> Loading Maf data...\n")
    if (is.character(Maf)) {
      if (file.exists(Maf)) {
        Maf <- data.table::fread(input = Maf)
      } else {
        stop("file ", Maf, " does not exist")
      }
    } else {
      if (inherits(Maf, "data.frame")) {
        data.table::setDT(Maf)
      } else {
        Stop("Unsupport Maf Format!")
      }
    }
  }


  if (verbose) cat("-> Checking data format of segmentation file...\n")
  #-- check Seg data
  seg_cols <- c("Sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")
  if (!all(seg_cols %in% colnames(Seg))) {
    stop("Columns ", paste(seg_cols, collapse = " "), " are should in file.")
  } else {
    Seg <- Seg[, c("Sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")]
  }

  Seg$Chromosome <- as.character(Seg$Chromosome)
  Seg$Chromosome <- gsub(pattern = "chr", replacement = "", Seg$Chromosome, ignore.case = TRUE)
  Seg$Chromosome <- gsub(pattern = "X", replacement = "23", Seg$Chromosome, ignore.case = TRUE)
  if (verbose) cat("-> Keeping only chr 1-23 for CNV data...\n")
  autosome <- as.character(seq(1, 23))
  Seg <- Seg[Chromosome %in% autosome, ]

  #-- check Maf data
  if (!is.null(Maf)) {
    maf_cols <- c(
      "Tumor_Sample_Barcode", "Hugo_Symbol", "Chromosome", "Start_position", "dbSNP_Val_Status",
      "t_ref_count", "t_alt_count"
    )
    maf_cols2 <- c(
      "Tumor_Sample_Barcode", "Hugo_Symbol", "Chromosome", "Start_Position", "dbSNP_Val_Status",
      "t_ref_count", "t_alt_count"
    )
    if (all(maf_cols %in% colnames(Maf))) {
      Maf <- Maf[, c(
        "Tumor_Sample_Barcode", "Hugo_Symbol", "Chromosome", "Start_position", "dbSNP_Val_Status",
        "t_ref_count", "t_alt_count"
      )]
    } else if (all(maf_cols2 %in% colnames(Maf))) {
      Maf <- Maf[, c(
        "Tumor_Sample_Barcode", "Hugo_Symbol", "Chromosome", "Start_Position", "dbSNP_Val_Status",
        "t_ref_count", "t_alt_count"
      )]
      colnames(Maf) <- maf_cols
    } else {
      stop("Necessary columns for Maf file: \n          ", paste(maf_cols, collapse = " "), "\n")
    }

    Maf$Chromosome <- as.character(Maf$Chromosome)
    Maf$Chromosome <- gsub(pattern = "chr", replacement = "", Maf$Chromosome, ignore.case = TRUE)
    Maf$Chromosome <- gsub(pattern = "X", replacement = "23", Maf$Chromosome, ignore.case = TRUE)
    if (verbose) cat("-> Keeping only chr 1-23 for Maf data...\n")
    Maf <- Maf[Chromosome %in% autosome, ]
  }



  #-- create temp files to generate input for RunAbsolute
  samples <- unique(Seg$Sample)
  SAMPLES <- samples
  seg_filepath <- vector(mode = "character", length = length(samples))
  maf_filepath <- vector(mode = "character", length = length(samples))

  if (verbose) cat("-> Creating temp directory ...\n")
  if (!dir.exists(temp.dir)) {
    dir.create(temp.dir, recursive = TRUE, showWarnings = FALSE)
  }

  if (verbose) cat("-> Spliting seg data of samples to different files...\n")
  for (i in seq_along(samples)) {
    if (verbose) cat("--> Processing sample ", samples[i], "...\n")
    seg <- Seg[Sample == samples[i], ]
    seg_filepath[i] <- file.path(temp.dir, paste0(samples[i], ".seg"))
    data.table::fwrite(x = seg, file = seg_filepath[i], sep = "\t")

    if (is.null(Maf)) {
      maf_filepath[i] <- NA_character_
    } else {
      maf <- Maf[Tumor_Sample_Barcode == samples[i], ]
      if (nrow(maf) == 0) {
        if (verbose) cat("---> This sample has not Maf data, skipping...\n")
        maf_filepath[i] <- NA_character_
      } else {
        if (verbose) cat("---> Filtering mutations which vaf<", min.mut.af, "...\n")
        maf <- maf[(t_alt_count / (t_ref_count + t_alt_count)) >= min.mut.af, ]
        if (verbose) cat("---> Filtering Maf which count<", min.no.mut, "...\n")
        if (nrow(maf) < min.no.mut) {
          if (verbose) cat("---> This sample has not Maf data (after filtering), skipping...\n")
          maf_filepath[i] <- NA_character_
        } else {
          if (verbose) cat("---> Outputing corresponding Maf file...\n")
          maf_filepath[i] <- file.path(temp.dir, paste0(samples[i], ".maf"))
          data.table::fwrite(x = maf, file = maf_filepath[i], sep = "\t")
          # write.csv(x = maf, file = maf_filepath[i], sep = "\t", quote = FALSE, row.names = FALSE)
        }
      }
    }
  }
  if (verbose) cat("-> Spliting seg data of samples done.\n")

  #-- match options
  platform <- match.arg(platform)
  copy_num_type <- match.arg(copy.num.type)

  cache.dir <- file.path(temp.dir, "cache")
  if (!dir.exists(cache.dir)) {
    dir.create(cache.dir, recursive = TRUE, showWarnings = FALSE)
  }

  # foreach loop
  nCores <- detectCores(logical = FALSE)
  if (nThread > nCores) {
    warning("Number of real physical core is ", nCores, " while user set ", nThread, "\nUse ", nCores, " Cores.")
    nThread <- nCores
  }

  if (Sys.info()[["sysname"]] == "Windows") {
    cl <- makeCluster(nThread)
    registerDoParallel(cl)
  } else {
    registerDoParallel(cores = nThread)
  }


  if (recover) {
    cat("-> recover mode is TRUE, checking samples have been called...\n")
    not_called <- c()
    maf_filepath <- c()
    seg_filepath <- c()
    for (i in seq_along(samples)) {
      if (!file.exists(file.path(cache.dir, paste0(samples[i], ".ABSOLUTE.RData")))) {
        cat("--> ", samples[i], "is not called.\n")
        not_called <- c(not_called, samples[i])
        if (!file.exists(file.path(temp.dir, paste0(samples[i], ".seg")))) {
          stop("Something wrong with splitted files.")
        }
        seg_filepath <- c(seg_filepath, file.path(temp.dir, paste0(samples[i], ".seg")))
        if (!file.exists(file.path(temp.dir, paste0(samples[i], ".maf")))) {
          maf_filepath <- c(maf_filepath, NA_character_)
        } else {
          maf_filepath <- c(maf_filepath, file.path(temp.dir, paste0(samples[i], ".maf")))
        }
      }
    }
    samples <- not_called

    if (length(not_called) != 0) {
      cat("-> ABSOLUTE calling for above samples will be recovered.\n")
    } else {
      cat("-> ABSOLUTE calling has been done. \n")
    }
  }

  if (length(samples) != 0) {
    if (verbose) cat("-> Running RunAbsolute...\n")

    if (nThread == 1) {
      for (i in seq_along(samples)) {
        if (is.na(maf_filepath[i])) {
          maf_fn <- NULL
        } else {
          maf_fn <- maf_filepath[i]
        }
        seg_fn <- seg_filepath[i]
        if (verbose) cat("--> Processing sample ", samples[i], "...\n")
        suppressWarnings(ABSOLUTE::RunAbsolute(
          seg.dat.fn = seg_fn, maf.fn = maf_fn, # output.fn.base = output.fn.base,
          sample.name = samples[i],
          sigma.p = sigma.p, max.sigma.h = max.sigma.h,
          min.ploidy = min.ploidy, max.ploidy = max.ploidy,
          primary.disease = primary.disease, platform = platform,
          results.dir = cache.dir, max.as.seg.count = max.as.seg.count,
          max.non.clonal = max.non.clonal, max.neg.genome = max.neg.genome,
          copy_num_type = copy_num_type,
          min.mut.af = min.mut.af, verbose = verbose
        ))
      }
    } else {
      foreach(i = seq_along(samples)) %dopar% {
        if (is.na(maf_filepath[i])) {
          maf_fn <- NULL
        } else {
          maf_fn <- maf_filepath[i]
        }
        seg_fn <- seg_filepath[i]
        if (verbose) cat("--> Processing sample ", samples[i], "...\n")
        suppressWarnings(ABSOLUTE::RunAbsolute(
          seg.dat.fn = seg_fn, maf.fn = maf_fn, # output.fn.base = output.fn.base,
          sample.name = samples[i],
          sigma.p = sigma.p, max.sigma.h = max.sigma.h,
          min.ploidy = min.ploidy, max.ploidy = max.ploidy,
          primary.disease = primary.disease, platform = platform,
          results.dir = cache.dir, max.as.seg.count = max.as.seg.count,
          max.non.clonal = max.non.clonal, max.neg.genome = max.neg.genome,
          copy_num_type = copy_num_type,
          min.mut.af = min.mut.af, verbose = verbose
        ))
      }
    }
    if (verbose) cat("-> RunAbsolute done. Retrieving results...\n")
  }

  ## This will make problem when iteration computation with same temp dir
  # absolute_files = file.path(cache.dir, grep("RData", dir(cache.dir), value = TRUE))
  absolute_files <- file.path(cache.dir, paste0(SAMPLES, ".ABSOLUTE.RData"))
  if (verbose) cat("-> Checking result files...\n")
  for (f in absolute_files) {
    if (!file.exists(f)) {
      warning("--> Result file ", f, " does not exist, drop it.", immediate. = TRUE)
      absolute_files <- setdiff(absolute_files, f)
    }
  }
  if (length(absolute_files) < 1) {
    stop("No result file to proceed.")
  }

  if (verbose) cat("-> Checked.\n")

  review.dir <- file.path(cache.dir, "review")
  if (dir.exists(review.dir)) {
    if (verbose) cat("-> Removed previous temp review result directory.\n")
    unlink(review.dir, recursive = TRUE)
  }

  if (verbose) cat("-> Running Absolute summarize...\n")
  suppressWarnings(ABSOLUTE::CreateReviewObject(
    obj.name = "DoAbsolute",
    absolute.files = absolute_files,
    indv.results.dir = review.dir,
    copy_num_type = copy_num_type,
    plot.modes = TRUE, verbose = verbose
  ))
  if (verbose) cat("\n-> Absolute summarize done. Prepare auto-reviewing...\n")

  # pp_call_fn = file.path(review.dir, grep("PP-calls_tab.txt", dir(review.dir), value = TRUE))
  # modes_fn = file.path(review.dir, grep("PP-modes.data.RData", dir(review.dir), value = TRUE))
  pp_call_fn <- file.path(review.dir, "DoAbsolute.PP-calls_tab.txt")
  modes_fn <- file.path(review.dir, "DoAbsolute.PP-modes.data.RData")
  suppressWarnings(ABSOLUTE::ExtractReviewedResults(
    reviewed.pp.calls.fn = pp_call_fn,
    analyst.id = "wsx",
    modes.fn = modes_fn,
    out.dir.base = review.dir,
    obj.name = "DoAbsolute",
    copy_num_type = copy_num_type,
    verbose = verbose
  ))
  if (verbose) cat("-> Absolute Auto-reviewing done.\n")

  reviewed.dir <- file.path(review.dir, "reviewed")

  if (verbose) cat("-> Outputing final results...\n")
  if (keepAllResult) {
    cat("--> Choose keeping all results...\n")
  } else {
    cat("--> Choose not keeping all results. Keepping only final results...\n")
  }


  seg.dir <- file.path(results.dir, "seg")
  maf.dir <- file.path(results.dir, "maf")
  dir.create(seg.dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(maf.dir, showWarnings = FALSE, recursive = TRUE)

  t <- file.copy(file.path(reviewed.dir, grep("DoAbsolute", dir(reviewed.dir), value = TRUE)), results.dir)
  files_seg <- file.path(
    file.path(reviewed.dir, "SEG_MAF"),
    grep("segtab.txt", dir(file.path(reviewed.dir, "SEG_MAF")), value = TRUE)
  )
  t <- file.copy(from = files_seg, to = seg.dir)
  files_maf <- file.path(
    file.path(reviewed.dir, "SEG_MAF"),
    grep("ABS_MAF.txt", dir(file.path(reviewed.dir, "SEG_MAF")), value = TRUE)
  )
  t <- file.copy(from = files_maf, to = maf.dir)

  if (keepAllResult) {
    call.dir <- file.path(results.dir, "summary")
    samples_before_call.dir <- file.path(results.dir, "sample_before_summary")
    samples_called.dir <- file.path(results.dir, "sample_final_called")
    dir.create(call.dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(samples_before_call.dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(samples_called.dir, showWarnings = FALSE, recursive = TRUE)

    files_call <- file.path(review.dir, grep("DoAbsolute", dir(review.dir), value = TRUE))
    t <- file.copy(from = files_call, to = call.dir)
    t <- file.copy(from = absolute_files, to = samples_before_call.dir)
    files_samplesCalled <- file.path(
      file.path(reviewed.dir, "samples"),
      dir(file.path(reviewed.dir, "samples"))
    )
    t <- file.copy(from = files_samplesCalled, samples_called.dir)
  }

  if (clean.temp) {
    unlink(temp.dir, recursive = TRUE, force = TRUE)
  }

  cat("-> Done.\n")
}


utils::globalVariables(
  c(
    "Chromosome",
    "Sample",
    "Stop",
    "Tumor_Sample_Barcode",
    "t_alt_count",
    "t_ref_count"
  )
)
