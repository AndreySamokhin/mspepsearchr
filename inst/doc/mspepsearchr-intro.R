## ----knitr_setup, include = FALSE---------------------------------------------
knitr::opts_chunk$set(
collapse = TRUE,
comment = "#>",
fig.align = 'center',
fig.height = 3.5,
fig.width = 4.5
)

## ----setup, include = FALSE---------------------------------------------------
library(mspepsearchr)
r_info <- paste0("`", "mspepsearchr_", packageVersion("mspepsearchr"), "`")

## ----cli_call_example, include = FALSE, eval = FALSE--------------------------
# # jobs <- .PrepareJobs(
# #   spectra = "test.msp",
# #   libraries = c("C:/NIST23/MSSEARCH/mainlib/", "C:/NIST23/MSSEARCH/replib/"),
# #   algorithm = "identity_normal",
# #   presearch = "default",
# #   valid_presearch_vals = "default",
# #   addl_columns = c("/OutRevMF", "/OutMW", "/OutCAS", "/OutChemForm", "/OutIK")
# # )
# # cat(paste(c("<path>\\MSPepSearch64.exe", jobs[[1L]]$cli_args),
# #           collapse = " "))

## ----identity_normal----------------------------------------------------------
alkanes    <- system.file("extdata", "spectra", "alkanes_ei_lr.msp",
                             package = "mspepsearchr")
eims_lib <- system.file("extdata", "libraries", "massbank_subset_ei_lr_with_ri",
                        package = "mspepsearchr")
hitlists   <- IdentitySearchEiNormal(alkanes, eims_lib)

## ----print_hitlist------------------------------------------------------------
col_names <- c("name", "mf", "rmf", "prob", "formula", "mw", "inchikey", "cas")
head(hitlists[[4L]][, col_names], 3L)

## ----select_hitlist-----------------------------------------------------------
unknown_names <- vapply(hitlists, attr, which = "unknown_name", character(1L))
idx <- which(unknown_names == "Dodecane")
head(hitlists[[idx]][, col_names], 3L)

## ----spectra_object, eval = FALSE---------------------------------------------
# ex_spectrum <- readRDS("data/ex_spectrum.rds")
# 
# data_origin <- Spectra::dataOrigin(ex_spectrum)
# rtime <- Spectra::rtime(ex_spectrum)
# ms_level <- Spectra::msLevel(ex_spectrum)
# precursor_mz <- Spectra::precursorMz(ex_spectrum)
# ion_mode <- switch(Spectra::polarity(ex_spectrum) + 1L, "NEGATIVE", "POSITIVE")
# mz <- Spectra::mz(ex_spectrum)
# intst <- Spectra::intensity(ex_spectrum)
# 
# fenamiphos <- lapply(seq_along(ex_spectrum), function(i) {
#   list(name = paste0(basename(data_origin[[i]]), ", rt=", rtime[[i]], "s"),
#        spectrum_type = paste0("MS", ms_level[[i]]),
#        precursormz = precursor_mz[[i]],
#        ion_mode = ion_mode[[i]],
#        mz = mz[[i]],
#        intst = intst[[i]])
# })

## ----fenamiphos---------------------------------------------------------------
msp_path <- system.file("extdata", "spectra", "fenamiphos_msms_hr.msp",
                        package = "mspepsearchr")
fenamiphos <- mssearchr::ReadMsp(msp_path)


## ----identity_msms------------------------------------------------------------
msms_lib <- system.file("extdata", "libraries", "massbank_subset_msms_hr",
                        package = "mspepsearchr")
hitlists <-
  IdentitySearchMsMs(fenamiphos, msms_lib,
                     precursor_ion_tol  = list(value = 0.05, utits = "mz"),
                     product_ions_tol = list(value = 0.05, utits = "mz"))
col_names <- c("name", "dot", "rdot", "formula", "prec_mz", "delta_mz")
head(hitlists[[1L]][, col_names], 3L)

## ----error_duplicated_flags, eval = FALSE-------------------------------------
# hitlists <- IdentitySearchEiNormal(spectra, ms_library,
#                                    addl_cli_args = "/HITS 10")
# #> Error in .PrepareJobs(spectra, libraries, algorithm = "identity_normal",  :
# #>   The following CLI flags are duplicated: /HITS

## ----tridecane----------------------------------------------------------------
tridecane <- list(
  list(name = "Tridecane",
       mz = c(53, 55, 56, 57, 69, 70, 71, 84, 85),
       intst = c(51, 314, 220, 999, 110, 126, 526, 54, 274),
       retention_index = 1300)
)
hitlists <- IdentitySearchEiNormal(tridecane, eims_lib)
col_names <- c("name", "mf", "rmf", "prob", "formula", "mw", "ri")
head(hitlists[[1L]][, col_names], 3L)

## ----addl_cli_args, warning = FALSE-------------------------------------------
hitlists <- IdentitySearchEiNormal(tridecane, eims_lib,
                                   addl_cli_args = "/RI nt15rIN")
head(hitlists[[1L]][, col_names], 3L)

## -----------------------------------------------------------------------------
old_par <- par(no.readonly = TRUE)
on.exit(par(old_par))
par(mar = c(4.1, 4.1, 1.1, 1.1), cex = 0.9)

df <- readRDS("data/speedup_vs_task_complexity.rds")
cols <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3")
representative_tasks <- c("100_nist_off",
                          "2594_nist_default",
                          "500_nist_default",
                          "100_nist_default")
idxs <- match(representative_tasks, df$info)
plot(df$time_one_thread, df$speedup,
     pch = 21L, cex = 1.5, bg = "gray90",
     log = "x", xaxt = "n", yaxt = "n",
     xlab = "Execution time with one thread, s",
     ylab = "Speedup with four threads")
points(df$time_one_thread[idxs], df$speedup[idxs],
       pch = 21L, cex = 1.5, bg = cols)
x_ticks <- c(1, 10, 100, 1000)
y_ticks <- c(1, 2, 3)
axis(1, at = x_ticks, labels = x_ticks)
axis(2, at = y_ticks, labels = y_ticks[y_ticks])
rug(x = c(0.5, 1.5, 2.5, 3.5), ticksize = -0.03, side = 2)
legend(5900, 0.5, legend = representative_tasks,
       xjust = 1, yjust = 0,
       pch = 21L, pt.cex = 1.5, pt.bg = cols)

## -----------------------------------------------------------------------------
old_par <- par(no.readonly = TRUE)
on.exit(par(old_par))
par(mar = c(4.1, 4.1, 1.1, 1.1), cex = 0.9)

df <- readRDS("data/speedup_vs_num_threads.rds")
n_bars <- length(unique(df$info))
x = matrix(df$speedup, nrow = n_bars, byrow = TRUE)
colnames(x) <- seq_len(ncol(x))
first_bar_x_pos <- 1.5 + (seq_len(ncol(x)) - 1) * (n_bars + 1)
bar_x_pos <- rep(first_bar_x_pos, n_bars) +
  rep(seq(0, by = 1, length.out = n_bars), each = 7L)
barplot(x,
        xlab = "Number of threads",
        ylab = "Speedup",
        ylim = c(0, max(x) + 0.2),
        col = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3"),
        beside = TRUE,
        legend = c("3.8 min", "1.5 min", "17.9 s", "3.4 s"),
        args.legend = list(x = 0.5, y = max(x) + 0.1, xjust = 0))
arrows(x0 = bar_x_pos, x1 = bar_x_pos,
       y0 = df$speedup - df$speedup_err, y1 = df$speedup + df$speedup_err,
       code = 3, angle = 90, length = 0.04) # lwd = 2
box(bty = "l")

