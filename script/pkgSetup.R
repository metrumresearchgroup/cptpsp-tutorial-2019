pkgs <- c("tidyverse", "mrgsolve", "cowplot")

repo <- "https://cran.microsoft.com/snapshot/2019-08-12"

pkgRoot <- "pkg"
pkgDir <- file.path(pkgRoot, "src", "contrib")
pkgDir <- normalizePath(pkgDir)
libDir <- "lib"

if(!dir.exists(pkgDir)) dir.create(pkgDir, recursive = TRUE)
if(!dir.exists(libDir)) dir.create(libDir)

.libPaths(libDir)

install.packages(pkgs, lib = libDir, destdir=pkgDir, repos=repo)