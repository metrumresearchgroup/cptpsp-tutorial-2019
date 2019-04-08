author <- c("--")

pkgs <- c("tidyverse", "qapply", "pmplots", "yspec", "mrggsave", 
          "mrgtable", "tidynm", "metrumrg", "fork", "review",
          "mrgsolve", "cowplot", "nloptr", "FME", "sensitivity")


pkgRoot <- "pkg"
pkgDir <- file.path(pkgRoot, "src", "contrib")
pkgDir <- normalizePath(pkgDir)
libDir <- "lib"

if(!dir.exists(pkgDir)) dir.create(pkgDir, recursive = TRUE)
if(!dir.exists(libDir)) dir.create(libDir)

.libPaths(libDir)

user <- Sys.info()["user"]

fromCRAN <- user %in% author | '*' %in% author

local_repos <- paste0("file://",file.path(getwd(),pkgRoot))
metrum_repos <- "https://metrumresearchgroup.github.io/r_validated/"
cran_repos <- "https://cran.rstudio.com/"
repos <- c(mrg = metrum_repos, cran = cran_repos, local = local_repos)


deps <- tools::package_dependencies(
  packages = pkgs,
  which = c("Depends", "Imports", "LinkingTo"),
  recursive = TRUE,
  db = available.packages(repos=repos[c("mrg", "cran")])
)

deps <- unlist(deps, use.names=FALSE)

pkgs <- unique(c(pkgs,deps))

base <- rownames(installed.packages(priority=c("base", "recommended")))

pkgs <- setdiff(pkgs,base)

installed <- row.names(installed.packages(libDir))

tools::write_PACKAGES(pkgDir)

if(file.exists(file.path(pkgDir,"PACKAGES"))){
  available <- available.packages(repos = repos["local"])[,"Package"]
} else{
  available <- NULL
  file.create(file.path(pkgDir,"PACKAGES"))
  tools::write_PACKAGES(pkgDir)
}

if(fromCRAN){
  
  newpkgs <- setdiff(pkgs, available)
  
  if(length(newpkgs) > 0){
    ## These packages are installed either from mrg or cran
    install.packages(newpkgs,
                     lib=libDir,
                     repos = repos[c("mrg", "cran")],
                     destdir=pkgDir,
                     type="source", 
                     INSTALL_opts="--no-multiarch")
    
    tools::write_PACKAGES(pkgDir)
  }
  
  ## If multiple authors qcing each other, a package could be available
  ## but uninstalled.  Install from local.
  uninstalled <- setdiff(pkgs, installed)
  
  if(length(uninstalled)>0){
    install.packages(uninstalled,
                     lib = libDir,
                     repos = repos["local"],
                     type = "source",
                     INSTALL_opts="--no-multiarch")
  }    
}


if(!fromCRAN){
  newpkgs <- setdiff(pkgs, installed)
  if(length(newpkgs)>0){
    install.packages(newpkgs,
                     lib = libDir,
                     repos = repos["local"],
                     type = "source",
                     INSTALL_opts="--no-multiarch")
    
  }
}

# Check for packages that were requested, but aren't installed at this point
not_found <- setdiff(pkgs, installed.packages(libDir)[,"Package"]) 
if(length(not_found) > 0) {
  warning("Some requested packages are not installed:", immediate. = TRUE) 
  cat(paste0(" - ",not_found), sep = "\n")
}



.ignore_libs <- function(root=getwd(),lib="lib", ci=FALSE) {
  
  if(!missing(root) & file.exists(root)) {
    lib <- file.path(root,"lib")
  }
  if(!file.exists(lib)) stop("Could not find lib directory")
  libs <- list.files(lib, full.names=FALSE)
  libs <- c(libs, "ignore.txt", "PACKAGES", "PACKAGES.gz")
  writeLines(con=file.path(lib,"ignore.txt"), libs)
  setwd(lib)
  system("svn propset svn:ignore -F ignore.txt .")
  setwd("..")
  if(ci) system("svn ci -m \"ignoring libs\" .")
}

.update_packages <- function(pkg = "pkg", lib = "lib") {
  repos <- paste0("file://", normalizePath(pkg,winslash = "/"))
  lib <- normalizePath(lib)
  avail <- available.packages(repos = repos, filters = "duplicates")
  n_pkg <- length(unique(avail[,"Package"]))
  message(paste0("updating from ", repos))
  message(paste0("available packages: ", n_pkg))
  message(paste0("lib directory: ", lib))
  update.packages(lib.loc = lib, repos = repos, ask = FALSE)
  return(invisible(list(repos = repos, lib = lib, avail  = avail)))
}