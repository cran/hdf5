.noGenerics <- TRUE

.onLoad <- function(lib, pkg)
{
  opath <-  Sys.getenv("PATH")
  libbin <- file.path(R.home(), "library/hdf5/libs")
  if(!exists("Sys.setenv")) Sys.setenv <- Sys.putenv
  Sys.setenv(PATH=paste(libbin, opath, sep=";"))
  library.dynam("hdf5", pkg, lib)
  Sys.setenv(PATH=opath)
}

.onUnload <- function (libpath)
{
   library.dynam.unload("hdf5", libpath)
}
