.noGenerics <- TRUE

.onLoad <- function(lib, pkg)
{
  opath <-  Sys.getenv("PATH")
  libbin <- file.path(R.home(), "library/hdf5/libs")
  Sys.putenv(PATH=paste(libbin, opath, sep=";"))
  library.dynam("hdf5", pkg, lib)
  Sys.putenv(PATH=opath)
}

.onUnload <- function (libpath)
{
   library.dynam.unload("hdf5", libpath)
}
