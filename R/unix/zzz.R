.noGenerics <- TRUE

.onLoad <- function (lib, pkg)
{
  library.dynam("hdf5", pkg, lib)
}

.onUnload <- function (libpath)
{
   library.dynam.unload("hdf5", libpath)
}
