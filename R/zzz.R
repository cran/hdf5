.noGenerics <- TRUE

.onUnload <- function (libpath)
{
   library.dynam.unload("hdf5", libpath)
}
