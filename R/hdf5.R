"hdf5load" <-  function (file, load = TRUE)
{
  call <- sys.call()
  .External("do_hdf5load", call, sys.frame(sys.parent()), file, load)
}


"hdf5save" <- function (fileout, ...)
{
  call <- sys.call()
  invisible(.External("do_hdf5save", call, sys.frame(sys.parent()),
                      fileout, ...))
}

"hdf5cleanup" <- function (fid)
{
  call <- sys.call()
  invisible(.External("do_hdf5cleanup", call, sys.frame(sys.parent()), fid))
}
