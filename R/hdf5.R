"hdf5load" <-  function (file, load = TRUE, verbosity = 0, tidy = FALSE)
{
  call <- sys.call()
  .External("do_hdf5load", call, sys.frame(sys.parent()), file, load,
            as.integer (verbosity), as.logical(tidy))
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
