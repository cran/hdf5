\name{hdf5}
\title{Interface to the HDF5 Library}
\alias{hdf5save}
\alias{hdf5load}
\alias{hdf5cleanup}  % added by Ph. Grosjean

% synopsis added with hdf5cleanup to avoid warning on undocumented function
% added by Ph. Grosjean (phgrosjean@sciviews.org) - 01/03/2002
\synopsis{
hdf5save(fileout, ...)
hdf5load(file, load = TRUE, verbosity = 0, tidy = FALSE)
hdf5cleanup(fid)
}
\usage{
hdf5save(fileout, \dots)
hdf5load(file, load = TRUE, verbosity = 0, tidy = FALSE)
}

\arguments{
  \item{fileout}{the name of the file in which the objects will be stored.}
  \item{\dots}{the names of the objects to be saved.}
  % This is the original instructions for "load", but they don't
  % describe what happens.
  %  \item{load}{a logical value.  If \code{FALSE}, a character vector
  %    containing the names of the objects in the file is returned.
  %    If \code{TRUE} (the default), the objects are returned as the
  %    components of a named list.}
  \item{file}{The name of an HDF5 file to be read in.}
  \item{load}{A logical value. If \code{FALSE}, the objects are returned
    as the components of a named list.
    If \code{TRUE} (the default), the objects are loaded as individual
    variables with their own names -- the function returns nothing in
    this case.}
    \item{verbosity}{An integer controlling the verbosity. With \code{verbosity = 0}
    (the default) the hdf5 file is loaded quietly. 
    With \code{verbosity = 1} names of groups and datasets are printed as
    they are encountered -- this can be reassuring
    if you have a large file which is taking a while to load. Greater 
    values of verbosity produce more messages which are probably only useful
    for debugging.}
  \item{tidy}{A logical value. If \code{FALSE} (the default), the names of the
    \R objects will be the same as the HDF5 groups and datasets. HDF5 names
    can contain various characters which are not permissible in \R variable
    names (e.g. +, -, \_, space etc.). If \code{tidy = TRUE}, then all these
    characters are converted
    to ".". The first character will also be converted to "." if it is a
    digit. This is clearly not foolproof as "foo+" and "foo-" will both be
    tidied to "foo.".}
}
\description{
  \code{hdf5save} and \code{hdf5load} provide an experimental interface
  to the NCSA HDF5 library.
}
\details{
  \code{hdf5save} writes a
  representation of \R objects to the specified file in a form which
  can be read by software which understands the HDF5 format.
  The objects can be read back from the file at a later
  date by using the function \code{hdf5load}.

  Not all \R types are supported and it probably doesn't make sense to
  put some of them into an HDF file (e.g. closures).  However, lists,
  strings, vectors, matrices, and higher-dimensional arrays work.
  Lists map to HDF groups. Vectors and matrices map to datasets.

  This capability is only available on machines which have the HDF5
  library, version 1.2 or higher (freely available from the reference below).
}
\author{Marcus G. Daniels \email{mgd@swarm.org},
        Hugh C. Pumphrey \email{hcp@met.ed.ac.uk},
        Philippe Grosjean \email{phgrosjean@sciviews.org}}
\references{\url{http://hdf.ncsa.uiuc.edu}}
\seealso{
  \code{\link[base]{save}}, \code{\link[base]{load}}
}
\examples{
(m <- cbind(A = 1, diag(4)))
ll <- list(a=1:10, b=letters[1:8]);
l2 <- list(C="c", l=ll); PP <- pi
hdf5save("ex1.hdf", "m","PP","ll","l2")
rm(m,PP,ll,l2)  # and reload them:
hdf5load("ex1.hdf")
m        # read from "ex1.hdf"; buglet: dimnames dropped
str(ll)  # buglet: "pairlist" instead of "list"
str(l2)
}
\keyword{file}
