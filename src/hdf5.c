/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1998-2003  Marcus G. Daniels <mgd@swarm.org>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */
#define USE_RINTERNALS
#include <Rversion.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rconfig.h>
#include <R_ext/RS.h>

#include <hdf5.h>

#define STRING2REF_CONV "string->ref"
#define REF2STRING_CONV "ref->string"
/* FIXME: should really use SEXP R_RowNamesSymbol : */
#define ROWNAMES "row.names"
#define MAX_ROWNAMES_COUNT 100000

#if R_VERSION < R_Version(2, 3, 0)
#define warningcall Rf_warningcall
void warningcall (SEXP call, const char *format, ...);
#define errorcall Rf_errorcall
void errorcall (SEXP call, const char *format, ...);
#endif

/* global variables */
int hdf5_global_verbosity = 3;
int hdf5_global_nametidy = 0;
int hdf5_global_attrcnt;

#if ((H5_VERS_MAJOR > 1) || (H5_VERS_MAJOR == 1 && H5_VERS_MINOR > 6) || (H5_VERS_MAJOR == 1 && H5_VERS_MINOR == 6 && H5_VERS_RELEASE >= 3))
#define count_size_t size_t
#else
#define count_size_t hsize_t
#endif

void
nametidy (char *name)
{
  /* This function operates on a string to remove any chars regarded
     by R as not suitable for a variable name i.e. the first char must
     be alpha or . and the rest can also be digits */
  unsigned i;

  if (!isalpha (name[0]) && name[0] !='.')
    name[0] = '.';

  for (i = 1; i < strlen (name); i++)
    if(!isalnum (name[i]) && name[i] !='.')
      name[i] = '.';
}

static herr_t
ref_string (hid_t sid, hid_t did, H5T_cdata_t *cdata,
            count_size_t count, size_t stride, size_t bkg_stride,
            void *buf, void *bkg,
            hid_t dset_xfer_plid)
{
  if (cdata->command == H5T_CONV_CONV)
    {
      char *srcbuf[count];
      char *destbuf = buf;
      char **recptr = srcbuf;
      size_t i;
      size_t maxlen = H5Tget_size (did);

      memcpy (srcbuf, buf, sizeof (srcbuf));

      for (i = 0; i < count; i++)
        {
          strncpy (destbuf, *recptr, maxlen);
          recptr++;
          destbuf += maxlen;
        }
    }
  return 0;
}

static herr_t
string_ref (hid_t sid, hid_t did, H5T_cdata_t *cdata,
	    count_size_t count, size_t stride, size_t bkg_stride,
	    void *buf, void *bkg,
	    hid_t dset_xfer_plid)
{
  if (cdata->command == H5T_CONV_CONV)
    {
      size_t size = H5Tget_size (sid);
      unsigned char srcbuf[size * count], *srcptr = srcbuf;
      unsigned char cstring[size+1];
      size_t i;

      /* `count' is the number of strings in the array. size is the
         number of chars in each string.  HDF5 string arrays are not
         ragged, so there is one value of size for all the elements of
         the dataset */
      if (hdf5_global_verbosity > 1)
	Rprintf ("in string_ref: count=%d, size=%d srcbf=%d\n",
                 (int) count,(int) size, sizeof(srcbuf));
      /* Copy the whole array of strings from buf to srcbuf */
      memcpy (srcbuf, buf, sizeof (srcbuf));
      /* now copy the strings into an R string array, one at a time */
      for (i = 0; i < count; i++)
        {
	  char **ptr = &(((char **) buf)[i]);
	  /* The next lines are necessary to ensure that each
	     element of the array is null-terminated.  HDF5 files can
	     contain arrays of Fortran strings which are space-padded,
	     not null-terminated */
	  strncpy (cstring, srcptr, size);
	  cstring[size]='\0';
	  *ptr = calloc (strlen (cstring) + 1, sizeof (char));
	  if (!*ptr)
	    abort ();
	  strcpy (*ptr, cstring);
	  srcptr += size;
	}
      if (hdf5_global_verbosity > 1)
	Rprintf ("leaving string_ref\n");
    }
  return 0;
}

struct permute_info {
  SEXP call;
  int writeflag;
  SEXPTYPE type;
  unsigned rank;
  hssize_t *dims;
  hssize_t *coord;
  hid_t dataset;
  hid_t memtid;
  hid_t space;
  long tmpcnt;
  void *buf;
  void *tmpbuf;
};

static void
permute (struct permute_info *pinfo, unsigned dimnum)
{
  /* With large datasets, a great deal too much time (sometimes a
   * quasi-infinite amount!) was spent in this function. This was
   * because it was reading the dataset one element at a time. HCP
   * changed this so that the data is read in one go into tmpbuf
   * _before_ permute is called. permute then rearranges it.*/
     
  hssize_t i;

  if (dimnum < pinfo->rank)
    {

      for (i = 0; i < pinfo->dims[dimnum]; i++)
	{
	  pinfo->coord[dimnum] = i;
	  permute (pinfo, dimnum + 1);
	}

    }
  else
    {
      unsigned offset, mult;
      
      offset = pinfo->coord[0];
      mult = 1;
      for (i = 1; i < pinfo->rank; i++)
	{
	  mult *= pinfo->dims[i - 1];
	  offset += pinfo->coord[i] * mult;
	}

      switch (pinfo->type)
	{
	case STRSXP: 
	  {
	    if (pinfo->writeflag)
              {
	        char **tmpaddr;

	        tmpaddr = &((char **) pinfo->tmpbuf)[pinfo->tmpcnt];
	        *tmpaddr = CHAR (STRING_ELT (pinfo->buf, offset));
              }
	    else
	      {
		char *ptr =
		  ((char **) pinfo->tmpbuf)[pinfo->tmpcnt];
		
		SET_STRING_ELT (pinfo->buf, offset, mkChar (ptr));
	      }
	  }
	  
	  break;
	case VECSXP:
          abort ();
	  break;
	  
	case REALSXP:
	  {
	    double *pointaddr, *tmpaddr;
	    
	    pointaddr = &((double *) pinfo->buf)[offset];
	    tmpaddr = &((double *) pinfo->tmpbuf)[pinfo->tmpcnt];
	    if (pinfo->writeflag)
	      *tmpaddr = *pointaddr;
	    else
	      *pointaddr = *tmpaddr;
	  }
	  break;
	case INTSXP:
	  { 
	    int *pointaddr, *tmpaddr;

	    pointaddr = &((int *) pinfo->buf)[offset];
	    tmpaddr = &((int * ) pinfo->tmpbuf)[pinfo->tmpcnt];
	    if (pinfo->writeflag)
	      *tmpaddr = *pointaddr;
	    else
	      *pointaddr = *tmpaddr;
	  }
	  break;
	case LGLSXP:
	  {
	    int *pointaddr;
            unsigned char *tmpaddr;

            pointaddr = &((int *) pinfo->buf)[offset];
            tmpaddr = &((unsigned char *) pinfo->tmpbuf)[pinfo->tmpcnt];

	    if (pinfo->writeflag)
	      *tmpaddr = *pointaddr;
            else
              *pointaddr = *tmpaddr;
	  }
	  break;
	default:
	  errorcall (pinfo->call, "No support for R type: %d", pinfo->type);
	  /*pointaddr = &offset; unreached; -Wall */
	}
      pinfo->tmpcnt++;
      
    }
  return ;
} /* end of permute */

static hid_t
make_sexp_ref_type (SEXP call)
{
  hid_t memtid;

  if ((memtid = H5Tcopy (H5T_STD_REF_OBJ)) < 0)
    errorcall (call, "Unable to copy H5T_STD_REF_OBJ");
#if 0
  if (H5Tset_size (memtid, sizeof (SEXPREC *)) < 0)
    errorcall (call, "unable to set size of reference type");
#endif
  return memtid;
}

static hid_t
get_string_type (SEXP call, SEXP vec)
{
  hid_t stid;
  unsigned vecpos;
  size_t maxstrlen = 0;

  for (vecpos = 0; vecpos < LENGTH (vec); vecpos++)
    {
      SEXP stritem = STRING_ELT (vec, vecpos);

      if (LENGTH (stritem) > maxstrlen)
	maxstrlen = LENGTH (stritem);
    }

  if ((stid = H5Tcopy (H5T_C_S1)) < 0)
    errorcall (call, "Cannot copy string type");

  if (H5Tset_size (stid, maxstrlen + 1) < 0)
    errorcall (call, "Cannot set size of string type");

  return stid;
}

static hid_t
make_boolean_type (SEXP call)
{
  hid_t tid;

  if ((tid = H5Tcopy (H5T_NATIVE_UINT)) < 0)
    errorcall (call, "Cannot copy unsigned integer type");
  if (H5Tset_precision (tid, 1) < 0)
    errorcall (call, "Cannot set precision of boolean type");
  if (H5Tset_size (tid, 1) < 0)
    errorcall (call, "Cannot set size of boolean type");
  return tid;
}

static void
load_scalar (SEXP call, hid_t dataset, hid_t space, SEXP obj)
{
  SEXPTYPE type = TYPEOF (obj);
  void *buf;
  hid_t tid, memtid;

  if ((tid = H5Dget_type (dataset)) < 0)
    errorcall (call, "Unable to get type for dataset");

  if (H5Tget_class (tid) == H5T_COMPOUND)
    errorcall (call, "Not equipped to load compound scalar");

  switch (type) {

  case STRSXP:
    {
      size_t size = H5Tget_size (tid);

      if (sizeof (char *) > size)
	size = sizeof (char *);

      memtid = make_sexp_ref_type (call);
      buf = R_chk_calloc (1, size * 2);
    }
    break;
  case REALSXP:
    memtid = H5T_NATIVE_DOUBLE;
    buf = REAL (obj);
    break;
  case INTSXP:
    memtid = H5T_NATIVE_INT;
    buf = INTEGER (obj);
    break;
  case LGLSXP:
    memtid = make_boolean_type (call);
    buf = LOGICAL (obj);
    break;
  default:
      errorcall (call, "Can't get type for R type: %d (IO)", type);
      /* unreached (-Wall): */ memtid = tid;  buf = &tid;
  }
  
  if (H5Dread (dataset,
               memtid,
               space,
               space,
               H5P_DEFAULT,
               buf) < 0)
    errorcall (call, "Unable to read dataset");

  if (type == STRSXP)
    {
      SET_STRING_ELT (obj, 0, mkChar (buf));
      R_chk_free (buf);
    }

  if (type == LGLSXP)
    if (H5Tclose (memtid) < 0)
      errorcall (call, "can't close boolean type");
}

static void
vector_io (SEXP call, int writeflag, hid_t dataset, hid_t space, SEXP obj)
{
  int rank = H5Sget_simple_extent_ndims (space);
  hsize_t dims[rank], maxdims[rank];
  /* WARNING: Don't make bufsize a hsize_t. If you do this,
     then bufsize=-1 ; if(bufsize>0) is (ludicrously) regarded as
     true. We therefore make bufsize a long and cast it to a hsize_t
     when we pass it to a HDF5 routine. There are known weirdnesses
     regarding hsize_t and gcc which may be a compiler bug and this is
     probably another example.  */
  long bufsize = 0;
  SEXPTYPE type = TYPEOF (obj);
  hid_t memtid, tid, plist;
  void *buf;
  void *tmpbuf = NULL;
  long n_elements,i;
  if (hdf5_global_verbosity > 3)
    Rprintf ("in vector_io: rank=%d\n",rank); 
  if ((tid = H5Dget_type (dataset)) < 0)
    errorcall (call, "Unable to get type for dataset");
  
  if (H5Sget_simple_extent_dims (space, dims, maxdims) < 0)
    errorcall (call, "Unable to get dimensions of space");

  /* calculate the total number of elements in the dataset */
  n_elements = 1;
  for (i = 0; i < rank; i++)
    {
      if (hdf5_global_verbosity > 3)
        Rprintf ("in vector_io:size %d = %d into n_elements..",i,dims[i]); 
      n_elements = n_elements * dims[i];
      if (hdf5_global_verbosity > 3)
        Rprintf ("....=%d\n ",n_elements);
    }

  switch (type) {
  case STRSXP:
    { 
      size_t size;

      memtid = make_sexp_ref_type (call);
      buf = obj;
      if (!writeflag)
	{
	  size = H5Tget_size (tid);
	  if (sizeof (char *) > size)
	    size = sizeof (char *);
	  bufsize = n_elements * size;

          /* this accomplishes two things:
	     1) it avoids buffer-too-small errors from HDF5 for small strings
	   
	     2) it ensures that the string set is read in one pass, and
	        the string vector fully populated (as opposed to several reads
	        that fill only half of it.
          */
	  bufsize *= 2;

	  tmpbuf = Calloc (bufsize, char);
	}
      else
	{
	  tmpbuf = Calloc (n_elements, char *);
	  bufsize = -1;
	}
	  
      if (!tmpbuf)
	abort ();
    }
    break;
  case REALSXP:
    memtid = H5T_NATIVE_DOUBLE;
    buf = REAL (obj);
    tmpbuf = R_alloc (n_elements, sizeof (double));
    bufsize = n_elements * sizeof (double);
    break;
  case INTSXP:
    memtid = H5T_NATIVE_INT;
    buf = INTEGER (obj);
    tmpbuf = R_alloc (n_elements, sizeof (int));
    bufsize = n_elements * sizeof (int);
    break;
  case LGLSXP:
    memtid = make_boolean_type (call);
    buf = LOGICAL (obj);
    tmpbuf = R_alloc (n_elements, sizeof (unsigned char));
    bufsize = n_elements * sizeof (unsigned char);
    break;
  default:
    errorcall (call, "Can't get type for R type: %d (IO)", type);
    /* unreached (-Wall): */ memtid = tid;  buf = &tid;
  }

  /* set type conversion buffer size in property list . In most circs, 
     you don't have to do this, but if you have an array with dims
     a x b x c x d x .... and b x c x d x ..... adds up to > 1MB
     then HDF5 will barf. */

  if (bufsize > 0)
    {
      if (hdf5_global_verbosity > 2)
        Rprintf ("Setting buffer size in plist\n");
      plist = H5Pcreate(H5P_DATASET_XFER);

    if (H5Pset_buffer (plist, ((hsize_t) bufsize), NULL, NULL) < 0)
      errorcall (call,"Unable to set buffer size in property list");
    }
  else
    {
      if (hdf5_global_verbosity > 2)
        Rprintf ("Using default transfer plist\n");
      plist = H5P_DEFAULT;
    }
  
  /* Read data into temporary buffer before permuting it */
  if (!writeflag)
    {
      if (hdf5_global_verbosity > 2)
        Rprintf ("About to read with bufsize = %d\n",bufsize);
      
      if (H5Dread (dataset, memtid, H5S_ALL, H5S_ALL, plist, tmpbuf) < 0)
        errorcall (call, "Unable to read dataset");
      if (hdf5_global_verbosity > 2)
        Rprintf (" Done read\n");
    }
  
  if (hdf5_global_verbosity > 2)
    Rprintf ("in vector_io: permuting\n");  
  
  {
    struct permute_info pinfo;
    hssize_t coord[rank];

    pinfo.call = call;
    pinfo.writeflag = writeflag;
    pinfo.type = type;
    pinfo.rank = rank;
    pinfo.coord = coord;
    pinfo.dims = dims;
    pinfo.dataset = dataset;
    pinfo.memtid = memtid;
    pinfo.space = space;
    pinfo.buf = buf;
    pinfo.tmpbuf = tmpbuf;
    pinfo.tmpcnt = 0;

    /* The grinding slowness happened inside this call to permute()
       and was caused quite specifically by calling H5Dread() with a
       dataset size of 1 for every element of the dataset. The fix I
       have applied is to have an o/p buffer tmpbuf as big as obj that
       we do the read/write to/from and use permute to transfer data
       from tmpbuf to buf. If we read directly into buf the
       indices get scrambled. Whatever permute does to them is
       necessary.*/

    permute (&pinfo, 0);

    if (writeflag)
      {
        if (hdf5_global_verbosity > 2)
          Rprintf ("About to write\n");
        if (H5Dwrite (dataset, memtid, H5S_ALL, H5S_ALL, plist, tmpbuf) < 0)
          errorcall (call, "Unable to write dataset");
        if (hdf5_global_verbosity > 2)
          Rprintf ("About to write\n");
      }
    else
      {
	if (type == STRSXP)
	  Free (tmpbuf);
      }
  }
  
  if (hdf5_global_verbosity > 2)
    Rprintf ("in vector_io: tidying\n");  
  if (bufsize > 0)
    {
      if (H5Pclose(plist) < 0)	
        errorcall (call, "Unable to close plist");
    }

  if (type == STRSXP || type == LGLSXP)
    {
      if (H5Tclose (memtid) < 0)
	errorcall (call, "Unable to close reference type");
    }
}

static void
hdf5_save_attributes (SEXP call, hid_t loc_id, SEXP val)
{
  SEXP l;

  for (l = ATTRIB (val); l != R_NilValue; l = CDR (l))
    {
      SEXP attr = CAR (l);
      SEXPTYPE type = TYPEOF (attr);
      const char *name = CHAR (PRINTNAME (TAG (l)));
      void *buf;
      hid_t tid, memtid;
      hid_t sid, aid;
      unsigned count = LENGTH (attr);
      /*SEXPREC *stringptrs[count];*/

      if (TAG (l) == R_RowNamesSymbol
	  || TAG (l) == R_ClassSymbol
	  || TAG (l) == R_NamesSymbol
	  || TAG (l) == R_DimNamesSymbol)
	continue;
      {
	hsize_t dims[1];

	dims[0] = count;

	if ((sid = H5Screate_simple (1, dims, NULL)) < 0)
	  errorcall (call,
		     "unable to create vector space for attribute `%s'", name);
      }

      if (type == STRSXP)
	{
	  unsigned i;
	  memtid = make_sexp_ref_type (call);
	  tid = get_string_type (call, attr);
	  buf = Calloc (count, const char *);
	  for (i = 0; i < count; i++)
	    ((const char **) buf)[i] = CHAR (STRING_ELT (attr, i));
	}
      else if (type == LGLSXP)
	{
	  memtid = make_boolean_type (call);
          tid = make_boolean_type (call);
	  buf = LOGICAL (attr);
	}
     else if (type == INTSXP)
	{
	  memtid = H5T_NATIVE_INT;
	  tid = H5T_NATIVE_INT;
	  buf = INTEGER (attr);
	}
      else if (type == REALSXP)
	{
	  memtid = H5T_NATIVE_DOUBLE;
	  tid = H5T_NATIVE_DOUBLE;
	  buf = REAL (attr);
	}
      else
	abort ();

      if ((aid = H5Acreate (loc_id, name, tid, sid, H5P_DEFAULT)) < 0)
	errorcall (call, "unable to create attribute `%s'", name);

      if (H5Awrite (aid, memtid, buf) < 0)
	errorcall (call, "unable to write attribute `%s'", name);

      if (H5Aclose (aid) < 0)
	errorcall (call, "unable to close attribute `%s'", name);

      if (type == STRSXP || type == LGLSXP)
	{
	  if (type == STRSXP)
	    Free (buf);
	  if (H5Tclose (memtid) < 0)
	    errorcall (call,
		       "unable to close string reference type `%s'",
		       name);
	  if (H5Tclose (tid) < 0)
	    errorcall (call, "unable to close output type `%s'", name);
	}
      if (H5Sclose (sid) < 0)
	errorcall (call, "unable to close space for attribute `%s'", name);
    }
}

static void
hdf5_write_vector (SEXP call, hid_t id, const char *symname, SEXP val)
{
  unsigned i, rank;
  SEXP dimvec;
  hid_t space, dataset;
  SEXPTYPE type = TYPEOF (val);
  hid_t tid;

  dimvec = getAttrib (val, R_DimSymbol);
  rank = (dimvec == R_NilValue) ? 1 : LENGTH (dimvec);

  {
    hsize_t dims[rank];

    if (rank > 1)
      for (i = 0; i < rank; i++)
	dims[i] = INTEGER (dimvec)[i];
    else
      dims[0] = length(val);

    if ((space = H5Screate_simple (rank, dims, NULL)) < 0)
      errorcall (call, "Unable to create file dataspace");

    if (type == STRSXP)
      tid = get_string_type (call, val);
    else if (type == LGLSXP)
      tid = make_boolean_type (call);
    else if (type == INTSXP)
      tid = H5T_NATIVE_INT;
    else if (type == REALSXP)
      tid = H5T_NATIVE_DOUBLE;
    else
      {
        errorcall (call, "Can't get type for R type: %d (Creating)", type);
        tid = H5T_NATIVE_INT;/*unreached; -Wall*/
      }

    if ((dataset = H5Dcreate (id,
			      symname,
			      tid,
			      space,
			      H5P_DEFAULT)) < 0)
      errorcall (call, "Unable to create dataset");

    vector_io (call, TRUE, dataset, space, val);
    hdf5_save_attributes (call, dataset, val);

    if (type == LGLSXP || type == STRSXP)
      if (H5Tclose (tid) < 0)
	errorcall (call, "Unable to close type");

    if (H5Dclose (dataset) < 0)
      errorcall (call, "Unable to close dataset");
    if (H5Sclose (space) < 0)
      errorcall (call, "Unable to close space");
  }
}

static void
hdf5_write_string (SEXP call, hid_t fid, const char *symname, const char *str)
{
  hid_t stringtype;
  hid_t dataset;
  hid_t dataspace;

  dataspace = H5Screate (H5S_SCALAR);

  stringtype = H5Tcopy (H5T_C_S1);
  H5Tset_size (stringtype, strlen (str) + 1);

  if ((dataset = H5Dcreate (fid,
			    symname,
			    stringtype,
			    dataspace,
			    H5P_DEFAULT)) < 0)
    errorcall (call, "Unable to create dataset");

  if (H5Dwrite (dataset,
		stringtype,
		H5S_ALL,
		H5S_ALL,
		H5P_DEFAULT,
		str) < 0)
    errorcall (call, "Unable to write dataset");

  H5Dclose (dataset);
  H5Sclose (dataspace);
  H5Tclose (stringtype);
}

static unsigned
align (unsigned offset, unsigned alignto)
{
#if 0
  unsigned mask = alignto - 1;

  if ((offset & mask) == 0)
    return offset;
  else
    return (offset + alignto) & ~mask;
#else
  return offset;
#endif
}

static void
create_rownames_dataset_attribute (SEXP call, hid_t dataset, SEXP rownames)
{
  hid_t stringtid = get_string_type (call, rownames);
  hid_t rtid = make_sexp_ref_type (call);
  hid_t rnattrib, rndataspace;
  unsigned rowcount = LENGTH (rownames);
  hsize_t dims[1];
  const char **buf = Calloc (rowcount, const char *);
  unsigned i;

  dims[0] = rowcount;

  for (i = 0; i < rowcount; i++)
    buf[i] = CHAR (STRING_ELT (rownames, i));

  if ((rndataspace = H5Screate_simple (1, dims, NULL)) < 0)
    errorcall (call, "Unable to create row names vector space");

  if ((rnattrib = H5Acreate (dataset, ROWNAMES,
			     stringtid, rndataspace, H5P_DEFAULT)) < 0)
    errorcall (call, "unable to create row names dataset");

  if (H5Awrite (rnattrib, rtid, buf) < 0)
    errorcall (call, "unable to write row names dataset");

  if (H5Aclose (rnattrib) < 0)
    errorcall (call, "unable to close row names dataset");

  if (H5Sclose (rndataspace) < 0)
    errorcall (call, "unable to close row names dataspace");

  if (H5Tclose (stringtid) < 0)
    errorcall (call, "unable to close row names string type");

  if (H5Tclose (rtid) < 0)
    errorcall (call, "unable to close reference type");

  Free (buf);
}

static void
hdf5_save_object (SEXP call, hid_t fid, const char *symname, SEXP val)
{
  if (isFrame (val))
    {
      unsigned colcount = length (val), pos;
      size_t offsets[colcount];
      hid_t hdftypes[colcount];
      hid_t hdfbooltype;
      size_t size = 0;

      hdfbooltype = make_boolean_type (call);

      for (pos = 0; pos < colcount; pos++)
	{
	  SEXPTYPE type = TYPEOF (VECTOR_ELT (val, pos));

	  switch (type)
	    {
	    case REALSXP:
	      hdftypes[pos] = H5T_NATIVE_DOUBLE;
	      break;
	    case INTSXP:
	      hdftypes[pos] = H5T_NATIVE_INT;
	      break;
	    case STRSXP:
	      hdftypes[pos] = get_string_type (call, VECTOR_ELT (val, pos));
	      break;
	    case LGLSXP:
	      hdftypes[pos] = make_boolean_type (call);
	      break;
	    default:
	      errorcall (call,
			 "No support for converting R type: %d to HDF5",
			 type);
	      break;
	    }
	  if (H5Tget_class (hdftypes[pos]) == H5T_STRING)
	    size = align (size, 8);
	  else
	    size = align (size, H5Tget_size (hdftypes[pos]));
	  offsets[pos] = size;
	  size += H5Tget_size (hdftypes[pos]);
	}
      size = align (size, H5Tget_size (hdftypes[0]));
      {
	hid_t ctid;
	hid_t dataset;
	hid_t dataspace;
	unsigned rowcount = length (VECTOR_ELT (val, 0));
	{
	  hsize_t dims[1];

	  dims[0] = rowcount;
	  if ((dataspace = H5Screate_simple (1, dims, NULL)) < 0)
	    errorcall (call, "Unable to create dataframe vector space");
	}

	{
	  SEXP colnames = getAttrib (val, R_NamesSymbol);

	  if ((ctid = H5Tcreate (H5T_COMPOUND, size)) < 0)
	    errorcall (call, "unable to create compound type");

	  for (pos = 0; pos < colcount; pos++)
	    if (H5Tinsert (ctid,
			   CHAR (STRING_ELT (colnames, pos)),
			   offsets[pos],
			   hdftypes[pos]) < 0)
	      errorcall (call, "unable to insert type into compound type");
	}
	if (H5Tpack (ctid) < 0)
	  errorcall (call, "Unable to pack type");

	if (H5Tlock (ctid) < 0)
	  errorcall (call, "Unable to lock type");

	if ((dataset = H5Dcreate (fid, symname,
				  ctid, dataspace, H5P_DEFAULT)) < 0)
	  errorcall (call, "unable to create dataframe dataset");

	{
	  unsigned ri;
	  unsigned char buf[rowcount][size];

	  for (ri = 0; ri < rowcount; ri++)
	    for (pos = 0; pos < colcount; pos++)
	      {
		SEXP item = VECTOR_ELT (val, pos);
		SEXPTYPE type = TYPEOF (item);
		void *ptr = &buf[ri][offsets[pos]];

		switch (type)
		  {
		  case REALSXP:
		    memcpy (ptr, &REAL (item)[ri], sizeof (double));
		    break;
		  case INTSXP:
		    memcpy (ptr, &INTEGER (item)[ri], sizeof (int));
		    break;
		  case STRSXP:
		    {
		      SEXP stritem = STRING_ELT (item, ri);

		      memset (ptr, 0, H5Tget_size (hdftypes[pos]));
		      strcpy ((char *)ptr, CHAR (stritem));
		    }
		    break;
		  case LGLSXP:
		    *(unsigned char *) ptr = LOGICAL (item)[ri];
		    break;
		  default:
		    abort ();
		  }
	      }
	  if (H5Dwrite (dataset,
			ctid,
			dataspace,
			dataspace,
			H5P_DEFAULT,
			buf) < 0)
	    errorcall (call, "Unable to write dataframe");
	}
	hdf5_save_attributes (call, dataset, val);
        if (rowcount < MAX_ROWNAMES_COUNT || MAX_ROWNAMES_COUNT == 0)
	  {
	    SEXP rownames = getAttrib (val, R_RowNamesSymbol);

	    if (rownames != R_NilValue)
	      create_rownames_dataset_attribute (call, dataset, rownames);
          }
	if (H5Dclose (dataset) < 0)
	  errorcall (call, "Cannot close dataset");
	if (H5Sclose (dataspace) < 0)
	  errorcall (call, "Cannot close dataspace");
      }

      for (pos = 0; pos < colcount; pos++)
	if (H5Tget_class (hdftypes[pos]) == H5T_STRING)
	  if (H5Tclose (hdftypes[pos]) < 0)
	    errorcall (call, "Cannot close string type");
      if (H5Tclose (hdfbooltype) < 0)
	errorcall (call, "Cannot close boolean type");

    }
  else if (isNull (val))
    {
      hid_t gid;

      if ((gid = H5Gcreate (fid, symname, 0)) < 0)
        errorcall (call, "unable to create group");
      
      if (H5Gclose (gid) < 0)
        errorcall (call, "unable to close group");
    }
  else
    {
      SEXPTYPE type = TYPEOF (val);

      switch (type)
	{
	case LGLSXP: case INTSXP: case REALSXP: case STRSXP:
	  hdf5_write_vector (call, fid, symname, val);
	  break;
	case LISTSXP: case VECSXP:
	  {
	    unsigned len = length (val);
	    hid_t gid;
	    unsigned pos;
	    char buf[(sizeof (pos) * 8 / 3 + 1) + 1];

	    if ((gid = H5Gcreate (fid, symname, len * 8)) < 0)
	      errorcall (call, "unable to create group");

	    if (type == LISTSXP)
	      {
		SEXP l;

		for (l = val, pos = 0; l != R_NilValue; l = CDR (l), pos++)
		  {
		    SEXP s = CAR (l);

		    if (!isNull (TAG (l)))
		      hdf5_save_object (call, gid,
				       CHAR (PRINTNAME (TAG (l))), s);
		    else
		      {
			sprintf (buf, "%u", pos);
			hdf5_save_object (call, gid, buf, s);
		      }
		  }
	      }
	    else
	      {
		for (pos = 0; pos < len; pos++)
		  {
		    SEXP s = VECTOR_ELT (val, pos);
		    SEXP names = getAttrib (val, R_NamesSymbol);

		    if (!isNull (names))
		      hdf5_save_object (call, gid,
				       CHAR (STRING_ELT (names, pos)),
				       s);
		    else
		      {
			sprintf (buf, "%u", pos);
			hdf5_save_object (call, gid, buf, s);
		      }
		  }
	      }
	    hdf5_save_attributes (call, gid, val);

	    if (H5Gclose (gid) < 0)
	      errorcall (call, "unable to close group");
	  }
	  break;
	case SYMSXP:
	  {
	    const char *pn = CHAR (PRINTNAME (val));

	    hdf5_write_string (call, fid, symname, pn);
	  }
	  break;
	default:
	  errorcall (call, "unhandled type: %d", type);
	  break;
	}
    }
}

void
setup_onexit (hid_t fid, SEXP env)
{
  eval (lang2 (install ("on.exit"),
               lang2 (install ("hdf5cleanup"),
                      ScalarInteger (fid))),
        env);
}

SEXP
do_hdf5cleanup (SEXP args)
{
  SEXP call, env;
  hid_t fid;

  args = CDR (args); call = CAR (args);
  args = CDR (args); env = CAR (args);
  args = CDR (args);

  if (TYPEOF (CAR (args)) != INTSXP)
    abort ();

  fid = (hid_t) INTEGER (CAR (args))[0];

  H5Tunregister (H5T_PERS_SOFT, STRING2REF_CONV, -1, -1, string_ref);
  H5Tunregister (H5T_PERS_SOFT, REF2STRING_CONV, -1, -1, ref_string);
  
  if (H5Fclose (fid) < 0)
    errorcall (call, "unable to close HDF file");
  return R_NilValue;
}

SEXP
do_hdf5save (SEXP args)
{
  int i, nobjs;
  const char *path; 
  hid_t fid;
  SEXP s, env, call, names, sym;

  args = CDR (args); call = CAR (args);
  args = CDR (args); env = CAR (args);
  args = CDR (args); 

  if (TYPEOF (CAR (args)) != STRSXP)
    errorcall (call, "first argument must be a pathname");

  path = CHAR (STRING_ELT (CAR (args), 0));

  H5dont_atexit ();

  if (H5Tregister (H5T_PERS_SOFT,
		   REF2STRING_CONV,
		   H5T_STD_REF_OBJ,
		   H5T_C_S1, ref_string) < 0)
    errorcall (call, "Unable to register ref->string converter");

  if ((fid = H5Fcreate (path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)) < 0)
    errorcall (call, "unable to create HDF file: %s", path);

  setup_onexit (fid, env);

  names = args = CDR (args);

  nobjs = length (names);
  if (nobjs < 1)
    errorcall (call, "no objects to save");

  for (i = 0; i < nobjs; i++)
    {
      if (TYPEOF (CAR (names)) != STRSXP)
        errorcall (call, "expecting a symbol name");
      {
        char *name = CHAR (STRING_ELT (CAR (names), 0));
        
        PROTECT (sym = install (name));
        PROTECT (s = findVar (sym, env));
        if (s == R_UnboundValue)
          errorcall (call, "symbol `%s' has no value", name);
        hdf5_save_object(call, fid, name, s);
        UNPROTECT(2);
        names = CDR (names);
      }
    }

  return names;
}

struct hdf5_iterate_info {
  SEXP call;
  void (*add) (struct hdf5_iterate_info *, const char *, SEXP);
  SEXP env;
  SEXP ret;
};

static void
add_to_list (struct hdf5_iterate_info *iinfo, const char *name, SEXP obj)
{
  SEXP pair;
  char newname[strlen(name)+1];

  strcpy (newname,name);
  if (hdf5_global_nametidy)
    {
      if (hdf5_global_verbosity > 1)
        Rprintf(" Tidying name %s ", newname);
      nametidy (newname);
    }
  PROTECT (pair = CONS (obj, CDR (iinfo->ret)));
  TAG (pair) = install ((char *) newname);
  SETCDR (iinfo->ret, pair);
  UNPROTECT (1);
  if (hdf5_global_verbosity > 1 && hdf5_global_nametidy) 
    Rprintf (".. to %s \n ",newname);
  if (hdf5_global_verbosity > 2)
    Rprintf ("Adding `%s' to list\n", newname);
}


static SEXP
collect (SEXP call, hid_t id, H5G_iterate_t iterate_func, SEXP env)
{
  struct hdf5_iterate_info iinfo;

  iinfo.call = call;
  iinfo.add = add_to_list;
  PROTECT (iinfo.ret = CONS (R_NilValue, R_NilValue));
  iinfo.env = env;

  if (H5Giterate (id, ".", NULL, iterate_func, &iinfo) < 0)
    errorcall (call, "unable to collect HDF group");

  UNPROTECT (1);
  return CDR (iinfo.ret);
}

static void
load_rownames_dataset_attribute (SEXP call, hid_t dataset, SEXP vec)
{
  hid_t rnattrib, rnspace, rntid;
  unsigned rowcount;
  SEXP rownames;
  H5E_auto_t errfunc;
  void *client_data;

  H5Eset_auto (NULL, NULL);
  rnattrib = H5Aopen_name (dataset, ROWNAMES);
  H5Eget_auto (&errfunc, &client_data);
  H5Eset_auto (errfunc, client_data);

  if (rnattrib < 0)
    {
      hid_t space;
      unsigned rank;

      if ((space = H5Dget_space (dataset)) < 0)
	errorcall (call, "unable to get dataset space");

      if (H5Sis_simple (space) != TRUE)
	errorcall (call, "space not simple");

      if ((rank = H5Sget_simple_extent_ndims (space)) < 0)
	errorcall (call, "unable to get space rank");

      if (rank != 1)
        errorcall (call, "rank not 1");

      {
	hsize_t dims[rank];
	hsize_t maxdims[rank];
        unsigned i;

	if (H5Sget_simple_extent_dims (space, dims, maxdims) < 0)
	  errorcall (call, "unable to get space extent");

        rowcount = dims[0];
        PROTECT (rownames = allocVector (STRSXP, rowcount));

        for (i = 0; i < rowcount; i++)
          {
            char buf[256];
 
            sprintf (buf, "%u", i);
            SET_STRING_ELT (rownames, i, mkChar (buf));
          }
        setAttrib (vec, R_RowNamesSymbol, rownames);
      }      
      UNPROTECT (1);

      if (H5Sclose (space) < 0)
	errorcall (call, "unable to close dataset space");

      return;
    }

  if ((rnspace = H5Aget_space (rnattrib)) < 0)
    errorcall (call, "could not get space for rownames attribute");

  if ((rntid = H5Aget_type (rnattrib)) < 0)
    errorcall (call, "could not get element type of rownames attribute");

  if (H5Sget_simple_extent_ndims (rnspace) != 1)
    errorcall (call, "rownames space should be of rank 1");

  {
    hsize_t dims[1], maxdims[1];

    if (H5Sget_simple_extent_dims (rnspace, dims, maxdims)< 0)
      errorcall (call, "can't get attribute space dims");
    rowcount = dims[0];
  }
  PROTECT (rownames = allocVector (STRSXP, rowcount));
  {
    size_t insize = H5Tget_size (rntid);
    size_t outsize = sizeof (const char *);
    size_t size = (insize > outsize ? insize: outsize);

    const char **strptrs = (const char **) R_chk_calloc (rowcount, size * 2);
    unsigned ri;
    hid_t rtid = make_sexp_ref_type (call);

    if (H5Aread (rnattrib, rtid, strptrs) < 0)
      errorcall (call, "can't read rownames");

    for (ri = 0; ri < rowcount; ri++)
      SET_STRING_ELT (rownames, ri, mkChar (strptrs [ri]));

    if (H5Tclose (rtid) < 0)
      errorcall (call, "can't close close reference type");
 
    R_chk_free (strptrs);
  }
  setAttrib (vec, R_RowNamesSymbol, rownames);
  UNPROTECT (1);

  if (H5Sclose (rnspace) < 0)
    errorcall (call, "unable to close row name space");
  
  if (H5Tclose (rntid) < 0)
    errorcall (call, "unable to close row name type");
    
  if (H5Aclose (rnattrib) < 0)
    errorcall (call, "unable to close row name attribute");
}

struct hdf5_attribute_info {
  SEXP call;
  SEXP obj;
  const char *name;
};

/*static*/ 
herr_t
hdf5_process_attribute (hid_t loc_id, const char *attrName, void *data)
{
  struct hdf5_attribute_info *ainfo = data;
  hid_t aid, sid, tid;
  H5T_class_t class;
  size_t tid_size;
  char newname[strlen(attrName) +1];

  if (strcmp (attrName, ROWNAMES) == 0)
    {
      if (hdf5_global_verbosity > 1)
        Rprintf ("Skipping attribute %s\n", attrName);
        
        return 0;
    }

  hdf5_global_attrcnt++;

  if (hdf5_global_verbosity > 1)
    Rprintf ("Processing attribute %d called %s\n",
	     hdf5_global_attrcnt,attrName);

  if ((aid = H5Aopen_name (loc_id, attrName)) < 0)
    errorcall (ainfo->call, "could not open attribute `%s'", attrName);

  if ((sid = H5Aget_space (aid)) < 0)
    errorcall (ainfo->call, "could not open space of attribute `%s'",
	       attrName);

  if ((tid = H5Aget_type (aid)) < 0)
    errorcall (ainfo->call, "could not get type of attribute `%s'", attrName);

  if ((tid_size = H5Tget_size (tid)) < 0)
    errorcall (ainfo->call, "could not get size of attribute `%s' tid",
	       attrName);

  if ((class = H5Tget_class (tid)) < 0)
    errorcall (ainfo->call, "could not get type class of attribute `%s'",
	       attrName);

  {
    int rank,r_rank;
    
    if ((rank = H5Sget_simple_extent_ndims (sid)) < 0)
      errorcall (ainfo->call,
		 "could not get rank of attribute space `%s'",
		 attrName);
    
    if (hdf5_global_verbosity > 1) 
      Rprintf ("attribute %s has rank %d \n",attrName,rank);    

    if (rank == 0)
      r_rank = 1;
    else
      r_rank = rank;

    {
      hsize_t dims[r_rank];

      if(rank == 1)
        {
          if (H5Sget_simple_extent_dims (sid, dims, NULL) < 0)
            errorcall (ainfo->call,
                       "could not get extent of attribute space `%s'",
                       attrName);
        }
      else
        {
          /* scalar, but we'll pretend it is a vector of length 1 */
          dims[0] = 1; 
          if (hdf5_global_verbosity > 2) 
            Rprintf ("Rank 0 attribute treated as rank 1 size 1\n");
        }
      if (r_rank == 1)
	{
	  unsigned count = dims[0];
	  SEXPTYPE type;
	  hid_t memtid;
	  SEXP vec;
	  void *buf;

	  switch (class)
	    {
	    case H5T_INTEGER:
	      if (tid_size == 1)
                {
                  memtid = make_boolean_type (ainfo->call);
                  type = LGLSXP;
                }
              else
                {
                  memtid = H5Tcopy (H5T_NATIVE_INT);
                  type = INTSXP;
                }
	      break;
	    case H5T_FLOAT:
	      type = REALSXP;
	      memtid = H5Tcopy (H5T_NATIVE_DOUBLE);
	      break;
	    case H5T_STRING:
	      type = STRSXP;
	      if (hdf5_global_verbosity > 2) 
		Rprintf ("Attribute is a string\n");
	      memtid = make_sexp_ref_type (ainfo->call);
	      break;
	    default:
	      warningcall (ainfo->call, "skipping attribute `%s' due to type",
			   attrName);
	      goto done;
	    }
	  PROTECT (vec = allocVector (type, count));
	  switch (class)
	    {
	    case H5T_INTEGER:
              if (tid_size == 1)
                buf = LOGICAL (vec);
              else
                buf = INTEGER (vec);
	      break;
	    case H5T_FLOAT:
	      buf = REAL (vec);
	      break;
	    case H5T_STRING:
              {
                const char **ptr;
		size_t insize = H5Tget_size (tid);
		size_t outsize = sizeof (const char *);
		size_t size = insize > outsize ? insize : outsize;
 
	        ptr = (const char **) R_chk_calloc (count, size * 2);

                buf = ptr;
              }
	      break;
	    default:
	      abort ();
	    }
	  if (H5Aread (aid, memtid, buf) < 0)
	    errorcall (ainfo->call, "unable to read attribute `%s'", attrName);
          if (class == H5T_STRING)
            {
              const char **ptr = buf;
              unsigned i;

              for (i = 0; i < count; i++)
                SET_STRING_ELT (vec, i, mkChar (ptr[i]));
              R_chk_free (buf);
            }
	  if (hdf5_global_verbosity > 2) 
	    Rprintf ("string length of new name =%d\n",strlen(attrName)+1);
	  strcpy(newname,attrName);
	  if (hdf5_global_nametidy)
            {
              if (hdf5_global_verbosity > 1) 
                Rprintf (" Tidying attribute name %s ",newname);
              nametidy (newname);
              if (hdf5_global_verbosity > 1) 
                Rprintf("....to %s\n",newname);
            }
	  
	  if (!isNull (ainfo->obj))
	    setAttrib (ainfo->obj, install ((char *) newname), vec);
	  
	  UNPROTECT (1);

	  if (H5Tclose (memtid) < 0)
	    errorcall (ainfo->call,
		       "unable to close reference type in attribute `%s'",
		       attrName);
	}
      else
	warningcall (ainfo->call, "skipping attribute `%s' due to rank",
		     attrName);
    }
  }
 done:
  if (H5Sclose (sid) < 0)
    errorcall (ainfo->call, "unable to close attribute `%s' space", attrName);
  if (H5Tclose (tid) < 0)
    errorcall (ainfo->call, "unable to close attribute `%s' type", attrName);
  if (H5Aclose (aid) < 0)
    errorcall (ainfo->call, "unable to close attribute `%s'", attrName);
  if (hdf5_global_verbosity > 1) 
    Rprintf ("Done processing attribute %s\n",attrName);
  if (hdf5_global_attrcnt > 100)
    {
      Rprintf ("WTF? More than 100 attributes? \n");
      return ((herr_t) 99);
      /* and here we sometimes get a seggie. Why?*/
    }
  else
    return 0;
}

static void
hdf5_load_attributes (SEXP call, hid_t id, SEXP obj, const char *name)
{
  unsigned *idx = NULL;
  struct hdf5_attribute_info ainfo;
  herr_t retval;

  ainfo.call = call;
  ainfo.obj = obj;
  ainfo.name = name;
  hdf5_global_attrcnt = 0;
  retval = H5Aiterate (id, idx, hdf5_process_attribute, &ainfo);
  if (retval < 0)
    errorcall (call, "unable to iterate over attributes");
}

static herr_t
hdf5_process_object (hid_t id, const char *name, void *client_data)
{
  struct hdf5_iterate_info *iinfo = client_data;
  H5G_stat_t statbuf;

  if (hdf5_global_verbosity > 0)
    Rprintf ("Processing object: %s ...",name);
  if (H5Gget_objinfo (id, name, 1, &statbuf) < 0)
    errorcall (iinfo->call, "Cannot query object `%s'", name);

  if (statbuf.type == H5G_GROUP)
    {
      SEXP l;
      hid_t gid = H5Gopen (id, name);

      if (hdf5_global_verbosity > 0)
        Rprintf ("... which is a Group \n");
      if (gid < 0)
	errorcall (iinfo->call, "unable to open group `%s'", name);

      PROTECT (l = PairToVectorList (collect (iinfo->call, gid, hdf5_process_object, iinfo->env)));

      if (hdf5_global_verbosity > 2)
        Rprintf ("Adding `%s'\n", name);
      iinfo->add (iinfo, name, l);

      hdf5_load_attributes (iinfo->call, gid, l, name);
      UNPROTECT (1);

      if (H5Gclose (gid) < 0)
	errorcall (iinfo->call, "unable to close group");
      if (hdf5_global_verbosity > 0)
        Rprintf ("... Done group %s \n",name);      
    }
  else if (statbuf.type == H5G_DATASET)
    {
      hid_t dataset, space, tid;
      int rank;
      SEXPTYPE type = NILSXP;

      if (hdf5_global_verbosity > 0)
        Rprintf ("... its a dataset...");

      if ((dataset = H5Dopen (id, name)) < 0)
	errorcall (iinfo->call, "unable to load dataset `%s'", name);

      if (hdf5_global_verbosity > 1)
        Rprintf ("Dataset has ID%d\n",dataset);

      if ((tid = H5Dget_type (dataset)) < 0)
	errorcall (iinfo->call, "unable to get dataset type");

      if (hdf5_global_verbosity > 1)
        Rprintf ("Dataset has tid %d\n",tid);

      switch (H5Tget_class (tid))
	{
	case H5T_INTEGER:
	  if (H5Tget_precision (tid) == 1)
	    type = LGLSXP;
	  else
	    type = INTSXP;
	  break;
	case H5T_FLOAT:
	  type = REALSXP;
	  break;
	case H5T_STRING:
	  type = STRSXP;
	  break;
	case H5T_COMPOUND:
	  type = VECSXP;
	  break;
	default:
	  errorcall (iinfo->call, "can't handle hdf type %d", tid);
	  break;
	}

      if ((space = H5Dget_space (dataset)) < 0)
	errorcall (iinfo->call, "unable to get dataset space");

      if (hdf5_global_verbosity > 1)
        Rprintf ("Dataset has space id %d\n", space);

      if (H5Sis_simple (space) != TRUE)
	errorcall (iinfo->call, "space not simple");

      if ((rank = H5Sget_simple_extent_ndims (space)) < 0)
	errorcall (iinfo->call, "unable to get space rank");

      if (hdf5_global_verbosity > 1)
        Rprintf ("Dataset has rank %d\n",rank);

      if (rank == 0)
        {
          SEXP obj;

          if (hdf5_global_verbosity > 2) 
            Rprintf ("Loading scalar\n");
          PROTECT (obj = allocVector (type, 1));
          load_scalar (iinfo->call, dataset, space, obj);
          iinfo->add (iinfo, name, obj);
          hdf5_load_attributes (iinfo->call, dataset, obj, name);
          UNPROTECT (1);
        }
      else
        {
          hsize_t dims[rank];
          hsize_t maxdims[rank];
          
          if (H5Sget_simple_extent_dims (space, dims, maxdims) < 0)
            errorcall (iinfo->call, "unable to get space extent");

          if (hdf5_global_verbosity > 1)
            {
              hsize_t irank;

              Rprintf ("Dataset has dims/maxdims:");
              for (irank = 0; irank < rank; irank++)
                Rprintf (" %d ",dims[irank]);
              Rprintf ("/");
              for (irank = 0;irank < rank; irank++)
                Rprintf (" %d ",maxdims[irank]);
              Rprintf ("\n");
            }
          if (type == VECSXP && rank == 1)
            {
              unsigned colcount = H5Tget_nmembers (tid), ci;
              SEXP vec;
              size_t size = H5Tget_size (tid);
              unsigned ri, rowcount = dims[0];
              char *buf = (char *) R_chk_calloc (rowcount, size * 2);
              hid_t rtid = make_sexp_ref_type (iinfo->call);
              SEXP names;

              if (hdf5_global_verbosity > 2) 
                {
                  Rprintf ("Dataset has type = VECSXP and rank 1\n");
                  Rprintf ("Reading...\n");
                }

              if (H5Dread (dataset, tid, space, space, H5P_DEFAULT,
                           buf) < 0)
                errorcall (iinfo->call, "can't read compound data vector");
              
              if (hdf5_global_verbosity > 2)
                Rprintf ("....done\n");

              PROTECT (vec = allocVector (VECSXP, colcount));
              PROTECT (names = allocVector (STRSXP, colcount));
              
              for (ci = 0; ci < colcount; ci++)
                {
                  hid_t ctid = H5Tget_member_type (tid, ci);
                  H5T_class_t class = H5Tget_class (ctid);
                  size_t csize = H5Tget_size (ctid);
                  size_t coffset = H5Tget_member_offset (tid, ci);
                  SEXPREC **rowptr = &VECTOR_ELT (vec, ci);
                  unsigned char itembuf[size]; /* for overrun */

#define ASSIGN(val) PROTECT (*rowptr = val)
                  
#define VECLOOP(vectype, vecref, dtid) \
  { \
    size_t dsize = H5Tget_size (dtid); \
    for (ri = 0; ri < rowcount; ri++) \
      { \
	memcpy (itembuf, buf + ri * size + coffset, csize); \
	if (H5Tconvert (ctid, dtid, 1, itembuf, NULL, H5P_DEFAULT) < 0) \
	  errorcall (iinfo->call, "type conversion failed"); \
	memcpy (&vecref (*rowptr)[ri], itembuf, dsize); \
      } \
  }

                  {
                    char *colname = H5Tget_member_name (tid, ci);
                    
                    if (colname)
                      SET_STRING_ELT (names, ci, mkChar (colname));
                    free (colname);
                  }
                  switch (class)
                    {
                    case H5T_INTEGER:
                      {
                        if (csize == 1)
                          {
                            ASSIGN (allocVector (LGLSXP, rowcount));
                            VECLOOP (LGLSXP, LOGICAL, H5T_NATIVE_INT);
                          }
                        else 
                          {
                            ASSIGN (allocVector (INTSXP, rowcount));
                            VECLOOP (INTSXP, INTEGER, H5T_NATIVE_INT);
                          }
                      }
                      break;
                    case H5T_FLOAT:
                      ASSIGN (allocVector (REALSXP, rowcount));
                      VECLOOP (REALSXP, REAL, H5T_NATIVE_DOUBLE);
                      break;
                    case H5T_STRING:
                      ASSIGN (allocVector (STRSXP, rowcount));
                      VECLOOP (STRSXP, STRING_PTR, rtid);
                      break;
                    default:
                      errorcall (iinfo->call, "can't handle hdf class %d",
                                 class);
                    }
                  if (H5Tclose (ctid) < 0)
                    errorcall (iinfo->call, "could not close member type");
                }
              UNPROTECT (colcount);

              if (rowcount < MAX_ROWNAMES_COUNT || MAX_ROWNAMES_COUNT == 0)
                load_rownames_dataset_attribute (iinfo->call, dataset, vec);
              else 
                {
                  SEXP rownames;

                  PROTECT (rownames = allocVector (STRSXP, 0));
                  setAttrib (vec, R_RowNamesSymbol, rownames);
                  UNPROTECT (1);
                }
              setAttrib (vec, R_NamesSymbol, names);
              UNPROTECT (1);
              setAttrib (vec, R_ClassSymbol, mkString ("data.frame"));
              iinfo->add (iinfo, name, vec);
              hdf5_load_attributes (iinfo->call, dataset, vec, name);
              UNPROTECT (1);

              if (H5Tclose (rtid) < 0)
                errorcall (iinfo->call, "could not close reference type");

	      R_chk_free (buf);
            }

          else
            {
              SEXP obj;
              
              if (rank == 1)
                {
                  if (hdf5_global_verbosity > 2) 
                    Rprintf ("Allocating vector with rank=%d dim=%d\n",
                             rank,dims[0]);
                  PROTECT (obj = allocVector (type, dims[0]));
                }
              else if (rank == 2)
                {
                  if (hdf5_global_verbosity > 2) 
                    Rprintf ("Allocating matrix with rank=%d dim=%d %d\n",
                             rank,dims[0],dims[1]);
                  PROTECT (obj = allocMatrix (type, dims[0], dims[1]));
                }
              else
                {
                  SEXP dimsvector;
                  hsize_t irank;

                  if (hdf5_global_verbosity > 2) 
                    Rprintf ("Allocating array with rank=%d dims:", rank);

                  PROTECT (dimsvector = allocVector (INTSXP, rank));
                  for (irank = 0; irank < rank; irank++)
                    {
                      if (hdf5_global_verbosity > 1)
                        Rprintf(" %d ", dims[irank]);
                      INTEGER(dimsvector)[irank] = dims[irank];
                    }
                  if (hdf5_global_verbosity > 2)
                    Rprintf ("\n");
                  if (hdf5_global_verbosity > 2)
                    Rprintf ("about to do actual allocation\n");
                  PROTECT (obj = allocArray (type, dimsvector));
                  if (hdf5_global_verbosity > 2)
                    Rprintf ("done allocation -- unprotecting\n");
                  UNPROTECT (2);
                  if (hdf5_global_verbosity > 2)
                    Rprintf ("done unprotecting -- re-protecting obj\n"); 
                  PROTECT (obj);
		  if (hdf5_global_verbosity > 2)
		    Rprintf ("done re-protecting\n");  

                }
	      if (hdf5_global_verbosity > 2)
		Rprintf ("calling vector_io. Hangs here with big datsets\n");  
              vector_io (iinfo->call, FALSE, dataset, space, obj);
	      if (hdf5_global_verbosity > 2)
		Rprintf ("Phew. Done it. calling iinfo->add\n");  
	      iinfo->add (iinfo, name, obj);
              if (hdf5_global_verbosity > 2)
                {
                  Rprintf ("Rank > 1 or not VECSXP\n");
                  Rprintf ("Calling  hdf5_load_attributes \n");
                }
              hdf5_load_attributes (iinfo->call, dataset, obj, name);
              if (hdf5_global_verbosity > 2) 
                Rprintf ("back from  hdf5_load_attributes \n");
              UNPROTECT (1);
            }
        }
      if (H5Sclose (space) < 0)
	errorcall (iinfo->call, "unable to close dataspace");
      if (H5Tclose (tid) < 0)
	errorcall (iinfo->call, "unable to close datatype");
      if (H5Dclose (dataset) < 0)
	errorcall (iinfo->call, "unable to close dataset");
      if (hdf5_global_verbosity > 0)
        Rprintf ("...Finished dataset \n");
    }
  else
    errorcall (iinfo->call, "no support for HDF object type: %d",
	       statbuf.type);
  return 0;
}

static void
add_to_symbol_table (struct hdf5_iterate_info *iinfo,
		     const char *name,
		     SEXP obj)
{
  char newname[strlen(name)+1];
  strcpy(newname,name);
  if(hdf5_global_nametidy)
    {
      if (hdf5_global_verbosity > 1)
        Rprintf("Tidying name %s ",newname);
      nametidy (newname);
    }
  setVar (install ((char *)newname), obj, iinfo->env);
  if (hdf5_global_verbosity > 1 && hdf5_global_nametidy)
    Rprintf (".. to %s \n ",newname);
}

SEXP do_hdf5load (SEXP args)
{
  const char *path;
  hid_t fid;
  int restore_syms;
  struct hdf5_iterate_info iinfo;
  SEXP call, env;

  args = CDR (args); call = CAR (args);
  args = CDR (args); env = CAR (args);
  args = CDR (args);

  if (!isValidString(CAR(args)))
    errorcall (call, "first argument must be a pathname\n");

  if (TYPEOF (CADR (args)) != LGLSXP)
    errorcall (call, "second argument must be a logical vector");

  path = CHAR (STRING_ELT (CAR(args), 0));
  restore_syms = INTEGER (CADR(args))[0];

  hdf5_global_verbosity = INTEGER (CADDR (args))[0];
  if (hdf5_global_verbosity > 2) 
    Rprintf ("hdf5_global_verbosity=%d load=%d\n",
             hdf5_global_verbosity, restore_syms);
  hdf5_global_nametidy  = INTEGER (CADDDR (args))[0];
  H5dont_atexit ();
  
  if ((fid = H5Fopen (path, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0)
    errorcall (call, "unable to open HDF file: %s", path);

  setup_onexit (fid, env);

  if (H5Tregister (H5T_PERS_SOFT,
		   STRING2REF_CONV,
		   H5T_C_S1,
		   H5T_STD_REF_OBJ, string_ref) < 0)
    errorcall (call, "Unable to register string->ref converter");

  iinfo.call = call;
  iinfo.add = restore_syms ? add_to_symbol_table : add_to_list;
  iinfo.env = env;
  PROTECT (iinfo.ret = CONS (R_NilValue, R_NilValue));

  if (H5Giterate (fid, "/", NULL, hdf5_process_object, &iinfo) < 0)
    errorcall (call, "unable to iterate over HDF file: %s", path);

  UNPROTECT (1);
  return CDR (iinfo.ret);
}
