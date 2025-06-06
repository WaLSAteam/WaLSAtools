/* Copyright (C) 1993-2003 Ghostgum Software Pty Ltd.  All rights reserved.

  This software is provided AS-IS with no warranty, either express or
  implied.

  This software is distributed under licence and may not be copied,
  modified or distributed except as expressly authorised under the terms
  of the licence contained in the file LICENCE in this distribution.

  For more information about licensing, please refer to
  http://www.ghostgum.com.au/ or contact Ghostsgum Software Pty Ltd, 
  218 Gallaghers Rd, Glen Waverley VIC 3150, AUSTRALIA, 
  Fax +61 3 9886 6616.
*/

/* $Id: cps.h,v 1.2 2003/01/12 06:32:44 ghostgum Exp $ */
/* Copying PostScript files */

int ps_extract(Doc *doc, LPCTSTR filename, PAGELIST *pagelist, int copies);
int ps_copy(GFile *outfile, GFile *infile, FILE_OFFSET begin, FILE_OFFSET end);
int ps_fgets(char *s, int n, GFile *f);

