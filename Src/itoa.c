
/* $Id$ */
/*
 * This work was produced at the University of California, Lawrence
 * Livermore National Laboratory (UC LLNL) under contract no.
 * W-7405-ENG-48 (Contract 48) between the U.S. Department of Energy
 * (DOE) and The Regents of the University of California (University)
 * for the operation of UC LLNL. Copyright is reserved to the University
 * for purposes of controlled dissemination, commercialization through
 * formal licensing, or other disposition under terms of Contract 48;
 * DOE policies, regulations and orders; and U.S. statutes. The rights
 * of the Federal Government are reserved under Contract 48 subject to
 * the restrictions agreed upon by the DOE and University as allowed
 * under DOE Acquisition Letter 97-1.
 *
 *
 * DISCLAIMER
 *
 * This work was prepared as an account of work sponsored by an agency
 * of the United States Government. Neither the United States Government
 * nor the University of California nor any of their employees, makes
 * any warranty, express or implied, or assumes any liability or
 * responsibility for the accuracy, completeness, or usefulness of any
 * information, apparatus, product, or process disclosed, or represents
 * that its use would not infringe privately-owned rights.  Reference
 * herein to any specific commercial products, process, or service by
 * trade name, trademark, manufacturer or otherwise does not necessarily
 * constitute or imply its endorsement, recommendation, or favoring by
 * the United States Government or the University of California. The
 * views and opinions of authors expressed herein do not necessarily
 * state or reflect those of the United States Government or the
 * University of California, and shall not be used for advertising or
 * product endorsement purposes.
 */

/*
 * itoa()	Convert integer digit(s) to characters and save them to s.
 *		This function is based on Brian and Dennis C's itoa function.
 */

#include	<stdio.h>
#include	<string.h>

void itoa(n, s)
int	n;
char	s[];
{
    int c, i, j, k, sign;

    if ((sign = n) < 0)		/* save the sign of n */
        n = -n;			/* make it to be a positive */
    i = 0;

    do
    {
        s[i++] = n % 10 + '0';	/* get next digit */
    }
    while ( ( n /= 10 ) > 0 );

    if (sign < 0 )
        s[i++] = '-';

    s[i] = '\0';

    for ( j = 0, k = strlen(s)-1; j < k; j++, k-- )
    {
        c = s[j];
        s[j] = s[k];
        s[k] = c;
    }
}
