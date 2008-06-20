
/* $Id$ */

/*
 * Error handling routines.
 *
 * The functions in this file are independent of any application
 * variables, and may be used with any C program.
 * Either of the names CLIENT or SERVER may be defined when compiling
 * this function.  If neither are defined, we assume CLIENT.
 *
 *
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
 * 
 */

#include	<stdio.h>
#include	<varargs.h>

#include	"systype.h"

#ifdef	CLIENT
#ifdef	SERVER
cant define both CLIENT and SERVER
#endif
#endif

#ifndef	CLIENT
#ifndef	SERVER
#define	CLIENT	1		/* default to client */
#endif
#endif

#ifndef	NULL
#define	NULL	((void *) 0)
#endif

char	*pname = NULL;

#ifdef	CLIENT			/* these all output to stderr */

/*
 * Fatal error.  Print a message and terminate.
 * Don't dump core and don't print the system's errno value.
 *
 *	err_quit(str, arg1, arg2, ...)
 *
 * The string "str" must specify the conversion specification for any args.
 */

/*VARARGS1*/
err_quit(va_alist)
va_dcl
{
	va_list		args;
	char		*fmt;

	va_start(args);
	if (pname != NULL)
		fprintf(stderr, "%s: ", pname);
	fmt = va_arg(args, char *);
	vfprintf(stderr, fmt, args);
	fputc('\n', stderr);
	va_end(args);

	exit(1);
}

/*
 * Fatal error related to a system call.  Print a message and terminate.
 * Don't dump core, but do print the system's errno value and its
 * associated message.
 *
 *	err_sys(str, arg1, arg2, ...)
 *
 * The string "str" must specify the conversion specification for any args.
 */

/*VARARGS1*/
err_sys(va_alist)
va_dcl
{
	va_list		args;
	char		*fmt;

	va_start(args);
	if (pname != NULL)
		fprintf(stderr, "%s: ", pname);
	fmt = va_arg(args, char *);
	vfprintf(stderr, fmt, args);
	va_end(args);

	my_perror();

	exit(1);
}

/*
 * Recoverable error.  Print a message, and return to caller.
 *
 *	err_ret(str, arg1, arg2, ...)
 *
 * The string "str" must specify the conversion specification for any args.
 */

/*VARARGS1*/
err_ret(va_alist)
va_dcl
{
	va_list		args;
	char		*fmt;

	va_start(args);
	if (pname != NULL)
		fprintf(stderr, "%s: ", pname);
	fmt = va_arg(args, char *);
	vfprintf(stderr, fmt, args);
	va_end(args);

	my_perror();

	fflush(stdout);
	fflush(stderr);

	return;
}

/*
 * Fatal error.  Print a message, dump core (for debugging) and terminate.
 *
 *	err_dump(str, arg1, arg2, ...)
 *
 * The string "str" must specify the conversion specification for any args.
 */

/*VARARGS1*/
err_dump(va_alist)
va_dcl
{
	va_list		args;
	char		*fmt;

	va_start(args);
	if (pname != NULL)
		fprintf(stderr, "%s: ", pname);
	fmt = va_arg(args, char *);
	vfprintf(stderr, fmt, args);
	va_end(args);

	my_perror();

	fflush(stdout);		/* abort doesn't flush stdio buffers */
	fflush(stderr);

	abort();		/* dump core and terminate */
	exit(1);		/* shouldn't get here */
}

/*
 * Print the UNIX errno value.
 */

my_perror()
{
	char	*sys_err_str();

	fprintf(stderr, " %s\n", sys_err_str());
}

#endif	/* CLIENT */

#ifdef	SERVER

#ifdef	BSD
/*
 * Under BSD, these server routines use the syslog(3) facility.
 * They don't append a newline, for example.
 */

#include	<syslog.h>

#else	/* not BSD */
/*
 * There really ought to be a better way to handle server logging
 * under System V.
 */

#define	syslog(a,b)	fprintf(stderr, "%s\n", (b))
#define	openlog(a,b,c)	fprintf(stderr, "%s\n", (a))

#endif	/* BSD */

char	emesgstr[255] = {0};	/* used by all server routines */

/*
 * Identify ourself, for syslog() messages.
 *
 * LOG_PID is an option that says prepend each message with our pid.
 * LOG_CONS is an option that says write to console if unable to send
 * the message to syslogd.
 * LOG_DAEMON is our facility.
 */

err_init(ident)
char	*ident;
{
	openlog(ident, (LOG_PID | LOG_CONS), LOG_DAEMON);
}

/*
 * Fatal error.  Print a message and terminate.
 * Don't print the system's errno value.
 *
 *	err_quit(str, arg1, arg2, ...)
 *
 * The string "str" must specify the conversion specification for any args.
 */

/*VARARGS1*/
err_quit(va_alist)
va_dcl
{
	va_list		args;
	char		*fmt;

	va_start(args);
	fmt = va_arg(args, char *);
	vsprintf(emesgstr, fmt, args);
	va_end(args);

	syslog(LOG_ERR, emesgstr);

	exit(1);
}

/*
 * Fatal error related to a system call.  Print a message and terminate.
 * Don't dump core, but do print the system's errno value and its
 * associated message.
 *
 *	err_sys(str, arg1, arg2, ...)
 *
 * The string "str" must specify the conversion specification for any args.
 */

/*VARARGS1*/
err_sys(va_alist)
va_dcl
{
	va_list		args;
	char		*fmt;

	va_start(args);
	fmt = va_arg(args, char *);
	vsprintf(emesgstr, fmt, args);
	va_end(args);

	my_perror();
	syslog(LOG_ERR, emesgstr);

	exit(1);
}

/*
 * Recoverable error.  Print a message, and return to caller.
 *
 *	err_ret(str, arg1, arg2, ...)
 *
 * The string "str" must specify the conversion specification for any args.
 */

/*VARARGS1*/
err_ret(va_alist)
va_dcl
{
	va_list		args;
	char		*fmt;

	va_start(args);
	fmt = va_arg(args, char *);
	vsprintf(emesgstr, fmt, args);
	va_end(args);

	my_perror();
	syslog(LOG_ERR, emesgstr);

	return;
}

/*
 * Fatal error.  Print a message, dump core (for debugging) and terminate.
 *
 *	err_dump(str, arg1, arg2, ...)
 *
 * The string "str" must specify the conversion specification for any args.
 */

/*VARARGS1*/
err_dump(va_alist)
va_dcl
{
	va_list		args;
	char		*fmt;

	va_start(args);
	fmt = va_arg(args, char *);
	vsprintf(emesgstr, fmt, args);
	va_end(args);

	my_perror();
	syslog(LOG_ERR, emesgstr);

	abort();		/* dump core and terminate */
	exit(1);		/* shouldn't get here */
}

/*
 * Print the UNIX errno value.
 * We just append it to the end of the emesgstr[] array.
 */

my_perror()
{
	register int	len;
	char		*sys_err_str();

	len = strlen(emesgstr);
	sprintf(emesgstr + len, " %s", sys_err_str());
}

#endif	/* SERVER */

			/* remainder is for both CLIENT and SERVER */
extern int	errno;		/* Unix error number */
extern int	sys_nerr;	/* # of error message strings in sys table */
extern char	*sys_errlist[];	/* the system error message table */

#ifdef	SYS5
int	t_errno;	/* in case caller is using TLI, these are "tentative
			   definitions"; else they're "definitions" */
int	t_nerr;
char	*t_errlist[1];
#endif


/*
 * Return a string containing some additional operating-system
 * dependent information.
 * Note that different versions of UNIX assign different meanings
 * to the same value of "errno" (compare errno's starting with 35
 * between System V and BSD, for example).  This means that if an error
 * condition is being sent to another UNIX system, we must interpret
 * the errno value on the system that generated the error, and not
 * just send the decimal value of errno to the other system.
 */

char *
sys_err_str()
{
	static char	msgstr[200];

	if (errno != 0) {
		if (errno > 0 && errno < sys_nerr)
			sprintf(msgstr, "(%s)", sys_errlist[errno]);
		else
			sprintf(msgstr, "(errno = %d)", errno);
	} else {
		msgstr[0] = '\0';
	}

#ifdef	SYS5
	if (t_errno != 0) {
		char	tmsgstr[100];

		if (t_errno > 0 && t_errno < sys_nerr)
			sprintf(tmsgstr, " (%s)", t_errlist[t_errno]);
		else
			sprintf(tmsgstr, ", (t_errno = %d)", t_errno);

		strcat(msgstr, tmsgstr);	/* catenate strings */
	}
#endif

	return(msgstr);
}
