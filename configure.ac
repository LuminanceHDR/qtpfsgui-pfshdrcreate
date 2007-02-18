dnl
dnl Initialization part
dnl
AC_INIT([pfshdrcreate], [0.1], [grota@users.sourceforge.net])
AM_INIT_AUTOMAKE
AM_CONFIG_HEADER(config.h)
AC_LANG(C++)			# perform tests in C++ language
AC_PROG_CXX
AC_PROG_CC        
AC_PROG_LIBTOOL


dnl
dnl Default setup
dnl
CXXFLAGS="-O3 -funroll-loops -fstrength-reduce -fschedule-insns2 -felide-constructors -frerun-loop-opt -fexceptions -fno-strict-aliasing -fexpensive-optimizations -ffast-math -pipe"

dnl
dnl Templates for autoheader defines
dnl
AH_TEMPLATE([DEBUG],
	    [Define as 1 when compiling in debug mode ])


dnl
dnl Commandline switches
dnl
AC_ARG_ENABLE(debug,
	AC_HELP_STRING([--enable-debug],
                       [Enable debug symbols and asserts]),
		 [if test "$enable_debug" = "yes" ; then
		     CXXFLAGS="-g"
		 fi],[enable_debug="no"])


dnl
dnl System detection and configuration
dnl

dnl PFS tools library - base for input sources
AC_MSG_CHECKING([for pfs library])
PKG_CHECK_MODULES(PFS, pfs >= 1.0)
AC_SUBST(PFS_CFLAGS)
AC_SUBST(PFS_LIBS)

dnl Debug version
if test "$enable_debug" = "yes" ; then
  AC_DEFINE([DEBUG], 1, [Compilation for debugging purposes.])
fi


dnl
dnl Autoheader defines
dnl
AH_BOTTOM([
/* Output stream for debug messages. */
#ifdef DEBUG
#define DEBUG_STR std::cerr
#else
#define DEBUG_STR if(1); else std::cerr
#endif 

/* Output stream for verbose mode messages. */
#define VERBOSE_STR if( verbose ) std::cerr << PROG_NAME << ": "
extern bool verbose; /* verbose should be declared for each standalone code */
])


dnl
dnl Finish
dnl
AC_CONFIG_FILES([ Makefile \
                  src/Makefile \
                ])

AC_OUTPUT

AC_MSG_RESULT([
pfshdrcreate is now configured. Type 
make 
to compile the sources and 
make install 
to install the package on your system.

  Source directory:     $srcdir
  Installation prefix:  $prefix
  C++ compiler:         $CXX $XTRA_CXXFLAGS $WARN_CXXFLAGS $CXXFLAGS
  LIBS:                 $LIBS

  Debug mode:           $enable_debug
])
