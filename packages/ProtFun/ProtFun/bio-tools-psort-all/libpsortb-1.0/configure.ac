#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.60)
AC_INIT([libpsortb], [1.0], [lairdm@sfu.ca])
AC_CONFIG_SRCDIR([svmloc/binding.cpp])
AC_CONFIG_HEADER([config.h])
AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_LIBTOOL

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [main])

# We need to set lib64 for 64-bit versions of Linux
libnn=lib
case "${host_os}" in
  linux*)
    ## Not all distros use this: some choose to march out of step  
    case "${host_cpu}" in
      x86_64|mips64|ppc64|powerpc64|sparc64|s390x)
        libnn=lib64
      ;;
    esac
    ;;
  solaris*)
    ## libnn=lib/sparcv9 ## on 64-bit only, but that's compiler-specific
    ;;  
esac
: ${LIBnn=$libnn}
AC_SUBST(LIBnn)
## take care not to  override the command-line setting
if test "${libdir}" = '${exec_prefix}/lib'; then
  libdir='${exec_prefix}/${LIBnn}'
fi


# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h limits.h memory.h netinet/in.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_HEADER_TIME
AC_STRUCT_TM

# Checks for library functions.
AC_FUNC_FSEEKO
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STAT
AC_FUNC_STRFTIME
AC_FUNC_VPRINTF
AC_FUNC_STRTOD
AC_CHECK_FUNCS([atexit floor memmove memset pow sqrt strchr strcspn strerror strpbrk strrchr strspn strstr strtol])

AC_CONFIG_FILES([Makefile
                 hmmer/Makefile
                 squid/Makefile
                 svmloc/Makefile
		 modhmm0.92b/Makefile
		 include/Makefile])
AC_OUTPUT
