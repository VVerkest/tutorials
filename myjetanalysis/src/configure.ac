AC_INIT(calojetrhoest,[1.00])

AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE

AC_PROG_CXX(CC g++)

LT_INIT([disable-static])


dnl   no point in suppressing warnings people should 
dnl   at least see them, so here we go for g++: -Wall
if test $ac_cv_prog_gxx = yes; then
   CXXFLAGS="$CXXFLAGS -Wall -Werror"
   # CXXFLAGS="$CXXFLAGS -Wall -Werror -ggdb3"
fi

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
