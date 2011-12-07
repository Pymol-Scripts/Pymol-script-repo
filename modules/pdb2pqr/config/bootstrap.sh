#! /bin/sh

# Basic configuration of autoconf et al, adapated from script by Mike Holst

rm -rf config.cache autom4te.cache

aclocal --verbose \
&& automake --verbose --gnu --add-missing --copy \
&& autoconf --verbose \
&& libtoolize --automake --copy --force

rm -rf config.cache autom4te.cache

