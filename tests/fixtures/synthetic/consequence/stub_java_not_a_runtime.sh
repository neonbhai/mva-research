#!/bin/sh
# mva_synthetic=true
#
# A stand-in for macOS's /usr/bin/java shim, used ONLY by
# tests/unit/test_consequence_adapter.py.
#
# macOS ships a `java` that is present, executable, and NOT a Java runtime: every
# check short of running it passes, and running it prints the line below and exits
# non-zero. This reproduces that exactly, so the adapter's JVM probe can be tested
# on a machine where a real JDK *is* installed.
#
# The wording is copied from the real shim; the adapter must diagnose this as "not
# a working Java runtime" rather than let it surface as a SnpEff failure.
echo "The operation couldn't be completed. Unable to locate a Java Runtime." >&2
echo "Please visit http://www.java.com for information on installing Java." >&2
exit 1
