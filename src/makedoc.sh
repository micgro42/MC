#!/bin/bash
if [ ! -f Doxyfile  ]
then
  doxygen -g
  sed --in-place=".bak" 's/^HAVE_DOT.*/HAVE_DOT = YES/' Doxyfile
  sed -i 's/^EXTRACT_PRIVATE.*/EXTRACT_PRIVATE = YES/' Doxyfile
  sed -i 's/^EXTRACT_STATIC.*/EXTRACT_STATIC = YES/' Doxyfile
  sed -i 's/^SOURCE_BROWSER.*/SOURCE_BROWSER = YES/' Doxyfile
  sed -i 's/^INLINE_SOURCES.*/INLINE_SOURCES = YES/' Doxyfile
  sed -i 's/^REFERENCED_BY_RELATION.*/REFERENCED_BY_RELATION = YES/' Doxyfile
  sed -i 's/^REFERENCES_RELATION.*/REFERENCES_RELATION = YES/' Doxyfile
  sed -i 's/^HTML_TIMESTAMP.*/HTML_TIMESTAMP = YES/' Doxyfile
  sed -i 's/^CALL_GRAPH.*/CALL_GRAPH = YES/' Doxyfile
  sed -i 's/^CALLER_GRAPHS.*/CALLER_GRAPH = YES/' Doxyfile
  sed -i 's/^DOT_MULTI_TARGETS.*/DOT_MULTI_TARGETS = YES/' Doxyfile
  sed -i 's/^WARN_LOGFILE.*/WARN_LOGFILE = doxywarnings.log/' Doxyfile
  sed -i 's/^EXCLUDE.*/EXCLUDE = Logging.cc Logging.hh/' Doxyfile
  sed -i 's/^OUTPUT_DIRECTORY.*/OUTPUT_DIRECTORY = ..\/doc/' Doxyfile
  sed -i 's/^PROJECT_NAME.*/PROJECT_NAME = "MC"/' Doxyfile
fi

doxygen Doxyfile

echo "This logfile was produced at `date`" >> doxywarnings.log

cat doxywarnings.log
