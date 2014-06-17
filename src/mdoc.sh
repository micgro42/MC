#!/bin/bash



if [ "`hostname`" = "Leawyn" ]; then
doxygen KonGrad.doxyfile
chmod a+r ../doc/html/*
scp -rpC ../doc/html/* micgro42@eiger.physik.hu-berlin.de:CP32014/KonGrad/trunk/doc/html
else
./makedoc.sh
chmod a+r ../doc/html/*
cp -rp ../doc/html/* $HOME/public_html/cp3/mc/
fi
cat doxywarnings.log
