#!/bin/bash
rm -rf ./doc
doxygen Doxyfile 
cd ./doc/html/
git init
git checkout -b gh-pages
git add .
git add origin git@github.com:maxhutch/nek.git
git remote add origin git@github.com:maxhutch/nek.git
git commit -m "Testing Doxygen"
git push origin gh-pages --force
