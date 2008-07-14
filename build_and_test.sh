#! /bin/bash

# Simple Bash script to automate the build and test process.

# build extensions and pure Python modules
python setup.py build

# make sure the tests use the current code
cp -rf build/lib/TOPKAPI/ tests/test_continuity/

# run tests
nosetests -v -w tests/test_continuity/