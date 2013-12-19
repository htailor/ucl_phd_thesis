#!/bin/sh

find . -name "core.*" -exec rm -rf {} \;
find . -name "*.*~" -exec rm -rf {} \;
find . -name "*.old" -exec rm -rf {} \;
find . -name "*.bak" -exec rm -rf {} \;
