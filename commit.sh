#!/bin/sh
message=$1

git add .; git commit -m ${message}; git push gkm main
