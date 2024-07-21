#!/bin/sh
message=$1
branch=$2
git add .; git commit -m ${message}; git push gkm ${branch}

