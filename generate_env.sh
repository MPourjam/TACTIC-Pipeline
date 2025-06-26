#!/bin/bash
git_tag=$(git describe --tags --abbrev=0)
grep -v '^GIT_TAG=' .env > .env.tmp
echo "GIT_TAG=$git_tag" >> .env.tmp
mv .env.tmp .env
