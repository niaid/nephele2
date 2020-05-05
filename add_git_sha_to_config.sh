#!/bin/bash
sha=$(git rev-parse --short HEAD)
long_sha=$(git log --format="%H" -n 1)

# rm SHA and LONG_SHA lines:
sed -i -r "/^SHA_SHORT[ ]?=/d" config.py
sed -i -r "/^SHA_LONG[ ]?=/d" config.py

printf "SHA_LONG = \"$long_sha\"\nSHA_SHORT = \"$sha\"\n" >> config.py 


