#!/bin/bash
cd $1
find . -type f | sed 's/.*\.//' | sort | uniq -c