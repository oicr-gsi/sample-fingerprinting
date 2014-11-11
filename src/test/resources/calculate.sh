#!/bin/bash
cd $1
ls | grep -v SWID | grep \. | grep -v total | grep -v finfiles | grep -v images | sed 's/.*\.//' | sort | uniq -c