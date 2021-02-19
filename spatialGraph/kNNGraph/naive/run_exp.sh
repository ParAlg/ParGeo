#!/bin/bash

# $1: file name
# $2: k
# $3: 1 if want to re compile

# if [ $3 -eq 1 ]
# then
#     make clean
#     GCILK=1 make -j
# fi

./knngraph -r 1 -k $2 -d /Users/sy/Desktop/MIT/OPTICS/small_datasets/$1.pbbs

# USEJEMALLOC=1  