#!/bin/sh
architecture="x86"

for choice in "$@" 
do
    case "${choice}" in
        "x86") architecture="x86";;
        "aarch") architecture="aarch";;
        *) echo "ERROR: Unknown architecture requested.\nDefault architecture (x86) selected.";;
    esac
done

cp lib/libmatplot_${architecture}.so lib/libmatplot.so;
echo "Selected ${architecture}.";
