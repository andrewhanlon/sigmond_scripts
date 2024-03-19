#!/usr/bin/env bash

infile=$1
margins=$2

stub=${infile%%.agr}
psfile=${stub}.ps
pdffile=${stub}.pdf

gracebat $infile -printfile $psfile
ps2pdf14 $psfile $pdffile
pdfcrop --margins $margins $pdffile $pdffile

rm -f $psfile

exit 0
