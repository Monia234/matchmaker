#!/bin/bash

AFRAMIDLIST="project/baharian_projects/HRS/data/AfricanAmericans/AfrAm.SubjID.list"

function SCCS {
    IBDDIR="project/baharian_projects/MergedData/phased/3_GERMLINE/cMcorrected" 
    BEDDIR="project/barakatt_projects/HRS/results/SCCS_notphased_20140923/outbed"

    for chr in $(seq $1 $2) ; do
        NAME="SCCS-$(date +%F-%T-chr$chr)"
        if python plot.py --ibd $IBDDIR/MERGED_chr${chr}.cM.IBD.match.gz --bed $BEDDIR --ids dummy \
            -o "plots/${NAME}.png" > "logs/${NAME}.log"
        then
            notify-send "SCCS IBD and local ancestry plots" "Finished generating 'plots/$NAME'."
        else
            notify-send "SCCS IBD and local ancestry plots" "Error generating 'plots/$NAME'."
        fi
    done
}

function HRS() {
    IBDDIR="project/baharian_projects/MergedData/phased/3_GERMLINE/cMcorrected" 
    BEDDIR="project/barakatt_projects/HRS/results/HRS_AFRAM_20140609/outbed"
    BEDDIRnotphased="project/barakatt_projects/HRS/results/HRS_AFRAM_notphased_20140921/outbed"
    NAME="HRS-$(date +%F-%T)"

    if python plot.py --ibd $IBDDIR/MERGED_chr1.cM.IBD.match.gz --bed $BEDDIRnotphased --ids $AFRAMIDLIST \
        -o "plots/${NAME}.png" > "logs/${NAME}.log"
    then
        notify-send "IBD and local ancestry plots" "Finished generating 'plots/$NAME'."
    else
        notify-send "IBD and local ancestry plots" "Error generating 'plots/$NAME'."
    fi
}

# HRS

if ! test -n "$1" ; then
    startchr=18
    endchr=18
else
    startchr=$1
    endchr=$2
fi

SCCS $startchr $endchr
