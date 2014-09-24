#!/bin/bash


AFRAMIDLIST="project/baharian_projects/HRS/data/AfricanAmericans/AfrAm.SubjID.list"

function SCCS() {
    IBDDIR="project/baharian_projects/MergedData/phased/3_GERMLINE/cMcorrected" 
    BEDDIR="project/barakatt_projects/HRS/results/SCCS_notphased_20140923/outbed"
    NAME="SCCS-$(date +%F-%T)"

        
    if python plot.py --ibd $IBDDIR/MERGED_chr1.cM.IBD.match.gz --bed $BEDDIR --ids dummy \
        -o "plots/${NAME}.png" > "logs/${NAME}.log"
    then
        notify-send "SCCS IBD and local ancestry plots" "Finished generating 'plots/$NAME'."
    else
        notify-send "SCCS IBD and local ancestry plots" "Error generating 'plots/$NAME'."
    fi
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

SCCS
