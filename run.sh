#!/bin/bash

function HRS() {
    NAME="$(date +%F-%T)"

    IBDDIR="project/baharian_projects/MergedData/phased/3_GERMLINE/cMcorrected" 
    BEDDIR="project/barakatt_projects/HRS/results/HRS_AFRAM_20140609/outbed"
    BEDDIRnotphased="project/barakatt_projects/HRS/results/HRS_AFRAM_notphased_20140921/outbed"
    IDLIST="project/baharian_projects/HRS/data/AfricanAmericans/AfrAm.SubjID.list"

    if python plot.py --ibd $IBDDIR/MERGED_chr1.cM.IBD.match.gz --bed $BEDDIRnotphased --ids $IDLIST \
        -o "plots/${NAME}.png" > "logs/${NAME}.log"
    then
        notify-send "IBD and local ancestry plots" "Finished generating 'plots/$NAME'."
    else
        notify-send "IBD and local ancestry plots" "Error generating 'plots/$NAME'."
    fi
}

HRS
