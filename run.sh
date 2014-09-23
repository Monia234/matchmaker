#!/bin/bash

NAME="$(date +%F-%T)"

python plot.py project "plots/$NAME" > "logs/${NAME}.log"
notify-send "IBD and local ancestry plots" "Finished generating 'plots/$NAME'."
