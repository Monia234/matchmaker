#!/bin/bash

NAME="$(date +%F-%T)"

python plot.py project "plots/$NAME" | tee "logs/${NAME}.log"
