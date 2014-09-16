#!/bin/bash

NAME="$(date +%F-%T)"

python plot.py project "plots/$NAME" > "logs/${NAME}.log"
