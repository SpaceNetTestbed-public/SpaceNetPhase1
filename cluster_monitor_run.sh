#!/bin/bash

# Set options
OUTPUT="output/resource/default_output2.txt"
INTERVAL=1
LOGCOUNT=30

python3 resource_monitor/cpu_ram_log.py -u $USER -o $OUTPUT --interval $INTERVAL --logcount $LOGCOUNT