#!/usr/bin/env python3
import sys

def log_to_status_file(msg):
    with open('/usr/local/bin/status.txt', 'w') as stat_f_h:
        stat_f_h.write(msg)

log_to_status_file(sys.argv[1])
