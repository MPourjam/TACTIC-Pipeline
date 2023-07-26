#!/usr/bin/env python3

"""script to help mounting network shares
"""

import os
import ntpath
import sys


def log_to_status_file(msg):
    with open('/usr/local/bin/status.txt', 'w') as stat_f_h:
        stat_f_h.write(msg)

log_to_status_file('Mounting drive containing files')
path = sys.argv[1]

drive, _ = ntpath.splitdrive(path)
drive_letter = drive.replace(':', '').lower()
mount_path = f'/mnt/{drive_letter}'

os.makedirs(mount_path, exist_ok=True)

cmd = f'mount -t drvfs "{drive}" "{mount_path}"'
os.system(cmd)
