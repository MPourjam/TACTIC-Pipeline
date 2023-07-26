#!/usr/bin/env python3

"""script to help unmounting network shares
"""

import os
import ntpath
import sys


def log_to_status_file(msg):
    with open('/usr/local/bin/status.txt', 'w') as stat_f_h:
        stat_f_h.write(msg)
log_to_status_file('Clean-up. Unmounting drive containing files')

path = sys.argv[1]

drive, _ = ntpath.splitdrive(path)
drive_letter = drive.replace(':', '').lower()
mount_path = f'/mnt/{drive_letter}'

cmd = f'umount "{mount_path}"'
os.system(cmd)

try:
    os.rmdir(mount_path)
except FileNotFoundError:
    print('Umount failed')
    exit(0)
except OSError:
    print('Umount failed. Busy.')
    exit(0)
