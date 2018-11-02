#!/usr/bin/env python

import subprocess


def run_cmd(cmd, logfile=None):
    if logfile is None:
        subprocess.call(cmd, shell=True)
    else:
        with open(logfile, 'a+') as f:
            f.write(cmd)
            subprocess.call(cmd, shell=True, stdout=f, stderr=f)
