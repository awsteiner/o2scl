#!/usr/bin/env python3

import sys
import subprocess

def main():
    # List of flags to filter out
    skip_flags = {'-fPIC', '-DPIC'}

    # Filter arguments
    filtered_args = [arg for arg in sys.argv[1:] if arg not in skip_flags]

    # Call nvcc with filtered arguments
    result = subprocess.run(['nvcc'] + filtered_args)
    sys.exit(result.returncode)

if __name__ == '__main__':
    main()
