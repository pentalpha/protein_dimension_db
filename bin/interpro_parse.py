#!/usr/bin/env python
import sys

from data.interpro_api.parsing import parse_interpro_raw

if __name__ == "__main__":
    output_file = sys.argv[1]
    input_files = sys.argv[2:]

    input_files = ",".join(input_files)
    parse_interpro_raw(input_files, output_file)
