#!/usr/bin/env python3
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(output_file, "w") as outf:
    outf.write("ids\tterms\n")
    with open(input_file, "r") as in_f:
        for line in in_f:
            outf.write(line)
