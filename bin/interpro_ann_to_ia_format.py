#!/usr/bin/env python3
"""
Q8Z0C0	IPR029039;IPR001279;IPR008254;IPR036866;IPR045761;IPR051285;IPR002563;IPR012349
Q6FKR3	IPR036395;IPR001083
P68470	IPR036058;IPR051597;IPR002350;IPR001239
Q20060	IPR027417;IPR036277;IPR003395;IPR024704;IPR010935
Q03LU2	IPR000836;IPR050137;IPR023050;IPR029057
A7HL08	IPR028989;IPR003728;IPR035956;IPR028998;IPR036847
P05693	IPR014710;IPR022379;IPR006044;IPR050253;IPR006045;IPR011051
A9KCC7	IPR018220;IPR027417;IPR042109;IPR042111;IPR033128;IPR042110;IPR001114
B3PYZ4	IPR036267;IPR013849;IPR011114;IPR012340;IPR010994;IPR060101;IPR000085
P46102	IPR000169;IPR025661;IPR025660;IPR013128;IPR000668;IPR039417;IPR038765;IPR013201
Q9C1X5	IPR020471;IPR036812;IPR023210;IPR018170
P56228	IPR010000

To:
EntryID term    aspect
target1 GO:0033058      BPO
target1 GO:0043056      BPO
target1 GO:0050879      BPO
target1 GO:0071965      BPO
target2 GO:0031987      BPO
target2 GO:0032501      BPO
"""

import sys


if __name__ == "__main__":
    ann_path = sys.argv[1]
    ia_path = sys.argv[2]

    output_stream = open(ia_path, "w")
    output_stream.write("EntryID\tterm\taspect\n")

    for ann_line in open(ann_path, "r"):
        parts = ann_line.strip().split("\t")
        if len(parts) < 2:
            continue
        entry_id = parts[0]
        terms = parts[1].split(";")
        for term in terms:
            output_stream.write(f"{entry_id}\t{term}\tinterpro\n")

    output_stream.close()
