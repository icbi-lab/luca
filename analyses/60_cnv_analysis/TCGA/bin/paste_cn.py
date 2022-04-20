#!/usr/bin/env python

# invoke with column nr to extract as first parameter followed by
# file with filenames. The files should all have the same number
# of rows. The 3rd arg is the cutoff of LOSS/GAIN

import sys
import os

col = int(sys.argv[1])
files = sys.argv[2]
cutoff = float(sys.argv[3])
f_list = open(files, "r")


res = {}
check_order = {}

file_nr = 0

for sample in f_list:
    file_name, sample_name = sample.strip().split()
    file_name = os.path.expanduser(file_name)
    for line_nr, line in enumerate(open(file_name)):
        if file_nr == 0:
            cols = line.strip().split('\t')
            #row_name = cols[0] + ":" + cols[1] + "-" + cols[2] "-" + cols[3] "-" + cols[4]
            #res.setdefault(line_nr, []).append(row_name)
            #res.setdefault(line_nr, []).append(cols[3])
            row_def = "\t".join(cols[0:col-1])

            res.setdefault(line_nr, []).append(row_def)
            check_order[line_nr] = cols[0]

        if line_nr == 0:
            res.setdefault(line_nr, []).append(sample_name)
        else:
            fields = line.strip().split('\t')
            if fields[0] != check_order[line_nr]:
                print("ERROR: wrong sort order: " + fields[0] + " <-> " + check_order[line_nr])
                sys.exit(1)
            if len(fields) >= (col):
                val = fields[col-1]
            else:
                val = ""
            res.setdefault(line_nr, []).append(val)
    file_nr += 1

for line_nr in sorted(res):
    if line_nr != 0:
        cn_n = [float(i) for i in res[line_nr][1:] if i != ""]
        n_ok = len(cn_n)
        if n_ok > 0:
            mean_cn = sum(cn_n) / n_ok
        else:
            mean_cn = ""

        cn_amp = [float(i) for i in res[line_nr][1:] if i != "" and float(i) > cutoff]
        n_amp = len(cn_amp)
        if n_amp > 0:
            mean_cn_amp = sum(cn_amp) / n_amp
        else:
            mean_cn_amp = ""

        cn_del = [float(i) for i in res[line_nr][1:] if i != "" and float(i) < -cutoff]
        n_del = len(cn_del)
        if n_del > 0:
            mean_cn_del = sum(cn_del) / n_del
        else:
            mean_cn_del = ""

        if n_ok > 0:
            f_amp = n_amp / n_ok
            f_del = n_del / n_ok
            n_no_change = n_ok - n_del - n_amp
            f_no_change = n_no_change / n_ok
        else:
            f_amp = 0
            f_del = 0
            f_no_change = 0
            n_no_change = 0


        res[line_nr].append(n_ok)
        res[line_nr].append(mean_cn)
        res[line_nr].append(mean_cn_amp)
        res[line_nr].append(mean_cn_del)
        res[line_nr].append(f_amp)
        res[line_nr].append(n_amp)
        res[line_nr].append(f_del)
        res[line_nr].append(n_del)
        res[line_nr].append(f_no_change)
        res[line_nr].append(n_no_change)

    else:
        res[line_nr].append("count_ok")
        res[line_nr].append("mean_cn")
        res[line_nr].append("mean_cn_amp")
        res[line_nr].append("mean_cn_del")
        res[line_nr].append("f_amp")
        res[line_nr].append("n_amp")
        res[line_nr].append("f_del")
        res[line_nr].append("n_del")
        res[line_nr].append("f_no_change")
        res[line_nr].append("n_no_change")

    new_line = [str(i) for i in res[line_nr]]
    # print('\t'.join(res[line_nr]))
    print('\t'.join(new_line))
