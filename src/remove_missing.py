#!/usr/bin/env python
import sys
path = sys.argv[1]

with open(path) as file:
    n_fields = len(file.readline().split("\t")) + 1
    empty_columns = set(range(n_fields))
    for line in file:
        fields = line.split("\t")
        ok_fields = set([i for i in empty_columns if fields[i] !="NA"])
        empty_columns = empty_columns.difference(ok_fields)
    keep_columns = set(range(n_fields)).difference(empty_columns)
    keep_columns.add(0)
    sys.stderr.write("%s empty columns (genes) will be removed post-normalization.\n" % len(empty_columns))

with open(path) as file:
    print "\t".join(field for i, field in enumerate(file.readline().split("\t")) \
                        if i in [j+1 for j in keep_columns])
    for line in file:
        print "\t".join(field for i,field in enumerate(line.split("\t")) if i in keep_columns)

