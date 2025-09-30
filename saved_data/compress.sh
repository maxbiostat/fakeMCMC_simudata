#!/usr/bin/bash
for i in *.csv
 do GZIP=-9 tar -czvf $i.tar.gz $i
done
