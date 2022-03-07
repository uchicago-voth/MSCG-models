#!/bin/bash

grep "Proxy mapping:" controlLog.out > map_grepped.dat

tr ' ' '\n' < map_grepped.dat > map_split.dat

sed '/^\s*$/d' < map_split.dat > map_trimmed.dat

tail -n +3 map_trimmed.dat > mapped.trans

