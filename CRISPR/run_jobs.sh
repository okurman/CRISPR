#!/usr/bin/env bash

for order in 2 3 4 5
do
    for cas in 1 4
    do
        python post_processing.py $order $cas &
    done;
    sleep 1
done;