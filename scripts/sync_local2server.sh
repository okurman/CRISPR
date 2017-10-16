#!/bin/bash

source_path=~/Projects/NewSystems/code/
dest_path=/Volumes/pan1/patternquest/Projects/NewSystems/

echo 'Starting: rsync '$source_path' '$dest_path

rsync -av --exclude ".git" --exclude "*__pycache__*" --exclude "*.pyc" $dest_path/code/ $dest_path/code_backup/
rsync -av --exclude ".git" --exclude "*__pycache__*" --exclude "*.pyc" $source_path $dest_path/code/

