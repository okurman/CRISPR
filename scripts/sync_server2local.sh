#!/bin/sh

usage ()
{
 echo 'Usage: sync_local2server.sh < code | data >' 
 exit
}

if [ "$#" -ne 1 ]
then 
 usage
fi

folder=$1
echo $folder

rsync -va --delete hudaiber@cbbdev11:~/Projects/NewSystems/$folder/ ~/Projects/NewSystems/$folder/  

