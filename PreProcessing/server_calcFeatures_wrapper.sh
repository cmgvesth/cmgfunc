#! /bin/sh

LD_LIBRARY_PATH=/home/cmgfunc/ProtFun/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH
exec python /home/cmgfunc/CMGfunc/feature_selection.py
