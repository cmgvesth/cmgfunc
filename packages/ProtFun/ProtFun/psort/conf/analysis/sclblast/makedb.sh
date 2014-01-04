#!/bin/sh

echo "Making Gram positive blast database"
cd grampos
formatdb -i sclblast -p T -o T

echo "Making Gram negative blast database"
cd ../gramneg
formatdb -i sclblast -p T -o T

echo "Making Gram archaea blast database"
cd ../archaea
formatdb -i sclblast -p T -o T

echo "Done."
cd ..
