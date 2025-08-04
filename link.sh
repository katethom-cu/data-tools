#!/usr/bin/env bash
install_dir=$HOME/.local/lib/util
mkdir -p $install_dir

for file in $(ls src/)
do
    ln -s $(pwd)/src/$file $install_dir/$file
done
