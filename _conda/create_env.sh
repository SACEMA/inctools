#!/bin/zsh
conda remove --name r-test-release --all
conda env create -f r-test-release.yml
conda remove --name r-test-oldrel --all
conda env create -f r-test-oldrel.yml
