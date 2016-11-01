#!/bin/bash

BASE_DIR="`dirname \"$0\"`/.."
cd $BASE_DIR

ROOT=`git rev-parse --show-toplevel`
cd ${ROOT}
VENV=${ROOT}/venv

echo "checking for virtualenv"
if [ ! -d ${VENV} ]; then
    cd ${ROOT}
    virtualenv --no-site-packages venv --python=python2.7
fi

echo "activate virtualenv"
cd ${VENV}
. ./bin/activate
if [ $? != 0 ]; then
    echo "failed to activate virtualenv at ${VENV}: ABORTING"
    exit 1
fi
cd -

echo "installing requirements"
cd ${ROOT}
pip install -U pip
pip install -r requirements.txt
if [ $? != 0 ]; then
    echo "pip install failed: ABORTING"
    exit 1
fi
cd -

echo "Deploy local complete"
