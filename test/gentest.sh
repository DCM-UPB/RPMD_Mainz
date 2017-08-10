#!/bin/bash
# First argument is expected to be the RPMD executable's full path 

TEST_PATH=$(pwd)

if [ "$1" == "" ]
  then
    RPMD_PATH=${TEST_PATH}/../qmd.x
  else
    RPMD_PATH=$1
fi

RUN_PATH=${TEST_PATH}/run/gentest.run
CMP_PATH=${TEST_PATH}/cmp/gentest.cmp

INPUT_NAME=input
PARAM_NAME=param_scott

rm -rv $RUN_PATH
mkdir -p $RUN_PATH

cp ${CMP_PATH}/$INPUT_NAME $RUN_PATH
cp ${CMP_PATH}/$PARAM_NAME $RUN_PATH

cd $RUN_PATH
$RPMD_PATH $INPUT_NAME $PARAM_NAME > output

cd ../
diff $CMP_PATH $RUN_PATH > gentest.diff
if [ "$(wc -c gentest.diff)" == "0 gentest.diff" ]
  then
    echo "Test finished succesfully!"
  else
    echo "Test produces unexpected output:"
    cat gentest.diff
fi

cd $TEST_PATH
