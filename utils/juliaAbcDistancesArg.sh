
# Runs ABC
# USAGE: bash juliaAbcDistancesArg.sh <dataset> <prefix> <nAccepted>

# Test that required args are supplied:
if [ -z "$1" ]
  then
    echo "Error: Must supply dataset." >&2
    exit 1
fi



if [ -z "$2" ]
  then
    echo "Error: Must supply prefic for outfile." >&2
    exit 1
fi

if [ -z "$3" ]
  then
    echo "Error: Must specify required number of accepted simulations." >&2
    exit 1
fi


# Check $1 is actually an existing file
if [ ! -e runJulia/"$1" ]
then 
    echo "Error: Data file $1 does not exist." >&2
    exit 1
fi

echo "Making tmp dir..."
tDir="$(mktemp -d -p ~/scratch)/"

echo "Tmp dir is $tDir"

echo "Copying dir..."
cp -r runJulia $tDir

echo "Changing to tmp dir..."
cd $tDir/runJulia

echo "Running julia..."
/home/hbecher/julia-1.8.5/bin/julia distsForEddie.jl $1 $2 $3

echo "Moving results..."
mv *abc ~/scratch/results/


echo "Removing tmp dir..."
rm -rf $tDir

echo "Done."

