# store the location of the home directory
homedir=$PWD

# update the submodules
./external/ExaHyPE-Engine/Submodules/updateSubmodules.sh

# the directory where we generate the ExaHyPE source code
gendir=spipack/ConservationEquations/ExaHyPE/GenerateSource

# the directory where we store the source for the ExaHyPE project
sourcedir=modules/ConservationEquations/src/ExaHyPE

Generate() {
  # create the source code for the generated project
  ./external/ExaHyPE-Engine/Toolkit/toolkit.sh spipack/ConservationEquations/ExaHyPE/Generate$1.exahype

  # make the generated project
  cd spipack/ConservationEquations/ExaHyPE/GenerateSource
  export MODE=Release
  export COMPILER=GNU
  make -j10

  cmk='@CMAKE_CURRENT_SOURCE_DIR@'
  sed -i 's,'"$homedir"','"$cmk"',' buildinfo.h
  mv buildinfo.h $homedir/cmake/exahype/$1/buildinfo.h.in

  # copy the generated source if it does not yet exist
  cp -n *.h ../$1
  cp -n *.cpp $homedir/$sourcedir/$1

  # read the files that we need to compile
  value=`cat cfiles.mk`

  # edit the files so that to use paths relative to the home dir
  value=$(echo $value | cut -c 10-)
  value=${value// /;}
  dir="$sourcedir/$1"
  value=${value//$gendir/$dir}
  value=${value//$homedir/'${CMAKE_CURRENT_SOURCE_DIR}'}

  # output the list of source files to cmake
  value="set(EXAHYPE_$1_SOURCE $value)"
  echo $value > $homedir/cmake/exahype/$1/ExaHyPEFiles.cmake

  # we don't actually need the build, so just clean it up
  make clean
  rm *

  # move back to the home dir
  cd $homedir
}

#Generate CompressibleEuler
Generate Boltzmann

echo ""
echo "Finished generating code using ExaHyPE!"
