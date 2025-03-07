# This script sets up and runs compilation of the Bering BGC ROMS executables

# Step 1: Copy necessary files from ROMS and App folders

myappfol="/gscratch/bumblereem/kearney/bering-Apps/Apps/Bering_BGC_variants"
myromssrc="/gscratch/bumblereem/kearney/roms-kate-ice"

cp ${myappfol}/build_roms.bash .
cp ${myappfol}/ana_psource.h .
cp ${myappfol}/bering_10k.h .

# Set project directory (will be used in place of value in build_roms.sh)

export MY_PROJECT_DIR="/gscratch/bumblereem/kearney/BGC_hindcasts_workdir"
export MY_ROMS_SRC="/gscratch/bumblereem/kearney/roms-kate-ice"

# Build ROMS in all 4 variants

# ./build_roms.bash -variant phys
./build_roms.bash -variant bestnpz

export MY_CPP_FLAGS="-DNLIM_ONEILL -DCOBALT_ALLBURIAL"
./build_roms.bash -variant cobalt
# ./build_roms.bash -variant banas








