#!/bin/sh
cd $PBS_O_WORKDIR

# e.g. 4
echo ${PBS_ARRAYID} 

# eg 58973[4].moby.cs.adelaide.edu.au
echo ${PBS_JOBID}

# eg KITTI-4
echo $PBS_JOBNAME


matlab -nodesktop -nosplash -r "runGridSearchCluster('$PBS_JOBNAME','$PBS_ARRAYID')"

