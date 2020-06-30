#export OMP_NUM_THREADS=24

#        pcout << "The structure of the input parameters should be : "
#              << (string)global::argv(0) << " tmp Peclet (e.g. 10) Da (e.g. 0.5) deltaRho (e.g. 1e-5) thickness (eg. 6 for real case) ...
# flag (1 for starting from steady state, 0 for starting from beginning) geometry.txt" << endl;;


#OUT_F=/cluster/work/petro/pandrea/kong_paolo/precipi

MY_EXE_1=singleSpecies_BGK_precipitation_pseudo_3D

#bsub -n 48 -W 24:00 -o Pe1Da01.txt -J Pe1Da01 -u pandrea  mpirun ./$MY_EXE_1 $OUT_F/test1 1 0.1 1e-2 ./myGeometry_P.txt   
#bsub -n 48 -W 24:00 -o Pe1Da01.txt -J Pe1Da01 -u pandrea  mpirun $MY_EXE_1 $OUT_F/test1 1 0.1 1e-2 ./myGeometry_P.txt

#OUT_F2=/cluster/work/petro/pandrea/kong_paolo/precipi

OUT_F2=/cluster/scratch/mahkami/3D_sim_Pe1

#OUT_F2=$ProjectFolder/3D_sim_Pel

#bsub -n 48 -W 24:00 -o Pe10Da01.txt -J Pe10Da01 -u pandrea  mpirun ./$MY_EXE_1 $OUT_F2/test2 10 0close all.1 1e-2 ./myGeometry_P.txt   
#bsub -n 48 -W 24:00 -o Pe10Da01.txt -J Pe10Da01 -u pandrea  mpirun $MY_EXE_1 $OUT_F2/test2 10 0.1 1e-2 ./myGeometry_P.txt


# Don't forget to check the geometry!!!!!!!!!
#Geometry_m is for the short boundary on left and right and Geometry_m2 is for the extended boundary on left and right


bsub -n 48 -W 20:00 -R "rusage[mem=2048]" -o Pe1Da0001.txt -J Pe1Da0001 -u mahkami  mpirun $MY_EXE_1 $OUT_F2/Pe1Da0001 1 0.001 1e-1 60 1 ./myGeometry_m2.txt

#bsub -n 24 -W 40:00 -R "rusage[mem=2048]" -o Pe1Da0005.txt -J Pe1Da0005 -u mahkami  mpirun $MY_EXE_1 $OUT_F2/Pe1Da0005 1 0.005 1e-1 60 1 ./myGeometry_m2.txt

bsub -n 48 -W 20:00 -R "rusage[mem=2048]" -o Pe1Da001.txt -J  Pe1Da001  -u mahkami  mpirun $MY_EXE_1 $OUT_F2/Pe1Da001  1 0.01 1e-1 60 1  ./myGeometry_m2.txt

#bsub -n 24 -W 40:00 -R "rusage[mem=2048]" -o Pe1Da005.txt -J  Pe1Da005  -u mahkami  mpirun $MY_EXE_1 $OUT_F2/Pe1Da005  1 0.05 1e-1 60 1  ./myGeometry_m2.txt

bsub -n 48 -W 20:00 -R "rusage[mem=2048]" -o Pe1Da01.txt -J   Pe1Da01   -u mahkami  mpirun $MY_EXE_1 $OUT_F2/Pe1Da01   1 0.1 1e-1 60 1   ./myGeometry_m2.txt
    
#bsub -n 24 -W 40:00 -R "rusage[mem=2048]" -o Pe1Da05.txt -J   Pe1Da05   -u mahkami  mpirun $MY_EXE_1 $OUT_F2/Pe1Da05   1 0.5 1e-1 60 1   ./myGeometry_m2.txt

bsub -n 48 -W 20:00 -R "rusage[mem=2048]" -o Pe1Da1.txt -J    Pe1Da1    -u mahkami  mpirun $MY_EXE_1 $OUT_F2/Pe1Da1    1 1 1e-1 60 1     ./myGeometry_m2.txt

#bsub -n 24 -W 40:00 -R "rusage[mem=2048]" -o Pe1Da1_5.txt -J  Pe1Da1_5  -u mahkami  mpirun $MY_EXE_1 $OUT_F2/Pe1Da1_5  1 1.5 1e-1 60 1   ./myGeometry_m2.txt
   
#bsub -n 24 -W 40:00 -R "rusage[mem=2048]" -o Pe1Da2.txt -J    Pe1Da2    -u mahkami  mpirun $MY_EXE_1 $OUT_F2/Pe1Da2    1 2 1e-1 60 1     ./myGeometry_m2.txt

bsub -n 48 -W 20:00 -R "rusage[mem=2048]" -o Pe1Da2_5.txt -J  Pe1Da2_5  -u mahkami  mpirun $MY_EXE_1 $OUT_F2/Pe1Da2_5  1 2.5 1e-1 60 1   ./myGeometry_m2.txt

#bsub -n 24 -W 40:00 -R "rusage[mem=2048]" -o Pe1Da3.txt -J    Pe1Da3    -u mahkami  mpirun $MY_EXE_1 $OUT_F2/Pe1Da3    1 3 1e-1 60 1     ./myGeometry_m2.txt

#bsub -n 24 -W 40:00 -R "rusage[mem=2048]" -o Pe1Da3_5.txt -J  Pe1Da3_5  -u mahkami  mpirun $MY_EXE_1 $OUT_F2/Pe1Da3_5  1 3.5 1e-1 60 1   ./myGeometry_m2.txt

bsub -n 48 -W 20:00 -R "rusage[mem=2048]" -o Pe1Da4.txt -J    Pe1Da4    -u mahkami  mpirun $MY_EXE_1 $OUT_F2/Pe1Da4    1 4 1e-1 60 1     ./myGeometry_m2.txt

bsub -n 48 -W 20:00 -R "rusage[mem=2048]" -o Pe1Da0.txt -J    Pe1Da0    -u mahkami  mpirun $MY_EXE_1 $OUT_F2/Pe1Da0    1 0 1e-1 60 1     ./myGeometry_m2.txt
