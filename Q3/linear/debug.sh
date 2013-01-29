runname=restart++
current_dir=/scratch/01523/meyerson/Q3/linear/data_re_eig64_.5
old_dir=/scratch/01523/meyerson/Q3/linear/data_eig64_.5
NP=16
rm $current_dir/BOUT.dmp*
cp $old_dir/BOUT.restart.* $current_dir

qsub -q development -pe 16way 16 -N $runname -cwd -V -m e -A BOUT++_startup -l h_rt=2:00:00 -j y ./subrun.sh 16 ./q3_simp $current_dir -re
