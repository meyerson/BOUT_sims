qsub -q development -pe 16way 16 -N q3_simp -cwd -V -m e -A BOUT++_startup -l h_rt=2:00:00 -j y ./runcase.sh

#qsub -q development -pe 16way 16 -N q3_simp_re_debug -cwd -V -m e -A TG-PHY090081 -l h_rt=2:00:00 -j y ./runcase_re.sh
