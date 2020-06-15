
#!/bin/bash

dirname=/home/mwu/MING_V9T/PhD_Pro/Test/Simulation/Node2Vec_BEELINE_Data_Res/
for file in $dirname/HSC*
do
    for drop in 0.25 0.5 0.75
      do
        filename=$(echo `basename "$file"`)
	if [[ $filename == *"PCA"* ]]; then
	break
	fi
        echo $filename':'$drop
        python /home/mwu/MING_V9T/PhD_Pro/PySCNet/model_compare.py $filename 0.1 1.0 $drop
        python /home/mwu/MING_V9T/PhD_Pro/PySCNet/model_compare.py $filename 1.0 0.1 $drop
      done
done
