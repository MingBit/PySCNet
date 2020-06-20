
#!/bin/bash

dirname=/home/mwu/MING_V9T/PhD_Pro/Test/Simulation/BoolODE_Data/BoolODE_Data_Repeat_Jun2020/
for file in $dirname/HSC_MEP_C*
do
    for drop in 0.75
      do
        filename=$(echo `basename "$file"`)
	if [[ $filename == *"PCA"* ]]; then
	break
	fi
        echo $filename':'$drop
        #python /home/mwu/MING_V9T/PhD_Pro/SCNode2Vec/Simulation_GRN.py $filename 0.1 1.0 $drop
        #python /home/mwu/MING_V9T/PhD_Pro/SCNode2Vec/Simulation_GRN.py $filename 1.0 0.1 $drop
      done
done
