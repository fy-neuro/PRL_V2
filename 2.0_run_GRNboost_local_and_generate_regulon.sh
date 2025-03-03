#!/bin/bash
# Activate the conda env#slurm options
#SBATCH -p intel-sc3        #选择多个分区用逗号隔开
#SBATCH -q normal                       #Qos只能选一个，否则会报错
#SBATCH -c 12                           #申请8个CPU核心
#SBATCH -J CellChat_FY	
#SBATCH --mem=100G
#SBATCH -o %j.log                       #%j表示实际运行时的作业号ironment
#SBATCH --mail-user=fanyi@westlake.edu.cn



# module load R/4.4.0

# Get the script's directory
script_dir="/storage/liuxiaodongLab/fanyi/FY/PRL/Analysis/20241218_SCENIC_GABAergic"

# List of sample directories (modify as needed)
sample_dirs=(
    "GABAergic1220"
)

# Path to the transcription factors list relative to the script directory
tf_list="/storage/liuxiaodongLab/fanyi/FY/PRL/Analysis/20241218_SCENIC_GABAergic/SCENIC/hg38/allTFs_hg38.txt"

# Number of workers and seed (modify if needed)
num_workers=15
seed=123

# Loop over each sample directory
for sample_dir in "${sample_dirs[@]}"
do
    echo "Processing sample: $sample_dir"
    
    # Define the paths
    int_dir="$script_dir/$sample_dir/int"
    loom_file="$int_dir/${sample_dir}_exprMat_filtered.loom"
    
    # Change to the int directory inside the sample directory
    cd "$int_dir" || { echo "Failed to navigate to $int_dir"; exit 1; }
    
    # --- Python Script Execution ---
    # Make the Python script executable
    # chmod +x "$script_dir/2.1_arboreto_with_multiprocessing.py"
    
    # Run GRNBoost2
    "$script_dir/2.1_arboreto_with_multiprocessing.py" \
        "$loom_file" \
        "$tf_list" \
        --method grnboost2 \
        --output adj.tsv \
        --num_workers $num_workers \
        --seed $seed
    
    # --- End of Python Script Section ---
    
    # Return to the sample directory
    cd "$script_dir/$sample_dir" || { echo "Failed to navigate to $script_dir/$sample_dir"; exit 1; }
    
    # Run the corrected R script to generate regulons
    Rscript "$script_dir/2.2_generate_regulon.R"
    
    # Check if the R script ran successfully
    if [ $? -eq 0 ]; then
        echo "Successfully processed sample: $sample_dir"
    else
        echo "Error processing sample: $sample_dir"
        # Optionally, you can exit the script or continue with the next sample
        # exit 1
    fi
    
    # Return to the script directory
    cd "$script_dir" || { echo "Failed to return to $script_dir"; exit 1; }
    
done

echo "All samples have been processed."