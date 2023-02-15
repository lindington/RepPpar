Yes_No ()
{
  # print question
  echo -n "This runs the NGSadmix replication wrapper for specified K values and beagle files, submitting each as a SLURM job. Ensure that you have backed up any previous outputs with matching names as they will be overwritten. Would you like to continue? (yes/no)."

  echo
  
  # read answer
  read YnAnswer

  # all to lower case
  YnAnswer=$(echo $YnAnswer | awk '{print tolower($0)}')

  # check and act on given answer
  case $YnAnswer in
    "yes")  Start_Install ;;
	"y") Start_Install ;;
	"no") exit 0 ;;
	"n") exit 0 ;;
    *)      echo "Please answer yes or no" ; Yes_No ;;
  esac
}

Start_Install ()
{

#Set the output directory.
OUTPUT_DIR="01.output"

if [ ! -d ${OUTPUT_DIR} ]; then

	mkdir ${OUTPUT_DIR}

fi

if [ ! -d 00.slurm_scripts ]; then

	mkdir 00.slurm_scripts
	
fi

#Set which K values you would like to run the analysis for.
for K in {2..10};

	do
		
		#Set the beagle files you would like to run the analysis on. Do not include the .beagle.gz extension.
		for beagle_file in "219_out_maxdepth"

			do
				sbatch -J K${K}_${beagle_file} -M biohpc_gen -p biohpc_gen_production -t 4-00:00:00 --mail-type=FAIL --mail-user=hagberg@lmu.de -o 00.slurm_scripts/${beagle_file}_K${K}_50reps.out -- wrapper_ngsAdmix.sh -likes ../03.beagle/01.output/${beagle_file}.beagle.gz -K ${K} -P 4 -o ${OUTPUT_DIR}/${beagle_file}_K${K}_50reps -minMaf 0.05 
			done

	done

echo "done..."

}

Yes_No
