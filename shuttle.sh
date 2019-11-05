PBS_fl=$1;

echo -e "[INFO] using PBS '$PBS_fl' for submitting.";
# Read line by line
while read line;do
# echo $line; 
# Get comma-separated parameters
IFS="," read METRIC METHOD NGENES <<< $line;
# Print information about the job submission
JOBNAME="${METRIC}_${METHOD}_${NGENES}";
echo -e "Submitting job '$JOBNAME'";
# Submit job
qsub ${PBS_fl} -v METRIC=$METRIC,METHOD=$METHOD,NGENES=$NGENES -N "$JOBNAME";
done < settings.csv

