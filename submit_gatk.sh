for EmbryoID in `ls /your/path/Data/Samples`
do
echo $EmbryoID
sbatch /your/path/NEWgatk.sh $EmbryoID
done
