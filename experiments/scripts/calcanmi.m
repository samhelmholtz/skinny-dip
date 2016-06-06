function calcanmi(truelabelsfile,predictedlabelsfile,outputFilename)

tlf=fopen(truelabelsfile,'r');
plf=fopen(predictedlabelsfile,'r');
fid=fopen(outputFilename,'w');

tline=fgets(tlf);
pline=fgets(plf);
while ischar(tline)
targetLabels=strread(tline);
predictedLabels=strread(pline);
targetLabels=targetLabels+(1-min(targetLabels)); % offset so that they start at 1 (ANMI calculation complains if there are nonpositive values)
predictedLabels=predictedLabels+(1-min(predictedLabels)); % offset so that they start at 1 (ANMI calculation complains if there are nonpositive values)
anmi=ANMI_analytical_11(predictedLabels,targetLabels);
numTargetClusters=size(unique(targetLabels),2);
numPredictedClusters=size(unique(predictedLabels), 2);
fprintf(fid,'%d,%d,%f\n',numTargetClusters, numPredictedClusters, anmi);
tline=fgets(tlf);
pline=fgets(plf);
end

fclose(fid);
fclose(plf);
fclose(tlf);
exit;
  
end
