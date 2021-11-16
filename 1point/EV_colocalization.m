%%%%find best threshold for free antibody fitting (1 point for intensity), and then deconvolution
close all
clear

pathnam='C:\Users\80933\Desktop\20210828_cd9\abcam\3\'; % input for data file path, only one step needs user input.

FirstFileName=fullfile(pathnam,'*.csv');
[fnam, pathnam, filterindex] = uigetfile(FirstFileName, 'Load the first track');
fullname1 = fullfile(pathnam, fnam);
SecondFileName=fullfile(pathnam,'*.csv')
[fnam, pathnam, filterindex] = uigetfile(SecondFileName, 'Load the second track');
fullname2 = fullfile(pathnam, fnam);
%data=data(1:600000);
%data=data(600001:1200000);
%data=data(1200001:1800000)'; %Uncomment if you want to crop the trajectory.
fileID1=fopen(fullname1);
data1=fscanf(fileID1,'%d');
%data1=data1(1:600000);
fileID2=fopen(fullname2);
data2=fscanf(fileID2,'%d');
%data2=data2(1:600000);

foldFit = 9:12; % the fold of mad for setting threshold in part 1
rsquareBest = 0;
foldBest = 0;
for i = 1 : length(foldFit)
    output(i,:) = colocalization1(foldFit(i), data1, data2);
    if (output(i,2)>rsquareBest)
        rsquareBest = output(i,2);
        foldBest = foldFit(i);
    end
   % output(i,:) = thresh_fit(foldFit(i), foldTest, foldDist, data);  
end
%g1 = 'foldFit, gof.rsquare, detection_eff(%), foldTest, foldDist, percent, frequency(Hz), baseline, mad\n';
g1 = 'fold, rsquare, detection_eff(%%)\n';
fprintf(g1);
disp(output);
disp([foldBest, rsquareBest]);
colocalization2(foldBest, data1, data2)

