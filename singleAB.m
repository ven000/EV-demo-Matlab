close all
clear


pathnam='C:\Users\80933\Desktop\20210828_cd9\abcam\3'; % input for data file path, only one step needs user input.
FileName=fullfile(pathnam,'*.csv');
[fnam, pathnam, filterindex] = uigetfile(FileName, 'Pick a csv-file');
fullname = fullfile(pathnam, fnam);
fileID=fopen(fullname);
data=fscanf(fileID,'%d');
%data=data(1:600000);
%data=data(600001:1200000);
%data=data(1200001:1800000)'; %Uncomment if you want to crop the trajectory.

foldFit = 8; % the fold of mad for setting threshold in part 1
foldDist = 10; % compared with the best fitting situation, the percentage of estimated false positives
% from real distribution (before fitting) at a specific intensity
% foldTest = foldFit; 
rsquareBest = 0;
foldBest = 0;
for i = 1 : length(foldFit)
    output(i,:) = thresh_fit(foldFit(i), foldFit(i), foldDist, data);
    if (output(i,2)>rsquareBest)
        rsquareBest = output(i,2);
        foldBest = foldFit(i);
    end
   % output(i,:) = thresh_fit(foldFit(i), foldTest, foldDist, data);  
end
%g1 = 'foldFit, gof.rsquare, detection_eff(%), foldTest, foldDist, percent, frequency(Hz), baseline, mad\n';
g1 = 'foldFit, gof.rsquare, detection_eff(%%), frequency(Hz), baseline, mad\n';
fprintf(g1);
disp(output);
disp([foldBest, rsquareBest]);

