%%% convert channel F dat file to csv file, then analyze with our matlab
%%% codes.
close all
clear

[filename path]=uigetfile('*.DAT');
cd(path);
fid = fopen(filename,'r');
Channel = fread(fid,'uint8');%how many data, 1D array
Channel = uint8(Channel);

fid = fopen(filename,'r','b');
channel2 = fread(fid,'uint32','b');%how many data, 1D array
% fclose('all');
interest = Channel(71:end);%for channel F,set the begining point as 71; for channel E, set as 72
data2=swapbytes(typecast(interest,'uint32'));
data2 = double(data2);
data3 = data2;
for i = (2:length(data2))
    data3(i,1) = data2(i,1)-data2(i-1,1);
end

figure(1)
plot(data3);

writematrix(data3,'channel F.csv');%change the name as you want, try not overwrite same file