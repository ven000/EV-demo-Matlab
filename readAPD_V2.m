%%% convert channel F dat file to csv file, then analyze with our matlab
%%% codes.
close all
clear

[filename path]=uigetfile('*.DAT');
cd(path);
fid = fopen(filename,'r');
channel = fread(fid,'uint8');%how many data, 1D array
channel = uint8(channel);

fid = fopen(filename,'r','b');
channel2 = fread(fid,'uint32','b');%how many data, 1D array
% fclose('all');
if (mod(length(channel)-70,4) == 0)
    interest = channel(71:end);%for channel F,set the begining point as 71; for channel E, set as 72
end

if (mod(length(channel)-71,4) == 0)
    interest = channel(72:end);%for channel F,set the begining point as 71; for channel E, set as 72
end

data2=swapbytes(typecast(interest,'uint32'));
data2 = double(data2);
data3 = data2;
for i = (2:length(data2))
    data3(i,1) = data2(i,1)-data2(i-1,1);
end
data3(1) = 0;
figure(1)
plot(data3);

writematrix(data3,'change name F.csv');%change the name as you want, try not overwrite same file