
%foldFit for r square, fold used for fitting, (10 or 11; a range from 8:12)
%foldTest for calcualting fitting area percentage start at the foldTest point
%foldDist for distribution area percentage, compared with whole fitting area
function f = thresh_fit(foldFit, foldTest, foldDist, data)
%Part1

baseline = mode(data);
mad1=mad(data,1);
%disp(mad1);
%disp(baseline);
threshold = baseline + foldFit * mad1;

threshold_upper=120; %use threshold_upper if you want to choose signals with upper threshold (eg. intensity < 200)
B(:,1)=data;
for t=4:length(data)-3
    %if B(t,1)>B(t-3,1) & B(t,1)>B(t-2,1) & B(t,1)>=B(t-1,1) & B(t,1)>B(t+1,1) & B(t,1)>B(t+2,1) & B(t,1)>B(t+3,1) & B(t,1)>threshold;   %peak finding
    % if you want to use threshold_upper
    if B(t,1)>B(t-3,1) & B(t,1)>B(t-2,1) & B(t,1)>=B(t-1,1) & B(t,1)>B(t+1,1) & B(t,1)>B(t+2,1) & B(t,1)>B(t+3,1) & B(t,1)>threshold & B(t,1)<threshold_upper;
    
        B(t,3)=B(t,1);
        B(t-3,3)=0;
        B(t-2,3)=0;
        B(t-1,3)=0;
        B(t+1,3)=0;
        B(t+2,3)=0;
        B(t+3,3)=0;
    else
        B(t,3)=0;
    end
end
index1 = find(B(:,3)>0); %Index of peaks
%disp(length(index1));
index_initial = index1;
j=1;
index0=zeros(1,length(index1));%save index of too-close peaks
window = 5; % windown for removing too-close peaks. value 5 is based on bin width of peaks (5 for single ab, 
%15 for beads with high intensity eg. 1000)
for i=1:(length(index1)-1)
    if index1(i+1)-index1(i)<=window 
        if data(index1(i+1))>data(index1(i))
            index0(j)=index1(i);
            index1(i)=0;
            j=j+1;
        else
            index0(j)=index1(i+1);
            index1(i+1)=0;
            j=j+1;
        end
    end
end
index0=index0(find(index0>0));
%disp(length(index0)); %uncomment, if you wanna check occurence of
%too-close peaks, especially for bright 200nm beads (shoulder peaks)
index1=index1(find(index1>0));

%2nd round of smoothing

j=1;
index00=zeros(1,length(index1));%save index of too-close peaks
for i=1:(length(index1)-1)
    if index1(i+1)-index1(i)<=window %based on bin width of peaks
        if data(index1(i+1))>data(index1(i))
            index00(j)=index1(i);
            index1(i)=0;
            j=j+1;
        else
            index00(j)=index1(i+1);
            index1(i+1)=0;
            j=j+1;
        end
    end
end
index00=index00(find(index00>0));
%disp(length(index00));
index1=index1(find(index1>0));


%disp(length(index1));
% index3 = peak(index1);
data1=zeros(length(data),1)+mode(data);
data1(index1)=B(index1,3); %Data1 is used to check baseline estimation and peak finding results.
figure(1)
plot(data);
hold on
plot(data1,'r'); %Original data and data1 are plotted together. Check peakfinding efficiency and adjust threshold if necessary.

sumi1=data(index1)-baseline;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
interval = 4 * mad1; %default value is mad1
bins = 0:interval:max(sumi1); %Default bin width is mad1, when event frequency is low, increase bin size to 4*mad1 
[bincounts,index] = histc(sumi1,bins);

figure(9)
bar(bins,bincounts);
title('Intensity histogram of 1 point')
xlabel('Intensity')
ylabel('Occurrence')

figure(2)

hist(sumi1,1500);
title('Intensity histogram of 1 point')
xlabel('Intensity')
ylabel('Occurrence')
%set(gca, 'XScale', 'log'); %uncomment if you want to plot X axis in log.

%fitting
irf=bincounts/max(bincounts); % normalize to 1
% puddle=zeros(1,200);
% irf=[irf puddle];
x=1:length(irf);
% g = fittype('a/(b*sqrt(2*pi))*exp(-(x-c)^2/2/b^2)');% normal expression.
% [f1,gof,output] = fit(x',irf,g,'Lower',[1 1 1],'upper',[50 10 50]);

g = fittype('a/(b*x*sqrt(2*pi))*exp(-(log(x)-c)^2/2/b^2)');% log normal expression. a: scale, b:standard deviation, c:mean of the variable's natural logarithm
[f1,gof,output] = fit(x',irf,g,'Lower',[1 0 1],'upper',[60 6 6]);
%     disp(gof.rsquare);

max1 = max(sumi1)/interval;

figure(11)
bar(x,irf);
hold on
plot(f1);
%%plot(f1,x,irf);
xlim([0 2*max1])

MyCoeffs = coeffvalues(f1);
a=MyCoeffs(1);
b=MyCoeffs(2);
c=MyCoeffs(3);

whole_area=a;

cutoff=foldTest * mad1/interval;% change it to the cutoff point
detection_eff=sum(f1(cutoff:max1))/whole_area*100;
fitCount = length(index_initial)/detection_eff;
freq = length(index_initial)/length(data)*10000;
%disp(freq);
%disp(fitCount/length(data)*10000);
%disp(detection_eff);

%Part2
%based on fitted distribution (normalized to 1 for y axis), set a thresholdTest
%calculate how much percentage of detected events (including false positives above the fitting
%excluding false negative below thresholdTest.
thresholdTest = baseline + foldDist * mad1;
C(:,1)=data;
for t=4:length(data)-3
    if C(t,1)>C(t-3,1) & C(t,1)>C(t-2,1) & C(t,1)>=C(t-1,1) & C(t,1)>C(t+1,1) & C(t,1)>C(t+2,1) & C(t,1)>C(t+3,1) & C(t,1)>thresholdTest;   %peak finding
        C(t,3)=C(t,1);
        C(t-3,3)=0;
        C(t-2,3)=0;
        C(t-1,3)=0;
        C(t+1,3)=0;
        C(t+2,3)=0;
        C(t+3,3)=0;
    else
        C(t,3)=0;
    end
end
index2 = find(C(:,3)>0); %Index of peaks
%disp(length(index2));

jj=1;
index20=zeros(1,length(index2));%save index of too-close peaks
for i=1:(length(index2)-1)
    if index2(i+1)-index2(i)<=5 %based on bin width of peaks
        if data(index2(i+1))>data(index2(i))
            index20(jj)=index2(i);
            index2(i)=0;
            jj=jj+1;
        else
            index20(jj)=index2(i+1);
            index2(i+1)=0;
            jj=jj+1;
        end
    end
end
index20=index20(find(index20>0));
%disp(length(index20));
index2=index2(find(index2>0));

%2nd round of smoothing

jj=1;
index20=zeros(1,length(index2));%save index of too-close peaks
for i=1:(length(index2)-1)
    if index2(i+1)-index2(i)<=5 %based on bin width of peaks
        if data(index2(i+1))>data(index2(i))
            index20(jj)=index2(i);
            index2(i)=0;
            jj=jj+1;
        else
            index20(jj)=index2(i+1);
            index2(i+1)=0;
            jj=jj+1;
        end
    end
end
index20=index20(find(index20>0));
%disp(length(index20));
index2=index2(find(index2>0));
%disp(length(index2));

sumi21=data(index2)-baseline;
%interval = mad1; %default value is 2
bins2 = 0:interval:max(sumi21); %Default bin number is 2, modify it to change bin size.
[bincounts2,index] = histc(sumi21,bins2);

% figure(19)
% bar(bins2,bincounts2);
% title('Intensity histogram of 1 point')
% xlabel('Intensity')
% ylabel('Occurrence')

irf2=bincounts2/max(bincounts); % normalize to 1
x=foldDist:length(irf2);
% figure(18)
% bar(x,irf2(x));

total=sum(irf2(x));
percent = total/whole_area;

%f = [foldFit, gof.rsquare, detection_eff, foldTest, foldDist, percent, freq, baseline, mad1];
f = [foldFit, gof.rsquare, detection_eff, freq, baseline, mad1];
end







