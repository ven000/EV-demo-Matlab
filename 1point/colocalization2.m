%%%%%%%%%%%%%%%function for EV_colocalization
function f = colocalization1(foldFit, data1, data2)
%part 1

baseline1 = mode(data1);
% 15 for 200nm beads, it should be adjusted based on samples
threshold1 = baseline1+15*mad(data1,1);
%threshold1 = baseline1+10*mad(data1,1); %Run flowhist first to determine threshold for peak finding. Input threshold for trajectory 1(membrane dye).
B(:,1)=data1;
for t=4:length(data1)-3
if B(t,1)>B(t-3,1) & B(t,1)>B(t-2,1) & B(t,1)>B(t-1,1) & B(t,1)>B(t+1,1) & B(t,1)>B(t+2,1) & B(t,1)>B(t+3,1) & B(t,1)>threshold1 %Peak finding loop. compare with nearby 6 points.  
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
index4 = find(B(:,3)>0); %index4 is for peak position in trajectory 1.

%data2=data2(1:600000);
figure(1)
plot(data2)
baseline2 = mode(data2);
threshold2 = baseline2+foldFit*mad(data2,1); %Run flowhist first to determine threshold for peak finding.
B(:,1)=data2;
for t=4:length(data2)-3
if B(t,1)>B(t-3,1) & B(t,1)>B(t-2,1) & B(t,1)>B(t-1,1) & B(t,1)>B(t+1,1) & B(t,1)>B(t+2,1) & B(t,1)>B(t+3,1) & B(t,1)>threshold2   %check  (channel B) peak height
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
index5 = find(B(:,3)>0);%index5 is for peak position in trajectory 2.
[acor, lag]=xcorr(data1,data2);% Run cross-correlation between the two raw data.
figure(2)
delay =lag(find(acor==max(acor)));%Use correlation peak delay time to determine transition delay time between the two channels.

delayindex = find(acor==max(acor));%x value of peak
% [acor, lag]=xcorr(data1-baseline1,data2-baseline2); %Subtract baseline before correlation analysis. Debug only.
plot(lag,acor);
xlim([-200 300])

corrbase = mean(acor(delayindex-80:delayindex-30));%baseline
meandelay = lag(delayindex-50:delayindex+50)*(acor(delayindex-50:delayindex+50)-corrbase)/sum(acor(delayindex-50:delayindex+50)-corrbase);%???
title('Crosscorrelation plot');
fprintf('most probable delay time is %4.2f bins \n', delay);
fprintf('averaged delay time is %4.2f ms, stop if inconsistent with figure 1\n', meandelay);

colocalize1_t = zeros(1,length(index4));%Empty array for colcolized event index in data1(membranedye).
colocalize2_t = zeros(1,length(index4));%Empty array for colcolized event index in data2(antibody).

delay=20;
left=20;
right=80; % defalut: delay=40;left=20; right=40;

% for i=1:length(index4)
%     for j=1:length(index5)
%          if index5(j)+delay+right>=index4(i)&&index4(i)>=index5(j)+delay-left % shift the two trajectories for colocalization.Comments and uncomments to adjust delay time range based on correlation plot.
%              colocalize1_t(i)=index4(i);
%              colocalize2_t(i)=index5(j);  
%              
%          end
%          
%     end
% end

for j=1:length(index5)
    for i=1:length(index4)
         if index5(j)+delay+right>=index4(i)&&index4(i)>=index5(j)+delay-left % shift the two trajectories for colocalization.Comments and uncomments to adjust delay time range based on correlation plot.
             colocalize1_t(j)=index4(i);
             colocalize2_t(j)=index5(j);
             
             break 
         end
         
    end
end

real_event_ab = colocalize2_t(find(colocalize2_t>0));%Sort out non-zero index for peaks in data2. a647
sumi2_ab=zeros(1,length(real_event_ab));
for i=1:length(real_event_ab)
sumi2_ab(1,i)=(data2(real_event_ab(i))-baseline2)+(data2(real_event_ab(i)-1)-baseline2)+(data2(real_event_ab(i)+1)-baseline2);%3 point sum intensity for peaks in data2.
end

real_event1_dye = colocalize1_t(find(colocalize1_t>0));%Sort out non-zero index for peaks in data1.
sumi1_dye=zeros(1,length(real_event1_dye));
delay1_t = zeros(1,length(real_event1_dye));
for i=1:length(real_event1_dye)
sumi1_dye(1,i)=(data1(real_event1_dye(i))-baseline1)+(data1(real_event1_dye(i)-1)-baseline1)+(data1(real_event1_dye(i)+1)-baseline1);%3 point sum intensity for peaks in data1.
delay1_t(1,i)=real_event1_dye(i)-real_event_ab(i);%anepps-a647 index
end

%figure(3)
% a=max(delay1_t);
% b=min(delay1_t);
% nbins=a-b;
%histogram(delay1_t,nbins);
%title('delay time destribution in a estimated range')
%xlabel('delay time (bins)');
%ylabel('counts');

delay_p=mode(delay1_t);%chose the most frequency delay time
flow_rate=10/(delay_p/10)*3;%channel ef is 10, channel eh is 20
fprintf('flow rate is %4.2f pL/s \n', flow_rate);
fprintf('most probable delay time is %4.2f bins \n', delay_p);
%part 2

colocalize1 = zeros(1,length(index4));%Empty array for colcolized event index in data1(membranedye).
colocalize2 = zeros(1,length(index4));%Empty array for colcolized event index in data2(antibody).


%%%%%%%index4 remove too close peaks in channel E for beads, eg. shoulder peak
%disp(length(index4));
window = 15; %windown for removing too-close peaks. value 15 for exosomes and beads with high intensity eg. 1000)
for i=1:(length(index4)-1)
    if index4(i+1)-index4(i)<=window %based on bin width of peaks
        if data1(index4(i+1))>=data1(index4(i))
            index4(i)=0;
        else
            index4(i+1)=0;
        end
    end
end

index4=index4(find(index4>0));
%disp(length(index4));

for i=1:(length(index4)-1)
    if index4(i+1)-index4(i)<=window %based on bin width of peaks
        if data1(index4(i+1))>=data1(index4(i))
            index4(i)=0;
        else
            index4(i+1)=0;
        end
    end
end

index4=index4(find(index4>0));
%disp(length(index4));
%disp(length(index5));
for j=1:(length(index5)-1)
    if index5(j+1)-index5(j)<=window %based on bin width of peaks
        if data2(index5(j+1))>=data2(index5(j))
            index5(j)=0;
        else
            index5(j+1)=0;
        end
    end
end

index5=index5(find(index5>0));

%disp(length(index5));

for j=1:(length(index5)-1)
    if index5(j+1)-index5(j)<=window %based on bin width of peaks
        if data2(index5(j+1))>=data2(index5(j))
            index5(j)=0;
        else
            index5(j+1)=0;
        end
    end
end

index5=index5(find(index5>0));

%disp(length(index5));


% default
% left=5;
% right=40;
% delay_p=-50;
left=15;
right=200;

a=1;%for smart colocalization, removing redundant countings

for j=1:length(index5)
    for i=a:length(index4)
         if index5(j)+delay_p+right>=index4(i)&&index4(i)>=index5(j)+delay_p-left % shift the two trajectories for colocalization.Comments and uncomments to adjust delay time range based on correlation plot.
             colocalize1(j)=index4(i);
             colocalize2(j)=index5(j);
             a=i+1;
             break
        
         end
    end
end

realevent = colocalize2(find(colocalize2>0));%Sort out non-zero index for peaks in data2.
sumi2=zeros(1,length(realevent));


for i=1:length(realevent)
sumi2(1,i)=(data2(realevent(i))-baseline2);%1 point intensity for peaks in data2.
end

realevent1 = colocalize1(find(colocalize1>0));%Sort out non-zero index for peaks in data1.
sumi1=zeros(1,length(realevent1));
delay_t=zeros(1,length(realevent));

for i=1:length(realevent1)
sumi1(1,i)=data1(realevent1(i))-baseline1;
delay_t(1,i)=realevent1(i)-realevent(i);
end

figure(21)
scatter(delay_t,sumi2);
title('channel H intensity 1-point vs delay time')
xlabel('delay time')
ylabel('intensity')
%sqrtI=sqrt(sumi1); %take square root of membrane dye intensity.

figure(3)
a=max(delay_t);
b=min(delay_t);
nbins=a-b;
histogram(delay_t,nbins);
title('exact delay time destribution')
xlabel('delay time (bins)');
ylabel('counts');

figure(80)
i=1:length(realevent);
scatter(i,sumi2);
title('channel H intensity vs event count')
xlabel('event index')
ylabel('intensity')

figure(81)
scatter(realevent,sumi2);
title('channel H intensity vs event count')
xlabel('time (bins)')
ylabel('intensity')

figure(22)
scatter(sumi1,sumi2);
title('channel H vs Anepps intensity ');
xlabel('Anepps intensity')
ylabel('channel H intensity')
set(gca,'XScale','log','YScale','log')
xlim([1 10000])
ylim([1 10000])

% figure(4)
% hist(sumi2,160); %Intensity histogram without compensation. Default bin number is 160, modify it to change bin size.%
% title('Intensity Histogram of Colocalized Antibody Events')


delaywidth = left+right;
fprintf('number of bins in the correlation peak used for colocalization = %4.2f\n', delaywidth);

colocalpercentage = length(realevent)/length(index4);%index4 is for peak position in trajectory 1. membrane dye
expresspercentage = (length(realevent)-length(index4)*length(index5)*delaywidth/length(data2))/length(index4)/(1-length(index4)*delaywidth/length(data2));%
randompercentage = (length(index5)-length(index4)*expresspercentage)*delaywidth/length(data2);

fprintf('colocalized percatge of all events in 488 channel is %4.2f\n', colocalpercentage);
fprintf('expected random colocalization percatge is %4.2f\n', randompercentage);
fprintf('percatge of vesicles with target protein is %4.2f\n', expresspercentage);

bins = 0:8:max(sumi2); %Default bin number is 50, modify it to change bin size.
[bincounts,ind] = histc(sumi2,bins);%3 point sum intensity for peaks in data2.
figure(4) %Intensity histogram without compensation.
plot(bins,bincounts);
xlim([0 500])
title('Intensity Histogram of Colocalized Antibody Events wo compensation')

figure(7)%Bin counts plot for Sumi1.
bins1 = 0:25:max(sumi1); %Default bin number is 50, modify it to change bin size.
[bincounts1,ind1] = histc(sumi1,bins1);
plot(bins1,bincounts1);
xlim([0 500])
pmatrix1=[bins1',bincounts1'/sum(bincounts1)]';
title('Intensity Histogram Bincounts Plot of Colocalized Membrane Dye Events');

data3=data2;
data3(realevent)=-10000;% Set colocalized events peak intensity to negative values. yifei original
sumi3=zeros(1,length(index5));

index_ab = setxor(index5,realevent);

for i=1:length(index_ab)
sumi3(1,i)=data3(index_ab(i))-baseline2;
end

k=find(sumi3>0);
fprintf('Total number for free antibody is %4.0f\n', length(k));

sumi4=zeros(1,length(k));% free dye (only positives)
for i=1:length(k)
sumi4(1,i)=sumi3(k(i));
end


[bincounts2,ind2] = histc(sumi3(find(sumi3>0)),bins);%Sort out postive peaks, which should be free AB signal.

%make a cutoff for free AB to remove tail problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bincounts2(find(bincounts2 == max(bincounts2))*4:end)=0;
figure(5) %Free AB events analysis.
plot(bins,bincounts2);
xlim([0 500])
pmatrix2=[bins',bincounts2'/sum(bincounts2)]';
averageAB=bins*bincounts2'/sum(bincounts2);%Free AB intensity average.
fprintf('Averaged Intensity of Free Antibody Events After Correction is %4.2f\n', averageAB);
title('Intensity Histogram Bincounts Plot of Free Antibody Events');
ABcon = length(find(sumi2>0))/length(data2);% Free AB concentration
fprintf('Antibody concentration is %4.2f\n', ABcon);


correctbincounts=round(bincounts-1*bincounts2*(length(realevent)*randompercentage/colocalpercentage/sum(bincounts2)));%length(realevent)*randompercentage/colocalpercentage=random free dye colocalization !ratio:randompercentage/colocalpercentage
bins_tp=bins';
correctbincounts_tp=correctbincounts';
totalcolev_tp=sumi2';
freedye_tp=sumi4';



figure(6)%Use free AB intensity histogtam and random colocalization percentage to do compensation.
plot(bins,correctbincounts);
title('Intensity Histogram Bincounts Plot of Colocalized Antibody Events After Correction');

correctbincounts_pos = correctbincounts(find(correctbincounts>0));
random_ev_col=sum(correctbincounts_pos);


averageI=bins*correctbincounts'/sum(correctbincounts);
fprintf('Averaged Intensity of Colocalized Antibody Events After Correction is %4.2f\n', averageI);
fprintf('Averaged Copy Number is %4.2f\n', averageI/averageAB);
% disp(fnam); 
% disp('total colocalized events');
% disp(length(realevent));
% disp('random ev colocalization, positive events after correction');
% disp(random_ev_col);
% disp('random free ab colocalization');
% disp(length(realevent)-random_ev_col);
% disp('random colocalization of free dyes based on probability of data1 and 2')
% disp(length(realevent)-sum(correctbincounts));
% disp('sum(correctbincounts), if its negative means free dye compensation too much')
% disp(sum(correctbincounts));
% disp('total free ab=uncolocalized ab + random colocalized ab');
% disp(length(realevent)-sum(correctbincounts)+length(k));

%part 3 

tic;
% bincounts3=round(bincounts-bincounts2*(length(realevent)*randompercentage/colocalpercentage/sum(bincounts2)));
irf=bincounts2/max(bincounts2); 
% puddle=zeros(1,200);
% irf=[irf puddle];
x=1:length(irf);
g = fittype('a/(b*x*sqrt(2*pi))*exp(-(log(x)-c)^2/2/b^2)');% log normal expression
[f1,gof,output] = fit(x',irf',g,'Lower',[0 0 0],'upper',[20 3 15]);
disp(gof);
figure(11)
%plot(f1);

bar(x,irf);
hold on
plot(f1,x,irf);

xlim([0.1 50]);
MyCoeffs = coeffvalues(f1);
a=MyCoeffs(1);
b=MyCoeffs(2);
c=MyCoeffs(3);
parti=zeros(30,length(irf));
f = [gof.rsquare, foldFit];

for i=1:30
    for j=1:length(irf)
parti(i,j)=a/(b*(j/i)*sqrt(2*pi))*exp(-(log(j/i)-c)^2/2/b^2);
    end
parti(i,:)=parti(i,:)/sum(parti(i,:));
end
coeffs=zeros(1,30);
coeffs(1)=0.5;
EVii=correctbincounts/max(correctbincounts);
i=1;
sumresidue=0;
while i<20000 %default as 20000;
    product=coeffs*parti;
    residue=EVii-product;
    figure(12)
    plot(residue);
    drawnow
    sumresidue(i)=sum(abs(residue));
    P=rand(1);
    Q=rand(1);
    coeffs1=coeffs;
    if Q<0.5
    coeffs1((round(1+P*29)))=coeffs1(round(1+P*29))+0.1;
    else
    coeffs1((round(1+P*29)))=coeffs1(round(1+P*29))-0.1;
    end
    product=coeffs1*parti;
    residue1=EVii-product;
    sumresidue1=sum(abs(residue1));
    if sumresidue1<min(sumresidue) && coeffs1((round(1+P*29)))>0
    coeffs=coeffs1;
    end
    i=i+1;
end

figure(13)
bar(coeffs);
copy=(1:length(coeffs))*coeffs'/sum(coeffs);
fprintf('copy number is %4.2f\n', copy);
toc;
end
