%% Load 

clear all; close all;
load('/Users/joshsalvi/Desktop/BullfrogPlasmaOsm.mat')

%% Input data

ind = 2;
raw_values = [194 195];
comment = 'F1/(H)';

% Calculate mean and error
data(ind).comments = comment;
data(ind).rawdata = raw_values;
data(ind).means = mean(data(ind).rawdata);
data(ind).length = length(data(ind).rawdata);
data(ind).std = std(data(ind).rawdata);
data(ind).sem = data(ind).std/sqrt(data(ind).length);



%% Calculate mean of means and propagated error


avgdata.comments = 'Averaged data';
avgdata.mean = mean(horzcat([data.means]));
avgdata.std = sqrt(sum(horzcat([data.std]).^2));
avgdata.sem = sqrt(sum(horzcat([data.sem]).^2));
avgdata.iqr = iqr(horzcat([data.means]));

%% Plot the data
clear rawdata rawindices
for j = 1:length(data)
    rawdata(1:data(j).length,j) = data(j).rawdata;
    for k = 1:data(j).length
        rawindices{k,j} = num2str(data(j).comments);
    end
end
rawdata(rawdata==0)=[];
sizeR = size(rawdata);
rawdata=reshape(rawdata,1,sizeR(1)*sizeR(2));
rawindices=reshape(rawindices(1:length(rawdata)),1,sizeR(1)*sizeR(2));
rawdata = [rawdata avgdata.mean];
rawindices = [rawindices avgdata.comments];

figure(1);subplot(2,1,1);setfiguredefaults()
boxplot(rawdata,rawindices,'notch','off');hold on;
errorbar(length(data)+1,avgdata.mean,1.5*avgdata.iqr,'r');
errorbar(length(data)+1,avgdata.mean,avgdata.iqr,'b');
for j = 1:length(data)
    scatter(j*ones(1,data(j).length),data(j).rawdata,'r');
    xlabels{j} = data(j).comments;
end
grid on;ylabel('Osmolality (mmol/kg)');title('Data:   Range');
xlabels{length(data)+1}=avgdata.comments;

figure(1);subplot(2,1,2);setfiguredefaults()
errorbar(1:(length(data)+1),horzcat([data.means avgdata.mean]),horzcat([data.sem avgdata.sem]),'ko');
set(gca,'Xtick',[1:length(data)+1],'XTickLabel',xlabels)
grid on;ylabel('Mean osmolality (mmol/kg)');title('Data:   Mean ± SEM');

%% Save

clear ans j k ind
save('/Users/joshsalvi/Desktop/BullfrogPlasmaOsm.mat')

    
