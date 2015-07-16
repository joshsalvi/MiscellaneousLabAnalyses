nm1=5;
nm2=5;
Fc=[linspace(-1,1,nm1) linspace(-1,1,nm2) linspace(1,3.4,nm1+nm2)];
kc=[2.6 5 0];
noiselevel=linspace(0,0.7,2);
t=linspace(0,1000,2e4);

for j = 1:length(Fc)
for k = 1:length(kc)
for l = 1:length(noiselevel)
    if j <= nm1+nm2
if l==1
[Xdet{j,k},Xsto{j,k,l}]=hbtoymodel(Fc(j),kc(k),noiselevel(l),0,0,t);
else
[~,Xsto{j,k,l}]=hbtoymodel(Fc(j),kc(k),noiselevel(l),0,0,t);
end
    else
        if l==1
[Xdet{j,k},Xsto{j,k,l}]=hbtoymodel(kc(k),Fc(j),noiselevel(l),0,0,t);
else
[~,Xsto{j,k,l}]=hbtoymodel(kc(k),Fc(j),noiselevel(l),0,0,t);
end
    end
end

end
end
disp('simulation complete')

for j = 1:length(Fc)
for k = 1:length(kc)
for l = 1:length(noiselevel)
Xsto{j,k,l}=Xsto{j,k,l}(1,round(length(Xsto{j,k,l})/2):end)./max(Xsto{j,k,l}(1,round(length(Xsto{j,k,l})/2):end));
end
end
end
for j = 1:length(Fc)
for k = 1:length(kc)
Xdet{j,k}=Xdet{j,k}(1,round(length(Xdet{j,k})/2):end);
if var(Xdet{j,k})>0
Xdet{j,k}=Xdet{j,k}./max(Xdet{j,k});
end
end
end
for j = 1:length(Fc)
for k = 1:length(kc)
Xdet{j,k}=Xdet{j,k}-mean(Xdet{j,k});
end
end
for j = 1:length(Fc)
for k = 1:length(kc)
for l = 1:length(noiselevel)
Xsto{j,k,l}=Xsto{j,k,l}-mean(Xsto{j,k,l});
end
end
end
disp('detrended data');

for j = 1:length(Fc)
for k = 1:length(kc)
for l = 1:length(noiselevel)
        if j <= nm1+nm2
if k == 1
index{j,k,l}='osc';
else
index{j,k,l}='mono';
end
        else
            if k == 3
            index{j,k,l}='osc';
            else
                index{j,k,l}='mono';
            end
        end
end
end
end
disp('classifiers named')

for j = 1:length(Fc)
for k = 1:length(kc)
for l = 1:length(noiselevel)
X=Xsto{j,k,l};
Xctr = mean(X);
Xupper = X(X>=Xctr) - Xctr;
Xlower = Xctr - X(X<=Xctr);
Xupper = [Xupper' -Xupper']';Xupper=reshape(Xupper,1,size(Xupper,1)*size(Xupper,2));
Xlower = [Xlower' -Xlower']';Xlower=reshape(Xlower,1,size(Xlower,1)*size(Xlower,2));
[khsymm(j,k,l),kpsymm(j,k,l),kKSstat(j,k,l)] = kstest2(Xupper,Xlower,10^-2);
[~,kpks(j,k,l),kpksstat(j,k,l)]=kstest(X);
NK=length(X);
SES=sqrt(6*NK*(NK-1)/((NK-2)*(NK+1)*(NK+3)));
SEK=2*SES*sqrt((NK^2-1)/((NK-3)*(NK+5)));
kpkurt(j,k,l)=(kurtosis(X,0)-3)/SEK;
kpK(j,k,l) = cdf('Normal',kpkurt(j,k,l),0,1);
[kpdip(j,k,l), kpdipp(j,k,l)]=HartigansDipSignifTest(X,10);
end
end
end
disp('completed statistics');

%%

sizeX=size(Xsto);
kpkurt=reshape(kpkurt,1,sizeX(1)*sizeX(2)*sizeX(3));
kKSstat=reshape(kKSstat,1,sizeX(1)*sizeX(2)*sizeX(3));
kpdipp=reshape(kpdipp,1,sizeX(1)*sizeX(2)*sizeX(3));
kpK=reshape(kpK,1,sizeX(1)*sizeX(2)*sizeX(3));
kpks=reshape(kpks,1,sizeX(1)*sizeX(2)*sizeX(3));
kpksstat=reshape(kpksstat,1,sizeX(1)*sizeX(2)*sizeX(3));
kpdip=reshape(kpdip,1,sizeX(1)*sizeX(2)*sizeX(3));
kpsymm=reshape(kpsymm,1,sizeX(1)*sizeX(2)*sizeX(3));
index=reshape(index,1,sizeX(1)*sizeX(2)*sizeX(3));
Xsto=reshape(Xsto,1,sizeX(1)*sizeX(2)*sizeX(3));

%%
kpvalues = [kpkurt;kKSstat;kpdipp;kpK;kpks;kpksstat;kpdip;kpsymm]';

svmstruct1=svmtrain(kpvalues,index,'ShowPlot',false,'kernel_function','linear');
svmstruct2=svmtrain(kpvalues,index,'ShowPlot',false,'kernel_function','quadratic');
svmstruct3=svmtrain(kpvalues,index,'ShowPlot',false,'kernel_function','polynomial','polyorder',2);
svmstruct4=svmtrain(kpvalues,index,'ShowPlot',false,'kernel_function','polynomial','polyorder',3);
svmstruct5=svmtrain(kpvalues,index,'ShowPlot',false,'kernel_function','mlp');
svmstruct6=svmtrain(kpvalues,index,'ShowPlot',false,'kernel_function','rbf');
%%
close all

n=550;

figure;
plot(Xsto{n});title(svmclassify(svmstruct1,[kpkurt(n) kKSstat(n) kpdipp(n) kpK(n) kpks(n) kpksstat(n) kpdip(n) kpsymm(n)]));
%figure;
%svmclassify(svmstruct1,[kpkurt(n) kKSstat(n) kpdipp(n) kpK(n) kpks(n) kpksstat(n) kpdip(n) kpsymm(n)],'ShowPlot',true);

%%
close all

Fe=2;ke=2;noiselevel=1;t=linspace(0,1000,1e4);


[~,xsto1]=hbtoymodel(Fe,ke,noiselevel,0,0,t);
xsto1=xsto1(1,round(length(xsto1)/2):end)./max(xsto1(1,round(length(xsto1)/2):end));
xsto1=xsto1-mean(xsto1);

X=xsto1;
Xctr = mean(X);
Xupper = X(X>=Xctr) - Xctr;
Xlower = Xctr - X(X<=Xctr);
Xupper = [Xupper' -Xupper']';Xupper=reshape(Xupper,1,size(Xupper,1)*size(Xupper,2));
Xlower = [Xlower' -Xlower']';Xlower=reshape(Xlower,1,size(Xlower,1)*size(Xlower,2));
[khsymm1,kpsymm1,kKSstat1] = kstest2(Xupper,Xlower,10^-2);
[~,kpks1,kpksstat1]=kstest(X);
NK=length(X);
SES=sqrt(6*NK*(NK-1)/((NK-2)*(NK+1)*(NK+3)));
SEK=2*SES*sqrt((NK^2-1)/((NK-3)*(NK+5)));
kpkurt1=(kurtosis(X,0)-3)/SEK;
kpK1 = cdf('Normal',kpkurt1,0,1);
[kpdip1, kpdipp1]=HartigansDipSignifTest(X,200);

figure;
plot(xsto1);title(svmclassify(svmstruct1,[kpkurt1 kKSstat1 kpdipp1 kpK1 kpks1 kpksstat1 kpdip1 kpsymm1]));
%figure;
%svmclassify(svmstruct1,[kpkurt(n) kKSstat(n) kpdipp(n) kpK(n) kpks(n) kpksstat(n) kpdip(n) kpsymm(n)],'ShowPlot',true);

%%
close all;
Xdet2=Xdet;Xsto2=Xsto;


Fc1=linspace(-3,3,21);
kc1=linspace(1,5,21);
noiselevel=[0.001 0.01 0.1];
t=linspace(0,1000,2e4);

for j = 1:length(Fc1)
for k = 1:length(kc1)
for l = 1:length(noiselevel)

[Xdet{j,k},Xsto{j,k,l}]=hbtoymodel(Fc1(j),kc1(k),noiselevel(l),0,0,t);

[~,Xsto{j,k,l}]=hbtoymodel(Fc1(j),kc1(k),noiselevel(l),0,0,t);

    end
end

end
disp('finished simulation.');

for j = 1:length(Fc1)
for k = 1:length(kc1)
for l = 1:length(noiselevel)
Xsto{j,k,l}=Xsto{j,k,l}(1,round(length(Xsto{j,k,l})/2):end)./max(Xsto{j,k,l}(1,round(length(Xsto{j,k,l})/2):end));
end
end
end
for j = 1:length(Fc1)
for k = 1:length(kc1)
Xdet{j,k}=Xdet{j,k}(1,round(length(Xdet{j,k})/2):end);
if var(Xdet{j,k})>0
Xdet{j,k}=Xdet{j,k}./max(Xdet{j,k});
end
end
end
for j = 1:length(Fc1)
for k = 1:length(kc1)
Xdet{j,k}=Xdet{j,k}-mean(Xdet{j,k});
end
end
for j = 1:length(Fc1)
for k = 1:length(kc1)
for l = 1:length(noiselevel)
Xsto{j,k,l}=Xsto{j,k,l}-mean(Xsto{j,k,l});
end
end
end

for j = 1:length(Fc1)
for k = 1:length(kc1)
for l = 1:length(noiselevel)
X=Xsto{j,k,l};
Xctr = mean(X);
Xupper = X(X>=Xctr) - Xctr;
Xlower = Xctr - X(X<=Xctr);
Xupper = [Xupper' -Xupper']';Xupper=reshape(Xupper,1,size(Xupper,1)*size(Xupper,2));
Xlower = [Xlower' -Xlower']';Xlower=reshape(Xlower,1,size(Xlower,1)*size(Xlower,2));
[khsymm2(j,k,l),kpsymm2(j,k,l),kKSstat2(j,k,l)] = kstest2(Xupper,Xlower,10^-2);
[~,kpks2(j,k,l),kpksstat2(j,k,l)]=kstest(X);
NK=length(X);
SES=sqrt(6*NK*(NK-1)/((NK-2)*(NK+1)*(NK+3)));
SEK=2*SES*sqrt((NK^2-1)/((NK-3)*(NK+5)));
kpkurt2(j,k,l)=(kurtosis(X,0)-3)/SEK;
kpK2(j,k,l) = cdf('Normal',kpkurt2(j,k,l),0,1);
[kpdip2(j,k,l), kpdipp2(j,k,l)]=HartigansDipSignifTest(X,10);
end
end
disp(num2str(j))
end



for j = 1:length(Fc1)
for k = 1:length(kc1)
for l = 1:length(noiselevel)
    indclass{j,k,l}=svmclassify(svmstruct1,[kpkurt2(j,k,l) kKSstat2(j,k,l) kpdipp2(j,k,l) kpK2(j,k,l) kpks2(j,k,l) kpksstat2(j,k,l) kpdip2(j,k,l) kpsymm2(j,k,l)]);
    if isempty(indclass{j,k,l}{1})==0 && indclass{j,k,l}{1}(1) == 'o'
        indclasssig(j,k,l)=1;
    else
        indclasssig(j,k,l)=0;
    end
end
end
end

figure(1);
subplot(1,3,1);pcolor(kc1,Fc1,indclasssig(:,:,1));title(num2str(noiselevel(1)));
subplot(1,3,2);pcolor(kc1,Fc1,indclasssig(:,:,2));title(num2str(noiselevel(2)));
subplot(1,3,3);pcolor(kc1,Fc1,indclasssig(:,:,3));title(num2str(noiselevel(3)));

for j = 1:length(Fc1)
for k = 1:length(kc1)
for l = 1:length(noiselevel)
    indclass{j,k,l}=svmclassify(svmstruct2,[kpkurt2(j,k,l) kKSstat2(j,k,l) kpdipp2(j,k,l) kpK2(j,k,l) kpks2(j,k,l) kpksstat2(j,k,l) kpdip2(j,k,l) kpsymm2(j,k,l)]);
    if isempty(indclass{j,k,l}{1})==0 && indclass{j,k,l}{1}(1) == 'o'
        indclasssig(j,k,l)=1;
    else
        indclasssig(j,k,l)=0;
    end
end
end
end

figure(2);
subplot(1,3,1);pcolor(kc1,Fc1,indclasssig(:,:,1));title(num2str(noiselevel(1)));
subplot(1,3,2);pcolor(kc1,Fc1,indclasssig(:,:,2));title(num2str(noiselevel(2)));
subplot(1,3,3);pcolor(kc1,Fc1,indclasssig(:,:,3));title(num2str(noiselevel(3)));

for j = 1:length(Fc1)
for k = 1:length(kc1)
for l = 1:length(noiselevel)
    indclass{j,k,l}=svmclassify(svmstruct3,[kpkurt2(j,k,l) kKSstat2(j,k,l) kpdipp2(j,k,l) kpK2(j,k,l) kpks2(j,k,l) kpksstat2(j,k,l) kpdip2(j,k,l) kpsymm2(j,k,l)]);
    if isempty(indclass{j,k,l}{1})==0 && indclass{j,k,l}{1}(1) == 'o'
        indclasssig(j,k,l)=1;
    else
        indclasssig(j,k,l)=0;
    end
end
end
end

figure(3);
subplot(1,3,1);pcolor(kc1,Fc1,indclasssig(:,:,1));title(num2str(noiselevel(1)));
subplot(1,3,2);pcolor(kc1,Fc1,indclasssig(:,:,2));title(num2str(noiselevel(2)));
subplot(1,3,3);pcolor(kc1,Fc1,indclasssig(:,:,3));title(num2str(noiselevel(3)));


for j = 1:length(Fc1)
for k = 1:length(kc1)
for l = 1:length(noiselevel)
    indclass{j,k,l}=svmclassify(svmstruct4,[kpkurt2(j,k,l) kKSstat2(j,k,l) kpdipp2(j,k,l) kpK2(j,k,l) kpks2(j,k,l) kpksstat2(j,k,l) kpdip2(j,k,l) kpsymm2(j,k,l)]);
    if isempty(indclass{j,k,l}{1})==0 && indclass{j,k,l}{1}(1) == 'o'
        indclasssig(j,k,l)=1;
    else
        indclasssig(j,k,l)=0;
    end
end
end
end

figure(4);
subplot(1,3,1);pcolor(kc1,Fc1,indclasssig(:,:,1));title(num2str(noiselevel(1)));
subplot(1,3,2);pcolor(kc1,Fc1,indclasssig(:,:,2));title(num2str(noiselevel(2)));
subplot(1,3,3);pcolor(kc1,Fc1,indclasssig(:,:,3));title(num2str(noiselevel(3)));

for j = 1:length(Fc1)
for k = 1:length(kc1)
for l = 1:length(noiselevel)
    indclass{j,k,l}=svmclassify(svmstruct5,[kpkurt2(j,k,l) kKSstat2(j,k,l) kpdipp2(j,k,l) kpK2(j,k,l) kpks2(j,k,l) kpksstat2(j,k,l) kpdip2(j,k,l) kpsymm2(j,k,l)]);
    if isempty(indclass{j,k,l}{1})==0 && indclass{j,k,l}{1}(1) == 'o'
        indclasssig(j,k,l)=1;
    else
        indclasssig(j,k,l)=0;
    end
end
end
end

figure(5);
subplot(1,3,1);pcolor(kc1,Fc1,indclasssig(:,:,1));title(num2str(noiselevel(1)));
subplot(1,3,2);pcolor(kc1,Fc1,indclasssig(:,:,2));title(num2str(noiselevel(2)));
subplot(1,3,3);pcolor(kc1,Fc1,indclasssig(:,:,3));title(num2str(noiselevel(3)));



for j = 1:length(Fc1)
for k = 1:length(kc1)
for l = 1:length(noiselevel)
    indclass{j,k,l}=svmclassify(svmstruct6,[kpkurt2(j,k,l) kKSstat2(j,k,l) kpdipp2(j,k,l) kpK2(j,k,l) kpks2(j,k,l) kpksstat2(j,k,l) kpdip2(j,k,l) kpsymm2(j,k,l)]);
    if isempty(indclass{j,k,l}{1})==0 && indclass{j,k,l}{1}(1) == 'o'
        indclasssig(j,k,l)=1;
    else
        indclasssig(j,k,l)=0;
    end
end
end
end

figure(6);
subplot(1,3,1);pcolor(kc1,Fc1,indclasssig(:,:,1));title(num2str(noiselevel(1)));
subplot(1,3,2);pcolor(kc1,Fc1,indclasssig(:,:,2));title(num2str(noiselevel(2)));
subplot(1,3,3);pcolor(kc1,Fc1,indclasssig(:,:,3));title(num2str(noiselevel(3)));

%%
close all
N=5;
figure;


for j = 1:length(Fc1)
for k = 1:length(kc1)
for l = 1:length(noiselevel)
    indclass{j,k,l}=svmclassify(svmstruct3,[kpkurt2(j,k,l) kKSstat2(j,k,l) kpdipp2(j,k,l) kpK2(j,k,l) kpks2(j,k,l) kpksstat2(j,k,l) kpdip2(j,k,l) kpsymm2(j,k,l)]);
    if isempty(indclass{j,k,l}{1})==0 && indclass{j,k,l}{1}(1) == 'o'
        indclasssig(j,k,l)=1;
    else
        indclasssig(j,k,l)=0;
    end
end
end
end

for n = 1:N
    
    for j = 1:length(Fc1)
for k = 1:length(kc1)
for l = 1:length(noiselevel)
indclass2{j,k,l}=indclass{j,k,l}{1};
end
end
end
sizeX=size(kpkurt2);
kpkurt3=reshape(kpkurt2,1,sizeX(1)*sizeX(2)*sizeX(3));
kKSstat3=reshape(kKSstat2,1,sizeX(1)*sizeX(2)*sizeX(3));
kpdipp3=reshape(kpdipp2,1,sizeX(1)*sizeX(2)*sizeX(3));
kpK3=reshape(kpK2,1,sizeX(1)*sizeX(2)*sizeX(3));
kpks3=reshape(kpks2,1,sizeX(1)*sizeX(2)*sizeX(3));
kpksstat3=reshape(kpksstat2,1,sizeX(1)*sizeX(2)*sizeX(3));
kpdip3=reshape(kpdip2,1,sizeX(1)*sizeX(2)*sizeX(3));
kpsymm3=reshape(kpsymm2,1,sizeX(1)*sizeX(2)*sizeX(3));
indclass3=reshape(indclass2,1,sizeX(1)*sizeX(2)*sizeX(3));
kpvalues3 = [kpkurt3;kKSstat3;kpdipp3;kpK3;kpks3;kpksstat3;kpdip3;kpsymm3]';

for j = 1:length(indclass3)
    if isempty(indclass3{j})==1
        indclass3{j}='mono';
    end
end
opts=statset('MaxIter',1e6);
if mod(n,2)==1
svmstructn=svmtrain(kpvalues3,indclass3','ShowPlot',false,'kernel_function','polynomial','polyorder',3,'method','SMO','options',opts);
else
    svmstructn=svmtrain(kpvalues3,indclass3','ShowPlot',false,'kernel_function','polynomial','polyorder',3,'method','LS','options',opts);
end


clear indclass indclass2 indclass3 indclasssig

for j = 1:length(Fc1)
for k = 1:length(kc1)
for l = 1:length(noiselevel)
    indclass{j,k,l}=svmclassify(svmstructn,[kpkurt2(j,k,l) kKSstat2(j,k,l) kpdipp2(j,k,l) kpK2(j,k,l) kpks2(j,k,l) kpksstat2(j,k,l) kpdip2(j,k,l) kpsymm2(j,k,l)]);
    if isempty(indclass{j,k,l}{1})==0 && indclass{j,k,l}{1}(1) == 'o'
        indclasssig(j,k,l)=1;
    else
        indclasssig(j,k,l)=0;
    end
end
end
end
figure(1);subplot(5,ceil(N/5),n);pcolor(kc1,Fc1,indclasssig(:,:,1));title(num2str(noiselevel(1)));
figure(2);subplot(5,ceil(N/5),n);pcolor(kc1,Fc1,indclasssig(:,:,2));title(num2str(noiselevel(2)));
figure(3);subplot(5,ceil(N/5),n);pcolor(kc1,Fc1,indclasssig(:,:,3));title(num2str(noiselevel(3)));
end
