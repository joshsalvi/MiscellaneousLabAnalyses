%%
clear;
[sig,Fs] = audioread('/Users/joshsalvi/Downloads/Richetti-John_17_Lewis-Carroll-Jabberwocky_100-Poems_2014.mp3');

N = 10;  % Number of times to divide signal
nn = 5;
N2 = 300;
choice1 = 1; % separate signal(1) or only use a snippet over and over defined by nn(2);

% define mu first
mu = [-1 -1 -1 -0.5 -0.25 -0.125 -0.07 -0.03 -0.015 -0.007 -0.003 -0.002 -0.001 linspace(0.001,0.01,14)];
N=length(mu);N2=length(mu);

warning off;
for j = 1:N
    sigp{j} = sig(1+(j-1)*length(sig)/N:j*length(sig)/N);
end

if choice1==2
    for j = 1:N2
        sigp{j} = sigp{nn};
    end
    N=N2;
end

%%
close all
%mu = linspace(-0.5,0.5,N);
fosc = 0.95;
Fshopf = 1;

%{
fo=790; % what is the frequency of oscillator?
fL = 0.99*fo;
fH = 1.01*fo;
[A,B,C,D] = butter(10,[fL fH]/(0.5*Fs));

d = designfilt('bandstopiir', ...
  'PassbandFrequency1',0.99*fL,'StopbandFrequency1',fL, ...
  'StopbandFrequency2',fH,'PassbandFrequency2',1.01*fH, ...
  'PassbandRipple1',0.5,'StopbandAttenuation',1, ...
  'PassbandRipple2', 0.5, ...
  'DesignMethod','butter','SampleRate', Fs);
%}
for j = 1:N
    [xo{j}, xi{j}] = vdpforced(mu(j),fosc,sigp{j},Fshopf);
    xovar(j) = var(xo{j}(1,:));
    close all
end

xovarmax = max(xovar);
for j = 1:N
    xomax=max(xo{j}(1,:));
    xors{j} = xo{j}.*(1e8);
end
%player = audioplayer(xo, Fs);play(player);
%}
%{
y=1;b=10;
for j = 1:b:N
    figure(1)
    ha=subplot(ceil(sqrt(N/b)),ceil(sqrt(N/b)),y);
    specgram(xo{j},linspace(1,5000,1e3),Fs,length(xo{j})/500); axis off     
    title(['mu = ' num2str(mu(j))]);  
    figure(2)
    subplot(ceil(sqrt(N/b)),ceil(sqrt(N/b)),y);
    [pxx,f]=pwelch(xo{j},[],[],[],Fs);
    plot(f,pxx);axis([min(f) 3000 min(pxx) 1e-3]);   axis off 
    title(['mu = ' num2str(mu(j))]);
    figure(3)
    subplot(ceil(sqrt(N/b)),ceil(sqrt(N/b)),y);
    plot(xo{j});axis off
    title(['mu = ' num2str(mu(j))]);
    y=y+1;
end
%}
%%
for j = 1:N
    player{j} = audioplayer(xo{j}(1,:), Fs);
end

%%
for j = 1:N
    player{j} = audioplayer(xors{j}(1,:).*0.01, Fs);
end


%% 
nplay = 2;
disp(['mu = ' num2str(mu(nplay))]);
play(player{nplay})
%%
stop(player{nplay});

%%
fo=800; % what is the frequency of oscillator?
fL = 0.98*fo;
fH = 1.02*fo;
[A,B,C,D] = butter(10,[fL fH]/(0.5*Fs));

d = designfilt('bandstopiir', ...
  'PassbandFrequency1',0.99*fL,'StopbandFrequency1',fL, ...
  'StopbandFrequency2',fH,'PassbandFrequency2',1.01*fH, ...
  'PassbandRipple1',0.5,'StopbandAttenuation',3, ...
  'PassbandRipple2', 0.5, ...
  'DesignMethod','butter','SampleRate', Fs);

%%
close all
xof = filter(d,xo(1,:));
xofmax = max(abs(xof(length(xof)/2:3*length(xof)/4)));
xof=xof./xofmax;
xof(abs(xof)>1)=[];

figure;
subplot(2,1,1);plot(xo(1,:));subplot(2,1,2);plot(xof(1,:));
figure;
subplot(2,1,1);pwelch(xo(1,:),[],[],[],Fs);subplot(2,1,2);pwelch(xof(1,:),[],[],[],Fs);

player = audioplayer(xof, Fs);play(player);
