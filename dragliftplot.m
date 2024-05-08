clear all;
close all;
data=load("force.dat");
time=data(:,1);
fx=data(:,2);
fy=data(:,3);
mom=data(:,4);
tl=length(time);
tls=floor(tl*5/8)+1;
time=time(tls:tl);
fx=fx(tls:tl);
fy=fy(tls:tl);
mom=mom(tls:tl);
t=time;
%% plot x force
plot(time,fx,'r-','LineWidth',0.75)
xlabel("time",'Interpreter','latex',FontSize=14);
ylabel("$C_D$",'Interpreter','latex',FontSize=14);
xlim([min(t) max(t)])
grid on;
grid minor;
ax=gca;
exportgraphics(ax,"CD.png","Resolution",900)

%% plot x force freq
[freq,Mag]=freqget(fx-mean(fx),t);
plot(freq,Mag,'r-','LineWidth',0.75)
xlabel("frequency $(1/T)$",'Interpreter','latex',FontSize=14);
ylabel("$|fft(C_D)|$",'Interpreter','latex',FontSize=14);
xlim([0 1]);
grid on;
grid minor;
ax=gca;
exportgraphics(ax,"CDfreq.png","Resolution",900)

%% plot y force
plot(time,fy,'b-','LineWidth',0.75)
xlabel("time",'Interpreter','latex',FontSize=14);
ylabel("$C_L$",'Interpreter','latex',FontSize=14);
xlim([min(t) max(t)])

grid on;
grid minor;
ax=gca;
exportgraphics(ax,"CL.emf","Resolution",900)

%% plot y force freq
[freq,Mag]=freqget(fy,t);
plot(freq,Mag,'b-','LineWidth',0.75)
xlim([0 1]);
xlabel("frequency $(1/T)$",'Interpreter','latex',FontSize=14);
ylabel("$|fft(C_L)|$",'Interpreter','latex',FontSize=14);
grid on;
grid minor;
ax=gca;
exportgraphics(ax,"CLfreq.emf","Resolution",900)
%% plot mom force
plot(time,mom,'k-','LineWidth',0.75)
xlim([min(t) max(t)])
xlabel("time",'Interpreter','latex',FontSize=14);
ylabel("$C_M$",'Interpreter','latex',FontSize=14);
grid on;
grid minor;
ax=gca;
exportgraphics(ax,"CM.emf","Resolution",900)

%% plot mom force freq
[freq,Mag]=freqget(mom,t);
plot(freq,Mag,'k-','LineWidth',0.75)
xlim([0 1])
xlabel("frequency $(1/T)$",'Interpreter','latex',FontSize=14);
ylabel("$|fft(C_M)|$",'Interpreter','latex',FontSize=14);
grid on;
grid minor;
ax=gca;
exportgraphics(ax,"CMfreq.emf","Resolution",900)
%% print values
cdp=mean(fx);
clp=rms(fy);
cmp=rms(mom);

%% fft function
function [freq,Mag]=freqget(X,t)
if(mod(length(t),2)==1) %%make data even
    t=t(2:end);
    X=X(2:end);
end
T = t(2)-t(1); % sampling period
Fs=1/T; % sampling frequency
t=t-t(1); % readjust time such that t(1)=0
L=length(t);
Y = fft(X);
freq=Fs/L*(0:L-1);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
Mag=P1(1:L/2);
freq=0:(Fs/L):(Fs/2-Fs/L);
end

