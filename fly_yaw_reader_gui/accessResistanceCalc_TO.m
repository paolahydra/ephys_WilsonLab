function [R_access, R_membrane] = accessResistanceCalc_TO(I, V, t, Fs)
%%% perform exponential fit to determine acceess resistance
%%% output access resistane in Mohms
%%% Tatsuo Okubo
%%% 2017/04/13

%%
V2 = V-mean(V); % demean (remove offset due to voltage clamp)
Thres = 0; % [V] threshold on the voltage

Onsets = find([Thres+eps; V2(1:end-1)]<Thres & V2>Thres); % list of pulse onset times [samples]

Offsets = find(V2>Thres & [V2(2:end); Thres+eps]<Thres); % list of pulse offset times [samples]
if Offsets(1)<Onsets(1) % started from offset
    Offsets = Offsets(2:end); % remove the first offset
end

if length(Onsets)>length(Offsets) % extra onset
    Onsets = Onsets(1:end-1); % make the number of onsets and offsets the same
end

H_mid = ceil((Onsets(1:end-1)+Offsets(1:end-1))/2); % middle point of the HIGH pulse [samples]
L_mid = ceil((Offsets(1:end-1)+Onsets(2:end))/2); % middle point of the LOW pulse [samples]

I_ss = nan(length(Onsets)-1,1);
I_peak = nan(length(Onsets)-1,1);
I_baseline = nan(length(Onsets)-1,1);
tau = nan(length(Onsets)-1,1);
R_m = nan(length(Onsets)-1,1);
R_s = nan(length(Onsets)-1,1);
for k=1:length(Onsets)-1 % for all pulses except for the last
    I_peak(k) = max(I(Onsets(k):Offsets(k))); % fast transient peak detection during high
    I_ss(k) = mean(I(Offsets(k)-round(Fs/1000):Offsets(k))); % use last 1 ms of HIGH as steady state
    I_baseline(k)= mean(I(L_mid(k):(Onsets(k+1)-1))); % use second half of LOW as baseline
%     y = I(Onsets(k)+1:Onsets(k)+round(Fs*0.002))-I_ss(k); % only use the first 2 ms to avoid fitting to noise
%     x = (0:length(y)-1)./Fs;
%     x = x'; % column vector
%     f = fit(x,y,'exp1'); % y = a*exp(bx)
%     tau(k) = (-1/f.b)*1e6; % time constant [us]
%     a(k) = f.a;
%     R_m(k) = f.a*5e-3/(I_ss(k)*(f.a+I_ss(k)));
%     R_s(k) = I_ss(k)*R_m(k)/f.a;
%     figure(1); clf;
%     hold on
%     plot(f,x,y)
%     xlabel('Time (s)','fontsize',10)
%     ylabel('I (A)','fontsize',10)
%     set(gca,'color','none')
end

I_diff_1 = mean(I_peak)-mean(I_baseline); % peak current
I_diff_2 = mean(I_ss)-mean(I_baseline); % steady-state current
R_access = (5e-3/(I_diff_1*1e-12))/1e6; % [Mohm] 5mV seal test
R_membrane = (5e-3/(I_diff_2*1e-12))/1e6; % [Mohm]
%Mean_tau = mean(tau); % [us]

% TO DO: plot the fit

%% plot the recorded data
figure(10); clf;
s(1) = subplot(211);
plot(1000*t,V,'color',[230/256 230/256 0/256])
ylabel('V (mV)')
box off

s(2) = subplot(212);
plot(1000*t,I,'color',[0 200/256 256/256])
hold on
line(xlim,[mean(I_peak) mean(I_peak)],'color','r','linestyle',':')
line(xlim,[mean(I_ss) mean(I_ss)],'color','g','linestyle',':')
line(xlim,[0 0],'color','k')
xlabel('Time (ms)')
ylabel('I (pA)')
box off
linkaxes(s,'x')
xlim([0 50])