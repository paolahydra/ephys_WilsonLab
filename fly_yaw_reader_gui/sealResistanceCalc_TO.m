function Rseal = sealResistanceCalc_TO(I, V, t)
%%% Given a voltage step (mV) and current trace (in pA) recorded with pipette in the bath 
%%% and seal test on (5 mV step), estimates seal resistance (in GOhm)
%%% Tatsuo Okubo
%%% 2017/04/03

%% plot the recorded data
figure(10); clf;
s(1) = subplot(211);
plot(1000*t,V,'color',[230/256 230/256 0/256]) 
ylabel('V (mV)')
box off

s(2) = subplot(212);
plot(1000*t,I,'color',[0 200/256 256/256])
xlabel('Time (ms)')
ylabel('I (pA)')
box off

linkaxes(s,'x')
xlim([0 100])

%% calculate seal resistance (same algorithm as pipette resistance)

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

% remove the transients
H_mid = ceil((Onsets(1:end-1)+Offsets(1:end-1))/2); % middle point of the HIGH pulse [samples]
L_mid = ceil((Offsets(1:end-1)+Onsets(2:end))/2); % middle point of the LOW pulse [samples]

IsHigh = false(size(I));
IsLow = false(size(I));
for k=1:length(Onsets)-1 % for all pulses except for the last
    IsHigh(H_mid(k):Offsets(k)) = true;
    IsLow (L_mid(k):(Onsets(k+1)-1)) = true;
end

I_High = I(IsHigh);
I_Low  = I(IsLow);
I_diff = mean(I_High)-mean(I_Low);
Rseal = 5e-3/(I_diff*1e-12); % [ohm] 5mV seal test

end