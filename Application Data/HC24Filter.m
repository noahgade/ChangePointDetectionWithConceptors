HC1 = bandpass(HC_LFP, [1, 4], 600);
HC1 = decimate(HC1, 600/4);

PFC1 = bandpass(PFC_LFP, [1, 4], 600);
PFC1 = decimate(PFC1, 600/4);

THAL1 = bandpass(THAL_LFP, [1, 4], 600);
THAL1 = decimate(THAL1, 600/4);

Time = decimate(timestamp, 600/4);

LFP1 = [Time.', HC1.', PFC1.', THAL1.'];
Sleep1 = nremSleep;

writematrix(LFP1, "LFP_Section1_DeltaBand.csv");
writematrix(Sleep1, "SleepTimes_Section1.csv");
