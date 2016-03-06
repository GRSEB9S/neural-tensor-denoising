%% sfn_data_script
clear
RCpath = 'Datasets/Motor/RC,2009-09-18,1-2,good-ss';
Spath = 'Datasets/Motor/Lara-20141105/Sstructs/';
%A = reshape(A,[n t 8 3]); %% for later
%% pre-process structs
bw = 10;
sd = 10;

[TL,AL,SL,T0L] = preS(Spath, bw, sd);
[TC,AC,T0C] = preRC(RCpath, bw, sd);

D1.L.T = TL;
D1.L.A = AL;
D1.L.S = SL;
D1.L.T0 = T0L;

D1.C.T = TC;
D1.C.A = AC;
D1.C.T0 = T0C;

D1.bw = bw;
D1.sd = sd;

%%
bw = 25;
sd = 25;

[TL,AL,SL,T0L] = preS(Spath, bw, sd);
[TC,AC,T0C] = preRC(RCpath, bw, sd);

D2.L.T = TL;
D2.L.A = AL;
D2.L.S = SL;
D2.L.T0 = T0L;

D2.C.T = TC;
D2.C.A = AC;
D2.C.T0 = T0C;

D2.bw = bw;
D2.sd = sd;

%%
D = D1; save('sfn_10bw_10sd','D','-v7.3')
%%
D = D2; save('sfn_25bw_10sd','D','-v7.3');