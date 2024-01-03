%from diba?
% 
species=['sugar','sugar1','sugar2','taro1','taro2','h2o','taro3','char','gco2','gch4','gc2h4','gco','gcoh2'];

ycoeff=[-1	0	0
        0.47 -1	0
        0.53 0	-1
        0	0.6	0
        0	0.4	0
        0	0.4	0.88
        0	0	0.13
        0	0	1.6
        0	0.9	1.3
        0	0	0.25
        0	0	0.1
        0	0	0.62
        0	0	0.73];

% check if TARO species end up in gas phase
g_index=[4? 5? 6 7];
s_index=[1 2 3 8 9 10 11 12 13];

istart=[1 2 3];

% add molecular weights
MW =[];

afac = [.8e10 .15e5 .2e2];
nfac = [0 0 0];

% convert to J/mol
ea = [26000 16000 20000];
