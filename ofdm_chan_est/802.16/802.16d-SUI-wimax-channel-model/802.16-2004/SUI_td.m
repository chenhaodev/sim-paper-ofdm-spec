%generate channel model per Stanford University Interim (SUI) model
%per IEEE802.16.3c-01/29r1
%
%use:
%[imp taps] = function SUI_model(SUI_indx, direct, fs );
% where:
%SUI_indx: index of SUI model allowed range 1...6
%direct: usage of directional antenna at CPE. 0=> omni 1=>directional
%fs: sampling rate in MHz
%imp: generated impulse response
%taps: non-zero taps of impulse response
%Tal Kaitz 6/3/01
function [out,imp] =SUI_model(dat,SUI_indx, direct, fs );

% check for input parameter range
if ~ismember(direct,[0,1])
error('direct must be in the range [0,1]' );
end;
if ~ismember(SUI_indx,[1:6])
error('direct must be in the range 1:6' );
end;
%model information matrix
SUI1_omni=[... %SU1 omni model
0 0.4 0.8; ... %tap delay in uSe
0 -15 -20;... %relative power
4 0 0 ]; %K factor per tap
SUI1_dir=[... %SU1 deirectional 30deg model
0 0.4 0.8; ...
0 -21 -32;...
16 0 0 ];
SUI2_omni=[... %SU2 omni model
0 0.5 1.0; ...
0 -12 -15;...
2 0 0 ];
SUI2_dir=[... %SU2 omni model
0 0.5 1.0; ...
0 -18 -27;...
8 0 0 ];
SUI3_omni=[... %SU3 omni model
0 0.5 1.0; ...
0 -5 -10;...
1 0 0 ];
SUI3_dir=[... %SU3 dir model
0 0.5 1.0; ...
0 -11 -22;...
3 0 0 ];
SUI4_omni=[... %SU4 omni model
0 2 4; ...
0 -4 -8;...
0 0 0 ];
SUI4_dir=[... %SU4 dir model
0 2 4; ...
0 -10 -20;...
0 0 0 ];
SUI5_omni=[... %SU5 omni model
0 4 10; ...
0 -5 -10;...
0 0 0 ];
SUI5_dir=[... %SU5 dir model
0 4 10; ...
0 -11 -22;...
0 0 0 ];
SUI6_omni=[... %SU6 omni model
0 14 20; ...
0 -10 -20;...
0 0 0 ];
SUI6_dir=[... %SU6 dir model
0 14 20; ...
0 -16 -26;...
0 0 0 ];
%overall SUI model
SUI={SUI1_omni SUI2_omni SUI3_omni SUI4_omni SUI5_omni SUI6_omni ;...
SUI1_dir SUI2_dir SUI3_dir SUI4_dir SUI5_dir SUI6_dir};
%select a specific model
SUI_mod=SUI{direct+1,SUI_indx};

%number of taps
n=size(SUI_mod,2);
%location of taps
locs=1+round(SUI_mod(1,:)*fs);
imp=zeros(1,max(locs));
%get taps
taps=gen_rice(10.^(SUI_mod(2,:)/10)... %power of taps
,SUI_mod(3,:),... %K factor of taps
1); %normalize the result
imp(locs)=taps;
out = filter(imp,1,dat);
