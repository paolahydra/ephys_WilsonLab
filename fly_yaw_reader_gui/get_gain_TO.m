function [ I_gain V_gain ] = get_gain_TO( X )
%%% Tatsuo Okubo
%%% 200B manual p.63

if X<0.85
    I_gain = 0.5;
    V_gain = nan;
elseif X<1.35
    I_gain = 0.1;
    V_gain = nan;
elseif X<1.85
    I_gain = 0.2;
    V_gain = nan;
elseif X<2.35
    I_gain = 0.5;
    V_gain = 0.5;
elseif X<2.85
    I_gain = 1;
    V_gain = 1;
elseif X<3.35
    I_gain = 2;
    V_gain = 2;
elseif X<3.85
    I_gain = 5;
    V_gain = 5;
elseif X<4.35
    I_gain = 10;
    V_gain = 10;
elseif X<4.85
    I_gain = 20;
    V_gain = 20;
elseif X<5.35
    I_gain = 50;
    V_gain = 50;
elseif X<5.85
    I_gain = 100;
    V_gain = 100;
elseif X<6.35
    I_gain = 200;
    V_gain = 200;
elseif X<6.85
    I_gain = 500;
    V_gain = 500;
end