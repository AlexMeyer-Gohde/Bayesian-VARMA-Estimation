
% clear all;
close all; home; format long g;  rng('shuffle');

[y, priorsARMA, proposalsARMA, states, settings, accepted, dropped, modelinfo, parameters] = doBayesian(0,0);