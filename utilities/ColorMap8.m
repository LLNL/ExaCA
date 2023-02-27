% Mtex color maps

% Needed to produce test ODF
%clc
%clear
%close all
%% set colomap with grayscale less than 1, colors to 4 MUD
mtexColorMap([0 0 0; .1 .1 .1 ;.2 .2 .2;.3 .3 .3;     .4 .4 .4 ;     .5 .5 .5;     .6 .6 .6;     .7 .7 .7;     .8 .8 .8;     .9 .9 .9;     0 0 1 ;...
    0 .1 .9;     0 .2 .8;     0 .3 .7;     0 .4 .6;     0 .5 .5;     0 .6 .4;     0 .7 .3;     0 .8 .2;     0 .9 .1;     0 1 0  ;...
    .1 1 0;     .2 1 0;     .3 1 0;     .4 1 0;     .5 1 0;     .6 1 0;     .7 1 0;     .8 1 0;     .9 1 0;     1 1 0      ;...
    1 .9 0;     1 .8 0 ;     1 .7 0;     1 .6 0;     1 .5 0;     1 .4 0;     1 .3 0;     1 .2 0;     1 .1 0;     1 0 0;...
    1 0 .025;   1 0 .05;     1 0 .075;   1 0 .1;     1 0 .125;   1 0 .15;   1 0 .175;   1 0 .2;    1 0 .225;    1 0 .25;...
    1 0 .275;   1 0 .30;     1 0 .325;   1 0 .35;    1 0 .375;   1 0 .4;    1 0 .425;   1 0 .45;   1 0 .475;    1 0 .5;...
    1 0 .525;   1 0 .55;     1 0 .575;   1 0 .60;    1 0 .625;   1 0 .65;   1 0 .675;   1 0 .70;   1 0 .725;    1 0 .75;...
    1 0 .775;   1 0 .80;     1 0 .825;   1 0 .85;    1 0 .875;   1 0 .90;   1 0 .925;   1 0 .95;   1 0 .975;    1 0 1 ])
cmap = colormap;
setMTEXpref('defaultColorMap',cmap);
setMTEXpref('FontSize',20);
%% test ODF
%{
savepath='.';

HW=20*degree;

cs = crystalSymmetry('m-3m');
ss = specimenSymmetry('orthorhombic');
Cube = orientation('Miller',[0 0 1],[1 0 0],cs,ss);
odf = unimodalODF(Cube,'halfwidth',HW,cs,ss);
figure; plot(odf,'phi2',[45]*degree,'projection','plain','minmax', 'off',cs,ss);CLim(gcm,[0, 4]);mtexColorbar;
export_fig(strcat(savepath,'/','ColorMapExample-phi2-45ODF.tiff'),'-r150') 

%}
%% reset to default Mtex colomaps
%setMTEXpref('defaultColorMap',WhiteJetColorpMap);