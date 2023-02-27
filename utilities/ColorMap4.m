% Mtex color maps

% Needed to produce test ODF
%clc
%clear
%close all
%% set colomap with grayscale less than 1, colors to 4 MUD
mtexColorMap([0 0 0; .1 .1 .1 ;.2 .2 .2;.3 .3 .3;     .4 .4 .4 ;     .5 .5 .5;     .6 .6 .6;     .7 .7 .7;     .8 .8 .8;     .9 .9 .9;     0 0 1 ;     0 .1 .9;     0 .2 .8;     0 .3 .7;     0 .4 .6;     0 .5 .5;     0 .6 .4;     0 .7 .3;     0 .8 .2;     0 .9 .1;     0 1 0  ;     .1 1 0;     .2 1 0;     .3 1 0;     .4 1 0;     .5 1 0;     .6 1 0;     .7 1 0;     .8 1 0;     .9 1 0;     1 1 0      ;     1 .9 0;     1 .8 0 ;     1 .7 0;     1 .6 0;     1 .5 0;     1 .4 0;     1 .3 0;     1 .2 0;     1 .1 0;     1 0 0 ])
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
export_fig(strcat(savepath,'/','ColorMapExample-phi2-45ODF.tiff'),'-r300') 
%export_fig(strcat('ColorMapExample-phi2-45ODF.tiff'),'-r300') 
%}

%% reset to default Mtex colomaps
%setMTEXpref('defaultColorMap',WhiteJetColorpMap);