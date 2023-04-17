% Copyright 2021-2023 Lawrence Livermore National Security, LLC, and other ExaCA Project Developers.
% See the top-level LICENSE file for details.
%
% SPDX-License-Identifier: MIT
%
% NIST-developed software is provided by NIST as a public service. You may use, copy, and distribute copies of the software in any medium, provided that you keep intact this entire notice. You may improve, modify, and create derivative works of the software or any portion of the software, and you may copy and distribute such modifications or works. Modified works should carry a notice stating that you changed the software and should note the date and nature of any such change. Please explicitly acknowledge the National Institute of Standards and Technology as the source of the software.

% NIST-developed software is expressly provided "AS IS." NIST MAKES NO WARRANTY OF ANY KIND, EXPRESS, IMPLIED, IN FACT, OR ARISING BY OPERATION OF LAW, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NON-INFRINGEMENT, AND DATA ACCURACY. NIST NEITHER REPRESENTS NOR WARRANTS THAT THE OPERATION OF THE SOFTWARE WILL BE UNINTERRUPTED OR ERROR-FREE, OR THAT ANY DEFECTS WILL BE CORRECTED. NIST DOES NOT WARRANT OR MAKE ANY REPRESENTATIONS REGARDING THE USE OF THE SOFTWARE OR THE RESULTS THEREOF, INCLUDING BUT NOT LIMITED TO THE CORRECTNESS, ACCURACY, RELIABILITY, OR USEFULNESS OF THE SOFTWARE.

% You are solely responsible for determining the appropriateness of using and distributing the software and you assume all risks associated with its use, including but not limited to the risks and costs of program errors, compliance with applicable laws, damage to or loss of data, programs or equipment, and the unavailability or interruption of operation. This software is not intended to be used in any situation where a failure could cause risk of injury or damage to property. The software developed by NIST employees is not subject to copyright protection within the United States.
%
%
% Mtex color maps
% Colormap developed by Adam Creuziger (NIST)
% Creates color map with a consistent range for comparison between data sets
% and publications. Developed to distinguish between values less than 1
% (grayscale) and values greater than 1 (color), as a value=1 is equal
% to random (uniform) texture intensity.
% Color range for colormap8.m extends the colormap4.m and are identical in the
% range of 0 to 4.

% Needed to produce test ODF
%clc
%clear
%close all
%% set colomap with grayscale 0 to 1, colors to 1 to 4 MUD    
mtexColorMap([0 0 0; .1 .1 .1 ;.2 .2 .2;.3 .3 .3;     .4 .4 .4 ;     .5 .5 .5;     .6 .6 .6;     .7 .7 .7;     .8 .8 .8;     .9 .9 .9;     0 0 1 ;     0 .1 .9;     0 .2 .8;     0 .3 .7;     0 .4 .6;     0 .5 .5;     0 .6 .4;     0 .7 .3;     0 .8 .2;     0 .9 .1;     0 1 0  ;     .1 1 0;     .2 1 0;     .3 1 0;     .4 1 0;     .5 1 0;     .6 1 0;     .7 1 0;     .8 1 0;     .9 1 0;     1 1 0      ;     1 .9 0;     1 .8 0 ;     1 .7 0;     1 .6 0;     1 .5 0;     1 .4 0;     1 .3 0;     1 .2 0;     1 .1 0;     1 0 0 ])
cmap = colormap;
setMTEXpref('defaultColorMap',cmap);
setMTEXpref('FontSize',20);
%% create an ODF to view or check the colormap
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
