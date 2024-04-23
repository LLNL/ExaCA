% Copyright Lawrence Livermore National Security, LLC, and other ExaCA Project Developers.
% See the top-level LICENSE file for details.
%
% SPDX-License-Identifier: MIT
%
% NIST-developed software is provided by NIST as a public service. You may use, copy, and distribute copies of the software in any medium, provided that you keep intact this entire notice. You may improve, modify, and create derivative works of the software or any portion of the software, and you may copy and distribute such modifications or works. Modified works should carry a notice stating that you changed the software and should note the date and nature of any such change. Please explicitly acknowledge the National Institute of Standards and Technology as the source of the software.

% NIST-developed software is expressly provided "AS IS." NIST MAKES NO WARRANTY OF ANY KIND, EXPRESS, IMPLIED, IN FACT, OR ARISING BY OPERATION OF LAW, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NON-INFRINGEMENT, AND DATA ACCURACY. NIST NEITHER REPRESENTS NOR WARRANTS THAT THE OPERATION OF THE SOFTWARE WILL BE UNINTERRUPTED OR ERROR-FREE, OR THAT ANY DEFECTS WILL BE CORRECTED. NIST DOES NOT WARRANT OR MAKE ANY REPRESENTATIONS REGARDING THE USE OF THE SOFTWARE OR THE RESULTS THEREOF, INCLUDING BUT NOT LIMITED TO THE CORRECTNESS, ACCURACY, RELIABILITY, OR USEFULNESS OF THE SOFTWARE.

% You are solely responsible for determining the appropriateness of using and distributing the software and you assume all risks associated with its use, including but not limited to the risks and costs of program errors, compliance with applicable laws, damage to or loss of data, programs or equipment, and the unavailability or interruption of operation. This software is not intended to be used in any situation where a failure could cause risk of injury or damage to property. The software developed by NIST employees is not subject to copyright protection within the United States.

function [] = PlotIPFColoredSection(MTEXFile)

    % Script created by Matt Rolchigo (ORNL) and Adam Creuziger (NIST)
    % Plots the cross-section specified by `MTEXFile` using the inverse pole figure colormap relative to the X, Y, and Z directions
    % crystal summetry
    CS = crystalSymmetry('m3m', [1 1 1], [90,90,90]*degree, 'X||a', 'Y||b*', 'Z||c*');

    % plotting convention
    setMTEXpref('xAxisDirection','east');
    setMTEXpref('yAxisDirection','north');

    % create an EBSD variable containing the data
    ebsd = EBSD.load(MTEXFile,CS,'interface','generic',...
      'ColumnNames', { 'phi1' 'Phi' 'phi2' 'Phase' 'x' 'y'}, 'Bunge');

    % Draw grain boundaries
    [grains,ebsd.grainId] = calcGrains(ebsd('indexed'),'threshold',0.001);

    % remove very small grains
    ebsd(grains(grains.grainSize<=6.25)) = [];
       
    % smooth the grains a bit
    grains = smooth(grains,5);
    ebsd = smooth(ebsd, halfQuadraticFilter, 'extrapolate');
    gB = grains.boundary;
  
    % and recompute grains
    [grains,ebsd.grainId] = calcGrains(ebsd('indexed'),'threshold',0.001);

    ipfKey = ipfColorKey(ebsd('1').orientations);

    Ext = strfind(MTEXFile,'.');
    SizeName = size(Ext);
    BaseFileName = extractBefore(MTEXFile,Ext(SizeName(2)));

    %plotting
    figure(1);
    ipfKey.inversePoleFigureDirection = vector3d.X;
    colors = ipfKey.orientation2color(grains('1').meanOrientation);
    plot(grains('1'),colors)
    hold on
    plot(gB,'lineWidth',1);
    hold off
    OutputFileName = strcat(BaseFileName,'_IPF-X.png');
    export_fig(OutputFileName,'-r150');

    figure(2);
    ipfKey.inversePoleFigureDirection = vector3d.Y;
    colors = ipfKey.orientation2color(grains('1').meanOrientation);
    plot(grains('1'),colors)
    hold on
    plot(gB,'lineWidth',1);
    hold off
    OutputFileName = strcat(BaseFileName,'_IPF-Y.png');
    export_fig(OutputFileName,'-r150');
    
    figure(3);
    ipfKey.inversePoleFigureDirection = vector3d.Z;
    colors = ipfKey.orientation2color(grains('1').meanOrientation);
    plot(grains('1'),colors)
    hold on
    plot(gB,'lineWidth',1);
    hold off
    OutputFileName = strcat(BaseFileName,'_IPF-Z.png');
    export_fig(OutputFileName,'-r150');

end
