% Copyright 2021-2024 Lawrence Livermore National Security, LLC, and other ExaCA Project Developers.
% See the top-level LICENSE file for details.
%
% SPDX-License-Identifier: MIT
%
% NIST-developed software is provided by NIST as a public service. You may use, copy, and distribute copies of the software in any medium, provided that you keep intact this entire notice. You may improve, modify, and create derivative works of the software or any portion of the software, and you may copy and distribute such modifications or works. Modified works should carry a notice stating that you changed the software and should note the date and nature of any such change. Please explicitly acknowledge the National Institute of Standards and Technology as the source of the software.

% NIST-developed software is expressly provided "AS IS." NIST MAKES NO WARRANTY OF ANY KIND, EXPRESS, IMPLIED, IN FACT, OR ARISING BY OPERATION OF LAW, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NON-INFRINGEMENT, AND DATA ACCURACY. NIST NEITHER REPRESENTS NOR WARRANTS THAT THE OPERATION OF THE SOFTWARE WILL BE UNINTERRUPTED OR ERROR-FREE, OR THAT ANY DEFECTS WILL BE CORRECTED. NIST DOES NOT WARRANT OR MAKE ANY REPRESENTATIONS REGARDING THE USE OF THE SOFTWARE OR THE RESULTS THEREOF, INCLUDING BUT NOT LIMITED TO THE CORRECTNESS, ACCURACY, RELIABILITY, OR USEFULNESS OF THE SOFTWARE.

% You are solely responsible for determining the appropriateness of using and distributing the software and you assume all risks associated with its use, including but not limited to the risks and costs of program errors, compliance with applicable laws, damage to or loss of data, programs or equipment, and the unavailability or interruption of operation. This software is not intended to be used in any situation where a failure could cause risk of injury or damage to property. The software developed by NIST employees is not subject to copyright protection within the United States.
%

function [] = PlotPoleFigure(MTEXFile, ColormapFile, ColormapUpperLimit)

    % Script created by Matt Rolchigo (ORNL) and Adam Creuziger (NIST)
    % Plots pole figure data, and inverse pole figure data for the file specified by `MTEXFile`, using the colormap specified by `ColormapFile`, with colormap bounds [0, `ColormapUpperLimit`]
    % crystal symmetry
    CS = crystalSymmetry('m3m', [1 1 1], [90,90,90]*degree, 'X||a', 'Y||b*', 'Z||c*');
    SS = specimenSymmetry('mmm');
    % plotting convention
    setMTEXpref('xAxisDirection','east');
    setMTEXpref('yAxisDirection','north');

    % set colorbar
    run(ColormapFile)

    % Create ODF
    odf_VPSClike = loadODF_generic(MTEXFile, 'Bunge', 'degree', 'cs', CS, ...
        'header', 5, 'ColumnNames', {'phi1' 'Phi' 'phi2' 'weights'}, 'density');

    % Create an Pole figure variable containing the data and plot
    h = [Miller(1,0,0,CS),Miller(1,1,0,CS),Miller(1,1,1,CS)];
    r = regularS2Grid('resolution',5*degree);
    pfs_VPSClike = calcPoleFigure(odf_VPSClike,h,r);

    figure(1); plot(pfs_VPSClike,'contourf');CLim(gcm,[0, ColormapUpperLimit]);mtexColorbar;
    Ext = strfind(MTEXFile,'.');
    SizeName = size(Ext);
    BaseFileName = extractBefore(MTEXFile,Ext(SizeName(2)));
    OutputFileName = strcat(BaseFileName,'_PoleFigure.png');
    export_fig(OutputFileName,'-r150');
  
    figure(2); plotIPDF(odf_VPSClike,[xvector,yvector,zvector],'contourf');CLim(gcm,[0, ColormapUpperLimit]);mtexColorbar;
    OutputFileName2 = strcat(BaseFileName,'_IPFDensity.png');
    export_fig(OutputFileName2,'-r150');
end
