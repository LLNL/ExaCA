% Copyright Lawrence Livermore National Security, LLC, and other ExaCA Project Developers.
% See the top-level LICENSE file for details.
%
% SPDX-License-Identifier: MIT
%
% NIST-developed software is provided by NIST as a public service. You may use, copy, and distribute copies of the software in any medium, provided that you keep intact this entire notice. You may improve, modify, and create derivative works of the software or any portion of the software, and you may copy and distribute such modifications or works. Modified works should carry a notice stating that you changed the software and should note the date and nature of any such change. Please explicitly acknowledge the National Institute of Standards and Technology as the source of the software.

% NIST-developed software is expressly provided "AS IS." NIST MAKES NO WARRANTY OF ANY KIND, EXPRESS, IMPLIED, IN FACT, OR ARISING BY OPERATION OF LAW, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NON-INFRINGEMENT, AND DATA ACCURACY. NIST NEITHER REPRESENTS NOR WARRANTS THAT THE OPERATION OF THE SOFTWARE WILL BE UNINTERRUPTED OR ERROR-FREE, OR THAT ANY DEFECTS WILL BE CORRECTED. NIST DOES NOT WARRANT OR MAKE ANY REPRESENTATIONS REGARDING THE USE OF THE SOFTWARE OR THE RESULTS THEREOF, INCLUDING BUT NOT LIMITED TO THE CORRECTNESS, ACCURACY, RELIABILITY, OR USEFULNESS OF THE SOFTWARE.

% You are solely responsible for determining the appropriateness of using and distributing the software and you assume all risks associated with its use, including but not limited to the risks and costs of program errors, compliance with applicable laws, damage to or loss of data, programs or equipment, and the unavailability or interruption of operation. This software is not intended to be used in any situation where a failure could cause risk of injury or damage to property. The software developed by NIST employees is not subject to copyright protection within the United States.
%

function [] = PlotPoleFigure(MTEXFile, ijk_vals, xyz_vals, ColormapFile, ColormapUpperLimit)

    % Script created by Matt Rolchigo (ORNL) and Adam Creuziger (NIST)
    % Plots pole figure data, and inverse pole figure data for the file specified by `MTEXFile`, for crystal planes specified by ijk_vals (for pole figures) and reference directions specified by xyz_vals (for inverse pole figures) using the colormap specified by `ColormapFile`, with colormap bounds [0, `ColormapUpperLimit`]. See `analysis/README.txt` for more details
    % crystal symmetry
    Ext = strfind(MTEXFile,'.');
    SizeName = size(Ext);
    BaseFileName = extractBefore(MTEXFile,Ext(SizeName(2)));

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

    num_pole_figs = size(ijk_vals,2);
    % Pole figure for each input set of planes
    for i=1:num_pole_figs
        % Create an Pole figure variable containing the data and plot
        ijk_pf = ijk_vals{i};
        i_pf = ijk_pf(1);
        j_pf = ijk_pf(2);
        k_pf = ijk_pf(3);
        h = Miller(i_pf,j_pf,k_pf,CS);
        r = regularS2Grid('resolution',5*degree);
        pfs_VPSClike = calcPoleFigure(odf_VPSClike,h,r);

        figure(i); plot(pfs_VPSClike,'contourf');
        hold on;
        clim(gca,[0, ColormapUpperLimit]);
        mtexColorbar;
        hold off;
        pf_plane_name = strcat(num2str(i_pf),num2str(j_pf),num2str(k_pf));
        OutputFileName = strcat(BaseFileName,'_',pf_plane_name,'_PoleFigure.png');
        export_fig(OutputFileName,'-r150');
    end
  
    num_xyz = size(xyz_vals,2);
    for i=1:num_xyz
        if ((strcmp(xyz_vals(i),'x')) || (strcmp(xyz_vals(i),'X')))
            ref_vector = xvector;
        elseif ((strcmp(xyz_vals(i),'y')) || (strcmp(xyz_vals(i),'Y')))
            ref_vector = yvector;
        elseif ((strcmp(xyz_vals(i),'z')) || (strcmp(xyz_vals(i),'Z')))
            ref_vector = zvector;
        else
            error('Error: Inputs for xyz_vals must be x, y, or z');
        end
        figure(i+num_pole_figs);
        hold on;
        plotIPDF(odf_VPSClike,ref_vector,'contourf');
        clim(gca,[0, ColormapUpperLimit]);
        mtexColorbar;
        hold off;
        OutputFileName2 = strcat(BaseFileName,'_IPFDensity.png');
        export_fig(OutputFileName2,'-r150');
    end
end
