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
