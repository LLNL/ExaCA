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
