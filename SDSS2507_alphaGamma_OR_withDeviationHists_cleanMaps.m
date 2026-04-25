function SDSS2507_paper_alphaGamma_OR_withDeviationHists_cleanMaps
%% ========================================================================
%  SDSS2507_paper_alphaGamma_OR_withDeviationHists_cleanMaps.m
%
%  Paper-style MTEX workflow for ferrite-austenite interphase OR.
%
%  Outputs per condition:
%    1) alpha/gamma misorientation histogram PNG + CSV
%    2) deviation-to-KS histogram PNG + CSV
%    3) deviation-to-NW histogram PNG + CSV
%    4) minimum-deviation histogram PNG + CSV
%    5) alpha/gamma segment table CSV
%    6) alpha/gamma OR summary CSV
%    7) CLEAN alpha/gamma OR map PNG (KS/NW only; Other hidden)
%
%  Combined output:
%    all_alphaGamma_OR_summary.csv
%
%  Important plotting choices:
%    - Histograms are normalized to NumberFraction
%    - Each histogram family uses the SAME y-axis maximum across samples
%    - OR maps are paper-style cleaned maps:
%         * Other not shown
%         * white specks suppressed by rasterized phase-map rendering
%         * internal striping suppressed by nearest-neighbor interpolation
%    - Other is still retained in the tables and summary outputs
%% ========================================================================

close all; clc;

%% ------------------------------------------------------------------------
% 0. Resolve paths and initialize MTEX
% -------------------------------------------------------------------------
thisFile = mfilename('fullpath');
if isempty(thisFile)
    error('Save this script as SDSS2507_paper_alphaGamma_OR_withDeviationHists_cleanMaps.m and run again.');
end

rootDir = fileparts(thisFile);
disp(['rootDir = ' rootDir]);

mtexDir = fullfile(rootDir,'mtex-6.1.0');
if exist(mtexDir,'dir')
    run(fullfile(mtexDir,'startup_mtex.m'));
else
    warning(['mtex-6.1.0 folder not found next to this script. ' ...
             'MTEX must already be on the MATLAB path.']);
end

setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');

oldDefaultFigureVisible = get(0,'DefaultFigureVisible');
set(0,'DefaultFigureVisible','off');

%% ------------------------------------------------------------------------
% 1. User parameters
% -------------------------------------------------------------------------

% ---------- input files ----------
files = { ...
    '01_AS.ang', ...
    '02_SR400.ang', ...
    '03_SR450.ang', ...
    '04_SR500.ang', ...
    '05_SR550.ang', ...
    '06_SA1100.ang'};

sampleLabels = { ...
    'AS', ...
    'SR400', ...
    'SR450', ...
    'SR500', ...
    'SR550', ...
    'SA1100'};

for i = 1:numel(files)
    assert(exist(fullfile(rootDir,files{i}),'file') == 2, ...
        'Missing file: %s', fullfile(rootDir,files{i}));
end

% ---------- import / reference-frame parameters ----------
pImport.edaxSetting = 'setting 2';
pImport.refConv     = 'convertEuler2SpatialReferenceFrame';
pImport.plotLikeOIM = true;

% ---------- display-only transform ----------
pDisplay.applyOIMLikeView = true;
pDisplay.flipX = false;
pDisplay.flipY = true;

% ---------- map colors ----------
pColor.ferrite    = [61 38 168] / 255;    % deep violet-purple
pColor.austenite  = [249 250 20] / 255;   % bright yellow
pColor.notIndexed = [1 1 1];              % white

% ---------- preprocessing ----------
pPrep.ciMin          = 0.10;
pPrep.grainTolDeg    = 5.0;
pPrep.minGrainPixels = 9;
pPrep.smoothN        = 4;

% ---------- OR classification ----------
pOR.tolDeg = 5.0;

% ---------- histogram bins ----------
pHist.misorientationEdgesDeg = 0:2:64;
pHist.deviationEdgesDeg      = 0:1:65;

% ---------- output ----------
outDir = fullfile(rootDir,'PAPER_alphaGamma_OR_withDeviationHists_cleanMaps');

if exist(outDir,'dir')
    oldPng = dir(fullfile(outDir,'*.png'));
    for k = 1:numel(oldPng)
        delete(fullfile(outDir,oldPng(k).name));
    end
    oldCsv = dir(fullfile(outDir,'*.csv'));
    for k = 1:numel(oldCsv)
        delete(fullfile(outDir,oldCsv(k).name));
    end
else
    mkdir(outDir);
end

TallSummary = table();

%% ------------------------------------------------------------------------
% 2. Main loop
% -------------------------------------------------------------------------
for iFile = 1:numel(files)

    dataFile    = files{iFile};
    sampleLabel = sampleLabels{iFile};

    fprintf('\n============================================================\n');
    fprintf('Processing %s (%s)\n', sampleLabel, dataFile);
    fprintf('============================================================\n');

    try
        %% ----------------------------------------------------------------
        % Load EBSD in self-consistent MTEX reference frame
        % -----------------------------------------------------------------
        fName = fullfile(rootDir,dataFile);
        ebsd = EBSD.load(fName, pImport.refConv, pImport.edaxSetting);

        if pImport.plotLikeOIM
            ebsd.how2plot.east = xvector;
            ebsd.how2plot.outOfScreen = -zvector;
        end

        %% ----------------------------------------------------------------
        % Resolve phase names
        % -----------------------------------------------------------------
        ferriteName   = resolvePhaseName(ebsd,'Ferrite');
        austeniteName = resolvePhaseName(ebsd,'Austenite');

        if isempty(ferriteName)
            error('Ferrite phase not found in %s.', dataFile);
        end
        if isempty(austeniteName)
            error('Austenite phase not found in %s.', dataFile);
        end

        %% ----------------------------------------------------------------
        % Keep Ferrite + Austenite + notIndexed only
        % -----------------------------------------------------------------
        keepMask = ~ebsd.isIndexed;
        keepMask = keepMask | ismember(ebsd.phaseId, unique(ebsd(ferriteName).phaseId));
        keepMask = keepMask | ismember(ebsd.phaseId, unique(ebsd(austeniteName).phaseId));
        ebsd = ebsd(keepMask);

        %% ----------------------------------------------------------------
        % CI thresholding
        % -----------------------------------------------------------------
        if isfield(ebsd.prop,'ci')
            lowCI = ebsd.isIndexed & (ebsd.prop.ci < pPrep.ciMin);
            ebsd(lowCI) = [];
        else
            warning('No CI field found in %s. CI threshold not applied.', dataFile);
        end

        %% ----------------------------------------------------------------
        % Grain reconstruction
        % -----------------------------------------------------------------
        [grains, ebsd('indexed').grainId] = calcGrains(ebsd('indexed'), ...
            'angle', pPrep.grainTolDeg * degree);

        smallGrains = grains(grains.numPixel < pPrep.minGrainPixels);
        if ~isempty(smallGrains)
            ebsd(smallGrains) = [];
        end

        [grains, ebsd('indexed').grainId] = calcGrains(ebsd('indexed'), ...
            'angle', pPrep.grainTolDeg * degree);

        grains = smooth(grains, pPrep.smoothN);

        %% ----------------------------------------------------------------
        % Ferrite-austenite boundaries only
        % -----------------------------------------------------------------
        aGrains = grains(austeniteName);
        fGrains = grains(ferriteName);
        gAF = grains.boundary(austeniteName, ferriteName);

        if isempty(gAF)
            warning('No ferrite-austenite boundaries found in %s.', sampleLabel);
            continue;
        end

        [aSegID, fSegID] = resolveBoundaryPhaseIDs(gAF, aGrains, fGrains);

        %% ----------------------------------------------------------------
        % alpha/gamma misorientation and OR classification
        % -----------------------------------------------------------------
        csA = ebsd(austeniteName).CS;
        csF = ebsd(ferriteName).CS;

        KS = orientation.KurdjumovSachs(csA, csF);
        NW = orientation.NishiyamaWassermann(csA, csF);

        mori   = gAF.misorientation;
        segLen = segLength(gAF);

        thetaAG = angle(mori) ./ degree;
        devKS   = angle(mori, KS) ./ degree;
        devNW   = angle(mori, NW) ./ degree;
        devMin  = min(devKS, devNW);

        [orClass, isKS, isNW, isOther] = classifyORNearest(devKS, devNW, pOR.tolDeg);

        %% ----------------------------------------------------------------
        % True boundary coordinates
        % -----------------------------------------------------------------
        [x1, y1, x2, y2] = extractBoundarySegmentsTrue(gAF);

        fprintf('True-boundary extraction successful for %s: %d alpha/gamma segments.\n', ...
            sampleLabel, numel(x1));

        %% ----------------------------------------------------------------
        % Segment table
        % -----------------------------------------------------------------
        Tseg = table( ...
            repmat(string(sampleLabel), length(gAF), 1), ...
            (1:length(gAF))', ...
            aSegID(:), fSegID(:), ...
            thetaAG(:), segLen(:), ...
            devKS(:), devNW(:), devMin(:), orClass(:), ...
            x1(:), y1(:), x2(:), y2(:), ...
            'VariableNames', {'Sample','SegmentID','AGrainID','FGrainID', ...
            'AlphaGammaMisorientation_deg','Length_um', ...
            'DevKS_deg','DevNW_deg','DevMin_deg','ORClass', ...
            'x1_um','y1_um','x2_um','y2_um'});

        writetable(Tseg, fullfile(outDir,[sampleLabel '_alphaGamma_segment_table.csv']));

        %% ----------------------------------------------------------------
        % OR summary table
        % -----------------------------------------------------------------
        Tsum = summarizeAlphaGammaOR(sampleLabel, thetaAG, segLen, isKS, isNW, isOther, ...
            devKS, devNW, devMin);
        writetable(Tsum, fullfile(outDir,[sampleLabel '_alphaGamma_OR_summary.csv']));
        TallSummary = [TallSummary; Tsum]; %#ok<AGROW>

        disp(' ');
        disp(['================ ' sampleLabel ' alpha/gamma OR summary ================']);
        disp(Tsum);

        %% ----------------------------------------------------------------
        % Histogram CSVs
        % -----------------------------------------------------------------
        writeHistogramCSV(thetaAG, pHist.misorientationEdgesDeg, ...
            fullfile(outDir,[sampleLabel '_alphaGamma_misorientation_histogram.csv']));

        writeHistogramCSV(devKS, pHist.deviationEdgesDeg, ...
            fullfile(outDir,[sampleLabel '_alphaGamma_devKS_histogram.csv']));

        writeHistogramCSV(devNW, pHist.deviationEdgesDeg, ...
            fullfile(outDir,[sampleLabel '_alphaGamma_devNW_histogram.csv']));

        writeHistogramCSV(devMin, pHist.deviationEdgesDeg, ...
            fullfile(outDir,[sampleLabel '_alphaGamma_devMin_histogram.csv']));

        %% ----------------------------------------------------------------
        % CLEAN paper-style OR map: KS / NW only, Other hidden
        % -----------------------------------------------------------------
        renderCleanORMapTrue(ebsd, x1, y1, x2, y2, isKS, isNW, isOther, ...
            ferriteName, austeniteName, outDir, ...
            [sampleLabel '_alphaGamma_OR_map_trueBoundary_clean'], ...
            [sampleLabel ' | alpha/gamma OR map (KS/NW only)'], ...
            pDisplay, pColor);

    catch ME
        warning('Processing failed for %s: %s', sampleLabel, ME.message);
    end
end

%% ------------------------------------------------------------------------
%% ------------------------------------------------------------------------
% 3. Render all histograms with family-specific common normalized y-axes
%    but force custom SR450 y-ticks exactly as requested
% -------------------------------------------------------------------------
globalY_mis  = getGlobalHistYMaxFromCSV(outDir, sampleLabels, '_alphaGamma_misorientation_histogram.csv');
globalY_ks   = getGlobalHistYMaxFromCSV(outDir, sampleLabels, '_alphaGamma_devKS_histogram.csv');
globalY_nw   = getGlobalHistYMaxFromCSV(outDir, sampleLabels, '_alphaGamma_devNW_histogram.csv');
globalY_min  = getGlobalHistYMaxFromCSV(outDir, sampleLabels, '_alphaGamma_devMin_histogram.csv');

for i = 1:numel(sampleLabels)
    s = sampleLabels{i};

    if strcmp(s,'SR450')
        % ---------------- misorientation ----------------
        renderHistogramFromCSV( ...
            fullfile(outDir,[s '_alphaGamma_misorientation_histogram.csv']), ...
            fullfile(outDir,[s '_alphaGamma_misorientation_histogram']), ...
            [s ' | \alpha/\gamma misorientation histogram'], ...
            '\alpha/\gamma misorientation angle (deg)', ...
            0.45, ...
            [0 0.09 0.18 0.27 0.36 0.45], ...
            {'0','0.09','0.18','0.27','0.36','0.45'});

        % ---------------- deviation from KS ----------------
        renderHistogramFromCSV( ...
            fullfile(outDir,[s '_alphaGamma_devKS_histogram.csv']), ...
            fullfile(outDir,[s '_alphaGamma_devKS_histogram']), ...
            [s ' | deviation from K-S'], ...
            'Deviation from K-S (deg)', ...
            0.25, ...
            [0 0.05 0.10 0.15 0.20 0.25], ...
            {'0','0.05','0.10','0.15','0.20','0.25'});

        % ---------------- deviation from NW ----------------
        renderHistogramFromCSV( ...
            fullfile(outDir,[s '_alphaGamma_devNW_histogram.csv']), ...
            fullfile(outDir,[s '_alphaGamma_devNW_histogram']), ...
            [s ' | deviation from N-W'], ...
            'Deviation from N-W (deg)', ...
            0.25, ...
            [0 0.05 0.10 0.15 0.20 0.25], ...
            {'0','0.05','0.10','0.15','0.20','0.25'});

        % ---------------- minimum deviation ----------------
        renderHistogramFromCSV( ...
            fullfile(outDir,[s '_alphaGamma_devMin_histogram.csv']), ...
            fullfile(outDir,[s '_alphaGamma_devMin_histogram']), ...
            [s ' | minimum deviation to nearest rational OR'], ...
            'Minimum deviation, min(\Delta\theta_{KS}, \Delta\theta_{NW}) (deg)', ...
            0.30, ...
            [0 0.06 0.12 0.18 0.24 0.30], ...
            {'0','0.06','0.12','0.18','0.24','0.30'});

    else
        renderHistogramFromCSV( ...
            fullfile(outDir,[s '_alphaGamma_misorientation_histogram.csv']), ...
            fullfile(outDir,[s '_alphaGamma_misorientation_histogram']), ...
            [s ' | \alpha/\gamma misorientation histogram'], ...
            '\alpha/\gamma misorientation angle (deg)', ...
            globalY_mis);

        renderHistogramFromCSV( ...
            fullfile(outDir,[s '_alphaGamma_devKS_histogram.csv']), ...
            fullfile(outDir,[s '_alphaGamma_devKS_histogram']), ...
            [s ' | deviation from K-S'], ...
            'Deviation from K-S (deg)', ...
            globalY_ks);

        renderHistogramFromCSV( ...
            fullfile(outDir,[s '_alphaGamma_devNW_histogram.csv']), ...
            fullfile(outDir,[s '_alphaGamma_devNW_histogram']), ...
            [s ' | deviation from N-W'], ...
            'Deviation from N-W (deg)', ...
            globalY_nw);

        renderHistogramFromCSV( ...
            fullfile(outDir,[s '_alphaGamma_devMin_histogram.csv']), ...
            fullfile(outDir,[s '_alphaGamma_devMin_histogram']), ...
            [s ' | minimum deviation to nearest rational OR'], ...
            'Minimum deviation, min(\Delta\theta_{KS}, \Delta\theta_{NW}) (deg)', ...
            globalY_min);
    end
end
%% ------------------------------------------------------------------------
% 4. Combined outputs
% -------------------------------------------------------------------------
writetable(TallSummary, fullfile(outDir,'all_alphaGamma_OR_summary.csv'));

disp(' ');
disp('Paper-style alpha/gamma OR outputs completed.');
disp(['All outputs are in: ' outDir]);

set(0,'DefaultFigureVisible',oldDefaultFigureVisible);

end

%% ========================================================================
% Helper functions
%% ========================================================================

function phaseName = resolvePhaseName(ebsd, phaseQuery)

names = cellstr(ebsd.mineralList(:));
names = names(~cellfun(@isempty,names));

phaseName = '';

if isempty(names)
    return;
end

idx = find(strcmpi(names, phaseQuery), 1, 'first');
if isempty(idx)
    idx = find(contains(lower(string(names)), lower(string(phaseQuery))), 1, 'first');
end

if ~isempty(idx)
    phaseName = names{idx};
end

end

function [aSegID, fSegID] = resolveBoundaryPhaseIDs(gAF, aGrains, fGrains)

g1 = gAF.grainId(:,1);
g2 = gAF.grainId(:,2);

aIDs = aGrains.id(:);
fIDs = fGrains.id(:);

isCol1A = ismember(g1, aIDs) & ismember(g2, fIDs);
isCol2A = ismember(g2, aIDs) & ismember(g1, fIDs);

if ~all(isCol1A | isCol2A)
    error('Could not resolve austenite/ferrite grain IDs on alpha/gamma boundaries.');
end

aSegID = zeros(length(gAF),1);
fSegID = zeros(length(gAF),1);

aSegID(isCol1A) = g1(isCol1A);
fSegID(isCol1A) = g2(isCol1A);

aSegID(isCol2A) = g2(isCol2A);
fSegID(isCol2A) = g1(isCol2A);

end

function [orClass, isKS, isNW, isOther] = classifyORNearest(devKS, devNW, tolDeg)

isKS = (devKS <= tolDeg) & (devKS < devNW);
isNW = (devNW <= tolDeg) & (devNW < devKS);

tieMask = (devKS <= tolDeg) & (devNW <= tolDeg) & (abs(devKS - devNW) < 1e-12);
isKS(tieMask) = true;

isOther = ~(isKS | isNW);

orClass = strings(size(devKS));
orClass(isKS)    = "KS";
orClass(isNW)    = "NW";
orClass(isOther) = "Other";

end

function [x1, y1, x2, y2] = extractBoundarySegmentsTrue(gB)

nSeg = length(gB);

if nSeg == 0
    x1 = zeros(0,1); y1 = zeros(0,1);
    x2 = zeros(0,1); y2 = zeros(0,1);
    return;
end

mp = gB.midPoint;
mx = double(mp.x(:));
my = double(mp.y(:));

dir = gB.direction;
dx = double(dir.x(:));
dy = double(dir.y(:));

L = double(segLength(gB));
L = L(:);

assert(numel(mx) == nSeg, 'midPoint.x size does not match number of segments.');
assert(numel(my) == nSeg, 'midPoint.y size does not match number of segments.');
assert(numel(dx) == nSeg, 'direction.x size does not match number of segments.');
assert(numel(dy) == nSeg, 'direction.y size does not match number of segments.');
assert(numel(L)  == nSeg, 'segLength size does not match number of segments.');

dnorm = sqrt(dx.^2 + dy.^2);
assert(all(dnorm > eps), 'At least one boundary segment has zero direction norm.');

dx = dx ./ dnorm;
dy = dy ./ dnorm;

halfL = 0.5 .* L;

x1 = mx - halfL .* dx;
y1 = my - halfL .* dy;
x2 = mx + halfL .* dx;
y2 = my + halfL .* dy;

assert(all(isfinite(x1)) && all(isfinite(y1)) && ...
       all(isfinite(x2)) && all(isfinite(y2)), ...
       'Non-finite endpoint coordinates produced.');

end

function T = summarizeAlphaGammaOR(sampleLabel, thetaAG, segLen, isKS, isNW, isOther, ...
    devKS, devNW, devMin)

nSeg = numel(thetaAG);
totLen = sum(segLen);

nKS = sum(isKS);
nNW = sum(isNW);
nOther = sum(isOther);

lenKS = sum(segLen(isKS));
lenNW = sum(segLen(isNW));
lenOther = sum(segLen(isOther));

T = table( ...
    string(sampleLabel), ...
    double(nSeg), double(totLen), ...
    double(median(thetaAG)), ...
    double(nKS), double(nNW), double(nOther), ...
    double(nKS/nSeg), double(nNW/nSeg), double(nOther/nSeg), ...
    double(lenKS/totLen), double(lenNW/totLen), double(lenOther/totLen), ...
    double(median(devKS)), double(median(devNW)), double(median(devMin)), ...
    'VariableNames', {'Sample','nAlphaGammaSegments','TotalAlphaGammaLen_um', ...
    'MedianAlphaGammaAngle_deg', ...
    'nKS','nNW','nOther', ...
    'CountFracKS','CountFracNW','CountFracOther', ...
    'LenFracKS','LenFracNW','LenFracOther', ...
    'MedianDevKS_deg','MedianDevNW_deg','MedianDevMin_deg'});

end

function writeHistogramCSV(valuesDeg, edgesDeg, outCsv)

counts = histcounts(valuesDeg, edgesDeg, 'Normalization', 'probability');

binLeft = edgesDeg(1:end-1).';
binRight = edgesDeg(2:end).';
binCenter = (binLeft + binRight) / 2;

Th = table(binLeft, binRight, binCenter, counts(:), ...
    'VariableNames', {'BinLeft_deg','BinRight_deg','BinCenter_deg','NumberFraction'});

writetable(Th, outCsv);

end

function globalYMax = getGlobalHistYMaxFromCSV(outDir, sampleLabels, suffix)

globalYMax = 0;

for i = 1:numel(sampleLabels)
    f = fullfile(outDir, [sampleLabels{i} suffix]);
    if ~exist(f,'file')
        continue;
    end

    T = readtable(f);
    if isempty(T)
        continue;
    end

    thisMax = max(T.NumberFraction);
    if thisMax > globalYMax
        globalYMax = thisMax;
    end
end

globalYMax = ceil((1.05 * globalYMax) / 0.05) * 0.05;

if globalYMax <= 0
    globalYMax = 0.10;
end

end
%==========================================================================
function renderHistogramFromCSV(inCsv, outStem, ttl, xlab, yMaxForced, yTickValsForced, yTickLabelsForced)

if nargin < 5 || isempty(yMaxForced)
    yMaxForced = [];
end
if nargin < 6
    yTickValsForced = [];
end
if nargin < 7
    yTickLabelsForced = {};
end

if ~exist(inCsv,'file')
    return;
end

T = readtable(inCsv);
if isempty(T)
    return;
end

if isempty(yMaxForced)
    globalYMax = max(T.NumberFraction);
    globalYMax = ceil((1.05 * globalYMax) / 0.05) * 0.05;
    if globalYMax <= 0
        globalYMax = 0.10;
    end
else
    globalYMax = yMaxForced;
end

if isempty(yTickValsForced)
    if globalYMax <= 0.30
        yStep = 0.05;
    else
        yStep = 0.10;
    end
    yTickVals = 0:yStep:globalYMax;
else
    yTickVals = yTickValsForced;
end

fig = figure('Visible','off', 'Color','w', 'Units','pixels', 'Position',[100 100 900 700]);
ax = axes('Parent',fig, 'Position',[0.13 0.14 0.82 0.78]);
hold(ax,'on');

bar(ax, T.BinCenter_deg, T.NumberFraction, 1.0, ...
    'FaceColor',[0.20 0.35 0.80], ...
    'EdgeColor','k', ...
    'LineWidth',1.8);

xlabel(ax, xlab, 'FontWeight','bold', 'FontSize',22, 'Interpreter','tex');
ylabel(ax, 'Number fraction', 'FontWeight','bold', 'FontSize',22);
title(ax, ttl, 'FontWeight','bold', 'FontSize',22, 'Interpreter','tex');

set(ax, ...
    'FontWeight','bold', ...
    'FontSize',18, ...
    'LineWidth',2.5, ...
    'Box','on', ...
    'Layer','top', ...
    'TickDir','in');

xlim(ax, [min(T.BinLeft_deg) max(T.BinRight_deg)]);
ylim(ax, [0 globalYMax]);
yticks(ax, yTickVals);

if ~isempty(yTickLabelsForced)
    yticklabels(ax, yTickLabelsForced);
else
    yticklabels(ax, makePrettyTickLabels(yTickVals));
end

grid(ax, 'on');
ax.GridAlpha = 0.25;

disp(['Saved histogram: ' outStem]);
disp('YTick values used:');
disp(get(ax,'YTick'));
disp('YTick labels used:');
disp(get(ax,'YTickLabel'));

saveFigureTrue(fig, outStem);
close(fig);

end
%==========================================================================
function labels = makePrettyTickLabels(vals)

labels = cell(size(vals));

for i = 1:numel(vals)
    v = vals(i);

    if abs(v) < 1e-12
        labels{i} = '0';
    elseif abs(v*10 - round(v*10)) < 1e-12
        labels{i} = sprintf('%.1f', v);   % 0.1, 0.2, 0.3
    else
        labels{i} = sprintf('%.2f', v);   % 0.05, 0.09, 0.18, etc.
    end
end

end

%===================================================================================================
function renderCleanORMapTrue(ebsd, x1, y1, x2, y2, isKS, isNW, isOther, ...
    ferriteName, austeniteName, outDir, fileStem, ttl, pDisplay, pColor) %#ok<INUSD>

[fig, ax] = createSidebarFigure();

% -------------------------------------------------------------------------
% Build true raster RGB phase map to eliminate white specks and striping
% -------------------------------------------------------------------------
[imgRGB, xlimData, ylimData] = buildRasterPhaseMapRGB( ...
    ebsd, ferriteName, austeniteName, pDisplay, pColor);

image(ax, 'XData', xlimData, 'YData', ylimData, 'CData', imgRGB);
set(ax, 'YDir', 'normal');
hold(ax, 'on');

% -------------------------------------------------------------------------
% Overlay ONLY KS and NW
% -------------------------------------------------------------------------
colKS = [0.85 0.33 0.10];
colNW = [0.00 0.45 0.74];

[x1p, y1p] = transformDisplayCoords(x1, y1, ebsd, pDisplay);
[x2p, y2p] = transformDisplayCoords(x2, y2, ebsd, pDisplay);

drawSegments(ax, x1p, y1p, x2p, y2p, isNW, colNW, 1.20);
drawSegments(ax, x1p, y1p, x2p, y2p, isKS, colKS, 1.20);

xlim(ax, xlimData);
ylim(ax, ylimData);
axis(ax, 'equal');
set(ax, 'YDir', 'normal');

ax.Visible = 'off';
ax.LineWidth = 1.2;
set(ax, 'Color', pColor.ferrite);

title(ax, ttl, 'FontWeight','bold', 'FontSize',14, 'Interpreter','none');

addORSidebarClean(fig, colKS, colNW);
addScaleBarSidebar(fig, 40);

saveFigureTrue(fig, fullfile(outDir,fileStem));
close(fig);

end

function [imgRGB, xlimData, ylimData] = buildRasterPhaseMapRGB( ...
    ebsd, ferriteName, austeniteName, pDisplay, pColor)

% Use indexed points only for phase reconstruction.
eIdx = ebsd('indexed');
if isempty(eIdx)
    error('No indexed EBSD points available for raster phase map rendering.');
end

% -------------------------------------------------------------------------
% Transform coordinates consistently with overlay geometry
% -------------------------------------------------------------------------
[xp, yp] = transformDisplayCoords(eIdx.x, eIdx.y, ebsd, pDisplay);

xmin = min(xp);
xmax = max(xp);
ymin = min(yp);
ymax = max(yp);

ux = unique(sort(xp));
uy = unique(sort(yp));

if numel(ux) > 1
    dx = median(diff(ux));
else
    dx = 1;
end

if numel(uy) > 1
    dy = median(diff(uy));
else
    dy = 1;
end

if ~isfinite(dx) || dx <= 0
    dx = 1;
end
if ~isfinite(dy) || dy <= 0
    dy = 1;
end

% Grid centers
gx = xmin:dx:xmax;
gy = ymin:dy:ymax;

Nx = numel(gx);
Ny = numel(gy);

% Binary phase label
aPhaseIds = unique(ebsd(austeniteName).phaseId);
isA = ismember(eIdx.phaseId, aPhaseIds);
phaseVal = double(isA(:));

% Deduplicate points before interpolation
XY = [xp(:), yp(:)];
[XYu, ~, ic] = unique(XY, 'rows', 'stable');
phaseValU = accumarray(ic, phaseVal, [], @mean);
phaseValU = double(phaseValU > 0.5);

% Nearest-neighbor interpolation across full raster
F = scatteredInterpolant(XYu(:,1), XYu(:,2), phaseValU, 'nearest', 'nearest');
[Xg, Yg] = meshgrid(gx, gy);
phaseGrid = F(Xg, Yg);

isAgrid = phaseGrid > 0.5;

% RGB image
imgRGB = zeros(Ny, Nx, 3);
for c = 1:3
    chan = pColor.ferrite(c) * ones(Ny, Nx);
    chan(isAgrid) = pColor.austenite(c);
    imgRGB(:,:,c) = chan;
end

% Use cell-edge limits
xlimData = [gx(1) - dx/2, gx(end) + dx/2];
ylimData = [gy(1) - dy/2, gy(end) + dy/2];

end

function [fig, ax] = createSidebarFigure()

fig = figure('Visible','off', 'Color','w', 'Units','pixels', 'Position',[100 100 1200 760]);
ax = axes('Parent',fig, 'Position',[0.06 0.08 0.66 0.84]);
hold(ax,'on');

end

function addORSidebarClean(fig, colKS, colNW)

x0 = 0.78;
y0 = 0.82;
dy = 0.065;
linew = 0.030;

annotation(fig,'textbox',[0.76 0.88 0.18 0.04], ...
    'String','\alpha/\gamma OR class', ...
    'LineStyle','none', ...
    'FontWeight','bold', ...
    'FontSize',12, ...
    'Interpreter','tex');

annotation(fig,'line',[x0 x0+linew],[y0 y0], ...
    'Color',colKS, 'LineWidth',2.8);
annotation(fig,'textbox',[x0+0.04 y0-0.015 0.12 0.03], ...
    'String','KS', ...
    'LineStyle','none', 'FontSize',10, 'FontWeight','bold');

annotation(fig,'line',[x0 x0+linew],[y0-dy y0-dy], ...
    'Color',colNW, 'LineWidth',2.8);
annotation(fig,'textbox',[x0+0.04 y0-dy-0.015 0.12 0.03], ...
    'String','NW', ...
    'LineStyle','none', 'FontSize',10, 'FontWeight','bold');

annotation(fig,'textbox',[0.76 0.56 0.18 0.04], ...
    'String','Base map = phase colors', ...
    'LineStyle','none', ...
    'FontWeight','bold', ...
    'FontSize',10);

end

function addScaleBarSidebar(fig, L)

annotation(fig,'textbox',[0.76 0.28 0.16 0.04], ...
    'String','Scale bar', ...
    'LineStyle','none', ...
    'FontWeight','bold', ...
    'FontSize',12);

annotation(fig,'line',[0.78 0.88],[0.24 0.24], ...
    'Color','k', 'LineWidth',3.0);
annotation(fig,'textbox',[0.79 0.19 0.12 0.03], ...
    'String',sprintf('%g \\mum', L), ...
    'LineStyle','none', ...
    'FontWeight','bold', ...
    'FontSize',11);

end

function [xp, yp] = transformDisplayCoords(x, y, ebsd, pDisplay)

xp = x(:);
yp = y(:);

xmin = min(ebsd.x);
xmax = max(ebsd.x);
ymin = min(ebsd.y);
ymax = max(ebsd.y);

if pDisplay.applyOIMLikeView
    if pDisplay.flipX
        xp = xmax - (xp - xmin);
    end
    if pDisplay.flipY
        yp = ymax - (yp - ymin);
    end
end

end

function drawSegments(ax, x1, y1, x2, y2, mask, colorRGB, lineWidth)

idx = find(mask);
if isempty(idx)
    return;
end

X = nan(3*numel(idx),1);
Y = nan(3*numel(idx),1);

X(1:3:end) = x1(idx);
X(2:3:end) = x2(idx);

Y(1:3:end) = y1(idx);
Y(2:3:end) = y2(idx);

line(ax, X, Y, 'Color', colorRGB, 'LineWidth', lineWidth);

end

function saveFigureTrue(fig, fileStem)

drawnow;
pause(0.2);

try
    exportgraphics(fig, [fileStem '.png'], 'Resolution', 500);
    return;
catch
end

try
    set(fig,'PaperPositionMode','auto');
    print(fig, [fileStem '.png'], '-dpng', '-r500');
    return;
catch
end

warning('Figure export failed for: %s', fileStem);

end