clear
clc
close all

%% 文件设置
scriptDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(scriptDir, '..', 'data');
outputDir = fullfile(scriptDir, '..', 'pics');
files = dir(fullfile(dataDir, '*.txt'));
files = files(~strcmpi({files.name}, '00.txt') & ...
    ~strcmpi({files.name}, 'table01.txt'));
[~, order] = sort({files.name});
files = files(order);
assert(~isempty(files), 'data文件夹中没有TXT波形文件。');
speedKmh = NaN(numel(files), 1);

%% 通道与单位设置
channelColumns = [1, 2];
scaleToG = [1, 1];

%% 时间段设置
segmentLengthSec = 5;
analysisRangeSec = [0, Inf];
manualWindowsSec = [];

%% 频谱与评价设置
analysisBandHz = [0.5, 40];
evaluationLimits = [2.5, 2.75];

%% 绘图设置
matlabBlue = [0, 0.4470, 0.7410];
matlabOrange = [0.8500, 0.3250, 0.0980];
fontSize = 12;
directionNames = {'垂向', '横向'};
subplotSizeCm = [25, 5];
cases = cell(numel(files), 1);
results = table;

%% 读取实时指标
realtime = readRealtimeIndices(fullfile(dataDir, 'table01.txt'));

for i = 1:numel(files)
    %% 读取波形
    [timeSec, accelerationG, deltat] = readWaveformTxt(fullfile(dataDir, files(i).name));
    accelerationG = accelerationG(:, channelColumns) .* scaleToG;
    [~, caseName] = fileparts(files(i).name);

    %% 生成计算区间
    windowsSec = makeWindows(timeSec, deltat, segmentLengthSec, ...
        analysisRangeSec, manualWindowsSec);
    segmentCount = size(windowsSec, 1);
    assert(segmentCount > 0, '没有有效计算段，请修改时间段设置。');
    actualStartSec = zeros(segmentCount, 1);
    actualEndSec = zeros(segmentCount, 1);
    sampleCounts = zeros(segmentCount, 1);
    frequencyStepHz = zeros(segmentCount, 1);
    sperlingIndex = zeros(segmentCount, 2);
    fprintf('%s：采样率%.8g Hz，有效计算段%d个。\n', ...
        files(i).name, 1/deltat, segmentCount);

    for s = 1:segmentCount
        %% 截取半开区间，避免边界重复
        selected = timeSec >= windowsSec(s, 1)-deltat*1e-8 & ...
            timeSec < windowsSec(s, 2)-deltat*1e-8;
        tSegment = timeSec(selected);
        xSegment = accelerationG(selected, :);
        N = size(xSegment, 1);
        assert(N >= 2, '第%d段不足两个采样点。', s);
        actualStartSec(s) = tSegment(1);
        actualEndSec(s) = tSegment(end) + deltat;
        sampleCounts(s) = N;
        frequencyStepHz(s) = 1 / (N*deltat);
        fprintf('  第%d段：%.8g~%.8g s，%d点，时长%.8g s。\n', ...
            s, actualStartSec(s), actualEndSec(s), N, N*deltat);

        %% 加Hann窗计算单边幅值谱，不补零
        [frequencyHz, amplitudeG] = amplitudeSpectrum(xSegment, deltat);
        inBand = frequencyHz > 0 & frequencyHz >= analysisBandHz(1) & ...
            frequencyHz <= analysisBandHz(2);
        assert(any(inBand), '第%d段在计算频带内没有频率点。', s);
        f = frequencyHz(inBand);

        %% 按频率修正并合成Sperling指标
        for c = 1:2
            F = frequencyCorrection(f, c);
            q = amplitudeG(inBand, c).^3 ./ f .* F;
            sperlingIndex(s, c) = 7.08 * sum(q)^0.1;
        end
    end

    %% 保存本组绘图数据
    cases{i} = struct('name', caseName, 'timeSec', timeSec, ...
        'accelerationG', accelerationG, ...
        'actualStartSec', actualStartSec, 'actualEndSec', actualEndSec, ...
        'sperlingIndex', sperlingIndex);

    %% 输出分段结果
    caseResults = table(repmat(string(files(i).name), segmentCount, 1), ...
        repmat(speedKmh(i), segmentCount, 1), (1:segmentCount).', ...
        windowsSec(:, 1), windowsSec(:, 2), actualStartSec, actualEndSec, ...
        sampleCounts*deltat, sampleCounts, frequencyStepHz, ...
        sperlingIndex(:, 1), sperlingIndex(:, 2), ...
        'VariableNames', {'File', 'Speed_km_h', 'Segment', ...
        'WindowStart_s', 'WindowEnd_s', 'ActualStart_s', 'ActualEnd_s', ...
        'Duration_s', 'SampleCount', 'FrequencyStep_Hz', 'W_vertical', 'W_lateral'});
    results = [results; caseResults];
end
disp(results)
writetable(results, fullfile(dataDir, 'sperling_results.csv'));

%% 汇总各组均值
caseNames = strings(numel(cases), 1);
meanValues = zeros(numel(cases), 4);
for i = 1:numel(cases)
    caseNames(i) = string(cases{i}.name);
    for c = 1:2
        selected = realtime.File == caseNames(i) & ...
            realtime.Direction == string(directionNames{c});
        meanValues(i, 2*c-1) = realtime.MeanW(selected);
        meanValues(i, 2*c) = mean(cases{i}.sperlingIndex(:, c));
    end
end
comparison = [table(caseNames, 'VariableNames', {'File'}), ...
    array2table(meanValues, 'VariableNames', ...
    {'W_vertical_realtime', 'W_vertical_postprocessed', ...
    'W_lateral_realtime', 'W_lateral_postprocessed'})];
disp(comparison)
writetable(comparison, fullfile(dataDir, 'sperling_comparison.csv'));

%% 绘制时域与分段指标汇总图
fig = drawSummary(cases, directionNames, '加速度波形图', 'time', ...
    [matlabBlue; matlabOrange], fontSize, evaluationLimits, subplotSizeCm);
savePng(fig, outputDir, 'time_waveforms');
fig = drawSummary(cases, directionNames, '平稳性指标', 'index', ...
    [matlabBlue; matlabOrange], fontSize, evaluationLimits, subplotSizeCm);
savePng(fig, outputDir, 'sperling_vs_time');

%% 绘制均值对比图
fig = drawComparison(caseNames, meanValues, directionNames, ...
    [matlabBlue; matlabOrange], fontSize, evaluationLimits);
savePng(fig, outputDir, 'sperling_comparison');


function realtime = readRealtimeIndices(inputFile)
lines = readlines(inputFile, 'Encoding', 'UTF-8');
lines = lines(strlength(strtrim(lines)) > 0);
caseNames = strings(numel(lines), 1);
directions = strings(numel(lines), 1);
meanW = zeros(numel(lines), 1);
for i = 1:numel(lines)
    fields = split(strtrim(lines(i)), char(9));
    caseNames(i) = fields(1);
    directions(i) = fields(2);
    meanW(i) = mean(str2double(fields(3:end)), 'omitnan');
end
realtime = table(caseNames, directions, meanW, ...
    'VariableNames', {'File', 'Direction', 'MeanW'});
end


function [timeSec, accelerationG, deltat] = readWaveformTxt(inputFile)
lines = readlines(inputFile);
settings = split(lines(startsWith(lines, 'delta t')), char(9));
deltat = str2double(settings(1, 2));
assert(all(str2double(settings(:, 2:3)) == deltat, 'all'), ...
    'TXT中的采样间隔必须一致。');
dataRows = ~cellfun('isempty', regexp(cellstr(lines), '^\d{4}/'));
data = split(lines(dataRows), char(9));
accelerationG = str2double(data(:, 2:3));
timeSec = (0:size(accelerationG, 1)-1).' * deltat;
end


function windowsSec = makeWindows(t, deltat, segmentLengthSec, ...
    analysisRangeSec, manualWindowsSec)
if ~isempty(manualWindowsSec)
    assert(size(manualWindowsSec, 2) == 2 && ...
        all(manualWindowsSec(:, 1) < manualWindowsSec(:, 2)) && ...
        all(manualWindowsSec(:, 1) >= t(1)) && ...
        all(manualWindowsSec(:, 2) <= t(end)+deltat+deltat*1e-8), ...
        '手动区间必须为[start,end]，且位于数据时间范围内。');
    rangesSec = sortrows(manualWindowsSec, 1);
    assert(all(rangesSec(2:end, 1) >= rangesSec(1:end-1, 2)), ...
        '手动区间不能重叠。');
else
    rangesSec = [max(t(1), analysisRangeSec(1)), ...
        min(t(end)+deltat, analysisRangeSec(2))];
end
segmentSamples = round(segmentLengthSec / deltat);
assert(segmentSamples >= 2, '分段时长不足两个采样点。');
windowsSec = zeros(0, 2);
for r = 1:size(rangesSec, 1)
    firstIndex = find(t >= rangesSec(r, 1)-deltat*1e-8, 1);
    lastIndex = find(t < rangesSec(r, 2)-deltat*1e-8, 1, 'last');
    assert(~isempty(firstIndex) && ~isempty(lastIndex) && firstIndex <= lastIndex, ...
        '有效分析区间内没有数据。');
    segmentCount = fix((lastIndex-firstIndex+1)*deltat / segmentLengthSec);
    % 按附录A：足够两段时分段并舍弃尾段，否则使用整段。
    if segmentCount > 1
        starts = firstIndex + (0:segmentCount-1).' * segmentSamples;
        ends = starts + segmentSamples;
    else
        starts = firstIndex;
        ends = lastIndex + 1;
    end
    windowsSec = [windowsSec; t(starts), t(ends-1)+deltat];
end
end


function [frequencyHz, amplitudeG] = amplitudeSpectrum(x, deltat)
N = size(x, 1);
x = x - mean(x, 1);
window = 0.5 - 0.5*cos(2*pi*(0:N-1).'/N); % 周期Hann窗
spectrum = fft(x .* window, [], 1);
positiveBins = floor(N/2) + 1;
frequencyHz = (0:positiveBins-1).' / (N*deltat);
amplitudeG = abs(spectrum(1:positiveBins, :)) / sum(window); % 补偿窗的相干增益
if mod(N, 2) == 0
    amplitudeG(2:end-1, :) = 2 * amplitudeG(2:end-1, :);
else
    amplitudeG(2:end, :) = 2 * amplitudeG(2:end, :);
end
end


function F = frequencyCorrection(f, direction)
F = zeros(size(f));
if direction == 1
    low = f >= 0.5 & f < 5.9;
    middle = f >= 5.9 & f < 20;
    F(low) = 0.325 * f(low).^2;
    F(middle) = 400 ./ f(middle).^2;
    F(f >= 20 & f <= 40) = 1;
else
    low = f >= 0.5 & f < 5.4;
    middle = f >= 5.4 & f < 26;
    F(low) = 0.8 * f(low).^2;
    F(middle) = 650 ./ f(middle).^2;
    F(f >= 26 & f <= 40) = 1;
end
end


function fig = drawSummary(cases, directionNames, heading, kind, ...
    colors, fontSize, evaluationLimits, subplotSizeCm)
%% 图像设置
figureSize = [2*subplotSizeCm(1), numel(cases)*subplotSizeCm(2)];
fig = figure;
set(fig, 'Units', 'centimeters', 'Position', [1, 1, figureSize], ...
    'PaperUnits', 'centimeters', 'PaperPosition', [0, 0, figureSize], ...
    'PaperSize', figureSize, 'Color', 'w')
layout = tiledlayout(fig, numel(cases), 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:numel(cases)
    item = cases{i};
    caseName = item.name;
    for c = 1:2
        color = colors(c, :);
        ax = nexttile(layout, 2*(i-1)+c);
        if strcmp(kind, 'time')
            plot(ax, item.timeSec, item.accelerationG(:, c), 'Color', color, 'LineWidth', 0.75)
            xLabel = '\fontname{SimSun}时间/\fontname{Times New Roman}s';
            yLabel = '\fontname{SimSun}加速度/\fontname{Times New Roman}\it g';
            limits = [item.timeSec(1), item.timeSec(end)];
        else
            %% 绘制分段指标
            starts = item.actualStartSec;
            ends = item.actualEndSec;
            values = item.sperlingIndex(:, c);
            x = reshape([starts, ends, NaN(size(starts))].', [], 1);
            y = reshape([values, values, NaN(size(values))].', [], 1);
            plot(ax, x, y, 'Color', color, 'LineWidth', 1.1)
            hold(ax, 'on')
            plot(ax, (starts+ends)/2, values, 'o', 'Color', color, ...
                'MarkerFaceColor', 'w', 'MarkerSize', 4)
            text(ax, (starts+ends)/2, values+0.08, compose('%.3f', values), ...
                'Color', color, 'FontName', 'Times New Roman', 'FontSize', fontSize, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                'Interpreter', 'none')
            for limit = evaluationLimits
                limitLine = yline(ax, limit, '--', num2str(limit), 'Color', [0.4, 0.4, 0.4], ...
                    'FontName', 'Times New Roman', 'FontSize', fontSize, ...
                    'LabelVerticalAlignment', 'top');
                if limit == 2.5
                    limitLine.LabelVerticalAlignment = 'bottom';
                end
            end
            xLabel = '\fontname{SimSun}时间/\fontname{Times New Roman}s';
            yLabel = '\fontname{Times New Roman}\it W';
            limits = [starts(1), ends(end)];
            ylim(ax, [0, 4])
        end
        styleAxes(ax, [directionNames{c} heading], caseName, fontSize);
        xlabel(ax, xLabel, 'Interpreter', 'tex', ...
            'FontName', 'SimSun', 'FontSize', fontSize)
        ylabel(ax, yLabel, 'Interpreter', 'tex', ...
            'FontName', 'SimSun', 'FontSize', fontSize)
        xlim(ax, limits)
    end
end
end


function fig = drawComparison(caseNames, meanValues, directionNames, colors, ...
    fontSize, evaluationLimits)
%% 图像设置
figureSize = [25, 28];
fig = figure;
set(fig, 'Units', 'centimeters', 'Position', [1, 1, figureSize], ...
    'PaperUnits', 'centimeters', 'PaperPosition', [0, 0, figureSize], ...
    'PaperSize', figureSize, 'Color', 'w')
layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
offsets = [-0.2, 0.2];
for c = 1:2
    ax = nexttile(layout);
    hold(ax, 'on')

    %% 每组依次排列实时、后处理均值
    pairColors = [0.55*colors(c, :)+0.45; colors(c, :)];
    values = meanValues(:, 2*c-1:2*c);
    bars = gobjects(2, 1);
    for j = 1:2
        y = (1:numel(caseNames)).' + offsets(j);
        bars(j) = barh(ax, y, values(:, j), 0.34, ...
            'FaceColor', pairColors(j, :), 'EdgeColor', 'none');
        text(ax, values(:, j)+0.03, y, compose('%.3f', values(:, j)), ...
            'FontName', 'Times New Roman', 'FontSize', fontSize, ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle')
    end
    set(ax, 'YDir', 'reverse', 'YTick', 1:numel(caseNames), ...
        'YTickLabel', caseNames, 'TickLabelInterpreter', 'none', ...
        'FontName', 'Times New Roman', 'FontSize', fontSize, ...
        'XGrid', 'on', 'YGrid', 'off', 'Box', 'on')
    xlim(ax, [0, 3])
    ylim(ax, [0.4, numel(caseNames)+1.2])

    %% 绘制评价界限
    for limit = evaluationLimits
        xline(ax, limit, '--', num2str(limit), ...
            'Color', [0.4, 0.4, 0.4], 'FontName', 'Times New Roman', ...
            'FontSize', fontSize, 'LabelOrientation', 'aligned', ...
            'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'right', ...
            'HandleVisibility', 'off');
    end
    xlabel(ax, '\fontname{SimSun}平均平稳性指标 \fontname{Times New Roman}\it W', ...
        'Interpreter', 'tex', 'FontSize', fontSize)
    ylabel(ax, '数据组', 'FontName', 'SimSun', 'FontSize', fontSize)
    title(ax, [directionNames{c} '平稳性指标均值对比'], ...
        'FontName', 'SimSun', 'FontSize', fontSize)
    legend(ax, bars, {'实时计算', '后处理计算'}, ...
        'FontName', 'SimSun', 'FontSize', fontSize, ...
        'Location', 'southeast', 'Box', 'on')
end
end


function savePng(fig, outputDir, imageName)
print(fig, fullfile(outputDir, sprintf('99_all_%s.png', imageName)), ...
    '-dpng', '-r300');
close(fig)
end


function styleAxes(ax, heading, caseName, fontSize)
grid(ax, 'on')
set(ax, 'FontName', 'Times New Roman', 'FontSize', fontSize)
title(ax, ['\fontname{SimSun}' heading ' \fontname{Times New Roman}- ' caseName], ...
    'Interpreter', 'tex', 'FontName', 'SimSun', 'FontSize', fontSize)
end
