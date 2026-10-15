function result = wvf_to_text(inputFile, exportCsv, outputDir)
% 默认批量转换相邻数据文件夹中的 WVF，输出到报告的 data 文件夹。
% 指定单个 WVF 时同时导出 CSV；第二个参数可控制 CSV 导出。
scriptDir = fileparts(mfilename('fullpath'));
if nargin == 0
    inputFile = fullfile(scriptDir, '..', '..', '数据 - 数采2');
end
if nargin < 2
    exportCsv = ~isfolder(inputFile);
end
if nargin < 3
    outputDir = fullfile(scriptDir, '..', 'data');
end

%% 遍历数据文件夹
if isfolder(inputFile)
    files = dir(fullfile(inputFile, '*.wvf'));
    assert(~isempty(files), '文件夹中没有 WVF 文件。');
    converted = cell(numel(files), 1);
    for i = 1:numel(files)
        converted{i} = wvf_to_text(fullfile(files(i).folder, files(i).name), exportCsv, outputDir);
    end
    result = vertcat(converted{:});
    fprintf('转换完成：共 %d 个 WVF 文件。\n', numel(files));
    return
end
[~, name] = fileparts(inputFile);
if ~isfolder(outputDir)
    mkdir(outputDir);
end
csvFile = '';
txtFile = fullfile(outputDir, [name '.txt']);
fprintf('\n文件：%s.wvf\n', name);

%% 读取文件头与记录索引
fid = fopen(inputFile, 'rb', 'ieee-be');
assert(fid >= 0, '无法打开 WVF 文件。');
inputCleanup = onCleanup(@() fclose(fid));
assert(strcmp(char(fread(fid, 4, '*uint8').'), 'DTLG'), ...
    '文件不是 LabVIEW DTLG 格式。');
version = fread(fid, 1, '*uint32');
recordCount = fread(fid, 1, 'uint32');
indexOffset = fread(fid, 1, 'uint32');
assert(version == uint32(hex2dec('08008000')), ...
    '当前脚本尚未验证此 DTLG 版本。');
fseek(fid, indexOffset, 'bof');
recordOffsets = fread(fid, recordCount, 'uint32');
fseek(fid, 0, 'eof');
fileSize = ftell(fid);
assert(numel(recordOffsets) == recordCount && ...
    all(diff(recordOffsets) > 0) && recordOffsets(end) < fileSize, ...
    '记录索引不完整。');

%% 逐条读取波形
records = cell(recordCount, 1);
for r = 1:recordCount
    fseek(fid, recordOffsets(r), 'bof');
    channelCount = fread(fid, 1, 'uint32');
    assert(channelCount == 2, '当前文件应包含两个通道。');
    channels = struct([]);
    for c = 1:channelCount
        channels(c).seconds = fread(fid, 1, '*int64');
        channels(c).fraction = fread(fid, 1, '*uint64');
        channels(c).dt = fread(fid, 1, 'double');
        sampleCount = fread(fid, 1, 'uint32');
        channels(c).Y = fread(fid, sampleCount, 'double');
        assert(numel(channels(c).Y) == sampleCount && channels(c).dt > 0, ...
            '波形数据不完整或采样间隔无效。');

        % 跳过波形内的错误簇，读取通道属性。
        errorStatus = fread(fid, 1, 'uint8');
        errorCode = fread(fid, 1, 'int32');
        errorSource = readString(fid);
        assert(errorStatus == 0 && errorCode == 0, ...
            '波形包含错误信息：%s', errorSource);
        [~, channels(c).attributes] = readVariant(fid);
    end
    assert(numel(channels(1).Y) == numel(channels(2).Y) && ...
        channels(1).dt == channels(2).dt, '两通道采样设置不一致。');
    if r < recordCount
        recordEnd = recordOffsets(r + 1);
    else
        recordEnd = fileSize;
    end
    assert(ftell(fid) == recordEnd, '第 %d 条记录解析边界不一致。', r);
    records{r} = channels;
end

%% 整理数据与相对时间
first = records{1};
sampleCounts = cellfun(@(ch) numel(ch(1).Y), records);
sampleTotal = sum(sampleCounts);
timeSec = zeros(sampleTotal, 1);
Y = zeros(sampleTotal, 2);
timestamps = strings(sampleTotal, 2);
row = 0;
for r = 1:recordCount
    channels = records{r};
    n = sampleCounts(r);
    rows = row + (1:n);
    offset = (0:n-1).' * channels(1).dt;
    timeSec(rows) = double(channels(1).seconds - first(1).seconds) ...
        + (double(channels(1).fraction) - double(first(1).fraction)) / 2^64 ...
        + offset;
    for c = 1:2
        Y(rows, c) = channels(c).Y;
        timestamps(rows, c) = formatTimestamp( ...
            channels(c).seconds, channels(c).fraction, offset);
    end
    row = row + n;
end
assert(all(diff(timeSec) > 0), '源文件有重复或倒序记录。');

%% 导出 CSV，保留双精度数值
if exportCsv
    csvFile = fullfile(outputDir, [name '.csv']);
    out = fopen(csvFile, 'wt');
    assert(out >= 0, '无法创建 CSV 文件。');
    csvCleanup = onCleanup(@() fclose(out));
    fprintf(out, 'Time_s,Y0,Y1\n');
    fprintf(out, '%.17g,%.17g,%.17g\n', [timeSec, Y].');
    clear csvCleanup
end

%% 导出 TXT，保持记录顺序
out = fopen(txtFile, 'wt');
assert(out >= 0, '无法创建 TXT 文件。');
txtCleanup = onCleanup(@() fclose(out));
row = 0;
for r = 1:recordCount
    channels = records{r};
    n = sampleCounts(r);
    fprintf(out, 'waveform\t[0]\t[1]\n');
    fprintf(out, 't0\t%s\t%s\n', ...
        char(timestamps(row + 1, 1)), char(timestamps(row + 1, 2)));
    fprintf(out, 'delta t\t%.8f\t%.8f\n\n', channels(1).dt, channels(2).dt);
    fprintf(out, 'time\tY[0]\tY[1]\n');
    for k = 1:n
        fprintf(out, '%s\t%.6E\t%.6E\n', ...
            char(timestamps(row + k, 1)), Y(row + k, 1), Y(row + k, 2));
    end
    fprintf(out, '\n');
    row = row + n;
end
clear txtCleanup

%% 输出解析结果
result = struct('records', {records}, 'timeSec', timeSec, 'Y', Y, ...
    'timestamps', timestamps, 'sampleCounts', sampleCounts, ...
    'csvFile', csvFile, 'txtFile', txtFile);
fprintf('记录数：%d；通道数：2；每通道采样数：%d\n', recordCount, sampleTotal);
fprintf('采样间隔：%.8f s；采样率：%.8g Hz\n', first(1).dt, 1/first(1).dt);
fprintf('通道：%s，%s；文件标注单位：%s\n', ...
    first(1).attributes.NI_ChannelName, first(2).attributes.NI_ChannelName, ...
    first(1).attributes.NI_UnitDescription);
if exportCsv
    fprintf('CSV：%s\n', csvFile);
end
fprintf('TXT：%s\n', txtFile);
end


function text = readString(fid)
n = fread(fid, 1, 'uint32');
assert(~isempty(n), '字符串长度字段不完整。');
bytes = fread(fid, n, '*uint8').';
assert(numel(bytes) == n, '字符串数据不完整。');
text = native2unicode(bytes, 'UTF-8');
end


function value = readIndex(fid)
value = fread(fid, 1, 'uint16');
if value >= 32768
    value = (value - 32768) * 65536 + fread(fid, 1, 'uint16');
end
end


function [value, attributes] = readVariant(fid)
variantVersion = fread(fid, 1, 'uint32');
assert(floor(variantVersion / 2^24) == 8, '不支持此属性版本。');
typeCount = fread(fid, 1, 'uint32');
types = zeros(typeCount, 1);
for i = 1:typeCount
    descriptorSize = fread(fid, 1, 'uint16');
    typeCode = fread(fid, 1, 'uint16');
    types(i) = mod(typeCode, 256);
    fseek(fid, descriptorSize - 4, 'cof');
end
usedCount = readIndex(fid);
assert(usedCount == 1, '当前属性应只包含一个根类型。');
typeIndex = readIndex(fid);
switch types(typeIndex + 1)
    case 0
        value = [];
    case 48
        value = readString(fid);
    otherwise
        error('尚未支持属性类型 0x%02X。', types(typeIndex + 1));
end
attributeCount = fread(fid, 1, 'uint32');
attributes = struct;
for i = 1:attributeCount
    name = readString(fid);
    [attributeValue, ~] = readVariant(fid);
    attributes.(matlab.lang.makeValidName(name)) = attributeValue;
end
end


function text = formatTimestamp(wholeSeconds, fraction, offset)
fractionalSeconds = double(fraction) / 2^64 + offset;
carry = floor(fractionalSeconds);
microseconds = floor((fractionalSeconds - carry) * 1e6);
epoch = datetime(1904, 1, 1, 'TimeZone', 'UTC');
clock = epoch + seconds(double(wholeSeconds)) + seconds(carry);
clock.TimeZone = 'Asia/Shanghai';
clock.Format = 'yyyy/M/d  HH:mm:ss';
text = string(clock) + "." + compose('%06d', microseconds);
end
