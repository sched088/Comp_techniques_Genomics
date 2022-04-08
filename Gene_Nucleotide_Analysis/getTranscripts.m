function transcripts = getTranscripts(annotationFileName, chr, x, y)
% #########################################################################
% Input:  annotationFileName - the human gene annotation file name
%         x - the start position in the human genome sequence
%         y - the end position in the long sequence
%         sequence(x:y) is a sample of the whole sequence
%         chr - chromosome name, e.g. 'chr1', 'chr2', ..., 'chr22', 'chrX',
%         'chrY'
% Output: a cell array transcripts that contains an array of forward
%         strands with start and end positions, and an array of reverse 
%         strands with start and end positions.
% #########################################################################

chromosome = ['''', chr, ''''];
% an array of (+) strands
ptranscripts = [];
% an array of (-) strands
ntranscripts = [];

fid = fopen(annotationFileName);
s = textscan(fid, '%s %s %s %s %s %s');
n = size(s{1});
for i = 1:n
    if strcmp(s{5}{i}, chromosome)
        startPosition = str2num(s{3}{i});
        endPosition = str2num(s{4}{i});        
        if (startPosition >= x && endPosition <= y)
            if strcmp(s{6}{i}, s{6}{2})
                if ~ismember([startPosition,endPosition],ptranscripts)
                    ptranscripts = [ptranscripts; [startPosition, endPosition]];
                end
            else
                if ~ismember([startPosition,endPosition],ntranscripts)
                    ntranscripts = [ntranscripts; [startPosition, endPosition]];
                end
            end
        end
    end
end

fclose(fid);
transcripts = {ptranscripts, ntranscripts, x};
end
