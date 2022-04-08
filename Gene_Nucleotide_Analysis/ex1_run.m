function ex1_run(FASTAfilename, annotationFilename, chr, x, y)
% #########################################################################
% EX1_RUN plots the locations of the genes together with the GC content
% over a sample from position x to y (inclusive) of the long sequence.
% Input:    FASTAfilename - the file name of the chromosome
%           annotationFilename - the human gene annotation file name
%           chr - chromosome name, e.g. 'chr1', 'chr2', ..., 'chr22', 
%           'chrX', 'chrY'
%           x - the start position in the human genome sequence
%           y - the end position in the long sequence
%           sequence(x:y) is a sample of the whole sequence
%           Condition: 0 < x < y <= sequence_length
% 
% Note that you should try sampling different positions.
% #########################################################################

    % When the sequence is too long, we cannot use method fastaread to read
    % the file. 
    % We have to read every blocksize of the fasta file.
    fidIn = fopen(FASTAfilename, 'r');
    header = fgetl(fidIn);
    [fullPath, filename, extension] = fileparts(FASTAfilename);
    mmFilename = [filename '.mm'];
    fidOut = fopen(mmFilename, 'w');
    newLine = sprintf('\n');
    blocksize = 2^20;
    while ~feof(fidIn)
        seq = fread(fidIn, blocksize, '*char')';
        seq = strrep(seq, newLine, '');
        seq = nt2int(seq);
        fwrite(fidOut, seq, 'uint8');
    end
    fclose(fidIn);
    fclose(fidOut);
    chromosome = memmapfile(mmFilename, 'format', 'uint8');

    % Take a sample from position x to position y from the long sequence
    sample = int2nt(chromosome.Data(x:y)');
    transcripts = getTranscripts(annotationFilename, chr, x, y);
    % Plot the sample
    % Here, the window length is set = 5000. You can try different
    % window lengths.
    ntdensity2(sample,transcripts,'window', 5000);
    clear chromosome;
    delete(mmFilename);
end
