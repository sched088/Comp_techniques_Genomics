function [density, highCG] = ntdensity2(seq, varargin)
%NTDENSITY2 plots the density of nucleotides in a sequence. If the sequence
%is longer than 500 Kb the output signal is decimated.
%
%   NTDENSITY2(SEQ) plots the density of nucleotides A,T,C,G in sequence
%   SEQ.
%
%   DENSITY = NTDENSITY2(SEQ) returns a structure of the density of
%   nucleotides A, C, G, and T.
%
%   NTDENSITY2(SEQ, TRANSCRIPTS) plot the locations of the genes together 
%   with the GC content
%
%   NTDENSITY2(...,'WINDOW',L) uses a window of length L for the density
%   calculation. The window length must be an odd integer >= 5 and has a
%   default value of length(SEQ)/20.
%
%   [DENSITY, HIGHCG] = NTDENSITY2(...,'CGTHRESHOLD',CGT) returns indices
%   for regions where the CG content of SEQ is greater than CGT. The
%   default value for CGT is 0.5.
%
%   Example:
%
%       % Create a random sequence and analyze the nucleotide density.
%       seq = randseq(240)
%       ntdensity2(seq)
%
%   See also BASECOUNT, CODONCOUNT, CPGISLAND, DIMERCOUNT, FILTER.

%   Copyright 2002-2005 The MathWorks, Inc.
%   $Revision: 1.5.6.6 $  $Date: 2005/06/09 21:56:53 $

% If the input is a structure then extract the Sequence data.
if isstruct(seq)
    try
        seq = seqfromstruct(seq);
    catch
        rethrow(lasterror);
    end
end
window = ceil(length(seq)/20);
cgthreshold = 0.5;
transcripts = {};
if nargin > 1
    transcripts = varargin{1};
    if nargin > 2
        if rem(nargin - 2,3) == 0
            error('Bioinfo:IncorrectNumberOfArguments',...
                'Incorrect number of arguments to %s.',mfilename);
        end
        okargs = {'window','cgthreshold'};
        for j=2:2:nargin-2
            pname = varargin{j};
            pval = varargin{j+1};
            k = strmatch(lower(pname), okargs); %#ok
            if isempty(k)
                error('Bioinfo:UnknownParameterName',...
                    'Unknown parameter name: %s.',pname);
            elseif length(k)>1
                error('Bioinfo:AmbiguousParameterName',...
                    'Ambiguous parameter name: %s.',pname);
            else
                switch(k)
                    case 1   % window
                        window = pval;
                    case 2  % gcthreshold
                        cgthreshold = pval;
                end
            end
        end
    end
end

if window < 5 % window has a minimum length of 5
    window = 5;
%elseif rem(window,2) == 0 % window must be odd
%    window = window + 1;
end

if ischar(seq)
    seq = nt2int(seq);
end

% Pad the end of the sequence to accomodate the window size
%seq = [seq zeros(1,window-1)];
len = length(seq);
sa = zeros(1,len);
sc = zeros(1,len);
sg = zeros(1,len);
st = zeros(1,len);

feqsa = zeros(1,len);
feqsc = zeros(1,len);
feqsg = zeros(1,len);
feqst = zeros(1,len);
IX = find(seq==1);
sa(IX) = 1;
IX = find(seq==2);
sc(IX) = 1;
IX = find(seq==3);
sg(IX) = 1;
IX = find(seq==4);
st(IX) = 1;
clear IX
% for i = 1:window
%     feqsa(i) = sum(sa(1:i))/i;
%     feqsc(i) = sum(sc(1:i))/i;
%     feqsg(i) = sum(sg(1:i))/i;
%     feqst(i) = sum(st(1:i))/i;
% end
feqsa(1:window) = cumsum(sa(1:window))./[1:window];
feqsc(1:window) = cumsum(sc(1:window))./[1:window];
feqsg(1:window) = cumsum(sg(1:window))./[1:window];
feqst(1:window) = cumsum(st(1:window))./[1:window];

sa = sa/window;
sc = sc/window;
sg = sg/window;
st = st/window;

for i = window+1:len
    feqsa(i) = feqsa(i-1)+sa(i)-sa(i-window);
    feqsc(i) = feqsc(i-1)+sc(i)-sc(i-window);
    feqsg(i) = feqsg(i-1)+sg(i)-sg(i-window);
    feqst(i) = feqst(i-1)+st(i)-st(i-window);
end

sa = feqsa;
sc = feqsc;
sg = feqsg;
st = feqst;

clear feqsa feqsc feqsg feqst;
sa = sa(50:end);
sc = sc(50:end);
sg = sg(50:end);
st = st(50:end);

if length(sa)>500000
    r=round(length(sa)/100000);
    sa=decimate(sa,r);
    sc=decimate(sc,r);
    sg=decimate(sg,r);
    st=decimate(st,r);
end

wsize=size(seq);
csize=size([sa' sc' sg' st']);
a=wsize(2)/csize(1);
if isempty(transcripts)
    p = 1;
else
    p = transcripts{3};
end

X=[p:a:(wsize(2)+p-1)];

% create two plots. One for all density one for AT and CG
subplot(2,1,1);
plot([X' X' X' X'],[sa' sc' sg' st']);
legend({'A','C','G','T'});
title('Nucleotide density')
subplot(2,1,2);
hold on
%xlim([0 10^8]);
ylim([0 1.5]);
plot([X' X'],[(sa+st)' (sc+sg)']);
if ~isempty(transcripts)
    plot(transcripts{1}', ones(size(transcripts{1}))'*1.3, '->');
    plot(transcripts{2}', ones(size(transcripts{2}))'*1.2, '-<');
end
hold off
legend({'A-T','C-G'});
title('A-T C-G density')

% Create outputs
if nargout > 0
    density.A = sa;
    density.C = sc;
    density.G = sg;
    density.T = st;
end

% with two outputs we output GC density and show the line threshold on the
% graph.

if nargout > 1
    highCG = find(sc+sg > cgthreshold);
    hold on;
    plot([1,length(seq)],[cgthreshold,cgthreshold],'k','linewidth',0.3,'linestyle','-.');
    hold off;
end



% Calculate nucleotide densities with a moving window
