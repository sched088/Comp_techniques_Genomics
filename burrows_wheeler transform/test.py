def suffixArray(s):
    ''' Given T return suffix array SA(T).  Uses "sorted"
        function for simplicity, which is probably very slow. '''
    satups = sorted([(s[i:], i) for i in range(len(s))])
    # print(satups)
    # print("SA: " + str(list(map(lambda x: x[1], satups))))
    return list(map(lambda x: x[1], satups))  # extract, return just offsets


def bwtFromSa(t, sa=None):
    ''' Given T, returns BWT(T) by way of the suffix array. '''
    bw = []
    dollarRow = None
    if sa is None:
        sa = suffixArray(t)
    print("SA: " + str(sa))
    for si in sa:
        print(si)
        if si == 0:
            dollarRow = len(bw)
            bw.append('$')
        else:
            bw.append(t[si - 1])
    return ''.join(bw), dollarRow  # return string-ized version of list bw


class FmCheckpoints(object):
    ''' Manages rank checkpoints and handles rank queries, which are
        O(1) time, with the checkpoints taking O(m) space, where m is
        length of text. '''

    def __init__(self, bw, cpIval=4):
        ''' Scan BWT, creating periodic checkpoints as we go '''
        self.cps = {}  # checkpoints
        self.cpIval = cpIval  # spacing between checkpoints
        tally = {}  # tally so far
        # Create an entry in tally dictionary and checkpoint map for
        # each distinct character in text
        for c in bw:
            if c not in tally:
                tally[c] = 0
                self.cps[c] = []
        # Now build the checkpoints
        for i, c in enumerate(bw):
            tally[c] += 1  # up to *and including*
            if i % cpIval == 0:
                for c in tally.keys():
                    self.cps[c].append(tally[c])

    def rank(self, bw, c, row):
        # print("BW: " + str(bw))
        ''' Return # c's there are in bw up to and including row '''
        if row < 0 or c not in self.cps:
            return 0
        i, nocc = row, 0
        # Always walk to left (up) when calculating rank
        while (i % self.cpIval) != 0:
            if bw[i] == c:
                nocc += 1
            i -= 1

        # print("nocc: " + str(nocc))
        # print("return: " + str(self.cps[c][i // self.cpIval]))
        return self.cps[c][i // self.cpIval] + nocc


st = 'teststring'
#     0123456789
cps = FmCheckpoints(st)
# should get back list of integers, where elt i gives
# # times 't' appears up to and including offset i
[cps.rank(st, 't', i) for i in range(10)]
# likewise for 'g'
[cps.rank(st, 'g', i) for i in range(10)]


class FmIndex():
    ''' O(m) size FM Index, where checkpoints and suffix array samples are
        spaced O(1) elements apart.  Queries like count() and range() are
        O(n) where n is the length of the query.  Finding all k
        occurrences of a length-n query string takes O(n + k) time.

        Note: The spacings in the suffix array sample and checkpoints can
        be chosen differently to achieve different bounds. '''

    @staticmethod
    def downsampleSuffixArray(sa, n=4):
        ''' Take only the suffix-array entries for every nth suffix.  Keep
            suffixes at offsets 0, n, 2n, etc with respect to the text.
            Return map from the rows to their suffix-array values. '''
        ssa = {}
        for i, suf in enumerate(sa):
            # We could use i % n instead of sa[i] % n, but we lose the
            # constant-time guarantee for resolutions
            if suf % n == 0:
                ssa[i] = suf
        return ssa

    def __init__(self, t, cpIval=4, ssaIval=4):

        if t[-1] != '$':
            t += '$'  # add dollar if not there already

        # Get BWT string and offset of $ within it
        sa = suffixArray(t)
        # print("2SA: " + str(sa))

        self.bwt, self.dollarRow = bwtFromSa(t, sa)
        # print('L/BWT: ' + str(self.bwt))
        # Get downsampled suffix array, taking every 1 out of 'ssaIval'
        # elements w/r/t T
        self.ssa = self.downsampleSuffixArray(sa, ssaIval)
        # print("2SSA: " + str(self.ssa))
        self.slen = len(self.bwt)
        # Make rank checkpoints
        self.cps = FmCheckpoints(self.bwt, cpIval)
        # Calculate # occurrences of each character
        tots = dict()
        for c in self.bwt:
            tots[c] = tots.get(c, 0) + 1
        # Calculate concise representation of first column
        self.first = {}
        totc = 0
        for c, count in sorted(tots.items()):
            self.first[c] = totc
            totc += count

    def count(self, c):
        ''' Return number of occurrences of characters < c '''
        if c not in self.first:
            # (Unusual) case where c does not occur in text
            for cc in sorted(self.first.iterkeys()):
                if c < cc: return self.first[cc]
            return self.first[cc]
        else:
            # print("count: " + str(self.first[c]))
            # print(self.first)
            return self.first[c]

    def range(self, p):
        ''' Return range of BWM rows having p as a prefix '''
        l, r = 0, self.slen - 1  # closed (inclusive) interval
        for i in range(len(p) - 1, -1, -1):  # from right to left
            # print(p[i])

            l = self.cps.rank(self.bwt, p[i], l - 1) + self.count(p[i])
            r = self.cps.rank(self.bwt, p[i], r) + self.count(p[i]) - 1
            # print(l, r + 1)

            if r < l:
                break
        # print(l, r+1)

        return l, r + 1

    def resolve(self, row):
        ''' Given BWM row, return its offset w/r/t T '''
        # print("r: " + str(row))
        def stepLeft(row):
            ''' Step left according to character in given BWT row '''
            c = self.bwt[row]
            # print("C: " + str(c))
            return self.cps.rank(self.bwt, c, row - 1) + self.count(c)

        nsteps = 0
        # print("SSA: " + str(self.ssa))
        while row not in self.ssa:
            row = stepLeft(row)
            nsteps += 1
            # print("row: " + str(row))
           # print("SSA[row]: " + str(self.ssa[row]))
        return self.ssa[row] + nsteps

    def hasSubstring(self, p):
        ''' Return true if and only if p is substring of indexed text '''
        l, r = self.range(p)
        return r > l

    def hasSuffix(self, p):
        ''' Return true if and only if p is suffix of indexed text '''
        l, r = self.range(p)
        off = self.resolve(l)
        return r > l and off + len(p) == self.slen - 1

    def occurrences(self, p):
        ''' Return offsets for all occurrences of p, in no particular order '''
        l, r = self.range(p)
        return [self.resolve(x) for x in range(l, r)]

p, t = "CAT", "TTGTGTGCATGTTGTTTCATCATTTAGAGATACATTGCGCTGCATCATGGTCA"
#              01234567890123456789012345678901234567890123456789012
# Occurrences:        *         *  *           *         *  *
fm = FmIndex(t)
print(sorted(fm.occurrences(p)))
