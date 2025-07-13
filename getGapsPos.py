# python3 ../scripts/getFullChainSize.py hg38ToPanTro5.over.chain.gz hg38.panTro5.full_size_chains.bed panTro5.hg38.full_size_chains.bed hg38.panTro5.gaps_chains.csv
import sys
import gzip

gap_sizes = open(sys.argv[2], 'w')
block_sizes = open(sys.argv[3], 'w')
chainid = ""

for line in gzip.open(sys.argv[1], "rt"):
    chain = False
    total_size = 0
    if not line.startswith("#") and not line.startswith("\n"):
        if line.startswith("chain"):
            line = line.strip().split(" ") # first is "chain" and then as defined on webpage
            chainid = line[-1]
            chr = line[2]
            # startChain = int(line[5])
            # endChain = line[6]
            # strand = line[4]
            gapEnd = int(line[5]) # initialize gapEnd with the start coordinate of the chain

        else:
            line = line.strip().split("\t") # ungapped / gapRef / gapQuery
            if(len(line) > 1): # skip last line that has only ungapped alignement size for last block
                gapStart = gapEnd + int(line[0]) # start after ungapped block
                if gapStart - gapEnd > 0: # ungapped block between end of last gap and start of next gap
                    block_sizes.write("\t".join((chr, str(gapEnd), str(gapStart), chainid)) + "\n")
                gapEnd = gapStart + int(line[1]) # end after start + size
                if gapEnd - gapStart > 0:
                    gap_sizes.write(chr + "\t" + str(gapStart) + "\t" + str(gapEnd) + "\t" + chainid + "\n")
                