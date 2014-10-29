ECS
===

Evolutionary conserved RNA structure prediction

 version 0.9beta
	 __  __          ______    _____  _____          _   _
	|  \/  |   /\   |  ____|  / ____|/ ____|   /\   | \ | |
	| \  / |  /  \  | |__    | (___ | |       /  \  |  \| |
	| |\/| | / /\ \ |  __|    \___ \| |      / /\ \ | . ` |
	| |  | |/ ____ \| |       ____) | |____ / ____ \| |\  |
	|_|  |_/_/    \_\_|      |_____/ \_____/_/    \_\_| \_|
	 SCAN MULTIPLE ALIGNMENTS FOR CONSERVED RNA STRUCTURES

Reads a set of maf files, realigns them with mafft, breaks them into windows,
calculates stats, scans with SISSIz or RNAz, outputs bed coordonates of high-confidence predictions
*** Known issues ***
-Using MAFFT realignment causes certain jobs to hang when finished; this is likely an executor service problem

Usage:     java -jar MafScanCcr.jar [options] -o output/directory -i input.maf (last parameter must be -i)
Options:
  -bs int       Block Size for splitting large MAF blocks (default 5000)
  -c  int       number of CPUs for calculations (default 4)
  -g  int       max gap percentage of sequences for 2D prediction (default 50)
  -mafft        Realign with mafft-ginsi (slower)
  -ml int       Max Length of MAF block for MAFFT realignment (default 10000)
  -s  int       step size (default 100)
  -v            verbose (messy but detailed) output
  -w  int       window size (default 200)
