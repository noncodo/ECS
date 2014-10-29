import java.util.*; import java.io.*; import java.math.* ; import java.lang.*;; 

public class ScanItFast implements Runnable {
	static boolean VERBOSE = false, FILTER_ID = true ;
	private boolean [] hasChars , keepMe, isAmbiguous, isNotUnique ;
	private String[][] MafTab ;
	private char[][] AlnTab;
	private static String Path, SSZBINARY, RNZBINARY;
	private int [] coordTab;
	private int goodSeqs, iterate, random, step, STEP, WINDOW, GAPS, SSZ_THRESHOLD, retainedColumns;
	private BufferedWriter WriteALN ;
	private BufferedReader ReadFile;
	private String [] OutAln, OutAlnRC, FilteredTab,	NameTab,  TempTab = new String [1] ;
	private String Line = "" ;
	private File Aln, AlnRC ;
	private double [] stats, chars, totalChars ;
	private double	mpi = 0, mpi2 = 0, var = 0, shanon = 0, uniqueComps = 0, uniqueSeqs , outCols ;
	private double [][]	pids, gaps ;
	
	ScanItFast(int step, String[][] MafTab, char[][] AlnTab, String Path, 
			 int [] coordTab , int WINDOW, int STEP, int GAPS, 
			 int SSZ_THRESHOLD, String SSZBINARY, String RNZBINARY, boolean VERBOSE, boolean FILTER_ID ) { 
		this.step = step ; 
		this.MafTab = MafTab ; 
		this.AlnTab = AlnTab ; 
		this.Path = Path ; 
		this.coordTab = coordTab ;
		this.WINDOW = WINDOW ; 
		this.STEP = STEP ; 
		this.GAPS = GAPS ; 
		this.SSZ_THRESHOLD = SSZ_THRESHOLD ;
		this.RNZBINARY = RNZBINARY ;
		this.SSZBINARY = SSZBINARY ;
		this.VERBOSE = VERBOSE ; 
		this.VERBOSE = FILTER_ID; // not implemented yet. 
	}
	
	public void run() {
        if (VERBOSE)
			System.out.println("- - -> Starting Scan") ;
		if (VERBOSE && AlnTab.length != MafTab.length ) {
			System.out.println(" #### Maf and AlnTab aren't same length" ) ;
			System.out.println( MafTab.length +" "+ AlnTab.length ) ;
		}
		isAmbiguous = new boolean [ AlnTab.length ];
		int startPos = Integer.parseInt(MafTab[0][2]);
		
		// remove identical rows or those with too many gaps & N's
		Set<String>	Uniques = new LinkedHashSet<String>(),
        UniquesWithGaps = new LinkedHashSet<String>(),
        UniqueNames = new LinkedHashSet<String>() ;
		for (int seq = 0 ; seq != AlnTab.length ; seq++ ) {
			String DegappedWindow = new String( AlnTab[seq], step, Math.min( WINDOW, AlnTab[seq].length - step )).replaceAll("[^ATCGUatcgu]", "" ).toUpperCase() ;
			// only retains non-identical unaligned sequences with at least one character
			if ( DegappedWindow.length() > 0 && Uniques.add( DegappedWindow ) ) {
				UniquesWithGaps.add(  new String( AlnTab[seq], step, Math.min( WINDOW, AlnTab[seq].length - step )).toUpperCase() ) ;
				UniqueNames.add( MafTab[seq][ 1 ] ) ;
			}
		}
		FilteredTab = new String [ UniquesWithGaps.size()]; // creating an Array from Hash
		FilteredTab = UniquesWithGaps.toArray(FilteredTab);
		NameTab = new String [ UniquesWithGaps.size() ] ;
		NameTab = UniqueNames.toArray( NameTab ) ;
		
		// first check for > 2 seqs
		goodSeqs = UniquesWithGaps.size()  ;
		if ( goodSeqs <= 3 ) {
			if (VERBOSE)
				System.out.println("-> Not Enough seqs in this window!") ;
			return;
		}
        
		// remove gappy sequences
		if (VERBOSE)
			System.out.println("- -> Gappy sequences") ;
		keepMe = new boolean [ FilteredTab.length ] ;
		for (int seq = 0 ; seq != FilteredTab.length ; seq++ ) {
			if ( FilteredTab[ seq ].replaceAll("[^ATCGUatcgu]", "" ).length() >= (int)((double)WINDOW*((double)GAPS/100))) {
				keepMe[seq] = true ;
			}
			else {
				keepMe[seq] = false ;
				goodSeqs-- ;
				//if (VERBOSE)
				//	System.out.println("  --> removed a GAPpy sequence form the alignment" ) ;
			}
		}
		if ( goodSeqs <= 3 ) {
			if (VERBOSE)
				System.out.println("-> Not Enough seqs in this window!") ;
			return;
		}
		// exit when human is shit
		if ( !keepMe[ 0 ] )
			return;

		// check for gap-only columns
		if (VERBOSE)
			System.out.println("- -> Gap only columns") ;
		retainedColumns = WINDOW;
		hasChars = new boolean [ FilteredTab[0].length() ];
        gapScan: for (int col = 0 ; col != FilteredTab[0].length() ; col++ ) {
            for (int seq = 0 ; seq != FilteredTab.length ; seq++ ) {
                if ( keepMe[ seq ] ) {
                    if ( FilteredTab[ seq ].charAt( col ) == 'A' || FilteredTab[ seq ].charAt( col ) == 'C'
                        || FilteredTab[ seq ].charAt( col ) == 'T' || FilteredTab[ seq ].charAt( col ) == 'G' ) {
                        hasChars[ col ] = true ;
                        continue gapScan ;
                    }
                }
            }
            if ( !hasChars[ col ] )
                retainedColumns-- ;
        //if ( !hasChars[ col ] && VERBOSE)
        //	System.out.println( "-> empty col!" );  
        }

		// prepare clustalw file
		if (VERBOSE)
			System.out.println("- -> preparing Clustal format") ; 
		OutAln = new String[ goodSeqs ];
		OutAlnRC = new String[ goodSeqs ]  ;  
        iterate = 0 ;
		for (int seq = 0 ; seq != FilteredTab.length  ; seq++ ) { //removed x < goodseqs
			if ( keepMe[ seq ] ) { 
				OutAln[ iterate ] = NameTab[ seq ].substring( 0, Math.min( NameTab[ seq ].length(), 20)) ;
				for (int i = 0 ; i != 25- Math.min(  NameTab[ seq ].length() , 20) ; i++ )
					OutAln[ iterate ] = OutAln[ iterate ] + " ";
				for (int i = 0 ; i != FilteredTab[ 0 ].length() ; i++ )
					if ( hasChars[ i ])
						OutAln[ iterate ] = OutAln[ iterate ] + FilteredTab[ seq ].charAt( i ) ;
				OutAln[ iterate ] = OutAln[ iterate ] + "\n" ;
				iterate++ ;
			}
		}
		
        //*********************************************************************
		//					calculate stats						*
		//*********************************************************************
		if (VERBOSE)
			System.out.println("- - -> calculating statistics") ;
		uniqueSeqs = goodSeqs;
		outCols = OutAln[0].length()-25 ; //change last variable if CLUSTAL properties changes
		stats = new double [6];
		chars = new double [5];
		totalChars = new double [5];
		pids = new double [ goodSeqs ][ goodSeqs ];
		gaps = new double [ goodSeqs ][ goodSeqs ]; // gaps (and potentially mismatches)
		isNotUnique = new boolean [ goodSeqs ] ;
		// calculate id matrix and mpi
		for ( int k = 25 ; k != OutAln[0].length()-1 ; k ++ ) { // -1 avoids line break char, 25 is clustal seq ID length
			lines :for ( int i = 0 ; i != goodSeqs ; i++ ) {
				// initiate gaps[] and pids[] to 0 ???????????????????????
				if ( isNotUnique[ i ] )
					continue lines ;
				for ( int j = i+1 ; j != goodSeqs ; j++ ) {
					try {
                        
                        if ( isNotUnique[ j ] )
                            continue;
                        else if ( OutAln[ i ].charAt( k ) == OutAln[ j ].charAt( k ) ) {
							// this DP matrix makes shit easy!
                            if ( OutAln[ i ].charAt( k ) == 'A' || OutAln[ i ].charAt( k ) == 'T'
                                || OutAln[ i ].charAt( k ) == 'C' || OutAln[ i ].charAt( k ) == 'G'
                                || OutAln[ i ].charAt( k ) == 'U') { // U is just in case alignments are RNA
                                pids[ i ][ j ]++ ;
                                pids[ j ][ i ]++ ;
                            }
                        }
                        // this ignores "-:-"
                        else if ( OutAln[ i ].charAt( k ) != OutAln[ j ].charAt( k ) ) {
                            pids[ j ][ i ]++ ; // mismatch
                            if ( OutAln[ i ].charAt( k ) == '-' || OutAln[ i ].charAt( k ) == '.')
                                gaps[i][j]++; // gap
                            if ( OutAln[ j ].charAt( k ) == '-' || OutAln[ j ].charAt( k ) == '.')
                                gaps[j][i]++; // gap
                        }
                    }catch (Exception E) {
                        E.printStackTrace();
                        System.err.println( "Caught Exception!\n");
                        System.err.println( i +" " +j +" " + k + " " + OutAln.length + " " + goodSeqs
                                           + " FilteredTab[i]=" + FilteredTab[i].length() + " FilteredTab[j]=" + FilteredTab[j].length());
                        System.err.println( "OutAln[i]=" + OutAln[i].length() + " OutAln[j]=" + OutAln[j].length()
                                           + "\n" + OutAln[i] + "\n" + OutAln[j]) ;
                    }
					// keep unique seqs ignoring gaps
					if ( k == OutAln[0].length()-2 ) {
						if (	pids[ j ][ i ]  - gaps[ i ][ j ] ==  pids[ i ][ j ]  ||
						    pids[ j ][ i ]  - gaps[ j ][ j ] ==  pids[ i ][ j ]  ) {
							//both sequences are identical without gaps
							//keep the longer one
							if ( gaps[ i ][ j ] > gaps[ j ][ i ] )
								isNotUnique[ i ] = true ;
							else 
								isNotUnique[ j ] = true ; // this should also consider identical seqs
						}
						else {
							uniqueComps++ ; 
							// old mean pairwise identity ( considers gaps ) 
							mpi = mpi + 100 * pids[ i ][ j ] / pids[ j ][ i ] ;  
							// classical average identity 
							mpi2 = mpi2 + 100 * pids[ i ][ j ] / Math.min( OutAln[ i ].replaceAll("[^ATCGU]","").length(), 
                                                                          OutAln[ j ].replaceAll("[^ATCGU]","").length() ) ;
						}
					}
				}
			}
		}
		// calculate gaps, GC, shanon entropy, and Reverse Complement
		for ( int k = 25 ; k != OutAln[0].length() ; k ++ ) { 
			chars = new double [5] ; 
			for ( int i = 0 ; i != goodSeqs ; i++ ) {
				if ( isNotUnique[ i ] ) { 
					if ( k == OutAln[0].length()-2 ) 
						uniqueSeqs-- ; 
					continue ; 
				}
				switch ( OutAln[i].charAt( k ) ) { 
					case 'A':
						chars[0]++ ;
						totalChars[0]++ ;
						OutAlnRC[ i ] = (k == 25 )? "T" : "T" + OutAlnRC[ i ] ; 
						break;
					case 'U':
						chars[1]++ ;
						totalChars[1]++ ;
						OutAlnRC[ i ] = (k == 25 )? "A" : "A" + OutAlnRC[ i ] ; 
						break;
					case 'T':
						chars[1]++ ;
						totalChars[1]++ ;
						OutAlnRC[ i ] = (k == 25 )? "A" : "A" + OutAlnRC[ i ] ; 
						break;
					case 'C':
						chars[2]++ ;
						totalChars[2]++ ;
						OutAlnRC[ i ] = (k == 25 )? "G" : "G" + OutAlnRC[ i ] ; 
						break;
					case 'G':
						chars[3]++ ;
						totalChars[3]++ ;
						OutAlnRC[ i ] = (k == 25 )? "C" : "C" + OutAlnRC[ i ] ; 
						break;
					case '\n':
						OutAlnRC[ i ] = OutAlnRC[ i ] + '\n' ;
						break;
					case 'N':
						chars[4]++ ;
						totalChars[4]++ ;
						OutAlnRC[ i ] = (k == 25 )? "N" : "N" + OutAlnRC[ i ] ; 
						break;
					default:
						chars[4]++ ; 
						totalChars[4]++ ;
						OutAlnRC[ i ] = (k == 25 )? "-" : "-" + OutAlnRC[ i ] ;
						break;
				}
			}
			for (int z = 0 ; z != 5 ; z++ ) 
				shanon = ( chars[z] == 0 )? shanon + 0 : shanon +  chars[z]/uniqueSeqs * ( Math.log( chars[z]/uniqueSeqs ) / Math.log( 2 )); 
		}
		//System.out.println( uniqueSeqs +"\t"+goodSeqs+"\t"+outCols+"\t"+totalChars[4]+"\t"+( outCols * goodSeqs));
		stats[0] = mpi / uniqueComps;																		// Mean Pairwise ID 
		stats[5] = mpi2 / uniqueComps;																    // classical MPI
		for (int seq1 = 0 ; seq1 != goodSeqs ; seq1++ ) 
			for (int seq2 = seq1 +1 ; seq2 != goodSeqs ; seq2++ ) 
				if ( !isNotUnique[ seq1 ] && !isNotUnique[ seq2 ] ) 
					var = var + (double) Math.pow( (100*pids[ seq1 ][ seq2 ]/ pids[ seq2 ][ seq1 ]) - stats[0] , 2) ; 
		stats[1] = var / uniqueComps ;																// Variance
		stats[2] = -1 * shanon / ((double)outCols) ;													    // Normalized Shanon entropy
		stats[3] = 100*(totalChars[2]+totalChars[3])/(totalChars[0]+totalChars[1]+totalChars[2]+totalChars[3]) ;	   // GC content
		stats[4] = 100 * totalChars[4] / (outCols * uniqueSeqs) ;										  // GAP content
		//System.out.println( stats[0]+"\t"+(Math.sqrt(stats[1]))+"\t"+stats[2]+"\t"+stats[3]+"\t"+stats[4]) ; 	 // print stats

		// save BED coords from MAF file
		if (VERBOSE) 
			System.out.println("- -> Calculating BED coords ") ; 
		String	BedFile = MafTab[0][1].substring( MafTab[0][1].lastIndexOf(".")+1) +"\t";
		int humanLength = WINDOW ; 
		if ( MafTab[0][4].equals("+") ){
			BedFile = BedFile + (startPos + coordTab[ step ] )+"\t"+(startPos + coordTab[step + WINDOW -1])+"\t" ;
			humanLength = (startPos + coordTab[step + WINDOW -1] -1) - (startPos + coordTab[ step ])+1 ;
		}
		else { // this should only occur in user specified cases
			BedFile = BedFile + (Integer.parseInt(MafTab[0][5]) - (startPos + coordTab[ step  + WINDOW ] )  )
			+"\t"+ (Integer.parseInt(MafTab[0][5]) - (startPos + coordTab[ step ] ) +1 )+"\t";
			
			humanLength =  (Integer.parseInt(MafTab[0][5]) - (startPos + coordTab[ step ] ) +1 ) 
			- (Integer.parseInt(MafTab[0][5]) - (startPos + coordTab[ step  + WINDOW ] )  )+1 ; 
		}		
		BedFile = BedFile	+ (int)uniqueSeqs+":"+((double)(int)(10*stats[0])/10)+":"      // MPI
			+ ((double)(int)(10*stats[5])/10)+":"					  // CLASSIC MPI
			+ ((double)(int)(10*stats[4])/10) +":"					 // GAPS 
			+ ((double)(int)(10*Math.sqrt(stats[1]))/10) +":"			// STDEV 
			+ ((double)(int)(100*stats[2])/100)  +":"			    // SHANON
			+ ((double)(int)(10*stats[3])/10) ;				   //      GC	
		if (VERBOSE) 
			System.out.println( "Pre SISSIz bed file: \n"+" "+BedFile ) ; 
		// select appropriate alignments
		for (int j = 0 ; j != WINDOW ; j++ ) 
			if ( hasChars[ j ] ) 
				outCols++ ; 
		int random = (int)((double)10000*Math.random()) ; 
		File Aln = new File( Path+"/"+BedFile.replaceAll("\t","_")+".aln."+ random ),	// 
			AlnRC = new File( Path+"/"+BedFile.replaceAll("\t","_")+"rc.aln." + random );  //  
		// v v v v v v v v    INCLUSION STATS     v v v v v v v v v v v v v
		if ( outCols > WINDOW / 2 && uniqueSeqs > 2 && humanLength >= Math.min( STEP, WINDOW/4 ) && stats[4] <= 75 && stats[0] > 0 ) {    
			// Write Sequences to ALN Format  
			try {
				BufferedWriter WriteClustal = new BufferedWriter(new FileWriter( Aln )),
				WriteClustalRC = new BufferedWriter(new FileWriter( AlnRC ));
				WriteClustal.write("CLUSTAL \n\n") ; 
				WriteClustalRC.write("CLUSTAL \n\n") ; 
				for ( int y = 0 ; y != goodSeqs ; y++ ) {
					if ( !isNotUnique[ y ] ) {
						WriteClustal.write( OutAln[ y ] ) ; 
						OutAlnRC[ y ] = OutAln[ y ].substring(0,25)+ OutAlnRC[ y ]  ;
						WriteClustalRC.write( OutAlnRC[ y ] ) ; 
					}
				}
				WriteClustal.close() ; 
				WriteClustalRC.close() ;
			} catch (IOException Fuck) {
				if (VERBOSE) 
					System.err.println("Arrgh... Couldn't write clustal file!");
				Fuck.printStackTrace();
				Aln.delete() ;
				AlnRC.delete() ;
				return; 
			}
		}
		else {
			if (VERBOSE) {
				System.out.println("---> rejected alignment" ) ; 
 	 	 		System.out.println("     outcols = "+outCols +"\tuniqueseqs = "+uniqueSeqs+"\thuman_len = "+humanLength+
							    "\tGAPS = "+stats[4]+"\n    PID = "+stats[0]);
				if ( stats[0] < 5 )
					System.out.println("-----> SUPER LOW PID"); 
			}
//			Aln.delete() ;
//			AlnRC.delete() ;
			return ; 
		}	
		String	FinalBedFile = "", 
				FinalBedFileRC = "", 
		Antisense = (MafTab[0][4].equals("+"))? "-" : "+" ; 
		//***************** 	SISSIz scan & parse		******************
		if ( stats[0] <= 85 ) {
			String [] SissizOutTab = new String [ 12 ] ; 
			try {
				SissizOutTab  = ScanSSZ( Path, BedFile , stats, 1, random) ;
				if ( SissizOutTab == null ) { // timeout
                  Aln.delete() ;
				}
			} catch (IOException Fuck) {						  
				Fuck.printStackTrace();
				System.err.println("ScanSSZ failed with ");
				for ( int y = 0 ; y != goodSeqs ; y++ ) {
					if ( !isNotUnique[ y ] ) {
						System.err.println( OutAln[ y ] ) ;
					}
				}
				Aln.delete() ;
			}
			// delete empty alignments
			if (  SissizOutTab == null || SissizOutTab[ 10 ] == null ) {
				Aln.delete() ;
			}
			else {
				// delete low scoring alignments
				if ( (SissizOutTab[ 1 ].equals("r") && Double.parseDouble( SissizOutTab[ 10 ] ) > -2.2)
                    || (SissizOutTab[ 1 ].equals("s") && Double.parseDouble( SissizOutTab[ 10 ] ) > -2.7) ){
					Aln.delete() ;
				}
				else {
                    FinalBedFile = BedFile +":"+SissizOutTab[ 1 ]+"_"+ (int)(Double.parseDouble(SissizOutTab[ 10 ])*-100)+"_"+ MafTab[0][4];
                    //write bed and rename alignment
                    System.out.println( FinalBedFile.replaceAll("_","\t")) ;
					File NewFile = new File( Path+"/"+FinalBedFile.replaceAll("\t","_")+".aln" ) ;
					int file_count = 0 ;
					while ( NewFile.exists() ) {
						file_count++;
						NewFile = new File( Path+"/"+FinalBedFile.replaceAll("\t","_")+".aln_"+ file_count) ;
					}
					boolean result = Aln.renameTo( NewFile );
				}
			}
			// * * * * * *  now for the RC  * * * * * * 
			try {
				SissizOutTab  = ScanSSZ( Path, BedFile+"rc" , stats, 1, random) ;
				if ( SissizOutTab == null ) { //
					AlnRC.delete() ;
				}
			} catch (IOException Fuck) {						  
				Fuck.printStackTrace();
				System.err.println("ScanSSZ failed in RC with ");
				for ( int y = 0 ; y != goodSeqs ; y++ ) {
					if ( !isNotUnique[ y ] ) {
						System.err.println( OutAln[ y ] ) ;
					}
				}
				AlnRC.delete() ;
			}
			if ( SissizOutTab == null || SissizOutTab[ 10 ] == null )  {
				AlnRC.delete() ; 
			}
			else {
				// delete low scoring alignments
				if ( (SissizOutTab[ 1 ].equals("r") && Double.parseDouble( SissizOutTab[ 10 ] ) > -2.2)
                    || (SissizOutTab[ 1 ].equals("s") && Double.parseDouble( SissizOutTab[ 10 ] ) > -2.7) ){
					AlnRC.delete() ;
				}
				else {
                    FinalBedFileRC = BedFile+":"+SissizOutTab[ 1 ]+"_"+ (int)(Double.parseDouble(SissizOutTab[ 10 ])*-100) +"_"+ Antisense ;
                    //write bedRC and rename alignment
                    System.out.println( FinalBedFileRC.replaceAll("_","\t")) ;
					File NewFile = new File( Path+"/"+FinalBedFileRC.replaceAll("\t","_")+".aln" ) ;
					int file_count = 0 ;
					while ( NewFile.exists() ) {
						file_count++;
						NewFile = new File( Path+"/"+FinalBedFileRC.replaceAll("\t","_")+".aln_"+ file_count) ;
					}
					boolean result = AlnRC.renameTo( NewFile );
				}
				return;
			}			
		}
		else {
			//***************** 	RNAz scan & parse		******************
			String [] RnazOutTab = new String [ 5 ] ; 
			try {
				RnazOutTab  = ScanRNZ( Path, BedFile , stats[ 0 ], 0, random) ;  
				if ( RnazOutTab == null ) { 
					Aln.delete() ;
				}
			} catch (IOException Fuck) {						  
				Fuck.printStackTrace();
				System.err.println("ScanRNAz failed with ");
				for ( int y = 0 ; y != goodSeqs ; y++ ) {
					if ( !isNotUnique[ y ] ) { // make this verbose? 
						System.err.println( OutAln[ y ] ) ;
					}
				}
				Aln.delete() ;
			}
			// delete when no RNAz output
			if ( RnazOutTab == null || RnazOutTab[4] == null ) {
				Aln.delete() ; 
			}
			//write bed and rename alignment
			else {
				// delete low scoring alignments
				if ( Double.parseDouble( RnazOutTab[ 4 ] ) < 0.32 )  {
					Aln.delete() ;
				}
				else {
                    FinalBedFile = BedFile +":z_"+ (int)(Double.parseDouble(RnazOutTab[4])*100) +"_"+ MafTab[0][4];
                    System.out.println( FinalBedFile.replaceAll("_","\t")) ;
                    File NewFile = new File( Path+"/"+FinalBedFile.replaceAll("\t","_")+".aln" ) ;
					int file_count = 0 ;
					while ( NewFile.exists() ) {
						file_count++;
						NewFile = new File( Path+"/"+FinalBedFile.replaceAll("\t","_")+".aln_"+ file_count) ;
					}
					boolean result = Aln.renameTo( NewFile );
				}
			}		
			// * * * * * *  now for the RC  * * * * * * 
			try {
				RnazOutTab  = ScanRNZ( Path, BedFile+"rc" , stats[ 0 ], 0, random) ;  
				if ( RnazOutTab == null ) { // timeout
					AlnRC.delete() ;
				}
			} catch (IOException Fuck) {						  
				Fuck.printStackTrace();
				System.err.println("ScanRNAz failed in RC with ");
				for ( int y = 0 ; y != goodSeqs ; y++ ) {
					if ( !isNotUnique[ y ] ) { // make this verbose? 
						System.err.println( OutAlnRC[ y ] ) ;
					}
				}
				AlnRC.delete() ;
			}
			// delete when no RNAz output
			if ( RnazOutTab == null || RnazOutTab[4] == null ) {
				AlnRC.delete() ; 
			}
			//write bed and rename alignment
			else {
				// delete low scoring alignments				
				if ( Double.parseDouble( RnazOutTab[ 4 ] ) < 0.32 )  {
					AlnRC.delete() ; 
				}
				else {
                    FinalBedFileRC = BedFile +":z_"+ (int)(Double.parseDouble(RnazOutTab[4])*100) +"_"+ Antisense ;
                    System.out.println( FinalBedFileRC.replaceAll("_","\t")) ;
					File NewFile = new File( Path+"/"+FinalBedFileRC.replaceAll("\t","_")+".aln" ) ;
					int file_count = 0 ;
					while ( NewFile.exists() ) {
						file_count++;
						NewFile = new File( Path+"/"+FinalBedFileRC.replaceAll("\t","_")+".aln_"+ file_count) ;
					}
					boolean result = AlnRC.renameTo( NewFile );
				}
				return ;
			}		
		}
		Aln.delete() ;
		AlnRC.delete() ;
	}
	//*********************************************************************
	//					SISSIz scan & parse						*
	//*********************************************************************
	// sissiz-di       cluster.109999_step.aln  8       150     0.8759  0.8542  0.0094  -13.88  -8.20   3.48    -1.63
	protected static String[] ScanSSZ( String Path, String BedFile, double [] stats,  int counter, int id ) throws IOException {
		//stats[0] Mean Pairwise ID 
		//stats[1] Variance
		//stats[2] Normalized Shanon entropy
		//stats[3] GC content
		//stats[4] GAP content
		//stats[5] classical MPI
		String SissizOutTab [] = new String [ 12 ] ; 
		String  Output = "", Error = "" ; 
		String Command = SSZBINARY ;  
		double threshold = 60 ; // threshold on statistic for SISSIz vs SISSIz-RIBOSUM selection
/*		switch ( counter ) {
			case 1:
				Command = Command + " -f 500 -p 0.05" ;
				break;
			case 2:
				Command = Command + " -j -f 500 -p 0.05" ;
				break;
			default:
				break;
		}
*/
		if ( stats[5] < threshold || stats[3] >= 70 ) {
            counter = 2 ;
			Command = Command +" -j "+Path+"/"+BedFile.replaceAll("\t","_")+".aln."+id; // RIBOSUM scoring
		}
		else {
			Command = Command +" "+Path+"/"+BedFile.replaceAll("\t","_")+".aln."+id;  // new scoring
		}
		try {
			long now = System.currentTimeMillis(); 
			long timeoutInMillis = 1000L * 300 ;						  // max 5 minutes  
			long finish = now + timeoutInMillis; 
			// launch initial SISSIz call
			Process Sissiz = Runtime.getRuntime().exec( Command );			
			BufferedReader SissizErr = new BufferedReader(new InputStreamReader(Sissiz.getErrorStream()));
			if (VERBOSE)
				System.out.println(counter+": Running "+Command);
			while ( isAlive( Sissiz ) ) {
				Thread.sleep( 100 );
				if ( System.currentTimeMillis() > finish ) {
					if (VERBOSE)
						System.out.println("SISSIz failed to run within time :-(") ; 					
					SissizErr.close();
					Sissiz.destroy();
					return null ; 
				}
				else {
					// avoid timeout
					if (SissizErr.ready()) {
						Error = SissizErr.readLine() ; 
/*						if ( Error.length() > 17 && Error.substring(0,17).equals( "WARNING: Negative") && counter <2 ) {
							if (VERBOSE) {
								System.out.println( Error ) ; 
								System.out.println("         Launching SISSIz with more flanking sites");
							}
							SissizErr.close(); 
							Sissiz.destroy();							
							// launch  SISSIz with more flanking sites 
							SissizOutTab = ScanSSZ( Path, BedFile, stats , (counter+1), id );
							return SissizOutTab; 
						}
*/					}
				}
			}
			SissizErr.close(); 
			// get Output if process didn't complete in recursion
			if (SissizOutTab[0] == null ) { 
				BufferedReader SissizOut = new BufferedReader(new InputStreamReader(Sissiz.getInputStream()));
				while ( (Output = SissizOut.readLine()) != null ) { 
					if ( Output.length() > 6 && Output.substring(0, 6).equals("sissiz")) {
						if (VERBOSE)
							System.out.println( Output ) ; 
						SissizOutTab = Output.split("\\s") ; 
						SissizOutTab[1] = ( counter == 1 )? "s" : "r";
					}				
				}
				SissizOut.close(); 
			}
			// rerun SISSIz if output is dodgy
			try {
/*>>>>>>*/      if ( Double.parseDouble( SissizOutTab[7]) == 0
						|| (Math.abs( Double.parseDouble( SissizOutTab[9] )) < 0.5 // variance of simulated alignments
						&& counter == 1 )) {
					if (VERBOSE)
						System.out.println("SISSIz gave dodgy output... retrying with RIBOSUM") ; 	
					SissizOutTab = new String [ 12 ] ;
					now = System.currentTimeMillis(); 
					finish = now + timeoutInMillis ;
					Sissiz = Runtime.getRuntime().exec( SSZBINARY+" -j "+Path+"/"+BedFile.replaceAll("\t","_")+".aln."+id); 
					if (VERBOSE)
						System.out.println("Running "+SSZBINARY+" -j "+Path+"/"+BedFile.replaceAll("\t","_")+".aln."+id);
					SissizErr = new BufferedReader(new InputStreamReader(Sissiz.getErrorStream()));
					while ( isAlive( Sissiz ) ) {
						Thread.sleep( 100 );  // do we need this? 
						if ( System.currentTimeMillis() > finish ) {
							if (VERBOSE)
								System.out.println("SISSIz failed to run within time") ; 					
							SissizErr.close(); 
							Sissiz.destroy();
							return null ; 
						}
						else {
							// avoid infinity loop
							if (SissizErr.ready()) {
								Error = SissizErr.readLine() ;
/*								if ( Error.length() > 17 && Error.substring(0,17).equals( "WARNING: Negative") && counter <2 ) {
									if (VERBOSE) {
										System.out.println( Error ) ; 
										System.out.println("         Adding extra flanking sites");
									}
									SissizErr.close(); 
									Sissiz.destroy();	
									// launching with more flanking sites
									SissizOutTab = ScanSSZ( Path, BedFile, stats, 4, id );
									break;
								}
*/							}
						}
					}
					SissizErr.close();
					if (SissizOutTab[0] == null ) { 
						BufferedReader SissizOut = new BufferedReader(new InputStreamReader(Sissiz.getInputStream()));
						while ( (Output = SissizOut.readLine()) != null ) { 
							if ( Output.length() > 6 && Output.substring(0, 6).equals("sissiz")) {
								if (VERBOSE)
									System.out.println( Output +"  ...DONE !!\n" ) ; 
								SissizOutTab = Output.split("\\s") ; 
								SissizOutTab[1] = "r";
							}
						}
						SissizOut.close(); 
					}
				}
				
            }
            catch (Exception Fuck) {
                System.err.println(" Null pointer after launching SSZ with counter = "+counter);
                System.err.println("  with file = "+BedFile.replaceAll("\t","_")+".aln");
                if (SissizOutTab[7] == null )
                    System.err.println("Null output data" );
                Fuck.printStackTrace();
            }
        }
    	catch (Exception err) {
			System.out.println(" Caught Error!\n ----> "+ Command +"\n  counter--> " +counter);
			System.err.println("!!!Caught Error!\n ----> "+ Command +"\n  counter--> " +counter);
			//String [] OutAln = new String[ goodSeqs ] ;
			err.printStackTrace();
			System.err.println("===============================" );
		}
		return SissizOutTab ;
	}
	//*********************************************************************
	//					RNAz scan & parse						*
	//*********************************************************************
	// RnazOutTab  = ScanRNZ( Path, BedFile , stats[ 0 ], 0) ;  
	protected static String[] ScanRNZ( String Path, String BedFile, double statistic, int counter, int id ) throws IOException {
		String RNAzOutTab [] = new String [ 5 ] ; 
		String Output = "", Error = "" ; 
		String Command = RNZBINARY + " -f -d " + Path+"/"+BedFile.replaceAll("\t","_")+".aln."+id;  
		try {
			long now = System.currentTimeMillis(); 
			long timeoutInMillis = 1000L * 300 ;						  // max 5 minutes  
			long finish = now + timeoutInMillis; 
			Process RNAz = Runtime.getRuntime().exec( Command );			// launch RNAz
			// add standard error?? 
			if (VERBOSE)
				System.out.println("Running "+Command);
			while ( isAlive( RNAz ) ) {
				Thread.sleep( 100 );
				if ( System.currentTimeMillis() > finish ) {
					if (VERBOSE)
						System.out.println("RNAz failed to run within time :-(") ; 					
					RNAz.destroy();
					return null ; 
				}
			}
			BufferedReader RNAzOut = new BufferedReader(new InputStreamReader(RNAz.getInputStream()));
            BufferedReader RNAzErr = new BufferedReader(new InputStreamReader(RNAz.getErrorStream()));
			if ( RNAzOut.ready() ) {
				while ( (Output = RNAzOut.readLine()) != null ) { 
					if ( Output.matches("^ Sequences: .+")) {
						RNAzOutTab[ 0 ] = Output.substring( Output.lastIndexOf(":")+1, Output.length());//.replaceAll("\\s","") ; 
					}
					else if ( Output.matches("^ Mean z-score: .+")) {
						RNAzOutTab[ 1 ] = Output.substring( Output.lastIndexOf(":")+1, Output.length());//.replaceAll("\\s","") ; 
					}
					else if ( Output.matches("^ Structure conservation index: .+")) {
						RNAzOutTab[ 2 ] = Output.substring( Output.lastIndexOf(":")+1, Output.length());//.replaceAll("\\s","") ; 
					} 
					else if ( Output.matches("^ SVM decision value: .+")) {
						RNAzOutTab[ 3 ] = Output.substring( Output.lastIndexOf(":")+1, Output.length());//.replaceAll("\\s","") ; 
					}  
					else if ( Output.matches("^ SVM RNA-class probability: .+")) {
						RNAzOutTab[ 4 ] = Output.substring( Output.lastIndexOf(":")+1, Output.length());//.replaceAll("\\s","") ; 
						break ; 
					}
				}
			}
			else {
				if (VERBOSE) {
					System.err.println(" RNAz failed to run :-(");
                    if (RNAzErr.ready()) {
                        String OhNo = RNAzErr.readLine();
                        while (OhNo != null ) {
                            System.err.println( OhNo ) ;
                            OhNo = RNAzErr.readLine();
                        }
                    }
                }
            }
			RNAzOut.close();
		}
		catch (Exception Sheeeit) {
			System.err.println( "Caught exception "+ Sheeeit + " while running RNAZ"); 
		}
		if (RNAzOutTab != null && VERBOSE){
			for (int t = 0 ; t != 5 ; t++ )
				System.out.print( RNAzOutTab[t] + "* *") ; 
			System.out.println(); 
		}
		return RNAzOutTab ;
	}
	//*********************************************************************
	//						Sample process						*
	//*********************************************************************
	private static boolean isAlive( Process p ) {
		try {
			p.exitValue();
			return false;
		} catch (IllegalThreadStateException e) {
			return true;
		}
	}
}
