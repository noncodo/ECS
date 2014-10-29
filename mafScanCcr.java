//
//  mafScan.java
//  - - reads a set of maf files, realigns them with mafft, breaks them into windows, calculates stats, 
//  - - scans with SISSIz, outputs bed coordonates of high-confidence predictions 
//  Created by Martin Smith on 27/04/11.
//  Copyright 2013 All rights reserved.
//
import java.util.*; import java.util.concurrent.* ; import java.io.*; import java.math.*; import java.lang.*;
public class mafScanCcr {

	static int	WINDOW = 200,
				STEP = 100,
				SSZ_THRESHOLD = -2,
				BLOCK_SIZE = 5000,
				MAFFT_LIMIT = 10000-(WINDOW+STEP), // this will not change when arguments are specified. Variable should be reinitiated for inputs
				GAPS = 50,
				NTHREDS = 4 ;
											
	static boolean FILTER_ID = true, //for removing identical sequences
				VERBOSE = false,
				WRITEALL = false,
				REALIGN = false ;  
	static String	FILENAME = "" ,
				OUT_PATH = "",
				SSZBINARY= "~/local/bin/SISSIz-0.3",
				RNZBINARY= "~/local/bin/RNAz-2.0pre";

	public static synchronized void main(String[] Args) throws IOException {
		// variables
		String [][] MafTab ; 
		String [] TempTab ; 
		String Temp = "" ; 
		char [][] AlnTab ; 
		// usage info
		if (Args.length == 0 ){
			System.out.println("\n\t\t\t  version 0.9beta \n"+
                               "\t __  __          ______    _____  _____          _   _\n"+
                               "\t|  \\/  |   /\\   |  ____|  / ____|/ ____|   /\\   | \\ | |\n"+
                               "\t| \\  / |  /  \\  | |__    | (___ | |       /  \\  |  \\| |\n"+
                               "\t| |\\/| | / /\\ \\ |  __|    \\___ \\| |      / /\\ \\ | . ` |\n"+
                               "\t| |  | |/ ____ \\| |       ____) | |____ / ____ \\| |\\  |\n"+
                               "\t|_|  |_/_/    \\_\\_|      |_____/ \\_____/_/    \\_\\_| \\_|\n"+
                               "\t SCAN MULTIPLE ALIGNMENTS FOR CONSERVED RNA STRUCTURES\n\n"+
                               "Reads a set of maf files, realigns them with mafft, breaks them into windows,\n"+
                               "calculates stats, scans with SISSIz or RNAz, outputs bed coordonates of high-confidence predictions\n"+
                               "*** Known issues ***\n"+
                               "-Using MAFFT realignment causes certain jobs to hang when finished; this is likely an executor service problem\n\n"+
							"Usage:     java -jar MafScanCcr.jar [options] -o output/directory -i input.maf (last parameter must be -i)\nOptions:\n"+
                            "  -bs int       Block Size for splitting large MAF blocks (default 5000)\n"+
                            "  -c  int       number of CPUs for calculations (default 4)\n"+
                          	//"  -f            Do not remove identical sequences (not recommended)\n"+
                            "  -g  int       max gap percentage of sequences for 2D prediction (default 50)\n"+ // for one or all ??
						    "  -mafft        Realign with mafft-ginsi (slower)\n"+
                            "  -ml int       Max Length of MAF block for MAFFT realignment (default 10000)\n"+
                            "  -s  int       step size (default 100)\n"+
                            "  -v            verbose (messy but detailed) output\n"+
                            // add a less verbose option, outputing all bed coordinates and scores? 
                            "  -w  int       window size (default 200)\n");
 			System.exit(0);
		}
		// get binary paths
		Process GetBinary = Runtime.getRuntime().exec("which SISSIz") ;
		BufferedReader ReadBin = new BufferedReader(new InputStreamReader(GetBinary.getInputStream()));
		if ( (SSZBINARY = ReadBin.readLine() ) == null ) {
			System.out.println("Please install SISSIz (version 2.0), and link it to your $PATH" );  
			System.exit(0); 
		}
		GetBinary = Runtime.getRuntime().exec("which RNAz") ;
		ReadBin = new BufferedReader(new InputStreamReader(GetBinary.getInputStream()));
		if ( (RNZBINARY = ReadBin.readLine() ) == null ) {
			System.out.println("Please install RNAz (version 2.0) and link it to your $PATH" );  
			System.exit(0); 
		}
		ReadBin.close() ; 
		// parse arguments
		for (int i = 0 ; i != Args.length ; i++ ) {
			// add filter_id flag < < < < < < < < < < < < < <
			if ( Args[ i ].equals("-w") ) {  // window size
				WINDOW = Integer.parseInt( Args[ i+1 ] ) ;
				i++ ; 
			}
			else if ( Args[ i ].equals("-c") ) {  // Threads
				NTHREDS = Integer.parseInt( Args[ i+1 ] ) ;
				i++ ;
			}
			else if ( Args[ i ].equals("-f") ) {  // Threads
				FILTER_ID = false ; 
				i++ ;
			}
			else if ( Args[ i ].equals("-g") ) {  // gap content
				GAPS = Integer.parseInt( Args[ i+1 ] ) ;
				i++ ; 
			}
			else if ( Args[ i ].equals("-o") ) { //output directory
				OUT_PATH= Args[ i+1 ] ;
				if ( !(new File(OUT_PATH)).isDirectory() )
					(new File(OUT_PATH)).mkdirs() ;	
				if (VERBOSE) 
					System.out.println("writing alignments to directory "+OUT_PATH) ; 
				i++ ; 
			}
			else if ( Args[ i ].equals("-s") ) { // step size
				STEP = Integer.parseInt( Args[ i+1 ] )	; 
				i++ ; 
			}
			else if ( Args[ i ].equals("-bs") ) { // step size
				BLOCK_SIZE = Integer.parseInt( Args[ i+1 ] )	; 
				i++ ; 
			}
			else if ( Args[ i ].equals("-ml") ) { // step size
				MAFFT_LIMIT = Integer.parseInt( Args[ i+1 ] )	; 
				i++ ; 
			}
			else if ( Args[ i ].equals("-v") ) { // verbose output
				VERBOSE = true; 
			}
			else if ( Args[ i ].equals("-mafft") ) { // realign
				REALIGN = true; 
			}
			else if ( Args[ i ].equals("-i") ) {			
				i++ ;
				FILENAME = Args[i].substring( Args[i].lastIndexOf("/")+1, Args[i].length()-4 ) ;
				//parse out individual alignment blocks from a multi maf file
				int lineCount= 0 ; 
				BufferedReader ReadFile = new BufferedReader(new FileReader( Args[ i ] ));
				String Line = "" ; 
				while ( (Line = ReadFile.readLine()) != null)   // count lines for array
					if ( Line.length() > 1 && Line.charAt(0) != '#' )
						lineCount++ ; 
					
				if (VERBOSE)
					System.out.println( "Read "+(lineCount-1)+" sequences from file "+ FILENAME);
				ReadFile.close() ; 
				ReadFile = new BufferedReader(new FileReader( Args[ i ] ));
				TempTab = new String [ lineCount ] ; // 7 maf columns
				// fill in array from file
				int newLineCount = 0; 
				/************************************************************************
				 ****   This stuff is messy, but avoids problems at last block       ****
				 ************************************************************************/
				readBlocks: while ( (Line = ReadFile.readLine()) != null || newLineCount == lineCount ) {
					if (  Line == null || Line.length() > 1 && Line.charAt(0) != '#' )
						newLineCount++ ; 
					if ( newLineCount > lineCount || Line.length() > 1 ) { // only lines with sequences
						try {
							if ( newLineCount <= lineCount && Line.length() > 1 && Line.substring(0,1).equals("s") ) {
								Temp = Temp + Line + "@" ; 
							}
							else if ( (Line != null && Line.length() > 1 && Line.substring(0,1).equals("a") ) || newLineCount > lineCount ) {
								if ( newLineCount == 1 ) 
									continue readBlocks ; 
								
								if ( Temp.split("@").length >= 3 ) { // at least 3 sequences
									TempTab = new String [ Temp.split("@").length ] ; 
									TempTab = Temp.split("@") ; 
									Temp = "" ; 
									MafTab = new String [ TempTab.length ] [ 7 ] ; 
									for ( int x = 0 ; x != TempTab.length ; x++ ) {
										MafTab[ x ] = TempTab [ x ].split("\\s+");
									}
									// make sure block is bigger than window, or default to smaller window   <<<<<    TO DO !!!!!!!!!!!!
									if ( MafTab[0][6].length() < WINDOW ) {
										if ( VERBOSE ) 
											System.out.println("Block is too small!!");
										continue readBlocks;
									}
									else if (!REALIGN)		// quick and dirty, just how we like it
									{
										AlnTab = new char [ MafTab.length ][] ;
										//lose the sequence information
										String [][] SmallMafTab = new String [ MafTab.length ][ 5 ] ;
										for (int r = 0 ; r != MafTab.length ; r++ ){
											AlnTab[r] = MafTab[ r ][ 6 ].toCharArray() ; 
											for (int c = 0 ; c != 5 ; c++ )
												SmallMafTab[ r ][ c ] = MafTab[ r ][ c ];  
										}
										SplitNfold( SmallMafTab, AlnTab );
									}
									else {
										// get proper bed coordinates
										// perhaps next version? 
										
										// check and split if too large
										if ( MafTab[0][6].length() > MAFFT_LIMIT + WINDOW ) {
											if (VERBOSE)
												System.out.println("--> Splitting large input alignment "+ MafTab[0][6].length()) ; 
											int from = 0, to = BLOCK_SIZE ; 
											String [][] MafPartTab = new String [ MafTab.length ][ 7 ] ;
											// copy split maf tab
											while ( from < MafTab[0][6].length() ) {
												for (int y = 0 ; y != MafTab.length ; y++ ) {
													for ( int z = 0 ; z != 7 ; z++ ) { 
														// deal with strands
														if ( z == 2 && y == 0 ) {
															MafPartTab[ 0 ][2] = (MafTab[ 0 ][ 4 ].charAt(0) == '+' )?
															""+(Integer.parseInt(MafTab[ 0 ][2]) + from ): 
															""+(Integer.parseInt(MafTab[ 0 ][2]) - from ); 	
														}
														// split maf sequence
														else if ( z == 6 ) {
															MafPartTab[ y ][ 6 ] = MafTab[ y ][ 6 ].substring( from, to ) ;
															// change alignment length annotation
															MafPartTab[ y ][ 3 ] = ""+MafPartTab[ y ][ 6 ].replaceAll( "[\\.-]", "").length() ; 
														}
														//ignore length
														else if ( z !=4 ) 
															MafPartTab[ y ][ z ] = MafTab[ y ][ z ] ; 
														
													}
												}				
												if (VERBOSE) {
													System.out.println( "-----> Realigning LARGE ("+ MafTab[0][6].length() +"nt) ALN of "+MafTab.length+" sequences");
													System.out.println("- - -> Block "+ from + " : " + to  ) ;
												}
												AlnTab = 	Realign( MafPartTab ) ;
												// make sure it's not empty!
												from = from + BLOCK_SIZE - ( WINDOW - STEP ) ; 
												to = Math.min( MafTab[0][6].length(), to + BLOCK_SIZE-( WINDOW - STEP )); 
												if (VERBOSE) 	
													System.out.println( "------> Scanning alignment");
												if ( AlnTab[0][0] == '@' ) 
													continue ; 
												//lose the sequence information
												String [][] SmallMafTab = new String [ MafTab.length ][ 5 ] ;
												for (int r = 0 ; r != MafTab.length ; r++ )
													for (int c = 0 ; c != 5 ; c++ )
														SmallMafTab[ r ][ c ] = MafTab[ r ][ c ];  
												SplitNfold( SmallMafTab, AlnTab );
											}
										}
										else {
											if (VERBOSE) 
												System.out.println( "----> Realigning small ALN of "+TempTab.length+" sequences");
											AlnTab = 	Realign( MafTab ) ; 
											if ( AlnTab[0][0] == '@' )
												continue ; 
											if (VERBOSE) 
												System.out.println("Free!!!!"); 
											
											if (AlnTab[0].length < 100 ) 
												continue ;
											if (VERBOSE) 
												System.out.println( "-------> Scanning alignment");
											//lose the sequence information
											String [][] SmallMafTab = new String [ MafTab.length ][ 5 ] ;
											for (int r = 0 ; r != MafTab.length ; r++ )
												for (int c = 0 ; c != 5 ; c++ )
													SmallMafTab[ r ][ c ] = MafTab[ r ][ c ];  
											SplitNfold( SmallMafTab, AlnTab );
										}
									}
									Temp = "" ; 								
								}
								else
									Temp = "" ; 
							}
						} 
						catch ( Exception Shit ){
							System.err.println("Problem parsing alignment blocks\n"+Shit);
							Shit.printStackTrace();
						}
					}
				}
				// Wait until all threads are finished

				ReadFile.close();
			}
		}
	}
	//*********************************************************************
	//						scan and parse windows				*
	//*********************************************************************
	private static synchronized boolean SplitNfold (String[][] MafTab, char[][] AlnTab) throws IOException {
        ExecutorService MultiThreads = Executors.newFixedThreadPool( NTHREDS );
        try {

		// check if dir exists
		// add Path Flag ? < < < < < < < < < < < < < < < < < < < < < < < < < < < < < 
		String Path = OUT_PATH+"/aln/"+MafTab[0][1].substring( MafTab[0][1].lastIndexOf(".")+1);
		if ( !(new File(Path)).isDirectory() )
			(new File(Path)).mkdirs() ;
		// parse human maf for bed file
		int [] coordTab = new int [ AlnTab[0].length ] ; 
		// may have to check exact coordonates
		int charCounter =  0 ; 
		for (int i = 0 ; i != coordTab.length ; i++ ) {
			if ( AlnTab[0][i] == 'A'  || AlnTab[0][i] =='C' || AlnTab[0][i] == 'T' || AlnTab[0][i] == 'G' || AlnTab[0][i] == 'N' ||
			    AlnTab[0][i] == 'a'  || AlnTab[0][i] =='c' || AlnTab[0][i] == 't' || AlnTab[0][i] == 'g' || AlnTab[0][i] == 'n' ) 
				charCounter++ ;
			else if ( AlnTab[0][i] != '-' && AlnTab[0][i] != '_' && AlnTab[0][i] != '.' ) 
				System.out.println ( "##################    unknown character : "+ AlnTab[0][i] ) ; 
			coordTab[ i ] = charCounter ; 
		}
        
        List<Future<Runnable>> futures = new ArrayList<Future<Runnable>>();

		//************    make and scan windows concurrently    **************
		for (int step = 0 ; step < AlnTab[0].length-WINDOW+STEP-1 ; step = step+STEP) {
			// avoid small step when last 2 windows are closer than STEP/2 
			if ( step != 0 && step > AlnTab[0].length-(WINDOW+STEP/2) )
				step = AlnTab[0].length-WINDOW ; 
			// remember 'N' residues and gap-only sequences
			//>>>>>>> add FILTER_ID <<<<<<<<
            Future f = MultiThreads.submit(new ScanItFast( step, MafTab, AlnTab, Path, coordTab, WINDOW, STEP, GAPS, SSZ_THRESHOLD, SSZBINARY, RNZBINARY, VERBOSE, FILTER_ID ));
            futures.add(f);
            //Runnable Worker = new ScanItFast( step, MafTab, AlnTab, Path, coordTab, WINDOW, STEP, GAPS, SSZ_THRESHOLD, SSZBINARY, RNZBINARY, VERBOSE );
			//MultiThreads.execute( Worker );
			// avoid infinite loop!
			if (step == (AlnTab[0].length-WINDOW) || AlnTab[0].length < WINDOW  ) 
				break  ;
		}
        for (Future<Runnable> f : futures)
        {
            f.get();
        }
        MultiThreads.shutdown();
            MultiThreads.awaitTermination(60*10L, TimeUnit.SECONDS);
        } catch (Exception Fuck) {
            System.err.println("MultiThreads took too long!  OOps!" );
            Fuck.printStackTrace();
        }
        return true; 
	}
		
	//*********************************************************************
	//						MAFFT Realign                                 *
	//*********************************************************************
	private static synchronized char[][] Realign( String[][] MafTab ) throws IOException {
        char [][] realignment = new char [ 1 ][ 1 ];
        try{
		// write to fasta format
        BufferedWriter Fasta = new BufferedWriter(new FileWriter( FILENAME+".fasta" ));
		for ( int x = 0 ; x != MafTab.length ; x++ ) {
			// make sure no empty lines are included
//			if ( MafTab[x][6].replaceAll("[\\.-]", "").length() != 0 ) {
				Fasta.write( ">"+ MafTab[x][1] + "\n" ) ; 
				Fasta.write( MafTab[x][6].replaceAll("[\\.-]", "") +"\n") ;
//				if ( VERBOSE)
//                   System.out.println( MafTab[x][1].substring(0,3)+ "\t"+MafTab[x][6].substring(0, Math.min(MafTab[x][6].length(), 100 ) ) ) ;
//			}
		}
		Fasta.close() ; 
		if ( MafTab.length < 3 ) {
			char [][] empty = new char [ 1 ][ 1 ] ; 
			empty[0][0] = '@' ; 
			if (VERBOSE)
				System.out.println("..Less than 3 seqs to align!");
			return empty ;
		}
		//find where mafft is installed
		Process GetBinary = Runtime.getRuntime().exec("which mafft-ginsi") ;
		BufferedReader ReadBin = new BufferedReader(new InputStreamReader(GetBinary.getInputStream()));
		String MafftBin = "" ; 
		if ( (MafftBin = ReadBin.readLine() ) == null ) {
			System.out.println("Please install MAFFT and link its binaries to your $PATH" );  
			System.exit(0);
		}
		ReadBin.close(); 
		Process Mafft = Runtime.getRuntime().exec(MafftBin+" --thread "+ ( NTHREDS )+" --quiet "+FILENAME+".fasta" ) ;

		BufferedReader ReAligned = new BufferedReader(new InputStreamReader(Mafft.getInputStream()));
		String Sequence = "", Line = "" ;
            
		String [] Gapless = new String [ MafTab.length ] ; 
		int seq = 0 ; 
  		ReAligned.readLine() ; 
		while ((Line = ReAligned.readLine()) != null) {
			if (Line.charAt(0) == '>'){
				seq++ ; 
				Gapless[ seq -1 ] = Sequence.toUpperCase() ; 
				Sequence = ReAligned.readLine() ;
			}
			else
				Sequence = Sequence + Line ; 
		}
		Gapless[ seq ] = Sequence.toUpperCase() ; 
            Mafft.waitFor() ;
        ReAligned.close() ;
            
		if (VERBOSE)
			System.out.println("Realigned "+ MafTab.length +" sequences");
        realignment = new char [ MafTab.length ][ Gapless[0].length() ];
		for ( int i = 0 ; i != MafTab.length ; i++ ) {
//null pointer !
			realignment[i] = Gapless[i].toCharArray();
		}
/*        if (VERBOSE && lineCount != MafTab.length ){
            System.out.println( "AlnTab = "+ lineCount+"; MafTab = "+ MafTab.length );
            for (int i = 0 ; i != lineCount ; i++ ){
                for (int j = 0 ; j != 100 ; j++ )
                    System.out.print( realignment[i][j]) ; 
                    System.out.println();
            }
        }
*/
		File Temp = new File( FILENAME+".fasta" ) ;
		Temp.delete() ; 

        }
        catch (InterruptedException Fuck) {
            System.err.println("Problem realigning...." );
            Fuck.printStackTrace();
        }
        return realignment ;

    }
 
    private static boolean isAlive( Process p ) {
		try {
			p.exitValue();
			return false;
		} catch (IllegalThreadStateException e) {
			return true;
		}
	}
    
}	


/*
public class SafeProgA {
	public static void main(String[] args) throws Exception {
		Runtime rt = Runtime.getRuntime();
		Process p = rt.exec("Whatever your exec does");
		StreamGobbler s1 = new StreamGobbler ("stdin", p.getInputStream ());
		StreamGobbler s2 = new StreamGobbler ("stderr", p.getErrorStream ());
		s1.start ();
		s2.start ();
		p.waitFor();
		System.out.println("Process Returned");
	}
}
*/