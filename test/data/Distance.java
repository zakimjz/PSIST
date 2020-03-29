// extract coordinates of Atoms

import java.io.*;
import java.util.*;

public class Distance
{

    public Distance()
    {
       
    }
    public static void main(String argv[])
    {
	if (argv != null) { 
	    int len = argv.length;
	    if (len != 4) {
		System.out.println("Usage: java pdblist <input dir> <output dir> <all/ca/backbone>");
		System.exit(1);
	    }
	}
	

	//String PDBL="PDB_CHAIN_N_LIST.txt";
	String PDBL= argv[0];
	System.out.println("list:"+PDBL); 
	int num_files=0;
	try
	    {
		File file=new File(PDBL);
		FileReader fis=new FileReader(file);
		BufferedReader buff=new BufferedReader(fis);
		String data;
		while((data=buff.readLine())!=null)
		    {
			String[] things=data.split(" ");
			System.out.println("PDB ID:"+things[0]+"Chain:"+things[1]);
			//File fil=new File("PDB_Select/"+things[0]+".pdb");
			File fil=new File(argv[1]+things[0]+".pdb");
			FileReader fir=new FileReader(fil);
			BufferedReader buf=new BufferedReader(fir);
			double[][] Coord=new double[15000][4];
			int pos=0;
			String line;
			int start_pos=1;
			boolean found_seq=false;
			int seq_len=0;
			while((line=buf.readLine())!=null)
			    {
				if(!found_seq&&line.startsWith("SEQRES"))
				    {
					String sl=line.substring(12,17);
					seq_len=(int)Double.parseDouble(sl);
					found_seq=true;
				    }
				if(line.startsWith("ENDMDL"))
				    break;
				if(line.startsWith("ATOM"))
				    {
				    
					String atom=line.substring(13,15);
					//System.out.println("line "+ atom); 
					if(atom.startsWith("CA")){
					    //System.out.println("CA " + atom);
					    Coord[pos][0]= 1;
					}
					
					else{
					    if(atom.startsWith("N") && atom.endsWith(" "))
						Coord[pos][0]= 2;
					    else{
						if(atom.startsWith("C") && atom.endsWith(" "))
						    Coord[pos][0]= 3;
						else
						    Coord[pos][0]= 4;
					    }
					}
					
					if(things[1].endsWith("_"))
					    {
						String x_c=line.substring(30,38);
						Coord[pos][1]=Double.parseDouble(x_c);
						String y_c=line.substring(38,46);
						Coord[pos][2]=Double.parseDouble(y_c);
						String z_c=line.substring(46,54);
						Coord[pos][3]=Double.parseDouble(z_c);
						if(pos==0)
						    {
							String st=line.substring(22,26);
							start_pos=(int)Double.parseDouble(st);
						    }
						pos++;
					    }
					else
					    {
						String chain_type=line.substring(21,22);
						if(chain_type.endsWith(things[1]))
						    {
							String x_c=line.substring(30,38);
							Coord[pos][1]=Double.parseDouble(x_c);
							String y_c=line.substring(38,46);
							Coord[pos][2]=Double.parseDouble(y_c);
							String z_c=line.substring(46,54);
							Coord[pos][3]=Double.parseDouble(z_c);
							if(pos==0)
							    {
								String st=line.substring(22,26);
								start_pos=(int)Double.parseDouble(st);
							    }
							pos++;
						    }
					    }
					
				    }
			    }
			System.out.println("SoM:"+pos+" seq_len:"+seq_len+" start:"+start_pos);
			buf.close();
			fir.close();
			//Spurt out the distance file
			File outer=new File(argv[2]+things[0]+things[1]+".coord");
			FileWriter fiw=new FileWriter(outer);
			BufferedWriter buw=new BufferedWriter(fiw);
			int end_l=seq_len-1;
			
			if (argv[3].startsWith("all")){
			
			    System.out.println("all atoms are written");
			    
			    int id =0;
			    for(int i=0;i<pos;i++)
				if(Coord[i][0] == 1.0)
				{
				  String output_c=""+ id + " " + Coord[i][0]+" "+Coord[i][1]+ " "+Coord[i][2]+ " "+Coord[i][3];
				  id++;
				  buw.write(output_c,0,output_c.length());
				  
				  buw.newLine();
				}

			    
			    for(int i=0;i<pos;i++)
				if(Coord[i][0] != 1.0)
				    {
					String output_c=""+ id + " " + Coord[i][0]+" "+Coord[i][1]+ " "+Coord[i][2]+ " "+Coord[i][3];
					id++;
					buw.write(output_c,0,output_c.length());
					
					buw.newLine();
				    }
			}

			if (argv[3].startsWith("backbone")){
			    
			    System.out.println("3 backbone atoms(CA, N, C) are written");
			    
			    int id =0;
			    for(int i=0;i<pos;i++)
				if(Coord[i][0] == 1.0 || Coord[i][0] == 2.0 || Coord[i][0] == 3.0)
				    {
					String output_c=""+ id + " " + Coord[i][0]+" "+Coord[i][1]+ " "+Coord[i][2]+ " "+Coord[i][3];
					id++;
					buw.write(output_c,0,output_c.length());
					
					buw.newLine();
				    }
			    
			}
			
			
			else{
			    System.out.println("only Ca atoms are written");
			    int id =0;
			    
			    for(int i=0;i<pos;i++)
				if(Coord[i][0] == 1.0)
				    {
					String output_c=""+ id + " "+Coord[i][1]+ " "+Coord[i][2]+ " "+Coord[i][3];
					id++;
					
					buw.write(output_c,0,output_c.length());
				    
					buw.newLine();
				    }
			}
			
			

			buw.close();
			fiw.close();
			num_files++;
		    }
		System.out.println("Total of "+num_files+" distance matrices calculated and written.");
		buff.close();
	    }
	catch(Exception e){System.out.println("Problemo"+e);}
    }
}
