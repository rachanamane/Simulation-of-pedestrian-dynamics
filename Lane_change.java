

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
// ------------------------------------------------------------------------------------------------------------------------
public class Lane_change 
{	
	
	static final double width = 2.5; // meter width of the corridor under study
	static final double ped_flow =500;  // pedestrian/hour (pedestrian flow)

	
	//---------------------------------------------------------------------------------------------------------------------------------
	static double Tsize = 0.37-(4*ped_flow/100000);  // length of the sides of triangular cell as a function of the pedestrian flow rate
	static double vmax=4;  // maximum a pedestrian can move is vmax cells in one step 
	static final int time_step = 1;		// variable to specify the increment in the time step 
	static int steps=0; 				// initiating the time of study
	static int total_steps= 10*60;			// total simulation time 
	static double spot = 45;			// point of analysis
	static int ped_count=0;				// initiating the counter for total pedestrians
	static int timeNext = 1;			// variable for time keeping during simulation
	static double bas_lc=-1.5;			// weight to the lane changing effect in basic desire of a pedestrian 
	static double bas_l=1.5;			// weight to the straight movement in basic desire of a pedestrian 
	static double inter_lc=3;			// weight to the lane changing effect in interference matrix of a pedestrian 
	static double inter_l=2;			// weight to the straight movement effect in interference matrix of a pedestrian 
	static int interference_count = -1;  // 4 * vmax^2 + 12*vmax + 6
	static int base_des_count = -1;  // vmax^2 + 2*vmax

	//---------------------------------------------------------------------------------------------------------------------------------

	// function is self explanatory
	
	public static double power(int a, int b)
		{
			return Math.pow(ped_flow*time_step/3600, a) * Math.exp(-ped_flow*time_step/3600)/b;
		}		
	
	//---------------------------------------------------------------------------------------------------------------------------------
	
	public static double [][][] bdesmat(double vmax, double [][][] bas_des_mat, int ped_count) throws FileNotFoundException{
	// Function to initiate the matrix with basic desire for each pedestrian
		PrintWriter basedesmat = new PrintWriter("base_mat.txt");  // Writing basic desire matrix in a file named "base_mat.txt"	
		
		basedesmat.println("This is base des mat");
		
		for (int i = 1; i <= ped_count; i++) 
			{
			int ii = -1;
				for (int v =0; v<=vmax;v++)
					{
						ii++;
						bas_des_mat[i][ii][0] = v;
						bas_des_mat[i][ii][1] = 0;
						bas_des_mat[i][ii][2] = 12+v*bas_l;
						basedesmat.println(bas_des_mat[i][ii][0]+"\t"+ bas_des_mat[i][ii][1] +"\t"+	bas_des_mat[i][ii][2]);
					}		
				for (int d = 1; d<=vmax;d++)
					{
						for (int v=d;v<=vmax;v++)
							{
								ii++;		
								bas_des_mat[i][ii][0] = v-d*0.5;
								bas_des_mat[i][ii][1] = d;  //*ydir; // ydir=1
								bas_des_mat[i][ii][2] =  12 + v*bas_l + d*bas_lc;
								basedesmat.println(bas_des_mat[i][ii][0]+"\t"+ bas_des_mat[i][ii][1] +"\t"+	bas_des_mat[i][ii][2]);
							}
					}		
				for (int d = 1; d<=vmax;d++)
					{
						for (int v=d;v<=vmax;v++)
							{
								ii++;
								bas_des_mat[i][ii][0] = v-d*0.5;
								bas_des_mat[i][ii][1] = -d;  //*ydir; // ydir=-1
								bas_des_mat[i][ii][2] =  12 + v*bas_l + d*bas_lc;
								basedesmat.println(bas_des_mat[i][ii][0]+"\t"+ bas_des_mat[i][ii][1] +"\t"+	bas_des_mat[i][ii][2]);
							}
					}
				base_des_count = ii;
			}
		
		basedesmat.close();
		return bas_des_mat;
		
	}	

	//----------------------------------------------------------------------------------------------------------------------------------
	
	public static double [][][] intmat(double vmax, double [][][] int_mat, int ped_count) throws FileNotFoundException{

		PrintWriter basedesmat = new PrintWriter("base_mat_int.txt");  // Writing interference matrix in a file named "base_mat.txt"	
		
		// Function to initiate the matrix with negative effect by each pedestrian on a nearby pedestrian
		basedesmat.println("This is base int mat");
		
	
		for (int ii = 1; ii <= ped_count; ii++) 
		{
			int i = -1;
			for (double v =-(vmax+1); v<=(vmax+1);v++)
				{
					i++;
					int_mat[ii][i][0] = v;
					int_mat[ii][i][1] = 0;
					int_mat[ii][i][2] = inter_l*(Math.abs(v))-12;	
					basedesmat.println(int_mat[ii][i][0]+"\t"+ int_mat[ii][i][1] +"\t"+	int_mat[ii][i][2]);				
				}	
					
			// for lane change = 0
			for (double d = 1; d<=(2*vmax+1);d++)
				{
					for (double v=-(vmax+1-d/2);v<=(vmax+1-d/2);v++)
						{
							i++;				
							int_mat[ii][i][0] = v;
							int_mat[ii][i][1] = d; //(*ydir;) // int ydir = 1;
							int_mat[ii][i][2] =  d*inter_lc+inter_l*(Math.abs(v))-20;
							basedesmat.println(int_mat[ii][i][0]+"\t"+ int_mat[ii][i][1] +"\t"+	int_mat[ii][i][2]);
						}
				}		
					
			// for lane change = 1	
			
			for (double d = 1; d<=(2*vmax+1);d++)
				{
					for (double v=-(vmax+1-d/2);v<=(vmax+1-d/2);v++)
						{
							i++;
							int_mat[ii][i][0] = v;
							int_mat[ii][i][1] = -d; //(*-ydir;) // int ydir = -1;
							int_mat[ii][i][2] = d*inter_lc+inter_l*(Math.abs(v))-20;
							basedesmat.println(int_mat[ii][i][0]+"\t"+ int_mat[ii][i][1] +"\t"+	int_mat[ii][i][2]);							
						}
												
				}
			interference_count = i;
		}
		basedesmat.close();
		return int_mat;							
}				

	//----------------------------------------------------------------------------------------------------------------------------------
	
	public static double [][][][] negeff(double [][][][] returnarray,int ped_count,int total_steps,double [][][] bas_des_mat, double [][][] int_mat,double [][][] main,int time,double vmax, double ywidth) throws FileNotFoundException{	
	// This function checks for the pedestrian in the vicinity of pedestrian under consideration. This function is performs its function for each pedestrian on the floor
		double [][][][] vicinityPed =new double [(2+total_steps)][(1+ped_count)][99][5];

		int i = 1;		
		int flag1 = 1;
		if(main[time][i][0] == 0){flag1 = 0;}
			while(flag1 == 1) 
				{
					int vicpedcount=0;			
					int j = 1; 
					int flag2 = 1;
					if(main[time][j][0] == 0){flag2 = 0;}
						while(flag2 == 1) 
							{
							
								if (j != i)
								{
									double relativeX = main[time][i][1] - main[time][j][1];
									double relativeY = main[time][i][2] - main[time][j][2];
									
									if(Math.abs(relativeX) + Math.abs(relativeY) / 2 <= 4 && Math.abs(relativeY) <= 7 &&  main[time][j][2]*main[time][j][2]+main[time][j][1]*main[time][j][1] != 0 && main[time][i][1]<main[time][j][1] )
									{
										vicpedcount++;
										for (int n = 0; n <= 4; n++) 
											{	
												vicinityPed[time][i][vicpedcount][n] = main[time][j][n];										
											}	
										double[][][] vic_intmat = new double [999][int_mat[j].length][3];
										for (int k = 0; k < interference_count; k++) 
											{ 
												for (int n = 0; n <= 2; n++) 
													{ 
														vic_intmat[vicpedcount][k][n]=int_mat[j][k][n];
													}
											}

										for (int k = 0; k <interference_count; k++) 
											{	
												vic_intmat[vicpedcount][k][0] += relativeX;		// shifting the X coordinates of Interference matrix as per the Pedestrian
												vic_intmat[vicpedcount][k][1] += relativeY;		// shifting the Y coordinates of Interference matrix as per the Pedestrian
												
												for (int l = 0; l < base_des_count; l++) 
													{	
														if (vic_intmat[vicpedcount][k][0] == bas_des_mat [i][l][0] && vic_intmat[vicpedcount][k][1] == bas_des_mat [i][l][1])
															{		
															bas_des_mat [i][l][2] = bas_des_mat [i][l][2] + vic_intmat[vicpedcount][k][2];
															}
													}
											}
									}
									
								}
								
								j++;
								if(main[time][i][0] == 0||j>ped_count){flag2 = 0;}
							
							}
				
					for (int k = 0; k < base_des_count; k++) 
					{								
						if ((bas_des_mat [i][k][1] + main[time][i][2]) > (ywidth - 1) || (bas_des_mat [i][k][1] + main[time][i][2]) < 1)
						{
							bas_des_mat [i][k][2] = -1000;								
						}				
				
					}
					
					i++;
					if(main[time][i][0] == 0||i>ped_count){flag1 = 0;}	
						
				}
							
		
		returnarray[0]=bas_des_mat;
		returnarray[1]=main;


		return returnarray;
	// Final output of the function is the updated basic desirability matrix.	
}

	//----------------------------------------------------------------------------------------------------------------------------------
	
	public static double [][][] updatePosition(double ywidth,int ped_count, double [][][] main, int timeNext, int time, double [][][] bas_des_mat) throws FileNotFoundException{
	// Based of the updated basic desire matrix, this function invokes logit model and accordingly update the position of each pedestrian
		
		PrintWriter averagespeed= new PrintWriter("Average_speed.txt");  // saving pedestrian data: ped_count, X, Y, random, V
		averagespeed.println( "Ped ID" + "\t"+ "Velocity");
		int i=1;			
		int flag1 = 1;
			if(main[time][i][0] == 0){flag1 = 0;}
				while(flag1 == 1) 
					{
						
						ArrayList<ArrayList<Double>> sorted = new ArrayList<ArrayList<Double>>();
						
						for (int j = 0; j <=base_des_count ; j++) 
							{
								ArrayList<Double> temp = new ArrayList<Double>();
								temp.add(bas_des_mat[i][j][0]);	
								temp.add(bas_des_mat[i][j][1]);
								temp.add(bas_des_mat[i][j][2]);
								sorted.add(temp);
								
							}
						Collections.shuffle(sorted);
						
						
						
						Collections.sort(sorted, new Comparator<ArrayList<Double>>()
											{
												public int compare(ArrayList<Double> a, ArrayList<Double> b)
													{
														return -Double.compare(a.get(2), b.get(2));						
													}
											}
										);
				
						
							double probabSum = Math.exp(sorted.get(0).get(2))+Math.exp(sorted.get(1).get(2))+Math.exp(sorted.get(2).get(2));				// LOGIT
			//				double probabSum = sorted.get(0).get(2)+sorted.get(1).get(2)+sorted.get(2).get(2);				
							
							if (sorted.get(0).get(2)*sorted.get(1).get(2)*sorted.get(2).get(2)>-1000000){
							
								double random =Math.random();
							
								if (random <= (  Math.exp(sorted.get(0).get(2)) / probabSum))										// Probability as per LOGIT
									{	
										int q=0;
										if(main[time][i][2] + sorted.get(q).get(1)>= 0 && main[time][i][2] + sorted.get(q).get(1)<6) 
										{	
											if ( main[time][i][1] < spot &&  main[time][i][1] + sorted.get(q).get(0) > spot )
												{
													averagespeed.println(main[timeNext][i][0] + "\t"+sorted.get(q).get(0)+Math.abs((sorted.get(q).get(1))/2));
												}
											
											main[timeNext][i][1] = main[time][i][1] + sorted.get(q).get(0);
											main[timeNext][i][2] = main[time][i][2] + sorted.get(q).get(1);
											main[timeNext][i][4] = Math.abs(sorted.get(q).get(0))+Math.abs((sorted.get(q).get(1))/2);		
										}	
										else if(main[time][i][2] + sorted.get(q).get(1)< 0) 
										{
											if ( main[time][i][1] < spot &&  main[time][i][1] + sorted.get(q).get(0) > spot )
											{
												averagespeed.println(main[timeNext][i][0] + "\t"+sorted.get(q).get(0)+Math.abs((sorted.get(q).get(1))/2));
											}
											
											main[timeNext][i][1] = main[time][i][1] + sorted.get(q).get(0);
											main[timeNext][i][2] =0;
											main[timeNext][i][4] = Math.abs(sorted.get(q).get(0))+Math.abs((0-main[timeNext][i][2])/2);	
										}
										else if (main[time][i][2] + sorted.get(q).get(1)>5) 
										{
											if ( main[time][i][1] < spot &&  main[time][i][1] + sorted.get(q).get(0) > spot )
											{
												averagespeed.println(main[timeNext][i][0] + "\t"+sorted.get(q).get(0)+Math.abs((sorted.get(q).get(1))/2));
											}
											
											main[timeNext][i][1] = main[time][i][1] + sorted.get(q).get(0);
											main[timeNext][i][4] = Math.abs(sorted.get(q).get(0))+Math.abs((5-main[timeNext][i][2])/2);	
											main[timeNext][i][2] = 5;
										}
										
										
											for (int j = 1; j < i; j++) 
												{
													if (main[timeNext][j][1] == main[timeNext][i][1] && main[timeNext][j][2] == main[timeNext][i][2] ) 
														{ 
														q=1;
														if(main[time][i][2] + sorted.get(q).get(1)>= 0 && main[time][i][2] + sorted.get(q).get(1)<6) 
														{
															main[timeNext][i][1] = main[time][i][1] + sorted.get(q).get(0);
															main[timeNext][i][2] = main[time][i][2] + sorted.get(q).get(1);
															main[timeNext][i][4] = Math.abs(sorted.get(q).get(0))+Math.abs(((sorted.get(q).get(1))/2));		
														}	
														else if(main[time][i][2] + sorted.get(q).get(1)< 0) 
														{
															main[timeNext][i][1] = main[time][i][1] + sorted.get(q).get(0);
															main[timeNext][i][2] =0;
															main[timeNext][i][4] = Math.abs(sorted.get(q).get(0))+Math.abs(((0-main[timeNext][i][2])/2));	
														}
														else if (main[time][i][2] + sorted.get(q).get(1)>5) 
														{
															main[timeNext][i][1] = main[time][i][1] + sorted.get(q).get(0);
															main[timeNext][i][4] = Math.abs(sorted.get(q).get(0))+Math.abs((5-main[timeNext][i][2])/2);	
															main[timeNext][i][2] = 5;
														}
															for (int j1 = 1; j1 < i; j1++) 
																{ 
																	if (main[timeNext][j1][1] == main[timeNext][i][1] && main[timeNext][j1][2] == main[timeNext][i][2] ) 
																		{
																		q=2;
																		if(main[time][i][2] + sorted.get(q).get(1)>= 0 && main[time][i][2] + sorted.get(q).get(1)<6) 
																		{
																			main[timeNext][i][1] = main[time][i][1] + sorted.get(q).get(0);
																			main[timeNext][i][2] = main[time][i][2] + sorted.get(q).get(1);
																			main[timeNext][i][4] = Math.abs(sorted.get(q).get(0))+Math.abs((sorted.get(q).get(1))/2);		
																		}	
																		else if(main[time][i][2] + sorted.get(q).get(1)< 0) 
																		{
																			main[timeNext][i][1] = main[time][i][1] + sorted.get(q).get(0);
																			main[timeNext][i][2] =0;
																			main[timeNext][i][4] = Math.abs(sorted.get(q).get(0))+Math.abs(((0-main[timeNext][i][2])/2));	
																		}
																		else if (main[time][i][2] + sorted.get(q).get(1)>5) 
																		{
																			main[timeNext][i][1] = main[time][i][1] + sorted.get(q).get(0);
																			main[timeNext][i][4] = Math.abs(sorted.get(q).get(0))+Math.abs((5-main[timeNext][i][2])/2);	
																			main[timeNext][i][2] = 5;
																		}
																			for (int j11 = 1; j11 < i; j11++) 
																			{
																				if (main[timeNext][j11][1] == main[timeNext][i][1] && main[timeNext][j11][2] == main[timeNext][i][2] ) 
																					{
																						main[timeNext][i][1] = main[time][i][1];
																						main[timeNext][i][2] = main[time][i][2];
																						main[timeNext][i][4] = 0;
																					
																					}
																			}
																			
																		}
																	
																}
																
															
														}
												}						
									}
								
								else if (random <= ((  Math.exp(sorted.get(0).get(2)) +  Math.exp(sorted.get(1).get(2))) / probabSum))
									{
									int q=1;
									if(main[time][i][2] + sorted.get(q).get(1)>= 0 && main[time][i][2] + sorted.get(q).get(1)<6) 
									{
										main[timeNext][i][1] = main[time][i][1] + sorted.get(q).get(0);
										main[timeNext][i][2] = main[time][i][2] + sorted.get(q).get(1);
										main[timeNext][i][4] =Math.abs( sorted.get(q).get(0))+Math.abs(((sorted.get(q).get(1))/2));		
									}	
									else if(main[time][i][2] + sorted.get(q).get(1)< 0) 
									{
										main[timeNext][i][1] = main[time][i][1] + sorted.get(q).get(0);
										main[timeNext][i][2] =0;
										main[timeNext][i][4] = Math.abs(sorted.get(q).get(0))+Math.abs((0-main[timeNext][i][2])/2);	
									}
									else if (main[time][i][2] + sorted.get(q).get(1)>5) 
									{
										main[timeNext][i][1] = main[time][i][1] + sorted.get(q).get(0);
										main[timeNext][i][4] = Math.abs(sorted.get(q).get(0))+Math.abs((5-main[timeNext][i][2])/2);	
										main[timeNext][i][2] = 5;
									}
										
										for (int j = 1; j < i; j++) 
										{
											if (main[timeNext][j][1] == main[timeNext][i][1] && main[timeNext][j][2] == main[timeNext][i][2] ) 
												{
												q=2;
												if(main[time][i][2] + sorted.get(q).get(1)>= 0 && main[time][i][2] + sorted.get(q).get(1)<6) 
												{
													main[timeNext][i][1] = main[time][i][1] + sorted.get(q).get(0);
													main[timeNext][i][2] = main[time][i][2] + sorted.get(q).get(1);
													main[timeNext][i][4] = Math.abs(sorted.get(q).get(0))+Math.abs((sorted.get(q).get(1))/2);		
												}	
												else if(main[time][i][2] + sorted.get(q).get(1)< 0) 
												{
													main[timeNext][i][1] = main[time][i][1] + sorted.get(q).get(0);
													main[timeNext][i][2] =0;
													main[timeNext][i][4] = sorted.get(q).get(0)+Math.abs((0-main[timeNext][i][2])/2);	
												}
												else if (main[time][i][2] + sorted.get(q).get(1)>5) 
												{
													main[timeNext][i][1] = main[time][i][1] + sorted.get(q).get(0);
													main[timeNext][i][4] = Math.abs(sorted.get(q).get(0))+Math.abs((5-main[timeNext][i][2])/2);	
													main[timeNext][i][2] = 5;
												}	
													
													for (int j1 = 1; j1 < i; j1++) 
														{
															if (main[timeNext][j1][1] == main[timeNext][i][1] && main[timeNext][j1][2] == main[timeNext][i][2] ) 
																{
																q=0;
																if(main[time][i][2] + sorted.get(q).get(1)>= 0 && main[time][i][2] + sorted.get(q).get(1)<6) 
																{
																	main[timeNext][i][1] = main[time][i][1] + sorted.get(q).get(0);
																	main[timeNext][i][2] = main[time][i][2] + sorted.get(q).get(1);
																	main[timeNext][i][4] =  Math.abs(sorted.get(q).get(0))+ Math.abs((sorted.get(q).get(1))/2);		
																}	
																else if(main[time][i][2] + sorted.get(q).get(1)< 0) 
																{
																	main[timeNext][i][1] = main[time][i][1] + sorted.get(q).get(0);
																	main[timeNext][i][2] =0;
																	main[timeNext][i][4] =  Math.abs(sorted.get(q).get(0))+Math.abs((0-main[timeNext][i][2])/2);	
																}
																else if (main[time][i][2] + sorted.get(q).get(1)>5) 
																{
																	main[timeNext][i][1] = main[time][i][1] + sorted.get(q).get(0);
																	main[timeNext][i][4] =  Math.abs(sorted.get(q).get(0))+ Math.abs((5-main[timeNext][i][2])/2);	
																	main[timeNext][i][2] = 5;
																}
																	for (int j11 = 1; j11 < i; j11++) 
																	{
																		if (main[timeNext][j11][1] == main[timeNext][i][1] && main[timeNext][j11][2] == main[timeNext][i][2] ) 
																			{
																				main[timeNext][i][1] = main[time][i][1];
																				main[timeNext][i][2] = main[time][i][2];
																				main[timeNext][i][4] = 0;
																			
																			}
																	}
																	
																}
															
														}
														
													
												}
										}
									}
								
								else
									{
										main[timeNext][i][1] = main[time][i][1] + sorted.get(2).get(0);
										main[timeNext][i][2] = main[time][i][2] + sorted.get(2).get(1);
										main[timeNext][i][4] = Math.abs(sorted.get(2).get(0))+Math.abs((sorted.get(2).get(1))/2);						
									
										for (int j = 1; j < i; j++) 
										{ 
											if (main[timeNext][j][1] == main[timeNext][i][1] && main[timeNext][j][2] == main[timeNext][i][2] ) 
												{
													main[timeNext][i][1] = main[time][i][1] + sorted.get(1).get(0);
													main[timeNext][i][2] = main[time][i][2] + sorted.get(1).get(1);
													main[timeNext][i][4] = Math.abs(sorted.get(1).get(0))+Math.abs((sorted.get(1).get(1))/2);		
													
													for (int j1 = 1; j1 < i; j1++) 
														{ 
															if (main[timeNext][j1][1] == main[timeNext][i][1] && main[timeNext][j1][2] == main[timeNext][i][2] ) 
																{
																	main[timeNext][i][1] = main[time][i][1] + sorted.get(0).get(0);
																	main[timeNext][i][2] = main[time][i][2] + sorted.get(0).get(1);
																	main[timeNext][i][4] = Math.abs(sorted.get(0).get(0))+Math.abs((sorted.get(0).get(1))/2);
																
																	for (int j11 = 1; j11 < i; j11++) 
																	{
																		if (main[timeNext][j11][1] == main[timeNext][i][1] && main[timeNext][j11][2] == main[timeNext][i][2] ) 
																			{
																				main[timeNext][i][1] = main[time][i][1];
																				main[timeNext][i][2] = main[time][i][2];
																				main[timeNext][i][4] = 0;
																			
																			}
																	}
																	
																}
															
														}
														
													
												}
										}
									
									
									
									
									}
						}
						i++;
						if(main[time][i][0] == 0||i>ped_count){flag1 = 0;}	
						
					}
			
		
		averagespeed.close();
		
		return main;				
		
}
	
	//-----------------------------------------------------------------------------------------------------------------------------------
	
	public static void main(String[] args) throws IOException {		
	    // Main function to execute all the function in a sequential manner and to write outputs in desired format
		
		System.out.println("       working ... ");
		
		
		double[] p = {power(0,1), power(1,1), power(2,2), power(3,6), power(4,24), power(5,120)};		// Considering poisons distribution for pedestrian flow 
		
		double ywidth = Math.floor(width / (Math.sin (Math.toRadians(60)) * Tsize));			// calculating the units available in Y-direction
		
		double [][][] main = new double [(total_steps+2)][999][5];		// 3-D array to store pedestrian information. [time step to enter][no. of possible node for movement][ped_ID,X,Y,rand,Velocity]
		
		ArrayList<Integer> y_pos = new ArrayList<Integer>();				// Array list with possible Y units available for movements										
		
		for (int i = 1; i < ywidth; i++) {
            y_pos.add((int) (i));
        }		//

		PrintWriter pedata = new PrintWriter("ped_data.txt");  //steps,ped_count,y_pos
		PrintWriter maindata= new PrintWriter("main_data.txt");  // saving pedestrian data: ped_count, X, Y, random, V
		
		maindata.println("Step" + "\t"+"ped_id" + "\t"+ "X" + "\t"+ "Y" + "\t"+ "\t"+ "\t"+ "v");
	// Following DO-WHILE loop generates pedestrians 
		do {
			double rand=Math.random();			// generating random number
			Collections.shuffle(y_pos);	
			if (rand>p[1])
			{}
			
			else if (rand>p[2])
				{
					ped_count++;
					pedata.println(steps +"\t"+ ped_count+"\t"+ y_pos.get(0));
					main[steps][ped_count][0]=ped_count;
					main[steps][ped_count][1]=(y_pos.get(0)% 2)*0.5;
					main[steps][ped_count][2]=y_pos.get(0);
					main[steps][ped_count][3]=Math.random();
					main[steps][ped_count][4]=vmax;
					
					maindata.println(steps+"\t"+main[steps][ped_count][0] +"\t" +"\t"+main[steps][ped_count][1] +"\t"+main[steps][ped_count][2] +"\t"+ main[steps][ped_count][4]);
				}				
			else if (rand>p[3]) 
				{
					for (int j=0;j<2;j++)
					{
						ped_count++;
						pedata.println(steps +"\t"+ ped_count+"\t"+ y_pos.get(j));
						main[steps][ped_count][0]=ped_count;
						main[steps][ped_count][1]=(y_pos.get(j)% 2)*0.5;
						main[steps][ped_count][2]=y_pos.get(j);
						main[steps][ped_count][3]=Math.random();
						main[steps][ped_count][4]=vmax;
					
						maindata.println(steps+"\t"+main[steps][ped_count][0]+"\t" +"\t"+main[steps][ped_count][1] +"\t"+main[steps][ped_count][2] +"\t"+ main[steps][ped_count][4]);
					}
				}				
			else if (rand>p[4]) 
				{
					for (int j=0;j<3;j++)
					{
						ped_count++;
						pedata.println(steps +"\t"+ ped_count+"\t"+ y_pos.get(j));
						main[steps][ped_count][0]=ped_count;
						main[steps][ped_count][1]=(y_pos.get(j)% 2)*0.5;
						main[steps][ped_count][2]=y_pos.get(j);
						main[steps][ped_count][3]=Math.random();
						main[steps][ped_count][4]=vmax;
					
						maindata.println(steps+"\t"+main[steps][ped_count][0] +"\t"+"\t"+main[steps][ped_count][1] +"\t"+main[steps][ped_count][2] +"\t"+ main[steps][ped_count][4]);
					}
				}
			else if (rand>p[5]) 
				{
					for (int j=0;j<4;j++)
					{
						ped_count++;
						pedata.println(steps +"\t"+ ped_count+"\t"+ y_pos.get(j));
						main[steps][ped_count][0]=ped_count;
						main[steps][ped_count][1]=(y_pos.get(j)% 2)*0.5;
						main[steps][ped_count][2]=y_pos.get(j);
						main[steps][ped_count][3]=Math.random();
						main[steps][ped_count][4]=vmax;
					
						maindata.println(steps+"\t"+main[steps][ped_count][0] +"\t"+"\t"+main[steps][ped_count][1] +"\t"+main[steps][ped_count][2] +"\t"+ main[steps][ped_count][4]);
					}
				}				
			else 
				{
					for (int j=0;j<5;j++)
					{
						ped_count++;
						pedata.println(steps +"\t"+ ped_count+"\t"+ y_pos.get(j));
						main[steps][ped_count][0]=ped_count;
						main[steps][ped_count][1]=(y_pos.get(j)% 2)*0.5;
						main[steps][ped_count][2]=y_pos.get(j);
						main[steps][ped_count][3]=Math.random();
						main[steps][ped_count][4]=vmax;
					
						maindata.println(steps+"\t"+main[steps][ped_count][0] +"\t"+"\t"+main[steps][ped_count][1] +"\t"+main[steps][ped_count][2] +"\t"+ main[steps][ped_count][4]);
					}
				}

			steps = steps + time_step ;
		} while(steps<total_steps);
		System.out.println("Total pedestrians = " + ped_count );
		pedata.close();
		maindata.close();
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------

			
		double [][][] int_mat = new double [ped_count+1][999][3];
		double [][] [] bas_des_mat  = new double [1+ped_count][98][3];
		double [][][][] returnarray = new double [2][][][];	
		int[] interesting_pedid={1,6,7,8,9,10,11,12,13,14,15,16};

        int countPedData = 12;
        PrintWriter[] finalAllPedData = new PrintWriter[countPedData]; // saving pedestrian data: step, X, Y
        for(int i=0;i<countPedData;i++) {
            // finalPedData[1] should write to file finalmain_data2.txt.
            finalAllPedData[i] = new PrintWriter("finalmain_data" + String.valueOf(i+1) + ".txt");
            finalAllPedData[i].println("Data for pedestrian " + interesting_pedid[i]);
            finalAllPedData[i].println();
            finalAllPedData[i].println("Step" +"\t"+ "X" +"\t"+ "Y"+"\t"+ "Velocity");
        }
		PrintWriter finalpeddata  =   new PrintWriter("finalmain_data.txt");
		
		finalpeddata.println(ped_count);
		finalpeddata.println("Data for all pedestrian "+"\t"+"lane width"+"\t"+ width +"\t"+ "Flow"+"\t"+ ped_flow+"\t"+"and size is"+"\t"+Tsize+"\t"+"duration"+"\t"+total_steps+"\t"+"sec");
		finalpeddata.println();
		finalpeddata.println("Step" +"\t"+ "X" +"\t"+ "Y"+"\t"+ "Velocity");
		
		timeNext=1;		
		steps=0;
		
		do {
            int_mat = intmat(vmax, int_mat, ped_count);
            bas_des_mat = bdesmat(vmax, bas_des_mat, ped_count);


            for (int j = 1; j <= ped_count; j++) {
                for (int j2 = 0; j2 < 5; j2++) {
                    main[timeNext][j][j2] += main[steps][j][j2];
                }
            }
            for (int i = 0; i < countPedData; i++) {
                finalAllPedData[i].println(steps + "\t" + main[timeNext][interesting_pedid[i]][1] + "\t" + main[timeNext][interesting_pedid[i]][2] + "\t" + main[timeNext][interesting_pedid[i]][4]);
            }

            for (int i1 = 1; i1 < ped_count; i1++) {
                finalpeddata.println(steps + "\t" + main[timeNext][i1][1] + "\t" + main[timeNext][i1][2] + "\t" + main[timeNext][i1][4]);
            }

            returnarray = negeff(returnarray, ped_count, total_steps, bas_des_mat, int_mat, main, steps, vmax, ywidth);

            bas_des_mat = returnarray[0];
            main = returnarray[1];

            if (ped_count >= 1) {
                if (steps >= 3 && steps <= 7) {
                    for (int j = 0; j <= base_des_count; j++) {
                        bas_des_mat[1][j][2] = -100;
                    }
                }
            }

            if (ped_count >= 10) {
                if (steps >= 10 && steps <= 12) {
                    for (int j = 0; j <= base_des_count; j++) {
                        bas_des_mat[10][j][2] = -100;
                    }
                }
            }

            if (ped_count >= 16) {
                if (steps >= 12) {
                    for (int j = 0; j <= base_des_count; j++) {
                        bas_des_mat[16][j][2] = -100;
                    }
                }
            }
            main = updatePosition(ywidth, ped_count, main, timeNext, steps, bas_des_mat);

            if ((total_steps - steps) % 10 == 0) {
                System.out.println(total_steps - steps);
            }


            steps++;
            timeNext++;
        } while(steps <= total_steps);

        for (int i=0;i<countPedData;i++) {
            finalAllPedData[i].close();
        }
        finalpeddata.close();

        System.out.println( "\n" + "     "  + "DONE");
        System.out.println( "\n" + "    " + "\\"+"(^.^)/");
    }

}

// That's all FOLKS  :)
