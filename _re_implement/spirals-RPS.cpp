#include <iostream>   //  image of matlab   // single spirals  copy right for Luo-Luo Jiang
#include <cmath>      
#include <cstdio>       
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <sstream>
using namespace std;
// warned: every arry have to be enlarge one 
FILE *fp2;

int main()
{  //------------------------------------------------------  //main program
//----------------------------------------
 int z1,z2,z3,z4,a,q,r,m;
 double f,f1,f2;
 long int size, mid_size, area,ix,iy,monte,i,j,i1,i2,j1,j2,cyc1,cyc2;
 long int add[10000],minus[10000];
 long int ceni1,ceni2,ceni3,ceni4,ceni5,ceni6,cenj1,cenj2,cenj3,cenj4,cenj5,cenj6,radius,radius2;
 double dou_mob[4],dou_d, dou_rep,dou_sel,dou_time,pai; 
 double dou_r1=0.0,dou_r2=0.0, dou_r3=0.0,dou_r4=0.0,dou_r5=0.0,dou_r6=0.0,dou_r7=0.0,dou_r8=0.0,dou_r9=0.0,dou_rat1=0.0,dou_rat2=0.0; 
 long int inset=0,seed=0,unit0=0,unit1=0,unit2=0,unit3=0,unit4=0;

    dou_d=0.00001;


    int*pmp1,*species[10001];  //++ 动态分配二维数组
    pmp1=(int*)new int [10001][10001]; //++
    for(int i=1;i<=10000;i++)          //++
	{   species[i]=pmp1+i*10000;   }  //++
  //------------------------------------------------------------
a=16807,q=127773,m=2147483647,r=2836;//这些都是随机数产生器的常量
     //++++++++++++++++++++++++++++++++++++++
     srand( (unsigned)time( NULL ) );//初始化随机数
     z3=rand();
     z1=rand();
     z4=rand();
     z2=z1+z3*97+z4*3;//初始化随机数
	//+++++++++++++++++++++++++++++++++++++++
//-------------------------------------------------

    size=256;   //******************************************size
    area=size*size;
    mid_size=size/2;
	radius=20;
    radius2=100;

    dou_rep=1.0;

    dou_sel=1.0;     //predation density
    pai=3.1415926;



 dou_time=0.0;
    ceni1=mid_size+radius2;
    cenj1=mid_size;
    ceni2=mid_size+(int)(cos((120.0/180.0)*pai)*radius2);
    cenj2=mid_size+(int)(sin((120.0/180.0)*pai)*radius2);
    ceni3=mid_size+(int)(cos((240.0/180.0)*pai)*radius2);
    cenj3=mid_size+(int)(sin((240.0/180.0)*pai)*radius2);
    


//  cyc2
//for (cyc2=1;cyc2<=3;cyc2++)
//{
   //  dou_d=0.0001+(cyc2-1)*0.0001;
//  cyc1
for (cyc1=1;cyc1<=10;cyc1++)
{
    for( i=1;i<=3;i++)
	{
     dou_mob[i]=dou_d*area*2.0;
	}
	for ( ix=1;ix<=size;ix++)
	{
		add[ix]=ix+1;
		minus[ix]=ix-1;
	}
	for (ix=1;ix<=size;ix++)
	{
		for (iy=1;iy<=size;iy++)
		{species[ix][iy]=0;}
	}
 //------------------------------------------
  

     //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	//   1
     unit1=ceni1-radius;
	 unit2=ceni1+radius;
	 unit3=cenj1-radius;
	 unit4=cenj1+radius;
     for (ix=unit1;ix<=unit2;ix++)
	 {   
	     for (iy=unit3;iy<=unit4;iy++)
		 {
			 dou_r1=fabs(ix-ceni1);
		     dou_r2=fabs(iy-cenj1);
             dou_r3=pow((double)radius,2);
             dou_r4=pow((double)dou_r1,2)+pow((double)dou_r2,2);
		     if(dou_r4<dou_r3)
			 { species[ix][iy]=1;}
                             
		 }
	 }
	 //  2
     unit1=ceni2-radius;
	 unit2=ceni2+radius;
	 unit3=cenj2-radius;
	 unit4=cenj2+radius;
     for (ix=unit1;ix<=unit2;ix++)
	 {
	     for (iy=unit3;iy<=unit4;iy++)
		 {
			 dou_r1=fabs(ix-ceni2);
		     dou_r2=fabs(iy-cenj2);
             dou_r3=pow((double)radius,2);
             dou_r4=pow((double)dou_r1,2)+pow((double)dou_r2,2);
		     if(dou_r4<dou_r3)
		           
			 {   species[ix][iy]=2;}
                      
		 }
	 }
	 // 3
     unit1=ceni3-radius;
	 unit2=ceni3+radius;
	 unit3=cenj3-radius;
	 unit4=cenj3+radius;
     for (ix=unit1;ix<=unit2;ix++)
	 {
	     for (iy=unit3;iy<=unit4;iy++)
		 {
			 dou_r1=fabs(ix-ceni3);
		     dou_r2=fabs(iy-cenj3);
             dou_r3=pow((double)radius,2);
             dou_r4=pow((double)dou_r1,2)+pow((double)dou_r2,2);
		     if(dou_r4<dou_r3)
			 {species[ix][iy]=3;}
		 }
	 }
	 





    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@







  inset=0;
  for (int step_mc=1;step_mc<=area;step_mc++)
  {                                               //step_mc    
 
	  if(dou_time>1000000)
	  { dou_time=0.0;
	    //++++++++++++++++++++++++++++++++++++++
        srand( (unsigned)time( NULL ) );//初始化随机数
        z3=rand();
        z1=rand();
        z4=rand();
        z2=z1+z3*97+z4*3;//初始化随机数
	    //+++++++++++++++++++++++++++++++++++++++
	  }

	  for ( monte=1;monte<=area;monte++)
	  {   //    one step of Monte Carlo
		   // +++++++++++++++++++++++++++++++++++++
           z2=a*(z2%q)-r*(z2/q);
           if(z2<0) z2=z2+m;
	       else z2=z2; 
	       f1=z2/2147483647.0;//[0-1) randon number
		   dou_time=dou_time+1.0;
           //++++++++++++++++++++++++++++++++++++++
           i1=(int)(f1*size)+1;
		   // +++++++++++++++++++++++++++++++++++++
           z2=a*(z2%q)-r*(z2/q);
           if(z2<0) z2=z2+m;
	       else z2=z2; 
	       f2=z2/2147483647.0;//[0-1) randon number
		   dou_time=dou_time+1.0;
           //++++++++++++++++++++++++++++++++++++++
		   //f2=rand()/RAND_MAX;
		   j1=(int)(f2*size)+1;
		   unit0=species[i1][j1];
           if( unit0>0)
		   {    //non vacation
                
			   // +++++++++++++++++++++++++++++++++++++

  label_1:     z2=a*(z2%q)-r*(z2/q);
               if(z2<0) z2=z2+m;
	           else z2=z2; 
	           f=z2/2147483647.0;//[0-1) randon number
			   dou_time=dou_time+1.0;
               //++++++++++++++++++++++++++++++++++++++
           
               seed=(int)(f*4.0)+1;
               if(seed==1)
			   {
                  i2=add[i1];
			      j2=j1;
			      if(i2>size||i2<1||j2>size||j2<1)
				  {goto label_1;}
			   }
			   else if(seed==2)
			   {
				   i2=i1;
				   j2=add[j1];
				   if(i2>size||i2<1||j2>size||j2<1)
				   {goto label_1;}
			   }
			   else if(seed==3)
			   {
				   i2=minus[i1];
				   j2=j1;
				   if(i2>size||i2<1||j2>size||j2<1)
				   {goto label_1;}
			   }
               else if(seed==4)
			   {
				   i2=i1;
				   j2=minus[j1];
				   if(i2>size||i2<1||j2>size||j2<1)
				   {goto label_1;}
			   }
               
               unit1=species[i1][j1];
			   unit2=species[i2][j2];
			   unit3=unit1-unit2;
			   dou_r1=dou_mob[unit1];
			   // +++++++++++++++++++++++++++++++++++++
               z2=a*(z2%q)-r*(z2/q);
               if(z2<0) z2=z2+m;
	           else z2=z2; 
	           f=z2/2147483647.0;//[0-1) randon number
			   dou_time=dou_time+1.0;
               //++++++++++++++++++++++++++++++++++++++
               dou_r2=f*(dou_rep+dou_sel+dou_r1);
			   dou_r3=dou_sel;
			   dou_r4=dou_sel+dou_rep;
               //   select
			   if(dou_r2<=dou_r3)
			   {      
				   if(unit2>0)
				   {
				       if(unit3==-1||unit3==2)
					   {species[i2][j2]=0;}
				   }
			   
                }
			   else if(dou_r2>dou_r3&&dou_r2<=dou_r4)
			   {
			       if(unit2==0)
				   {
				   species[i2][j2]=unit1;
				   }
			   }
			   else if(dou_r2>dou_r4)
			   {
                    species[i2][j2]=unit1;
                    species[i1][j1]=unit2;
			   }
		  }    //end non vacation
	  }   //    end one step of Monte Carlo




}    //---------------------------------------------------   //  end   step_mc

         char filename[20];   
   	     sprintf(filename,"%lf_%d.txt",dou_d,cyc1); 
//	     cout<<filename<<endl;
        
         fp2=fopen(filename,"w");
          for (ix=1;ix<=size;ix++)
	        {
		        for (iy=1;iy<=size;iy++)
		        {
                    unit1=species[ix][iy];
                    fprintf(fp2,"%d,%d,%d\n",ix,iy,unit1); 
                }
            }
         


}       // end of cyc1
//}       // end of cyc2

  delete[]pmp1;      //++  //动态分配二维数组
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  return 0;
}

