//有限元法平面问题程序
//#include"conio.h"
#include<stdio.h>
#include<math.h>


//#define  NZ 14
//#define NPF 12
#define nl 30   //板长方向划分个数
#define nb 10    //板高方向划分个数
#define  NE   (nl*nb)
#define  NJ   ((nl+1)*(nb+1))
#define  NJ2  (2*(nl+1)*(nb+1))
#define  DD   ((nb+3)*2) 
#define NPF   nl   //上边界均布荷载与板长方向划分个数一样
#define  NZ  ((nb+1)*2)    //支撑数



double  E00=2E10;    //定义板最下边界弹性模量
double  Ex=0;      //   E0=E00+Ex*X+Ey*Y+Exy*x*y  在X，Y方向增大系数
double  Ey=0;
double  Exy=0;
int LXM=0;
double junbu=-100000.0;
//double E0;
double MU=0.167;
double LOU=0.0;
double TE=1.0;
double ll=9.0;    //板长
double bb=3.0;    //板高
//int nl=12;    //板长方向划分个数
//int nb=6;    //板高方向划分个数


double E0[NE];                      //变阶的弹性模量         在mesh()中定义数组数值
double AJZ[NJ+1][3];
int JM[NE+1][5];                   //节点码  在mesh中定义 //={{0,0,0,0,0},{0,2,3,4,1},{0,3,5,6,4}};
double pf[NPF+1][7]; //={{0,0,0,0,0,0,0},{0,-100000.0,0.75,6.0,7.0,14.0,1.0},{0,-100000.0,0.75,12.0,14.0,21.0,1.0},{0,-100000.0,0.75,18.0,21.0,28.0,1.0},{0,-100000.0,0.75,24.0,28.0,35.0,1.0},{0,-100000.0,0.75,30.0,35.0,42.0,1.0},{0,-100000.0,0.75,36.0,42.0,49.0,1.0},{0,-100000.0,0.75,42.0,49.0,56.0,1.0},{0,-100000.0,0.75,48.0,56.0,63.0,1.0},{0,-100000.0,0.75,54.0,63.0,70.0,1.0},{0,-100000.0,0.75,60.0,70.0,77.0,1.0},{0,-100000.0,0.75,66.0,77.0,84.0,1.0},{0,-100000.0,0.75,72.0,84.0,91.0,1.0}};    //{0,荷载值，荷载位置、长度，单元，节点码1，节点码2，荷载类型}
int NZC[NZ+1];//={0,1,3,5,7,9,11,13,15,17,19,21,23,25,626};
//double PJ[NPJ+1][2+1]={{0,0,0}};
double AE,KZ[NJ2+1][DD+1],P[NJ2+1],S[3+1][8+1],KE[8+1][8+1];
double f0[3],SSS[3+1][NJ+1];                 //SSS[]存放每个节点的应力和
int IE,JE,ME;
int TJJD[NJ+1];
/*
void mesh()
{   
	
	double ccl;  //板长方向单元长度
	double ccb;   //板高方向单元长度

	int i,j,k;
	 
	ccl=ll/nl;    //单元x方向长度
	ccb=bb/nb;     //单元y方向长度
	for(i=0;i<=NJ;i++)
	{	
		for(j=0;j<=3;j++)
		AJZ[i][j]=0.0;
	}
	for(i=1;i<=(nl+1);i++)      //从左下角数  第i列 第k个单元  x坐标
	
	{
		j=(i-1)*ccl;
		for(k=1;k<=nb+1;k++)
		{
			AJZ[(nb+1)*(i-1)+k][1]=j;

		}
				
		
	}

}
*/



//DUGD()  函数：《功能：可以计算面积AE，生成S=D*B矩阵，KE矩阵》         函数声明 
void mesh();

void DUGD(int,int,double,double);
void gdnl(int);


//主函数
void main()

{
	int mm,a1,b1,NJ1,k,IM,IN,jn,m,i,j,z,J0,ii,jj,h,dh,E,l,zl,dl,hz;
	int jdh;     //取出节点号存放到jdh
	double c,SIG1,SIG2,SIG3,PYL,RYL,MAYL,MIYL,CETA,AAXX,BBYY;
	double WY[6+1],YL[3+1];
	FILE *fp;
    	
	system("mode con cols=200000  lines=600");     //输出DOS窗口口之函数   列数和行数
	system("color F0");                           //用 system("color 0A"); 其中color后面的0是背景色代号，A是前景色代号。各颜色代码如下：0=黑色 1=蓝色 2=绿色 3=湖蓝色 4=红色 5=紫色 6=黄色 7=白色 8=灰色 9=淡蓝色 A=淡绿色 B=淡浅绿色 C=淡红色 D=淡紫色 E=淡黄色 F=亮白色
//	window(10,10,40,11);
//	textcolor(BLACK);
//	textbackground(WHITE);
//	fp=fopen("shuju.txt","w");
//	FILE* fp2=fopen("shuju.txt","r"); 
	
     
	fp=fopen("shuju.txt","w");
	mesh();                                       //自动划分单元

	for(i=1;i<=(NJ);i++)                           //输出节点坐标
	{
//		printf("\n第%d个节点坐标: \t\t %f\t%f\t",i, AJZ[i][1], AJZ[i][2]);
		fprintf(fp,"\n第%-2d个节点坐标: \t\t %-4f\t%-4f\t",i, AJZ[i][1], AJZ[i][2]);
	}
	for(i=1;i<=(NE);i++)                           //输出单元节点码
 		
	{
//		printf("\n第%d个单元编码:\t%d\t%d\t%d\t%d\t",i, JM[i][1], JM[i][2], JM[i][3],JM[i][4]);
		fprintf(fp,"\n第%d个单元编码:\t%-4d\t%-4d\t%-4d\t%-4d\t",i, JM[i][1], JM[i][2], JM[i][3],JM[i][4]);
		
	
	}
	for(i=1;i<=NE;i++)                            //输出输出变阶后各单元E0
	{
//		printf("\n第%d个单元弹性模量E: \t\t%f\t",i,E0[i]);
		fprintf(fp,"\n第%-2d个单元弹性模量: \t\t %-4f",i,E0[i]);
	}

//	printf("\n\n\n\n\n");
	fprintf(fp,"\n\n");

	if(LXM!=0)
	{
		for(i=1;i<=NE;i++)
		{
			E0[i]=E0[i]/(1.0-MU*MU);
		}
		MU=MU/(1.0-MU);

	}

	//调用函数DUGD(),生成KZ矩阵
	for(i=0;i<=NJ2;i++)
	{
		for(j=0;j<=DD;j++)
		{
			KZ[i][j]=0.0;
		}
	}																//KZ清零
	
	for(E=1;E<=NE;E++)
	{
		DUGD(E,3,0,0);
			//输出第E个ke矩阵
			fprintf(fp,"\n第%d个单元的KE矩阵\n",E);
			for(i=1;i<=8;i++)
			{
				for(j=1;j<=8;j++)
					fprintf(fp,"%-18.6f\t",KE[i][j]);
				fprintf(fp,"\n");
			}
			
		for(i=1;i<=4;i++)
		{
			for(ii=1;ii<=2;ii++)
			{
				h=2*(i-1)+ii;
				dh=2*(JM[E][i]-1)+ii;
				for(j=1;j<=4;j++)
				{
					for(jj=1;jj<=2;jj++)
					{
						l=2*(j-1)+jj;
						zl=2*(JM[E][j]-1)+jj;
						dl=zl-dh+1;
						if(dl>0)
							KZ[dh][dl]=KZ[dh][dl]+KE[h][l];
					}
				}
			}
		}
	}
	//矩阵输出
    fprintf(fp,"\n\nKZ矩阵如下\n");
	for(i=1;i<=NJ2;i++)
	{
		for(j=1;j<=DD;j++)
		fprintf(fp,"%-18.6f\t",KZ[i][j]);
		fprintf(fp,"\n");

	}
	fprintf(fp,"\n");

	//***********形成P矩阵*****************
	for(i=1;i<=NJ2;i++)
		P[i]=0.0;
	if(NPF>0)
	{
		for(i=1;i<=NPF;i++)
		{
			hz=i;
			gdnl(hz);
//			e=(int)pf[hz][3];

//			zb(e);
/*			for(j=1;j<=6;j++)
			{
				pe[j]=0.0;
				for(k=1;k<=6;k++)
				{
					pe[j]=pe[j]-t[k][j]*f0[k];    //用的是转秩矩阵     zb()函数出来的是直接矩阵   在这先让列不动行动  相当于乘的转秩
				}
			}   */
//			a1=JM[e][1];
//			b1=JM[e][2];
			a1=(int)pf[hz][4];
	        b1=(int)pf[hz][5];
		
			P[2*a1-1]=P[2*a1-1]+f0[1];
			P[2*a1]=P[2*a1]+f0[2];
			P[2*b1-1]=P[2*b1-1]+f0[1];
			P[2*b1]=P[2*b1]+f0[2];

		}
	}
	//***************************输出p矩阵**************************
//	 printf("\n\nP矩阵如下\n");
	 fprintf(fp,"\n\nP矩阵如下\n");
	for(i=1;i<=NJ2;i++)
	{
//		printf("%d  %f\t",i,P[i]);
//		printf("\n");
		fprintf(fp,"\n%-4d  %-18.6f\t",i,P[i]);
		fprintf(fp,"\n");

	}

/*	if(NJ>0)
	{

		for(i=1;i<=NPJ;i++)
		{
			j=(int)PJ[i][2];
			P[j]=PJ[i][1]+p[j];

		}
	}*/
/*	if(LOU>0)
	{
		for(E=1;E<=NE;E++)
		{
			DUGD(E,1);
			PE=-LOU*(AE)*TE/3;
			P[2*IE]=P[2*IE]+PE;
			P[2*JE]=P[2*JE]+PE;
			P[2*ME]=P[2*ME]+PE;

		}
	}        */

	//*******引入边界条件*********
	for(i=1;i<=NZ;i++)
	{

		z=NZC[i];
		KZ[z][1]=1.0;
		for(j=2;j<=DD;j++)
			KZ[z][j]=0.0;
		if(z!=1)
		{
			if(z>DD)
				J0=DD;
			else J0=z;
			for(j=2;j<=J0;j++)
				KZ[z-j+1][j]=0.0;

		}
		P[z]=0.0;
	}
	//矩阵输出
    fprintf(fp,"\n\nKZ矩阵如下(变化后  对角线置一之后)\n");
	for(i=1;i<=NJ2;i++)
	{
		for(j=1;j<=DD;j++)
		//	printf("%f\t",KZ[i][j]);
		//printf("\n");
		fprintf(fp,"%-18.6f\t",KZ[i][j]);
		fprintf(fp,"\n");

	}
	//printf("\n");
	fprintf(fp,"\n");

	//高斯消元法
	NJ1=NJ2-1;
	for(k=1;k<=NJ1;k++)
	{
		if(NJ2>k+DD-1)
			IM=k+DD-1;
		else IM=NJ2;

		IN=k+1;
		for(i=IN;i<=IM;i++)
		{
			l=i-k+1;
			c=KZ[k][l]/KZ[k][1];
			jn=DD-l+1;
			for(j=1;j<=jn;j++)
			{
				m=j+i-k;
				KZ[i][j]=KZ[i][j]-c*KZ[k][m];         //有问题

			}
			P[i]=P[i]-c*P[k];

		}

	}

	//矩阵回代
	P[NJ2]=P[NJ2]/KZ[NJ2][1];
	for(i=NJ1;i>=1;i--)
	{

		if(DD>NJ2-i+1)J0=NJ2-i+1;
		else J0=DD;
		for(j=2;j<=J0;j++)
		{
			h=j+i-1;
			P[i]=P[i]-KZ[i][j]*P[h];
		}
		
		
			P[i]=(P[i]/KZ[i][1]);
		
	}
	fprintf(fp,"\n输出位移矩阵\n");
	fprintf(fp,"位移号---------位移\n");
	for(i=1;i<=NJ2;i++)
	{
	fprintf(fp,"%-4d---------      %-18.6f\n",i,P[i]);
	}
//	printf("\n\n\n\n\n");
//	printf("节点位移如下：\n");
//	printf("JD                            U                              V\n");
	
	fprintf(fp,"\n\n\n\n\n");
	fprintf(fp,"节点位移如下：\n");
	
	fprintf(fp,"JD                            U                              V\n");
	for(i=1;i<=NJ;i++)
	{
//		printf("%d                          %-9.15f                      %-9.15f\n",i,P[2*i-1],P[2*i]);

		fprintf(fp,"%-4d                          %-18.6f                      %-18.6f\n",i,P[2*i-1],P[2*i]);

	}

for(i=1;i<=3;i++)
{	
	for(j=1;j<=NJ+1;j++)
	
       SSS[i][j]=0.0;
}
			


	//**********求单元应力和主应力*************
	for(E=1;E<=NE;E++)
	{
		for(mm=1;mm<=4;mm++)
		{
				if(mm==1)
				{
					AAXX=-1.0;
		   			BBYY=-1.0;
				}
				if(mm==2)
				{
					AAXX=1.0;
					BBYY=-1.0;
				}
				if(mm==3)
				{
					AAXX=1.0;
					BBYY=1.0;
				}
				if(mm==4)
				{
					AAXX=-1.0;
					BBYY=1.0;
				}
				/*if(mm==5)
				{
					
					AAXX=0.0;          //求中间点应力
					BBYY=0.0;

				}*/
				

				DUGD(E,2,AAXX,BBYY);

				
				for(i=1;i<=4;i++)
				{
					for(j=1;j<=2;j++)
					{
						h=2*(i-1)+j;
						dh=2*(JM[E][i]-1)+j;     //每个单元8个位移拿出来存入wy数组
						WY[h]=P[dh];
					}
				}
				for(i=1;i<=3;i++)
				{
					YL[i]=0;
					for(j=1;j<=8;j++)
						YL[i]=YL[i]+S[i][j]*WY[j];
				}
				SIG1=YL[1];
				SIG2=YL[2];
				SIG3=YL[3];
				PYL=(SIG1+SIG2)/2;
				RYL=sqrt(pow((SIG1-SIG2)/2.0,2)+pow(SIG3,2));
				MAYL=PYL+RYL;
				MIYL=PYL-RYL;
				if(SIG2==MIYL)
					CETA=0;
				else
					CETA=90-57.29578*atan2(SIG3,(SIG2-MIYL));
				fprintf(fp,"\n\n\n\n\n");
				fprintf(fp,"%-4d单元应力如下\n",E);
					if(mm==1)
				{
					//AAXX=-1.0;
					//BBYY=-1.0;
					    fprintf(fp,"\n应力矩阵为");
						for(i=1;i<=3;i++)
				
						{
								
								fprintf(fp,"\n%-18.6f    %-18.6f   %-18.6f    %-18.6f    %-18.6f    %-18.6f    %-18.6f    %-18.6f",S[i][1],S[i][2],S[i][3],S[i][4],S[i][5],S[i][6],S[i][7],S[i][8]);
						}
					fprintf(fp,"\n单元左下角应力i节点\n");
					jdh=JM[E][1];
					
					SSS[1][jdh]=SIG1+SSS[1][jdh];
					SSS[2][jdh]=SIG2+SSS[2][jdh];
					SSS[3][jdh]=SIG3+SSS[3][jdh];

				}
				if(mm==2)
				{
					//AAXX=1.0;
					//BBYY=-1.0;
					fprintf(fp,"\n应力矩阵为");
					for(i=1;i<=3;i++)
				
						{
								
								fprintf(fp,"\n%-18.6f    %-18.6f   %-18.6f    %-18.6f    %-18.6f    %-18.6f    %-18.6f    %-18.6f",S[i][1],S[i][2],S[i][3],S[i][4],S[i][5],S[i][6],S[i][7],S[i][8]);
						}
					fprintf(fp,"\n单元右下角应力j节点\n");
					jdh=JM[E][2];
					SSS[1][jdh]=SIG1+SSS[1][jdh];
					SSS[2][jdh]=SIG2+SSS[2][jdh];
					SSS[3][jdh]=SIG3+SSS[3][jdh];



				}
				if(mm==3)
				{
					//AAXX=1.0;
					//BBYY=1.0;
					fprintf(fp,"\n应力矩阵为");

					for(i=1;i<=3;i++)
				
						{
								fprintf(fp,"\n%-18.6f    %-18.6f   %-18.6f    %-18.6f    %-18.6f    %-18.6f    %-18.6f    %-18.6f",S[i][1],S[i][2],S[i][3],S[i][4],S[i][5],S[i][6],S[i][7],S[i][8]);
						}
					fprintf(fp,"\n单元右上角应力m节点\n");
					jdh=JM[E][3];
					SSS[1][jdh]=SIG1+SSS[1][jdh];
					SSS[2][jdh]=SIG2+SSS[2][jdh];
					SSS[3][jdh]=SIG3+SSS[3][jdh];


				}
				if(mm==4)
				{
					//AAXX=-1.0;
					//BBYY=1.0;
					fprintf(fp,"\n应力矩阵为");
					for(i=1;i<=3;i++)
				
						{
							
								fprintf(fp,"\n%-18.6f    %-18.6f   %-18.6f    %-18.6f    %-18.6f    %-18.6f    %-18.6f    %-18.6f",S[i][1],S[i][2],S[i][3],S[i][4],S[i][5],S[i][6],S[i][7],S[i][8]);
						}
					fprintf(fp,"\n单元左上角p节点应力\n");
					jdh=JM[E][4];
					SSS[1][jdh]=SIG1+SSS[1][jdh];
					SSS[2][jdh]=SIG2+SSS[2][jdh];
					SSS[3][jdh]=SIG3+SSS[3][jdh];


				}
			/*	if(mm==5)
				{
					
					printf("单元中间点应力AAXX=%f,\t,BBYY=%f\n",AAXX,BBYY);
					
				

				}*/

				
				fprintf(fp,"E=%-4d\n",E);
			
				fprintf(fp,"sx=%-9.6f   sy=%-9.6f  tou=%-9.6f\n",SIG1,SIG2,SIG3);
				fprintf(fp,"s1=%-9.6f   s2=%-9.6f  theta=%-9.6f\n",MAYL,MIYL,CETA);

		}				
	}
	//*****************绕点平均法求各节点应力********************//
	for(i=1;i<=NJ;i++)
		TJJD[i]=0;
	for(E=1;E<=NE;E++)           //对单元循环
	{
		for(i=1;i<=4;i++)
		
		{
			  jdh=JM[E][i];                    //对每个单元中的节点循环累计统计节点在单元出现的次数
			  TJJD[jdh]=TJJD[jdh]+1;

		}

	}
	for(i=1;i<=NJ;i++)
	{
		fprintf(fp,"\n节点个数为TJJD[%-4d]=%-4d\n",i,TJJD[i]);
	}
	for(i=1;i<=NJ;i++)
	{
		
		SSS[1][i]=SSS[1][i]/TJJD[i];
		SSS[2][i]=SSS[2][i]/TJJD[i];
		SSS[3][i]=SSS[3][i]/TJJD[i];

		fprintf(fp,"绕点平均法    节点%-4d      sx应力为%-18.6f     sy应力为%-18.6f       tou应力为%-18.6f\n",i,SSS[1][i],SSS[2][i],SSS[3][i]);
	}

	fclose(fp);
	
	
}

//DUGD()函数；《功能：可以计算面积EA，生成S=D*B矩阵，KE矩阵》
void DUGD(int E,int ASK,double AX,double AY )
{
	

	int i,j;
	double U1,U2,ALF,BAT,HT,R1,H1,XIJ,YIJ,XJM,YJM,A2,B2,A1,B1;
	double SC,BBB,AAA;
	IE=JM[E][1];
	JE=JM[E][2];
	ME=JM[E][3];
//	PM=JM[E][4];
	
	XIJ=fabs(AJZ[JE][1]-AJZ[IE][1]);
	YIJ=fabs(AJZ[JE][2]-AJZ[IE][2]);
	XJM=fabs(AJZ[ME][1]-AJZ[JE][1]);
	YJM=fabs(AJZ[ME][2]-AJZ[JE][2]);
//	CM=AJZ[JE][1]-AJZ[IE][1];
//	BM=AJZ[IE][2]-AJZ[JE][2];
//	CJ=AJZ[IE][1]-AJZ[ME][1];
//	BJ=AJZ[ME][2]-AJZ[IE][2];
	A2=sqrt(XIJ*XIJ+XJM*XJM);
	B2=sqrt(YIJ*YIJ+YJM*YJM);
	AE=A2*B2;
	A1=A2/2;
	B1=B2/2;
 //对KE进行分块赋值

	H1=1/(1-MU*MU);
	R1=(1-MU)/2;
	HT=H1*E0[E]*TE;           //BT为厚度   
	U1=0.125*(1+MU)*HT;          
	U2=0.125*(1-3*MU)*HT;
	ALF=A1/B1/3.0;
	BAT=B1/A1/3.0;
//	V1=(BAT+R1*ALF)*HT;
	for(i=1;i<=8;i++)
	{
		for(j=1;j<=8;j++)
		{
			{
				KE[i][j]=0.0;

			}
		}

	}
	

	



	if(ASK==3)
	{
		
	KE[1][1]=(BAT+R1*ALF)*HT;
	KE[2][1]=U1;
	KE[2][2]=(ALF+R1*BAT)*HT;
	KE[3][1]=(-BAT+R1*ALF/2)*HT;
	KE[3][2]=U2;
	KE[3][3]=(BAT+R1*ALF)*HT;
	KE[4][1]=-U2;
	KE[4][2]=(ALF/2-R1*BAT)*HT;
	KE[4][3]=-U1;
	KE[4][4]=(ALF+R1*BAT)*HT;
	KE[5][1]=HT*(-BAT-R1*ALF)/2;
	KE[5][2]=-U1;
	KE[5][3]=(BAT/2-R1*ALF)*HT;
	KE[5][4]=U2;
	KE[5][5]=(BAT+R1*ALF)*HT;
	KE[6][1]=-U1;
	KE[6][2]=(-ALF-R1*BAT)/2*HT;
	KE[6][3]=-U2;
	KE[6][4]=(-ALF+(R1*BAT)/2)*HT;
	KE[6][5]=U1;	
	KE[6][6]=(ALF+R1*BAT)*HT;
	KE[7][1]=(BAT/2-R1*ALF)*HT;
	KE[7][2]=-U2;
	KE[7][3]=HT*(-BAT-R1*ALF)/2;
	KE[7][4]=U1;
	KE[7][5]=(-BAT+R1*ALF/2)*HT;
	KE[7][6]=U2;
	KE[7][7]=(BAT+R1*ALF)*HT;
	KE[8][1]=U2;
	KE[8][2]=(-ALF+R1*BAT/2)*HT;
 	KE[8][3]=U1;
	KE[8][4]=((-ALF-R1*BAT)/2)*HT;
	KE[8][5]=-U2;
	KE[8][6]=(ALF/2-R1*BAT)*HT;
	KE[8][7]=-U1;
	KE[8][8]=(ALF+R1*BAT)*HT;
	
	KE[1][2]=KE[2][1];
	KE[1][3]=KE[3][1];
	KE[1][4]=KE[4][1];
	KE[1][5]=KE[5][1];
	KE[1][6]=KE[6][1];
	KE[1][7]=KE[7][1];
	KE[1][8]=KE[8][1];
	KE[2][3]=KE[3][2];
	KE[2][4]=KE[4][2];
	KE[2][5]=KE[5][2];
	KE[2][6]=KE[6][2];
	KE[2][7]=KE[7][2];
	KE[2][8]=KE[8][2];
	KE[3][4]=KE[4][3];
	KE[3][5]=KE[5][3];
	KE[3][6]=KE[6][3];
	KE[3][7]=KE[7][3];
	KE[3][8]=KE[8][3];
	KE[4][5]=KE[5][4];
	KE[4][6]=KE[6][4];
	KE[4][7]=KE[7][4];
	KE[4][8]=KE[8][4];
	KE[5][6]=KE[6][5];
	KE[5][7]=KE[7][5];
	KE[5][8]=KE[8][5];
	KE[6][7]=KE[7][6];
	KE[6][8]=KE[8][6];
	KE[7][8]=KE[8][7];
	
	}
if(ASK==2)
{
	SC=E0[E]/4/A1/B1/(1-MU*MU);
	AAA=AX*A1;
	BBB=AY*B1;
	
	S[1][1]=SC*(-(B1-BBB));
	S[1][2]=SC*(-MU*(A1-AAA));
	S[1][3]=SC*(B1-BBB);
	S[1][4]=SC*(-MU*(A1+AAA));
	S[1][5]=SC*(B1+BBB);
	S[1][6]=SC*(MU*(A1+AAA));
	S[1][7]=SC*(-(B1+BBB));
	S[1][8]=SC*(MU*(A1-AAA));
	S[2][1]=SC*(-MU*(B1-BBB));
	S[2][2]=SC*(-(A1-AAA));
	S[2][3]=SC*(MU*(B1-BBB));
	S[2][4]=SC*(-(A1+AAA));
	S[2][5]=SC*(MU*(B1+BBB));
	S[2][6]=SC*(A1+AAA);
	S[2][7]=SC*(-MU*(B1+BBB));
	S[2][8]=SC*(A1-AAA);
	S[3][1]=SC*(-(1-MU)/2*(A1-AAA));
	S[3][2]=SC*(-(1-MU)/2*(B1-BBB));
	S[3][3]=SC*(-(1-MU)/2*(A1+AAA));
	S[3][4]=SC*((1-MU)/2*(B1-BBB));
	S[3][5]=SC*((1-MU)/2*(A1+AAA));
	S[3][6]=SC*((1-MU)/2*(B1+BBB));
	S[3][7]=SC*((1-MU)/2*(A1-AAA));
	S[3][8]=SC*((1-MU)/2*(B1+BBB));


	


}


	/*	for(i=1;i<=3;i++)
			for(j=1;j<=6;j++)
				B[i][j]=0.0;
			B[1][1]=(-BJ-BM)/(2.0*AE);
			B[1][3]=BJ/(2*AE);
			B[1][5]=BM/(2*AE);
			B[2][2]=(-CJ-CM)/(2*AE);
			B[2][4]=CJ/(2*AE);
			B[2][6]=CM/(2*AE);
			B[3][1]=B[2][2];
			B[3][2]=B[1][1];
			B[3][3]=B[2][4];
			B[3][4]=B[1][3];
			B[3][5]=B[2][6];
			B[3][6]=B[1][5];
			D[1][1]=E0/(1-MU*MU);
			D[1][2]=E0*MU/(1-MU*MU);
			D[1][3]=0;
			D[2][1]=D[1][2];
			D[2][2]=D[1][1];
			D[2][3]=0;
			D[3][1]=0;
			D[3][2]=0;
			D[3][3]=E0/(2*(1+MU)); */

			/*for(i=1;i<=3;i++)
			{
				for(j=1;j<=6;j++)
				{
					S[i][j]=0.0;
				for(k=1;k<=3;k++)
					S[i][j]=S[i][j]+D[i][k]*B[k][j];
				}

			}*/
//			if(ASK>2)
//			{
					/*	for(i=1;i<=8;i++)
						{
										for(j=1;j<=8;j++)
										{
											KE[i][j]=0.0;
											for(k=1;k<=3;k++)
												KE[i][j]=KE[i][j]+S[k][i]*B[k][j]*AE*TE;
										}
						}*/
				
//			}
	
	//输出第E个ke矩阵
/*	printf("\n第%d个KE矩阵\n",E);
	for(i=1;i<=6;i++)
	{
		for(j=1;j<=6;j++)
			printf("%f\t",KE[i][j]);
		printf("\n");
	}*/
}
//gdnl()函数：<功能：将非节点荷载下的杆端力计算出来存入f0[]>
//***********************************************************

void gdnl(int hz)
{

	int ind,e,JIEDIAN1,JIEDIAN2;
	double g,c;
	 g=pf[hz][1];
	 c=pf[hz][2];
	 e=(int)pf[hz][3];
	 ind=(int)pf[hz][6];
	 JIEDIAN1=(int)pf[hz][4];        //求节点码
	 JIEDIAN2=(int)pf[hz][5];
	 
//	 l0=gc[e];
//	 d=l0-c;

	 if(ind==1)
	 
	 {
	 
		 if((AJZ[JIEDIAN1][1]-AJZ[JIEDIAN2][1])!=0)     //如果不等于零  两点横向排列
		 {
			 f0[2]=fabs(AJZ[JIEDIAN1][1]-AJZ[JIEDIAN2][1])*g/2;         //横向的两点排列   f0值为Y方向力
			 f0[1]=0.0;


		 }
		 if((AJZ[JIEDIAN1][1]-AJZ[JIEDIAN2][1])==0)                 //如果不等于零  两点纵向排列
		 
		 {
			 f0[1]=fabs(AJZ[JIEDIAN2][2]-AJZ[JIEDIAN2][2])*g/2;      //这时为x方向节点荷载
			f0[2]=0.0;
		 }
	 }
																/*		 f0[1]=0.0;
																	 f0[2]=-(g*c*(2-2*c*c/(l0*l0)+(c*c*c)/(l0*l0*l0)))/2;
																	 f0[3]=-(g*c*c)*(6-8*c/l0+3*c*c/(l0*l0))/12;
																	 f0[4]=0.0;
																	 f0[5]=-g*c-f0[2];
																	 f0[6]=(g*c*c*c)*(4-3*c/l0)/(12*l0);


																 }
																 else 
																 {
																	 if(ind==2)
																	 {
																		 f0[1]=0.0;
																		 f0[2]=(-(g*d*d)*(l0+2*c))/(l0*l0*l0);
																		 f0[3]=-(g*c*d*d)/(l0*l0);
																		 f0[4]=0.0;
																		 f0[5]=(-(g*c*c)*(l0+2*d))/(l0*l0*l0);     //不知道书上是不是错了？
																		 f0[6]=(g*c*c*d)/(l0*l0);

																	 }
																	 else
																	 {
																		 f0[1]=-(g*d/l0);
																		 f0[2]=0.0;
																		 f0[3]=0.0;
																		 f0[4]=-g*c/l0;
																		 f0[5]=0.0;
																		 f0[6]=0.0;
																	 }
																  }        */
}



///自动划分网格函数
void mesh()
{   
	
	double ccl;  //板长方向单元长度
	double ccb;   //板高方向单元长度

	int i,j,k;
	 
	ccl=ll/nl;    //单元x方向长度
	ccb=bb/nb;     //单元y方向长度
	for(i=0;i<=NJ;i++)       //节点坐标数组清零
	{	
		for(j=0;j<=3;j++)
		{
			AJZ[i][j]=0.0;
		}
	}
	for(i=1;i<=NE;i++)          //弹性模量变阶数组清零
	{
		E0[i]=0.0;
	} 
	
	for(i=1;i<=(nl+1);i++)  
	{
		                                    //从左下角数  第i列 第k个单元  x坐标
		for(k=1;k<=nb+1;k++)
		{
			AJZ[(nb+1)*(i-1)+k][1]=(i-1)*ccl;   //求节点坐标x


		}
				
		
	}
	for(i=1;i<=nb+1;i++)
	{
		for(k=1;k<=nl+1;k++)
		{
			AJZ[(nb+1)*(k-1)+i][2]=(i-1)*ccb;            //求节点坐标y
		}
	}
	for(i=1;i<=nl;i++)                              //求单元节点码
	{
		for(j=1;j<=nb;j++)
		{
			JM[(i-1)*nb+j][1]=(i-1)*(nb+1)+j;
			JM[(i-1)*nb+j][2]=i*(nb+1)+j;
			JM[(i-1)*nb+j][3]=i*(nb+1)+j+1;
			JM[(i-1)*nb+j][4]=(i-1)*(nb+1)+j+1;

		}
		
	}
	for(i=1;i<=nl;i++)                              //求各单元变阶E0
	{
		for(j=1;j<=nb;j++)
		{
			E0[(i-1)*nb+j]=E00+(i*ccl-ccl/2)*Ex+(j*ccb-ccb/2)*Ey+(i*ccl-ccl/2)*(j*ccb-ccb/2)*Ex*Ey;   //   E0=E00+Ex*X+Ey*Y+Exy*x*y  在X，Y方向增大系数

		}
		
	}
	for(i=1;i<=NPF;i++)
	{
		for(j=1;j<=6;j++)
		{
			pf[i][j]=0.0;
		}
	}
	for(i=1;i<=NPF;i++)        //{0,荷载值，荷载位置、长度，单元，节点码1，节点码2，荷载类型}
	
	{
		pf[i][1]=junbu;  //均布荷载大小
		pf[i][2]=ccl;   //每个单元x方向长度
		pf[i][3]=i*nb;   //单元号
		pf[i][4]=(nb+1)*i; //节点码1
		pf[i][5]=(nb+1)*(i+1); //节点码2
		pf[i][6]=1.0;
	}
	for(i=1;i<=NZ;i++)
	{
		NZC[i]=i;
	}
//	NZC[i]=((nb+1)*(nl+1)-nb)*2;


}

