#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<conio.h>
#define  MAX  50

typedef   struct   wuli{
           float   d[MAX];
           char   name[50];
           int       LEN;
           float   ccha[MAX];       /*残差数组*/
           float   avg;                 /*data的平均值*/
           double   sx;                 /*标准偏差Sx*/

}wulidata;
  wulidata   *InputData();
  void   average(wulidata   *wl);
   void   YCZhi(wulidata   *wl);
  void   CanCha(wulidata   *wl);
  void   BZPianCha(wulidata   *wl);

  void   output(wulidata   *wl);
  void range(wulidata   *wl);


  /*----------------------------------------------------------*/
  void   line()
  {
                  int   i;
                  printf("\n");
                  for(i=0;i<74;i++)
					printf("=");
                  printf("\n");
  }
  /*-------------------------------------------------------*/
  wulidata   *InputData()
  {

                  int i=0,k;
                  float da;
                  char Z=0;
                  wulidata   *wl;
                  wl=(wulidata   *)malloc(sizeof(wulidata));
                  printf("请为你要处理的数据组命名：");
                  scanf("%s",wl->name);
                  printf("\n下面请你输入数据%s具体数值，数据不能超过50个\n",wl->name);
                  printf("当name='#'时输入结束\n");
                  do{
                        printf("%s%d=",wl->name,i+1);
                        scanf("%f",&da);
                        wl->d[i]=da;
						i++;
						if(getchar()=='#') break;
                  }while(wl->d[i-1]!=0.0&&i<MAX);
                  wl->LEN=i-1;
                  do{
                            printf("你输入的数据如下：\n");
                            for(i=0;i<wl->LEN;i++)
                                            printf("%s%d=%f\t",wl->name,i+1,wl->d[i]);
                                    printf("\n你是否要作出修改(Y/N)?");
                            while(   getchar()!='\n');
                                    Z=getchar();

                            if(   Z=='y'||Z=='Y'){
                                              printf("你须要修改哪一个元素，请输入其标号i=(1~%d)\n",wl->LEN);
                                                      while(   getchar()!='\n');
                              scanf("%d",&k);
                                              printf("\n%s%d=",wl->name,k);
                                              scanf("%f",&(wl->d[k-1]));
                      }
                            else   if(Z=='n'||Z=='N')
                                              printf("OK!下面开始计算。\n");

                  }while(Z!='N'&&Z!='n');
                  return(wl);
  }
  /*--------------------------------------------------------------------*/
  void   average(wulidata   *wl)
  {
                  float   ad,sum=0;
                  int   i;
                  for(i=0;i<wl->LEN;i++)
                  {
                                  sum=sum+(wl->d[i]);
                  }
                  ad=sum/(wl->LEN);
                  wl->avg=ad;
  }
  /*-------------------------------------------------------------------*/
  void   CanCha(wulidata   *wl)
  {
                  int   i;
                  for(i=0;i<wl->LEN;i++)
                    wl->ccha[i]=(wl->d[i])-(wl->avg);
  }
    /*---------------------------------------------------------------*/
  void   YCZhi(wulidata   *wl)/*检查并剔除异常值*/
  {
                  int   i,j;
				  float g,YCZhi;
				  double temp,CCha;

                  printf("下面开始检查并提出异常值！\n");
				  do{
						printf("当前共有%d个数，数据如下：\n",wl->LEN);
						for(i=0;i<wl->LEN;i++)
							printf("%s%d=%f\t",wl->name,i+1,wl->d[i]);
						j=-1;
						CCha=0.0;
						printf("\n请输入g的值\ng=");
						scanf("%f",&g);
						for(i=0;i<wl->LEN;i++)
						{
							temp=fabs((wl->d[i])-(wl->avg));
							if((temp>g*(wl->sx))&&(temp>CCha))
							{
								YCZhi=wl->d[i];
								CCha=temp;
								j=i;
							}
						}
						if(j>=0){
							printf("找到异常值为%s%d=%f，将它剔除。\n",wl->name,(j+1),wl->d[j]);
							for(i=j;i<wl->LEN-1;i++)
								wl->d[i]=wl->d[i+1];
							wl->LEN--;
						}
						else
							printf("本次未找到异常数据，数据中异常数据已剔除完毕！\n");
				  }while(j>=0);
				  printf("当前共有%d个数，数据如下：\n",wl->LEN);
				  for(i=0;i<wl->LEN;i++)
					printf("%s%d=%f\t",wl->name,i+1,wl->d[i]);
  }
  /*---------------------------------------------------------------*/
  void   BZPianCha(wulidata   *wl)/*标准偏差*/
  {
                  double   sum;
                  int   i;
                  sum=0.0;
                  for(i=0;i<wl->LEN;i++)
                      sum=sum+pow(wl->d[i],2);
				  sum=sum-wl->LEN*pow(wl->avg,2);
                  wl->sx=sqrt(sum/(wl->LEN-1));
  }
  /*--------------------------------------------------------------*/
  void leijinxwc(wulidata   *wl)/*判断累进性误差*/
	 {
	  double M,sum1,sum2,temp;
  int i;
  sum1=sum2=0.0;
  temp=wl->ccha[0];
 for (i=1;i<=wl->LEN;i++)
	 if (temp<wl->ccha[i])
		 temp=wl->ccha[i];
  if(wl->LEN%2==0)            /*数据为偶数个时*/
	 {for(i=0;i<(wl->LEN/2);i++)
	  {sum1=sum1+ wl->ccha[i];}
	  for (i=(wl->LEN/2);i<wl->LEN;i++)
	  {sum2=sum2+wl->ccha[i];}
  M=fabs(sum1-sum2);

     if(M>temp)
	  printf("存在累进性误差\n");
  else printf("不存在累进性误差\n");}
 else                      /*数据为奇数个时*/
 { for(i=0;i<(((wl->LEN)-1)/2);i++)
	  {sum1=sum1+ wl->ccha[i];}
	  for (i=((wl->LEN)+1)/2;i<wl->LEN;i++)
	  {sum2=sum2+wl->ccha[i];}
  M=fabs(sum1-sum2);
   if(M>temp)
	  printf("存在累进性误差\n");
  else printf("不存在累进性误差\n");}
  }
  /*---------------------------------------------------------------*/
  void zhouqixwc(wulidata   *wl)/*判断周期性误差*/
  {	double sum=0;
  int i;
  for(i=0;i<wl->LEN-1;i++)
	  sum=sum+(wl->ccha[i])*(wl->ccha[i+1]);
  if(fabs(sum)>(sqrt(wl->LEN-1))*pow(wl->sx,2))
	  printf("存在周期性误差\n");
  else printf("不存在周期性误差\n");
  }

  /*--------------------------------------------------------------------*/
   void range(wulidata   *wl)/*给出95%置信区间*/
   {
	   double u1,u2,x,ta;
	   printf("\n请输入ta的值\nta=");
		scanf("%f",&ta);
		ta=2.262;
		x=(ta)*(wl->sx/sqrt(wl->LEN));
		u1=wl->avg-x;
		u2=wl->avg+x;
		printf("平均值为%.5f\n平均值的标准偏差为%.5f\n",wl->avg,wl->sx/sqrt(wl->LEN));
		printf("计算得到所求数值的范围应取%.5f~%.5f。\n",u1,u2);
		 }


  void   output(wulidata   *wl)
  {
                  int   i;
                  printf("\n");
                  line();
                  printf("你输入的数据如下：\n");
                  for(i=0;i<wl->LEN;i++)
                  {printf("%s%d=%f\t",wl->name,i+1,wl->d[i]);}
          printf("\n");
                  printf("\n\t数据%s的平均值(A)%s=%f",wl->name,wl->name,wl->avg);
                  line();
                  printf("数据的残差如下:\n");
                  for(i=0;i<wl->LEN;i++)
                  {printf("△%s%d=%s%d-(A)%s=%f\t\t",wl->name,i+1,wl->name,i+1,wl->name,wl->ccha[i]);}
          line();
          printf("求得标准偏差Sx\n");
                  printf("Sx=%f",wl->sx);
                  printf("\n");


  }
  /*================================================================*/
 int main()
  {
            wulidata   *Hua=NULL;

            Hua=InputData();

            average(Hua);
            BZPianCha(Hua);/*标准偏差*/
			YCZhi(Hua);
              average(Hua);
			BZPianCha(Hua);
			CanCha(Hua);
            leijinxwc(Hua);
            zhouqixwc(Hua);
            output(Hua);
			 range(Hua);
            getch();
			return 0;
  }
