#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<conio.h>
#define  MAX  50

typedef   struct   wuli{
           float   d[MAX];
           char   name[50];
           int       LEN;
           float   ccha[MAX];       /*�в�����*/
           float   avg;                 /*data��ƽ��ֵ*/
           double   sx;                 /*��׼ƫ��Sx*/

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
                  printf("��Ϊ��Ҫ�����������������");
                  scanf("%s",wl->name);
                  printf("\n����������������%s������ֵ�����ݲ��ܳ���50��\n",wl->name);
                  printf("��name='#'ʱ�������\n");
                  do{
                        printf("%s%d=",wl->name,i+1);
                        scanf("%f",&da);
                        wl->d[i]=da;
						i++;
						if(getchar()=='#') break;
                  }while(wl->d[i-1]!=0.0&&i<MAX);
                  wl->LEN=i-1;
                  do{
                            printf("��������������£�\n");
                            for(i=0;i<wl->LEN;i++)
                                            printf("%s%d=%f\t",wl->name,i+1,wl->d[i]);
                                    printf("\n���Ƿ�Ҫ�����޸�(Y/N)?");
                            while(   getchar()!='\n');
                                    Z=getchar();

                            if(   Z=='y'||Z=='Y'){
                                              printf("����Ҫ�޸���һ��Ԫ�أ�����������i=(1~%d)\n",wl->LEN);
                                                      while(   getchar()!='\n');
                              scanf("%d",&k);
                                              printf("\n%s%d=",wl->name,k);
                                              scanf("%f",&(wl->d[k-1]));
                      }
                            else   if(Z=='n'||Z=='N')
                                              printf("OK!���濪ʼ���㡣\n");

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
  void   YCZhi(wulidata   *wl)/*��鲢�޳��쳣ֵ*/
  {
                  int   i,j;
				  float g,YCZhi;
				  double temp,CCha;

                  printf("���濪ʼ��鲢����쳣ֵ��\n");
				  do{
						printf("��ǰ����%d�������������£�\n",wl->LEN);
						for(i=0;i<wl->LEN;i++)
							printf("%s%d=%f\t",wl->name,i+1,wl->d[i]);
						j=-1;
						CCha=0.0;
						printf("\n������g��ֵ\ng=");
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
							printf("�ҵ��쳣ֵΪ%s%d=%f�������޳���\n",wl->name,(j+1),wl->d[j]);
							for(i=j;i<wl->LEN-1;i++)
								wl->d[i]=wl->d[i+1];
							wl->LEN--;
						}
						else
							printf("����δ�ҵ��쳣���ݣ��������쳣�������޳���ϣ�\n");
				  }while(j>=0);
				  printf("��ǰ����%d�������������£�\n",wl->LEN);
				  for(i=0;i<wl->LEN;i++)
					printf("%s%d=%f\t",wl->name,i+1,wl->d[i]);
  }
  /*---------------------------------------------------------------*/
  void   BZPianCha(wulidata   *wl)/*��׼ƫ��*/
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
  void leijinxwc(wulidata   *wl)/*�ж��۽������*/
	 {
	  double M,sum1,sum2,temp;
  int i;
  sum1=sum2=0.0;
  temp=wl->ccha[0];
 for (i=1;i<=wl->LEN;i++)
	 if (temp<wl->ccha[i])
		 temp=wl->ccha[i];
  if(wl->LEN%2==0)            /*����Ϊż����ʱ*/
	 {for(i=0;i<(wl->LEN/2);i++)
	  {sum1=sum1+ wl->ccha[i];}
	  for (i=(wl->LEN/2);i<wl->LEN;i++)
	  {sum2=sum2+wl->ccha[i];}
  M=fabs(sum1-sum2);

     if(M>temp)
	  printf("�����۽������\n");
  else printf("�������۽������\n");}
 else                      /*����Ϊ������ʱ*/
 { for(i=0;i<(((wl->LEN)-1)/2);i++)
	  {sum1=sum1+ wl->ccha[i];}
	  for (i=((wl->LEN)+1)/2;i<wl->LEN;i++)
	  {sum2=sum2+wl->ccha[i];}
  M=fabs(sum1-sum2);
   if(M>temp)
	  printf("�����۽������\n");
  else printf("�������۽������\n");}
  }
  /*---------------------------------------------------------------*/
  void zhouqixwc(wulidata   *wl)/*�ж����������*/
  {	double sum=0;
  int i;
  for(i=0;i<wl->LEN-1;i++)
	  sum=sum+(wl->ccha[i])*(wl->ccha[i+1]);
  if(fabs(sum)>(sqrt(wl->LEN-1))*pow(wl->sx,2))
	  printf("�������������\n");
  else printf("���������������\n");
  }

  /*--------------------------------------------------------------------*/
   void range(wulidata   *wl)/*����95%��������*/
   {
	   double u1,u2,x,ta;
	   printf("\n������ta��ֵ\nta=");
		scanf("%f",&ta);
		ta=2.262;
		x=(ta)*(wl->sx/sqrt(wl->LEN));
		u1=wl->avg-x;
		u2=wl->avg+x;
		printf("ƽ��ֵΪ%.5f\nƽ��ֵ�ı�׼ƫ��Ϊ%.5f\n",wl->avg,wl->sx/sqrt(wl->LEN));
		printf("����õ�������ֵ�ķ�ΧӦȡ%.5f~%.5f��\n",u1,u2);
		 }


  void   output(wulidata   *wl)
  {
                  int   i;
                  printf("\n");
                  line();
                  printf("��������������£�\n");
                  for(i=0;i<wl->LEN;i++)
                  {printf("%s%d=%f\t",wl->name,i+1,wl->d[i]);}
          printf("\n");
                  printf("\n\t����%s��ƽ��ֵ(A)%s=%f",wl->name,wl->name,wl->avg);
                  line();
                  printf("���ݵĲв�����:\n");
                  for(i=0;i<wl->LEN;i++)
                  {printf("��%s%d=%s%d-(A)%s=%f\t\t",wl->name,i+1,wl->name,i+1,wl->name,wl->ccha[i]);}
          line();
          printf("��ñ�׼ƫ��Sx\n");
                  printf("Sx=%f",wl->sx);
                  printf("\n");


  }
  /*================================================================*/
 int main()
  {
            wulidata   *Hua=NULL;

            Hua=InputData();

            average(Hua);
            BZPianCha(Hua);/*��׼ƫ��*/
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
