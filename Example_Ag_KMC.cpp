#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <fstream>
using namespace std;
const int xdim=33;                                                               //X dimension////////////////////////////////////////////////////////////Editable//////////////////////////////
const int ydim=33;                                                               //Y dimension////////////////////////////////////////////////////////////Editable//////////////////////////////
int coverspecies;                                                               //The number of total species(need input)
int freespecies; //The number of free component
int site[6*xdim*ydim+1];                                                        //The component of sites(1-totalspecies)
int sitetype[6*xdim*ydim+1];                                                    //The type of sites(1-A Top, 2-Fcc, 3-Hcp, 4-Bridge1, 5-Bridge2, 6-Bridge3)
int treeposition[6*xdim*ydim+1];                                                //The position of a site in the tree
double treerate[18*xdim*ydim+1];                                                //Total rate of one tree site
int effectreactnum[6*xdim*ydim+1]; //The number of effective reactions on one site
int effectreactname[6*xdim*ydim+1][71]; //The name of effective reactions on one site
double effectreactrate[6*xdim*ydim+1][71]; //The rate of effective reactions
int changednum; //The number of site need to be refreshed
int changedsite[10000];                                                         //The tree site number of changed sites
int changedchk[6*xdim*ydim+1]; //The list of sites that are alreaday checked
int position[xdim+1][ydim+1][7];                                                //The site number of one point
int coordinatesx[6*xdim*ydim+1]; //The X coordinate of one site
int coordinatesy[6*xdim*ydim+1]; //The Y coordinate of one site
int sn[6*xdim*ydim+1][151]; //The site number of one site's neighbor site
int refreshneighbornum[5]; //The number of one site's neighbor that need to be refreshed
int refreshneighbor[5][151]; //The sites that need to be refreshed near one site
double componentcoverage[51]; //The coverage of different speices
double freequantity[51]; //The quantity of free component
int reactnum[5][51]; //The number of reactions on sitetype X with component Y
int reactsum[5][51][101]; //The sum of reactions happened
double reactrate[5][51][101]; //The rate of reactions on sitetype X with component Y of certain reaction type
int reactreplacenum[5][51][101]; //The number of basic site that need to be replaced
int reactreplacepos[5][51][101][151]; //The sites that need to be replaced
int reactreplaceval[5][51][101][151]; //The value that need to be replaced
int reactfreenum[5][51][101]; //The number of free components that need to be changed
int reactinter[5][51][101];                                                     //The number of react interaction
int reactfreename[5][51][101][51]; //The free components that need to be changed
double reactfreeval[5][51][101][51]; //THe value change of free components
int reactlimitnum[5][51][101]; //The number of limits of one specific reaction
int reactlimitpos[5][51][101][151]; //The limit position of a specific raction
int reactlimitval[5][51][101][151]; //The value of these limit position
int totalsite,totaltreesite,treelevel,treetop; //Number of total sites and total tree sites
double coverage; //Initial coverage
int covercomponent,coversite; //Characters of initial coverage
int selectedsite,selectedtreesite; //Main Step
int step;
double eFactor,ro2ad,ro2dis,ro2des,ro2hyd,rohyd,rho2dis,rhodes,po,ea,temp,rohyd2,ke,ksf,ppo,rhoads,rho2hyd,rho2des,rh2o2dis,rh2o2des,kb,h,xOH,deltObinding,deltO2binding,deltOHbinding,deltOOHbinding,rho2hyd2,ro2hyd2,Eao2dis,Eao2hyd,Eao2hyd2,Eaohyd,Eaho2dis,Eaho2hor,Eaho2hyd2,rho2hor,rho2dis2,eIntercept; //parameters
double GAg6H2O,GH2O,GH2,GO2,GH2O2,GadO,GadO2,GadOH,GadOOH,Gslab,GH,GOH,GHO2,dGO2ads,dGO2dis,dGOpro,dGO2pro,dGOOHdes,dGOOHdis,dGOHdes,ActEO2dis,ActEO2pro,ActEOpro,ActEO2Hdis,ActEO2Hpro,ActHO2HOr;
double BEPO2dis,BEPO2Hdis,BEPO2Hpro,BEPO2pro,BEPOpro,BEPHO2HOr,sGadO,sGadO2,sGadOH,sGadOOH,BEPohohr,Actohohr,rohohr,Eaohohr,switch_Surf,switch_Elec,rHtrans,rOHoxidation;
int    switch_Ea;
////////////////////////////////////////////////////////////////////////////////
int pnt(int xx, int yy, int pp)
{
    xx=xx+xdim;
    yy=yy+ydim;
    xx=xx%xdim;
    yy=yy%ydim;
    if (xx==0) xx=xdim;
    if (yy==0) yy=ydim;
    return (position[xx][yy][pp]);
};
////////////////////////////////////////////////////////////////////////////////
void initialization()                                                           //Initialization begin
{
    int i,j,k,l,x,y,mmm;
    for (i=0;i<=4;i++)for(j=0;j<=50;j++) reactnum[i][j]=0;
    for (i=0;i<=4;i++)for(j=0;j<=50;j++)for(k=0;k<=100;k++) reactrate[i][j][k]=0;
    for (i=0;i<=4;i++)for(j=0;j<=50;j++)for(k=0;k<=100;k++) reactinter[i][j][k]=0;
    for (i=0;i<=4;i++)for(j=0;j<=50;j++)for(k=0;k<=100;k++) reactsum[i][j][k]=0;
    for (i=0;i<=4;i++)for(j=0;j<=50;j++)for(k=0;k<=100;k++) reactreplacenum[i][j][k]=0;
    for (i=0;i<=4;i++)for(j=0;j<=50;j++)for(k=0;k<=100;k++)for(l=0;l<=150;l++) reactreplacepos[i][j][k][l]=0;
    for (i=0;i<=4;i++)for(j=0;j<=50;j++)for(k=0;k<=100;k++)for(l=0;l<=150;l++) reactreplaceval[i][j][k][l]=0;
    for (i=0;i<=4;i++)for(j=0;j<=50;j++)for(k=0;k<=100;k++) reactfreenum[i][j][k]=0;
    for (i=0;i<=4;i++)for(j=0;j<=50;j++)for(k=0;k<=100;k++)for(l=0;l<=50;l++) reactfreename[i][j][k][l]=0;
    for (i=0;i<=4;i++)for(j=0;j<=50;j++)for(k=0;k<=100;k++)for(l=0;l<=50;l++) reactfreeval[i][j][k][l]=0;
    for (i=0;i<=4;i++)for(j=0;j<=50;j++)for(k=0;k<=100;k++) reactlimitnum[i][j][k]=0;
    for (i=0;i<=4;i++)for(j=0;j<=50;j++)for(k=0;k<=100;k++)for(l=0;l<=150;l++) reactlimitpos[i][j][k][l]=0;
    for (i=0;i<=4;i++)for(j=0;j<=50;j++)for(k=0;k<=100;k++)for(l=0;l<=150;l++) reactlimitval[i][j][k][l]=0;
    
    
    
    for (i=0;i<=50;i++) freequantity[i]=0;
    for (i=0;i<=4;i++) refreshneighbornum[i]=0;
    for (i=0;i<=4;i++)for(j=0;j<=150;j++) refreshneighbor[i][j]=0;
    
    
    totalsite=xdim*ydim*6;
    
    for (i=0;i<=totalsite;i++)
        changedchk[i]=0;
    
    k=1;                                                                       //Initialization of lattice
    for(i=1;i<=xdim;i++)
        for(j=1;j<=ydim;j++)
            for(l=1;l<=6;l++)
            {coordinatesx[k]=i;
                coordinatesy[k]=j;
                sitetype[k]=l;
                position[i][j][l]=k;
                k=k+1;
            };
    ////////////////////////////////////////////////////////////////////////////////
    //Initialization of neighbors
    for(mmm=1;mmm<=totalsite;mmm++)
        sn[mmm][0]=mmm;
    for(mmm=1;mmm<=totalsite;mmm++)
    { x=coordinatesx[mmm];
        y=coordinatesy[mmm];
        if (sitetype[mmm]==1){sn[mmm][1]=pnt(x,y+1,1); sn[mmm][2]=pnt(x+1,y,1);  sn[mmm][3]=pnt(x+1,y-1,1); sn[mmm][4]=pnt(x,y-1,1); sn[mmm][5]=pnt(x-1,y,1);  sn[mmm][6]=pnt(x-1,y+1,1);
            sn[mmm][7]=pnt(x,y,2) ; sn[mmm][8]=pnt(x,y-1,3);  sn[mmm][9]=pnt(x,y-1,2);   sn[mmm][10]=pnt(x-1,y-1,3); sn[mmm][11]=pnt(x-1,y,2);  sn[mmm][12]=pnt(x-1,y,3);
            sn[mmm][13]=pnt(x,y,3); sn[mmm][14]=pnt(x+1,y-1,2);  sn[mmm][15]=pnt(x,y-2,3); sn[mmm][16]=pnt(x-1,y-1,2); sn[mmm][17]=pnt(x-2,y,3);  sn[mmm][18]=pnt(x-1,y+1,2);
            sn[mmm][19]=pnt(x-1,y+1,3); sn[mmm][20]=pnt(x,y+1,2);  sn[mmm][21]=pnt(x+1,y,2); sn[mmm][22]=pnt(x+1,y-1,3); sn[mmm][23]=pnt(x+1,y-2,3);  sn[mmm][24]=pnt(x+1,y-2,2); sn[mmm][25]=pnt(x,y-2,2); sn[mmm][26]=pnt(x-1,y-2,3);  sn[mmm][27]=pnt(x-2,y-1,3); sn[mmm][28]=pnt(x-2,y,2); sn[mmm][29]=pnt(x-2,y+1,2);  sn[mmm][30]=pnt(x-2,y+1,3);
            sn[mmm][31]=pnt(x,y+2,1);   sn[mmm][32]=pnt(x+1,y+1,1);  sn[mmm][33]=pnt(x+2,y,1); sn[mmm][34]=pnt(x+2,y-1,1); sn[mmm][35]=pnt(x+2,y-2,1);  sn[mmm][36]=pnt(x+1,y-2,1); sn[mmm][37]=pnt(x,y-2,1); sn[mmm][38]=pnt(x-1,y-1,1);  sn[mmm][39]=pnt(x-2,y,1);  sn[mmm][40]=pnt(x-2,y+1,1);sn[mmm][41]=pnt(x-2,y+2,1);  sn[mmm][42]=pnt(x-1,y+2,1);
            sn[mmm][43]=pnt(x,y,5); sn[mmm][44]=pnt(x,y,4);  sn[mmm][45]=pnt(x+1,y-1,6);sn[mmm][46]=pnt(x,y-1,5); sn[mmm][47]=pnt(x-1,y,4);  sn[mmm][48]=pnt(x,y,6);
            sn[mmm][49]=pnt(x+1,y,6);   sn[mmm][50]=pnt(x+1,y-1,5);  sn[mmm][51]=pnt(x,y-1,4); sn[mmm][52]=pnt(x,y-1,6); sn[mmm][53]=pnt(x-1,y,5);    sn[mmm][54]=pnt(x-1,y+1,4);
            sn[mmm][55]=pnt(x,y+1,6);   sn[mmm][56]=pnt(x,y+1,4);    sn[mmm][57]=pnt(x+1,y,5); sn[mmm][58]=pnt(x+2,y-1,6); sn[mmm][59]=pnt(x+1,y-1,4);  sn[mmm][60]=pnt(x+1,y-2,5); sn[mmm][61]=pnt(x+1,y-2,6);sn[mmm][62]=pnt(x-1,y-1,4);  sn[mmm][63]=pnt(x-1,y-1,5); sn[mmm][64]=pnt(x-1,y,6);  sn[mmm][65]=pnt(x-2,y+1,4);  sn[mmm][66]=pnt(x-1,y+1,5);
            sn[mmm][67]=pnt(x,y+1,5);   sn[mmm][68]=pnt(x+1,y,4);    sn[mmm][69]=pnt(x+2,y-2,6);sn[mmm][70]=pnt(x,y-2,5); sn[mmm][71]=pnt(x-2,y,4);    sn[mmm][72]=pnt(x-1,y+1,6);
            sn[mmm][73]=pnt(x-1,y+2,4); sn[mmm][74]=pnt(x+1,y+1,6);  sn[mmm][75]=pnt(x+2,y,6); sn[mmm][76]=pnt(x+2,y-1,5); sn[mmm][77]=pnt(x+2,y-2,5);  sn[mmm][78]=pnt(x+1,y-2,4); sn[mmm][79]=pnt(x,y-2,4);  sn[mmm][80]=pnt(x,y-2,6);    sn[mmm][81]=pnt(x-1,y-1,6); sn[mmm][82]=pnt(x-2,y,5);  sn[mmm][83]=pnt(x-2,y+1,5);  sn[mmm][84]=pnt(x-2,y+2,4);
        };
        
        if (sitetype[mmm]==2){sn[mmm][1]=pnt(x+1,y,1);  sn[mmm][2]=pnt(x,y,1); sn[mmm][3]=pnt(x,y+1,1);
            sn[mmm][4]=pnt(x,y,3); sn[mmm][5]=pnt(x,y-1,3); sn[mmm][6]=pnt(x-1,y,3);
            sn[mmm][7]=pnt(x,y+1,2); sn[mmm][8]=pnt(x+1,y,2); sn[mmm][9]=pnt(x+1,y-1,2);  sn[mmm][10]=pnt(x,y-1,2); sn[mmm][11]=pnt(x-1,y,2); sn[mmm][12]=pnt(x-1,y+1,2);
            sn[mmm][13]=pnt(x+1,y+1,1); sn[mmm][14]=pnt(x+1,y-1,1); sn[mmm][15]=pnt(x-1,y+1,1);
            sn[mmm][16]=pnt(x+1,y-1,3); sn[mmm][17]=pnt(x-1,y-1,3); sn[mmm][18]=pnt(x-1,y+1,3);
            sn[mmm][19]=pnt(x,y+2,1); sn[mmm][20]=pnt(x+2,y,1); sn[mmm][21]=pnt(x+2,y-1,1); sn[mmm][22]=pnt(x,y-1,1); sn[mmm][23]=pnt(x-1,y,1); sn[mmm][24]=pnt(x-1,y+2,1);
            sn[mmm][25]=pnt(x,y+1,3); sn[mmm][26]=pnt(x+1,y,3); sn[mmm][27]=pnt(x+1,y-2,3); sn[mmm][28]=pnt(x,y-2,3); sn[mmm][29]=pnt(x-2,y,3); sn[mmm][30]=pnt(x-2,y+1,3);
            sn[mmm][31]=pnt(x+1,y+1,2); sn[mmm][32]=pnt(x+2,y-1,2); sn[mmm][33]=pnt(x+1,y-2,2);
            sn[mmm][34]=pnt(x-1,y-1,2); sn[mmm][35]=pnt(x-2,y+1,2); sn[mmm][36]=pnt(x-1,y+2,2);
            sn[mmm][37]=pnt(x,y+2,2); sn[mmm][38]=pnt(x+2,y,2); sn[mmm][39]=pnt(x+2,y-2,2); sn[mmm][40]=pnt(x,y-2,2); sn[mmm][41]=pnt(x-2,y,2); sn[mmm][42]=pnt(x-2,y+2,2);
            sn[mmm][43]=pnt(x-1,y+2,3); sn[mmm][44]=pnt(x+2,y-1,3); sn[mmm][45]=pnt(x+2,y-2,3); sn[mmm][46]=pnt(x-1,y-2,3); sn[mmm][47]=pnt(x-2,y-1,3); sn[mmm][48]=pnt(x-2,y+2,3);
            sn[mmm][49]=pnt(x,y+3,1); sn[mmm][50]=pnt(x+1,y+2,1); sn[mmm][51]=pnt(x+2,y+1,1); sn[mmm][52]=pnt(x+3,y,1); sn[mmm][53]=pnt(x+3,y-1,1); sn[mmm][54]=pnt(x+3,y-2,1); sn[mmm][55]=pnt(x+2,y-2,1); sn[mmm][56]=pnt(x+1,y-2,1); sn[mmm][57]=pnt(x,y-2,1); sn[mmm][58]=pnt(x-1,y-1,1); sn[mmm][59]=pnt(x-2,y,1); sn[mmm][60]=pnt(x-2,y+1,1); sn[mmm][61]=pnt(x-2,y+2,1); sn[mmm][62]=pnt(x-2,y+3,1); sn[mmm][63]=pnt(x-1,y+3,1);
            sn[mmm][64]=pnt(x+1,y,6); sn[mmm][65]=pnt(x,y,4);     sn[mmm][66]=pnt(x,y,5);
            sn[mmm][67]=pnt(x,y+1,4); sn[mmm][68]=pnt(x+1,y,5); sn[mmm][69]=pnt(x+1,y-1,5); sn[mmm][70]=pnt(x+1,y-1,6); sn[mmm][71]=pnt(x,y,6); sn[mmm][72]=pnt(x-1,y+1,4);
            sn[mmm][73]=pnt(x,y+1,5); sn[mmm][74]=pnt(x+1,y,4); sn[mmm][75]=pnt(x+2,y-1,6); sn[mmm][76]=pnt(x,y-1,5); sn[mmm][77]=pnt(x-1,y,4); sn[mmm][78]=pnt(x,y+1,6);
            sn[mmm][79]=pnt(x+1,y+1,6); sn[mmm][80]=pnt(x+2,y,6); sn[mmm][81]=pnt(x+1,y-1,4); sn[mmm][82]=pnt(x,y-1,4); sn[mmm][83]=pnt(x-1,y,5); sn[mmm][84]=pnt(x-1,y+1,5);
            sn[mmm][85]=pnt(x+2,y-1,5); sn[mmm][86]=pnt(x,y-1,6); sn[mmm][87]=pnt(x-1,y+2,4);
            sn[mmm][88]=pnt(x,y+2,6); sn[mmm][89]=pnt(x,y+2,5); sn[mmm][90]=pnt(x,y+2,4);   sn[mmm][91]=pnt(x+1,y+1,5); sn[mmm][92]=pnt(x+1,y+1,4); sn[mmm][93]=pnt(x+2,y,5);  sn[mmm][94]=pnt(x+2,y,4); sn[mmm][95]=pnt(x+3,y-1,6); sn[mmm][96]=pnt(x+2,y-1,4); sn[mmm][97]=pnt(x+3,y-2,6); sn[mmm][98]=pnt(x+2,y-2,5); sn[mmm][99]=pnt(x+2,y-2,6); sn[mmm][100]=pnt(x+1,y-2,5); sn[mmm][101]=pnt(x+1,y-2,6);sn[mmm][102]=pnt(x,y-2,5); sn[mmm][103]=pnt(x-1,y-1,4); sn[mmm][104]=pnt(x-1,y-1,5); sn[mmm][105]=pnt(x-2,y,4); sn[mmm][106]=pnt(x-1,y,6); sn[mmm][107]=pnt(x-2,y+1,4); sn[mmm][108]=pnt(x-1,y+1,6); sn[mmm][109]=pnt(x-2,y+2,4); sn[mmm][110]=pnt(x-1,y+2,6); sn[mmm][111]=pnt(x-1,y+2,5);
            sn[mmm][112]=pnt(x-1,y+3,4);sn[mmm][113]=pnt(x+1,y+2,6);sn[mmm][114]=pnt(x+2,y+1,6);
            sn[mmm][115]=pnt(x+3,y,6); sn[mmm][116]=pnt(x+3,y-1,5);sn[mmm][117]=pnt(x+3,y-2,5);sn[mmm][118]=pnt(x+2,y-2,4);sn[mmm][119]=pnt(x+1,y-2,4);sn[mmm][120]=pnt(x,y-2,4); sn[mmm][121]=pnt(x,y-2,6); sn[mmm][122]=pnt(x-1,y-1,6);sn[mmm][123]=pnt(x-2,y,5); sn[mmm][124]=pnt(x-2,y+1,5); sn[mmm][125]=pnt(x-2,y+2,5);sn[mmm][126]=pnt(x-2,y+3,4);
        };
        if (sitetype[mmm]==3){sn[mmm][1]=pnt(x+1,y+1,1);  sn[mmm][2]=pnt(x+1,y,1); sn[mmm][3]=pnt(x,y+1,1);
            sn[mmm][4]=pnt(x+1,y,2); sn[mmm][5]=pnt(x,y,2); sn[mmm][6]=pnt(x,y+1,2);
            sn[mmm][7]=pnt(x,y+1,3); sn[mmm][8]=pnt(x+1,y,3); sn[mmm][9]=pnt(x+1,y-1,3); sn[mmm][10]=pnt(x,y-1,3); sn[mmm][11]=pnt(x-1,y,3); sn[mmm][12]=pnt(x-1,y+1,3);
            sn[mmm][13]=pnt(x+2,y,1); sn[mmm][14]=pnt(x,y,1); sn[mmm][15]=pnt(x,y+2,1);
            sn[mmm][16]=pnt(x+1,y+1,2); sn[mmm][17]=pnt(x+1,y-1,2); sn[mmm][18]=pnt(x-1,y+1,2);
            sn[mmm][19]=pnt(x+1,y+2,1); sn[mmm][20]=pnt(x+2,y+1,1); sn[mmm][21]=pnt(x+2,y-1,1); sn[mmm][22]=pnt(x+1,y-1,1); sn[mmm][23]=pnt(x-1,y+1,1); sn[mmm][24]=pnt(x-1,y+2,1);
            sn[mmm][25]=pnt(x,y+2,2); sn[mmm][26]=pnt(x+2,y,2); sn[mmm][27]=pnt(x+2,y-1,2); sn[mmm][28]=pnt(x,y-1,2); sn[mmm][29]=pnt(x-1,y,2); sn[mmm][30]=pnt(x-1,y+2,2);
            sn[mmm][31]=pnt(x+1,y+1,3); sn[mmm][32]=pnt(x+2,y-1,3); sn[mmm][33]=pnt(x+1,y-2,3); sn[mmm][34]=pnt(x-1,y-1,3); sn[mmm][35]=pnt(x-2,y+1,3); sn[mmm][36]=pnt(x-1,y+2,3);
            sn[mmm][37]=pnt(x,y+2,3); sn[mmm][38]=pnt(x+2,y,3); sn[mmm][39]=pnt(x+2,y-2,3); sn[mmm][40]=pnt(x,y-2,3); sn[mmm][41]=pnt(x-2,y,3); sn[mmm][42]=pnt(x-2,y+2,3);
            sn[mmm][43]=pnt(x+1,y+2,2); sn[mmm][44]=pnt(x+2,y+1,2); sn[mmm][45]=pnt(x+2,y-2,2); sn[mmm][46]=pnt(x+1,y-2,2); sn[mmm][47]=pnt(x-2,y+1,2); sn[mmm][48]=pnt(x-2,y+2,2);
            sn[mmm][49]=pnt(x,y+3,1); sn[mmm][50]=pnt(x+1,y+3,1); sn[mmm][51]=pnt(x+2,y+2,1); sn[mmm][52]=pnt(x+3,y+1,1); sn[mmm][53]=pnt(x+3,y,1); sn[mmm][54]=pnt(x+3,y-1,1); sn[mmm][55]=pnt(x+3,y-2,1); sn[mmm][56]=pnt(x+2,y-2,1); sn[mmm][57]=pnt(x+1,y-2,1); sn[mmm][58]=pnt(x,y-1,1); sn[mmm][59]=pnt(x-1,y,1); sn[mmm][60]=pnt(x-2,y+1,1); sn[mmm][61]=pnt(x-2,y+2,1); sn[mmm][62]=pnt(x-2,y+3,1); sn[mmm][63]=pnt(x-1,y+3,1);
            sn[mmm][64]=pnt(x+1,y,5); sn[mmm][65]=pnt(x+1,y,6); sn[mmm][66]=pnt(x,y+1,4);
            sn[mmm][67]=pnt(x+1,y+1,6); sn[mmm][68]=pnt(x+2,y,6); sn[mmm][69]=pnt(x+1,y,4); sn[mmm][70]=pnt(x,y,4); sn[mmm][71]=pnt(x,y,5); sn[mmm][72]=pnt(x,y+1,5);
            sn[mmm][73]=pnt(x+1,y+1,5); sn[mmm][74]=pnt(x+1,y+1,4); sn[mmm][75]=pnt(x+2,y-1,6); sn[mmm][76]=pnt(x+1,y-1,5); sn[mmm][77]=pnt(x-1,y+1,4); sn[mmm][78]=pnt(x,y+1,6);
            sn[mmm][79]=pnt(x,y+2,4); sn[mmm][80]=pnt(x+2,y,5); sn[mmm][81]=pnt(x+2,y-1,5); sn[mmm][82]=pnt(x+1,y-1,6); sn[mmm][83]=pnt(x,y,6); sn[mmm][84]=pnt(x-1,y+2,4);
            sn[mmm][85]=pnt(x+2,y+1,6); sn[mmm][86]=pnt(x+1,y-1,4); sn[mmm][87]=pnt(x-1,y+1,5);
            sn[mmm][88]=pnt(x,y+2,5); sn[mmm][89]=pnt(x+1,y+2,6); sn[mmm][90]=pnt(x+1,y+2,5); sn[mmm][91]=pnt(x+1,y+2,4); sn[mmm][92]=pnt(x+2,y+1,5); sn[mmm][93]=pnt(x+2,y+1,4); sn[mmm][94]=pnt(x+3,y,6); sn[mmm][95]=pnt(x+2,y,4); sn[mmm][96]=pnt(x+3,y-1,6); sn[mmm][97]=pnt(x+2,y-1,4); sn[mmm][98]=pnt(x+3,y-2,6); sn[mmm][99]=pnt(x+2,y-2,5); sn[mmm][100]=pnt(x+2,y-2,6); sn[mmm][101]=pnt(x+1,y-2,5); sn[mmm][102]=pnt(x,y-1,4); sn[mmm][103]=pnt(x,y-1,5); sn[mmm][104]=pnt(x-1,y,4); sn[mmm][105]=pnt(x-1,y,5); sn[mmm][106]=pnt(x-2,y+1,4); sn[mmm][107]=pnt(x-1,y+1,6); sn[mmm][108]=pnt(x-2,y+2,4); sn[mmm][109]=pnt(x-1,y+2,6); sn[mmm][110]=pnt(x-1,y+2,5); sn[mmm][111]=pnt(x,y+2,6);
            sn[mmm][112]=pnt(x,y+3,4); sn[mmm][113]=pnt(x+2,y+2,6);sn[mmm][114]=pnt(x+3,y+1,6);sn[mmm][115]=pnt(x+3,y,5); sn[mmm][116]=pnt(x+3,y-1,5);sn[mmm][117]=pnt(x+3,y-2,5);sn[mmm][118]=pnt(x+2,y-2,4);sn[mmm][119]=pnt(x+1,y-2,4); sn[mmm][120]=pnt(x+1,y-2,6); sn[mmm][121]=pnt(x,y-1,6); sn[mmm][122]=pnt(x-1,y,6); sn[mmm][123]=pnt(x-2,y+1,5); sn[mmm][124]=pnt(x-2,y+2,5); sn[mmm][125]=pnt(x-2,y+3,4); sn[mmm][126]=pnt(x-1,y+3,4);
        };
        if (sitetype[mmm]==4){sn[mmm][1]=pnt(x+1,y,1); sn[mmm][2]=pnt(x,y,1);
            sn[mmm][3]=pnt(x,y-1,3); sn[mmm][4]=pnt(x,y,2);
            sn[mmm][5]=pnt(x+1,y-1,1); sn[mmm][6]=pnt(x,y+1,1);
            sn[mmm][7]=pnt(x,y,3); sn[mmm][8]=pnt(x+1,y-1,2); sn[mmm][9]=pnt(x,y-1,2); sn[mmm][10]=pnt(x-1,y,3);
            sn[mmm][11]=pnt(x+1,y,2); sn[mmm][12]=pnt(x+1,y-1,3); sn[mmm][13]=pnt(x-1,y-1,3); sn[mmm][14]=pnt(x-1,y,2);
            sn[mmm][15]=pnt(x+1,y+1,1); sn[mmm][16]=pnt(x+2,y-1,1); sn[mmm][17]=pnt(x,y-1,1); sn[mmm][18]=pnt(x-1,y+1,1);
            sn[mmm][19]=pnt(x,y+1,2); sn[mmm][20]=pnt(x+1,y-2,3); sn[mmm][21]=pnt(x,y-2,3); sn[mmm][22]=pnt(x-1,y+1,2);
            sn[mmm][23]=pnt(x+2,y,1); sn[mmm][24]=pnt(x-1,y,1);
            sn[mmm][25]=pnt(x+1,y,3); sn[mmm][26]=pnt(x+2,y-1,2); sn[mmm][27]=pnt(x-1,y-1,2); sn[mmm][28]=pnt(x-2,y,3);
            sn[mmm][29]=pnt(x+1,y-2,2); sn[mmm][30]=pnt(x-1,y+1,3);
            sn[mmm][31]=pnt(x,y+1,3); sn[mmm][32]=pnt(x+2,y-2,2); sn[mmm][33]=pnt(x,y-2,2); sn[mmm][34]=pnt(x-2,y+1,3);
            sn[mmm][35]=pnt(x+1,y+1,2); sn[mmm][36]=pnt(x+2,y-2,3); sn[mmm][37]=pnt(x-1,y-2,3); sn[mmm][38]=pnt(x-2,y+1,2);
            sn[mmm][39]=pnt(x+2,y,2); sn[mmm][40]=pnt(x+2,y-1,3); sn[mmm][41]=pnt(x-2,y-1,3); sn[mmm][42]=pnt(x-2,y,2);
            sn[mmm][43]=pnt(x,y+2,1); sn[mmm][44]=pnt(x+1,y+2,1); sn[mmm][45]=pnt(x+2,y+1,1); sn[mmm][46]=pnt(x+3,y,1); sn[mmm][47]=pnt(x+3,y-1,1); sn[mmm][48]=pnt(x+3,y-2,1); sn[mmm][49]=pnt(x+2,y-2,1); sn[mmm][50]=pnt(x+1,y-2,1); sn[mmm][51]=pnt(x,y-2,1); sn[mmm][52]=pnt(x-1,y-1,1); sn[mmm][53]=pnt(x-2,y,1); sn[mmm][54]=pnt(x-2,y+1,1); sn[mmm][55]=pnt(x-2,y+2,1); sn[mmm][56]=pnt(x-1,y+2,1);
            sn[mmm][57]=pnt(x+1,y,6); sn[mmm][58]=pnt(x+1,y-1,5); sn[mmm][59]=pnt(x+1,y-1,6); sn[mmm][60]=pnt(x,y,5);
            sn[mmm][61]=pnt(x+1,y,5); sn[mmm][62]=pnt(x+2,y-1,6); sn[mmm][63]=pnt(x,y-1,5); sn[mmm][64]=pnt(x,y,6);
            sn[mmm][65]=pnt(x+1,y,4); sn[mmm][66]=pnt(x-1,y,4);
            sn[mmm][67]=pnt(x,y+1,4); sn[mmm][68]=pnt(x+1,y-1,4); sn[mmm][69]=pnt(x,y-1,4); sn[mmm][70]=pnt(x-1,y+1,4);
            sn[mmm][71]=pnt(x,y+1,5); sn[mmm][72]=pnt(x+2,y-2,6); sn[mmm][73]=pnt(x+1,y-2,5); sn[mmm][74]=pnt(x,y+1,6);
            sn[mmm][75]=pnt(x+2,y,6); sn[mmm][76]=pnt(x+2,y-1,5); sn[mmm][77]=pnt(x,y-1,6); sn[mmm][78]=pnt(x-1,y,5);
            sn[mmm][79]=pnt(x+1,y+1,6); sn[mmm][80]=pnt(x+2,y-2,5); sn[mmm][81]=pnt(x+1,y-2,6); sn[mmm][82]=pnt(x-1,y+1,5);
            sn[mmm][83]=pnt(x+1,y+1,4); sn[mmm][84]=pnt(x+2,y-1,4); sn[mmm][85]=pnt(x-1,y-1,4); sn[mmm][86]=pnt(x-2,y+1,4);
            sn[mmm][87]=pnt(x+1,y+1,5); sn[mmm][88]=pnt(x+3,y-2,6); sn[mmm][89]=pnt(x,y-2,5); sn[mmm][90]=pnt(x-1,y+1,6);
            sn[mmm][91]=pnt(x+2,y,5); sn[mmm][92]=pnt(x+3,y-1,6); sn[mmm][93]=pnt(x-1,y-1,5); sn[mmm][94]=pnt(x-1,y,6);
            sn[mmm][95]=pnt(x+2,y,4); sn[mmm][96]=pnt(x-2,y,4);
            sn[mmm][97]=pnt(x,y+2,4); sn[mmm][98]=pnt(x+2,y+1,6); sn[mmm][99]=pnt(x+3,y,6); sn[mmm][100]=pnt(x+3,y-1,5); sn[mmm][101]=pnt(x+3,y-2,5); sn[mmm][102]=pnt(x+2,y-2,4); sn[mmm][103]=pnt(x+1,y-2,4); sn[mmm][104]=pnt(x,y-2,4); sn[mmm][105]=pnt(x,y-2,6); sn[mmm][106]=pnt(x-1,y-1,6); sn[mmm][107]=pnt(x-2,y,5); sn[mmm][108]=pnt(x-2,y+1,5); sn[mmm][109]=pnt(x-2,y+2,4); sn[mmm][110]=pnt(x-1,y+2,4);
        };
        if (sitetype[mmm]==5){sn[mmm][1]=pnt(x,y,1); sn[mmm][2]=pnt(x,y+1,1);
            sn[mmm][3]=pnt(x-1,y,3); sn[mmm][4]=pnt(x,y,2);
            sn[mmm][5]=pnt(x-1,y+1,1); sn[mmm][6]=pnt(x+1,y,1);
            sn[mmm][7]=pnt(x,y-1,3); sn[mmm][8]=pnt(x-1,y,2); sn[mmm][9]=pnt(x-1,y+1,2); sn[mmm][10]=pnt(x,y,3);
            sn[mmm][11]=pnt(x,y-1,2); sn[mmm][12]=pnt(x-1,y-1,3); sn[mmm][13]=pnt(x-1,y+1,3); sn[mmm][14]=pnt(x,y+1,2);
            sn[mmm][15]=pnt(x+1,y-1,1); sn[mmm][16]=pnt(x-1,y,1); sn[mmm][17]=pnt(x-1,y+2,1); sn[mmm][18]=pnt(x+1,y+1,1);
            sn[mmm][19]=pnt(x+1,y-1,2); sn[mmm][20]=pnt(x-2,y,3); sn[mmm][21]=pnt(x-2,y+1,3); sn[mmm][22]=pnt(x+1,y,2);
            sn[mmm][23]=pnt(x,y-1,1); sn[mmm][24]=pnt(x,y+2,1);
            sn[mmm][25]=pnt(x,y-2,3); sn[mmm][26]=pnt(x-1,y-1,2); sn[mmm][27]=pnt(x-1,y+2,2); sn[mmm][28]=pnt(x,y+1,3);
            sn[mmm][29]=pnt(x-2,y+1,2); sn[mmm][30]=pnt(x+1,y-1,3);
            sn[mmm][31]=pnt(x+1,y-2,3); sn[mmm][32]=pnt(x-2,y,2); sn[mmm][33]=pnt(x-2,y+2,2); sn[mmm][34]=pnt(x+1,y,3);
            sn[mmm][35]=pnt(x+1,y-2,2); sn[mmm][36]=pnt(x-2,y-1,3); sn[mmm][37]=pnt(x-2,y+2,3); sn[mmm][38]=pnt(x+1,y+1,2);
            sn[mmm][39]=pnt(x,y-2,2); sn[mmm][40]=pnt(x-1,y-2,3); sn[mmm][41]=pnt(x-1,y+2,3); sn[mmm][42]=pnt(x,y+2,2);
            sn[mmm][43]=pnt(x+2,y-1,1); sn[mmm][44]=pnt(x+2,y-2,1); sn[mmm][45]=pnt(x+1,y-2,1); sn[mmm][46]=pnt(x,y-2,1); sn[mmm][47]=pnt(x-1,y-1,1); sn[mmm][48]=pnt(x-2,y,1); sn[mmm][49]=pnt(x-2,y+1,1); sn[mmm][50]=pnt(x-2,y+2,1); sn[mmm][51]=pnt(x-2,y+3,1); sn[mmm][52]=pnt(x-1,y+3,1); sn[mmm][53]=pnt(x,y+3,1); sn[mmm][54]=pnt(x+1,y+2,1); sn[mmm][55]=pnt(x+2,y+1,1); sn[mmm][56]=pnt(x+2,y,1);
            sn[mmm][57]=pnt(x,y,4); sn[mmm][58]=pnt(x,y,6); sn[mmm][59]=pnt(x-1,y+1,4); sn[mmm][60]=pnt(x+1,y,6);
            sn[mmm][61]=pnt(x+1,y-1,6); sn[mmm][62]=pnt(x-1,y,4); sn[mmm][63]=pnt(x,y+1,6); sn[mmm][64]=pnt(x,y+1,4);
            sn[mmm][65]=pnt(x,y-1,5); sn[mmm][66]=pnt(x,y+1,5);
            sn[mmm][67]=pnt(x+1,y-1,5); sn[mmm][68]=pnt(x-1,y,5); sn[mmm][69]=pnt(x-1,y+1,5); sn[mmm][70]=pnt(x+1,y,5);
            sn[mmm][71]=pnt(x+2,y-1,6); sn[mmm][72]=pnt(x-2,y+1,4); sn[mmm][73]=pnt(x-1,y+1,6); sn[mmm][74]=pnt(x+1,y,4);
            sn[mmm][75]=pnt(x,y-1,4); sn[mmm][76]=pnt(x,y-1,6); sn[mmm][77]=pnt(x-1,y+2,4); sn[mmm][78]=pnt(x+1,y+1,6);
            sn[mmm][79]=pnt(x+1,y-1,4); sn[mmm][80]=pnt(x-1,y,6); sn[mmm][81]=pnt(x-2,y+2,4); sn[mmm][82]=pnt(x+2,y,6);
            sn[mmm][83]=pnt(x+1,y-2,5); sn[mmm][84]=pnt(x-1,y-1,5); sn[mmm][85]=pnt(x-1,y+2,5); sn[mmm][86]=pnt(x+1,y+1,5);
            sn[mmm][87]=pnt(x+2,y-2,6); sn[mmm][88]=pnt(x-2,y,4); sn[mmm][89]=pnt(x-1,y+2,6); sn[mmm][90]=pnt(x+1,y+1,4);
            sn[mmm][91]=pnt(x+1,y-2,6); sn[mmm][92]=pnt(x-1,y-1,4); sn[mmm][93]=pnt(x,y+2,6); sn[mmm][94]=pnt(x,y+2,4);
            sn[mmm][95]=pnt(x,y-2,5); sn[mmm][96]=pnt(x,y+2,5);
            sn[mmm][97]=pnt(x+2,y-2,5); sn[mmm][98]=pnt(x+1,y-2,4); sn[mmm][99]=pnt(x,y-2,4); sn[mmm][100]=pnt(x,y-2,6); sn[mmm][101]=pnt(x-1,y-1,6); sn[mmm][102]=pnt(x-2,y,5); sn[mmm][103]=pnt(x-2,y+1,5); sn[mmm][104]=pnt(x-2,y+2,5); sn[mmm][105]=pnt(x-2,y+3,4); sn[mmm][106]=pnt(x-1,y+3,4); sn[mmm][107]=pnt(x+1,y+2,6); sn[mmm][108]=pnt(x+2,y+1,6); sn[mmm][109]=pnt(x+2,y,5); sn[mmm][110]=pnt(x+2,y-1,5);
        };
        if (sitetype[mmm]==6){sn[mmm][1]=pnt(x-1,y+1,1); sn[mmm][2]=pnt(x,y,1);
            sn[mmm][3]=pnt(x-1,y,3); sn[mmm][4]=pnt(x-1,y,2);
            sn[mmm][5]=pnt(x,y+1,1); sn[mmm][6]=pnt(x-1,y,1);
            sn[mmm][7]=pnt(x-2,y,3); sn[mmm][8]=pnt(x-1,y+1,2); sn[mmm][9]=pnt(x,y,2); sn[mmm][10]=pnt(x-1,y-1,3);
            sn[mmm][11]=pnt(x-2,y+1,2); sn[mmm][12]=pnt(x-2,y+1,3); sn[mmm][13]=pnt(x,y-1,3); sn[mmm][14]=pnt(x,y-1,2);
            sn[mmm][15]=pnt(x-2,y+1,1); sn[mmm][16]=pnt(x-1,y+2,1); sn[mmm][17]=pnt(x+1,y,1); sn[mmm][18]=pnt(x,y-1,1);
            sn[mmm][19]=pnt(x-2,y,2); sn[mmm][20]=pnt(x-1,y+1,3); sn[mmm][21]=pnt(x,y,3); sn[mmm][22]=pnt(x-1,y-1,2);
            sn[mmm][23]=pnt(x-2,y+2,1); sn[mmm][24]=pnt(x+1,y-1,1);
            sn[mmm][25]=pnt(x-3,y+1,3); sn[mmm][26]=pnt(x-2,y+2,2); sn[mmm][27]=pnt(x+1,y-1,2); sn[mmm][28]=pnt(x,y-2,3);
            sn[mmm][29]=pnt(x,y+1,2); sn[mmm][30]=pnt(x-2,y-1,3);
            sn[mmm][31]=pnt(x-3,y,3); sn[mmm][32]=pnt(x-1,y+2,2); sn[mmm][33]=pnt(x+1,y,2); sn[mmm][34]=pnt(x-1,y-2,3);
            sn[mmm][35]=pnt(x-3,y+1,2); sn[mmm][36]=pnt(x-2,y+2,3); sn[mmm][37]=pnt(x+1,y-1,3); sn[mmm][38]=pnt(x,y-2,2);
            sn[mmm][39]=pnt(x-3,y+2,2); sn[mmm][40]=pnt(x-3,y+2,3); sn[mmm][41]=pnt(x+1,y-2,3); sn[mmm][42]=pnt(x+1,y-2,2);
            sn[mmm][43]=pnt(x-2,y,1); sn[mmm][44]=pnt(x-3,y+1,1); sn[mmm][45]=pnt(x-3,y+2,1); sn[mmm][46]=pnt(x-3,y+3,1); sn[mmm][47]=pnt(x-2,y+3,1); sn[mmm][48]=pnt(x-1,y+3,1); sn[mmm][49]=pnt(x,y+2,1); sn[mmm][50]=pnt(x+1,y+1,1); sn[mmm][51]=pnt(x+2,y,1); sn[mmm][52]=pnt(x+2,y-1,1); sn[mmm][53]=pnt(x+2,y-2,1); sn[mmm][54]=pnt(x+1,y-2,1); sn[mmm][55]=pnt(x,y-2,1); sn[mmm][56]=pnt(x-1,y-1,1);
            sn[mmm][57]=pnt(x-1,y,5); sn[mmm][58]=pnt(x-1,y+1,4); sn[mmm][59]=pnt(x,y,5); sn[mmm][60]=pnt(x-1,y,4);
            sn[mmm][61]=pnt(x-2,y+1,4); sn[mmm][62]=pnt(x-1,y+1,5); sn[mmm][63]=pnt(x,y,4); sn[mmm][64]=pnt(x,y-1,5);
            sn[mmm][65]=pnt(x-1,y+1,6); sn[mmm][66]=pnt(x+1,y-1,6);
            sn[mmm][67]=pnt(x-1,y,6); sn[mmm][68]=pnt(x,y+1,6); sn[mmm][69]=pnt(x+1,y,6); sn[mmm][70]=pnt(x,y-1,6);
            sn[mmm][71]=pnt(x-2,y,4); sn[mmm][72]=pnt(x,y+1,5); sn[mmm][73]=pnt(x,y+1,4); sn[mmm][74]=pnt(x-1,y-1,5);
            sn[mmm][75]=pnt(x-2,y+1,5); sn[mmm][76]=pnt(x-2,y+2,4); sn[mmm][77]=pnt(x+1,y-1,5); sn[mmm][78]=pnt(x,y-1,4);
            sn[mmm][79]=pnt(x-2,y,5); sn[mmm][80]=pnt(x-1,y+2,4); sn[mmm][81]=pnt(x+1,y,5); sn[mmm][82]=pnt(x-1,y-1,4);
            sn[mmm][83]=pnt(x-2,y+1,6); sn[mmm][84]=pnt(x-1,y+2,6); sn[mmm][85]=pnt(x+2,y-1,6); sn[mmm][86]=pnt(x+1,y-2,6);
            sn[mmm][87]=pnt(x-3,y+1,4); sn[mmm][88]=pnt(x-1,y+2,5); sn[mmm][89]=pnt(x+1,y,4); sn[mmm][90]=pnt(x,y-2,5);
            sn[mmm][91]=pnt(x-3,y+2,4); sn[mmm][92]=pnt(x-2,y+2,5); sn[mmm][93]=pnt(x+1,y-1,4); sn[mmm][94]=pnt(x+1,y-2,5);
            sn[mmm][95]=pnt(x-2,y+2,6); sn[mmm][96]=pnt(x+2,y-2,6);
            sn[mmm][97]=pnt(x-2,y,6); sn[mmm][98]=pnt(x-3,y+1,5); sn[mmm][99]=pnt(x-3,y+2,5); sn[mmm][100]=pnt(x-3,y+3,4); sn[mmm][101]=pnt(x-2,y+3,4); sn[mmm][102]=pnt(x,y+2,6); sn[mmm][103]=pnt(x+1,y+1,6); sn[mmm][104]=pnt(x+2,y,6); sn[mmm][105]=pnt(x+2,y-1,5); sn[mmm][106]=pnt(x+2,y-2,5); sn[mmm][107]=pnt(x+1,y-2,4); sn[mmm][108]=pnt(x,y-2,4); sn[mmm][109]=pnt(x,y-2,6); sn[mmm][110]=pnt(x-1,y-1,6);
        };
        
        
        
        
    };                                                                       //End of neighbors
    for(mmm=1;mmm<=totalsite;mmm++)
        if (sitetype[mmm]>4) sitetype[mmm]=4;
};                                                                             // Initialization end
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Editable//////////////////////////
void originaldata()                                                             //Input original data
{
    coverspecies=11;
    freespecies=4;
    freequantity[1]=100000000000;    //O2
    freequantity[2]=0;                  //OH
    freequantity[3]=0;                  //OOH
    freequantity[4]=0;                  //H2O2
    int i,j,k,h,jl,jm,jo,jp;
    kb=1.38e-23;
    h=6.626e-34;
    temp=298.15;
    po=ppo;
    ksf=10000000000000;//kb*temp/h;
    ke=1000000000;
    eFactor=0.45;
    eIntercept=0.235;
    
    double deltX;
    deltX=0;
    sGadO=2.046*deltX;sGadO2=1.5807*deltX;sGadOH=deltX; sGadOOH=1.0874*deltX;
    GH2O=-14.34;GH2=-6.88;GO2=-10.06;GH2O2=-18.28;
    
    GadO=-171.79+sGadO;GadO2=-176.77+sGadO2;GadOH=-176.56+sGadOH;GadOOH=-180.58+sGadOOH;Gslab=-166.14;
    GH=-4.21; GOH=-10.13;GHO2=-14.07;GAg6H2O=-180.34;
    
    dGO2ads=GadO2-GO2-Gslab;dGO2dis=2*GadO-Gslab-GadO2;dGO2pro=GadOOH-GadO2+GOH-GH2O;dGOpro=GadOH-GadO+GOH-GH2O;
    dGOOHdes=Gslab+GHO2-GadOOH; dGOOHdis=GadOH+GadO-Gslab-GadOOH;dGOHdes=GOH+GAg6H2O-GadOH-GH2O;
    
    switch_Ea=0;//0-->Act    1-->BEP
    switch_Surf=0*9999;
    switch_Elec=0*9999;
    
    
    BEPO2dis=0.2163*(2*sGadO-sGadO2);BEPO2Hdis=0.5954*(sGadOH+sGadO-sGadOOH);BEPO2Hpro=-0.03*(3*sGadOH-sGadOOH);BEPO2pro=0.1971*(sGadOOH+sGadOH-sGadO2);BEPOpro=0.35*(2*sGadOH-sGadO);BEPHO2HOr=0.7667*(sGadOOH+sGadOH-sGadO2);BEPohohr=0.5791*(sGadO-2*sGadOH);
    BEPO2dis=BEPO2dis*switch_Ea;BEPO2Hdis=BEPO2Hdis*switch_Ea;BEPO2Hpro=BEPO2Hpro*switch_Ea;BEPO2pro=BEPO2pro*switch_Ea;BEPOpro=BEPOpro*switch_Ea;BEPHO2HOr=BEPHO2HOr*switch_Ea;BEPohohr=BEPohohr*switch_Ea;
    ActEO2dis=0.75+BEPO2dis; ActEO2Hdis=0.38+BEPO2Hdis; ActEO2Hpro=0.29+BEPO2Hpro; ActEO2pro=0.25+BEPO2pro;ActEOpro=0.05+BEPOpro;ActHO2HOr=0.28+BEPHO2HOr;Actohohr=0.9+BEPohohr;
    
    
    ro2ad=420;
    ea=(-1*dGO2ads)*96486.9;                                            ro2des=ksf*exp(((-1)*ea)/(8.314*temp));     //ro2des=0;
    ea=(0)*96486.9;                                                     ro2dis=ksf*exp(((-1)*ea)/(8.314*temp));     //ro2dis=0;
    ea=((dGO2pro+1*po)*eFactor+eIntercept)*96486.9;                     ro2hyd=ke*exp(((-1)*ea)/(8.314*temp));      //ro2hyd=0;
    ea=(0)*96486.9;                                                     ro2hyd2=ksf*exp(((-1)*ea)/(8.314*temp));    //ro2hyd2=0;
    ea=((dGOpro+1*po)*eFactor+eIntercept)*96486.9;                      rohyd=ke*exp(((-1)*ea)/(8.314*temp));       //rohyd=0;
    ea=(0)*96486.9;                                                     rohyd2=ksf*exp(((-1)*ea)/(8.314*temp));     //rohyd2=0;
    ea=(0)*96486.9;                                                     rho2dis=ksf*exp(((-1)*ea)/(8.314*temp));    //rho2dis=0;
    ea=(0)*96486.9;                                                     rho2hyd2=ksf*exp(((-1)*ea)/(8.314*temp));   //rho2hyd2=0;
    ea=((dGOOHdes+1*po)*eFactor+eIntercept)*96486.9;                    rho2des=ke*exp(((-1)*ea)/(8.314*temp));     //rho2des=0;
    ea=((dGOHdes+1*po)*eFactor+eIntercept)*96486.9;                     rhodes=ke*exp(((-1)*ea)/(8.314*temp));
    ea=((-1*dGOHdes-1*po)*eFactor+eIntercept)*96486.9;                  rhoads=ke*exp(((-1)*ea)/(8.314*temp));  //rhoads=0;
    ea=0;                                                               rho2hor=ksf*exp(((-1)*ea)/(8.314*temp));    //rho2hor=0;
    ea=((dGOOHdis+dGOHdes+1*po)*eFactor+eIntercept)*96486.9;            rho2dis2=ke*exp(((-1)*ea)/(8.314*temp));
    ea=0;                                                               rohohr=ksf*exp(((-1)*ea)/(8.314*temp));
    ea=0.09*96486.9;                                                    rHtrans=ksf*exp(((-1)*ea)/(8.314*temp));
    ea=((-dGOpro-1*po)*eFactor+eIntercept)*96486.9;                     rOHoxidation=ke*exp(((-1)*ea)/(8.314*temp));
    
    
    
    refreshneighbornum[4]=111;
    for (i=1;i<=111;i++) refreshneighbor[4][i]=i-1;
    /*refreshneighbor[4][1]=1;refreshneighbor[4][2]=2;refreshneighbor[4][3]=5;refreshneighbor[4][4]=6;refreshneighbor[4][5]=15;
     refreshneighbor[4][6]=16;refreshneighbor[4][7]=17;refreshneighbor[4][8]=18;refreshneighbor[4][9]=23;refreshneighbor[4][10]=24;
     refreshneighbor[4][11]=57;refreshneighbor[4][12]=58;refreshneighbor[4][13]=59;refreshneighbor[4][14]=60;refreshneighbor[4][15]=61;
     refreshneighbor[4][16]=62;refreshneighbor[4][17]=63;refreshneighbor[4][18]=64;refreshneighbor[4][19]=65;refreshneighbor[4][20]=66;
     refreshneighbor[4][21]=67;refreshneighbor[4][22]=68;refreshneighbor[4][23]=69;refreshneighbor[4][24]=70;refreshneighbor[4][25]=75;
     refreshneighbor[4][26]=76;refreshneighbor[4][27]=77;refreshneighbor[4][28]=78;refreshneighbor[4][29]=0;*/
    
    refreshneighbornum[1]=85;
    for (i=1;i<=85;i++) refreshneighbor[1][i]=i-1;
    /*refreshneighbor[1][1]=1;refreshneighbor[1][2]=2;refreshneighbor[1][3]=3;refreshneighbor[1][4]=4;refreshneighbor[1][5]=5;refreshneighbor[1][6]=6;
     refreshneighbor[1][7]=43;refreshneighbor[1][8]=44;refreshneighbor[1][9]=45;refreshneighbor[1][10]=46;refreshneighbor[1][11]=47;refreshneighbor[1][12]=48;
     refreshneighbor[1][13]=49;refreshneighbor[1][14]=50;refreshneighbor[1][15]=51;refreshneighbor[1][16]=52;refreshneighbor[1][17]=53;refreshneighbor[1][18]=54;
     refreshneighbor[1][19]=55;refreshneighbor[1][20]=56;refreshneighbor[1][21]=57;refreshneighbor[1][22]=58;refreshneighbor[1][23]=59;refreshneighbor[1][24]=60;
     refreshneighbor[1][25]=61;refreshneighbor[1][26]=62;refreshneighbor[1][27]=63;refreshneighbor[1][28]=64;refreshneighbor[1][29]=65;refreshneighbor[1][30]=66;
     refreshneighbor[1][31]=67;refreshneighbor[1][32]=68;refreshneighbor[1][33]=69;refreshneighbor[1][34]=70;refreshneighbor[1][35]=71;refreshneighbor[1][36]=72;
     refreshneighbor[1][37]=0;*/
    
    
    reactnum[4][0]=1;
    
    reactrate[4][0][1]=ro2ad;
    reactinter[4][0][1]=7;
    reactlimitnum[4][0][1]=6;
    reactlimitpos[4][0][1][1]=1;reactlimitpos[4][0][1][2]=2;reactlimitpos[4][0][1][3]=15;reactlimitpos[4][0][1][4]=16;
    reactlimitpos[4][0][1][5]=17;reactlimitpos[4][0][1][6]=18;reactlimitpos[4][0][1][7]=5;reactlimitpos[4][0][1][8]=6;reactlimitpos[4][0][1][9]=23;reactlimitpos[4][0][1][10]=24;
    reactlimitval[4][0][1][1]=0;reactlimitval[4][0][1][2]=0;reactlimitval[4][0][1][3]=0;reactlimitval[4][0][1][4]=0;
    reactlimitval[4][0][1][5]=0;reactlimitval[4][0][1][6]=0;reactlimitval[4][0][1][7]=0;reactlimitval[4][0][1][8]=0;reactlimitval[4][0][1][9]=0;reactlimitval[4][0][1][10]=0;
    reactreplacenum[4][0][1]=3;
    reactreplacepos[4][0][1][1]=0;reactreplacepos[4][0][1][2]=1;reactreplacepos[4][0][1][3]=2;
    reactreplaceval[4][0][1][1]=1;reactreplaceval[4][0][1][2]=5;reactreplaceval[4][0][1][3]=5;
    reactfreenum[4][0][1]=1;
    reactfreename[4][0][1][1]=1;
    reactfreeval[4][0][1][1]=-1;
    
    reactnum[4][1]=8;
    
    reactrate[4][1][1]=ro2dis;
    reactinter[4][1][1]=238;
    reactreplacenum[4][1][1]=3;
    reactreplacepos[4][1][1][1]=0;reactreplacepos[4][1][1][2]=1;reactreplacepos[4][1][1][3]=2;
    reactreplaceval[4][1][1][1]=0;reactreplaceval[4][1][1][2]=2;reactreplaceval[4][1][1][3]=2;
    
    reactrate[4][1][2]=ro2hyd;
    reactinter[4][1][2]=58;
    reactreplacenum[4][1][2]=3;
    reactreplacepos[4][1][2][1]=0;reactreplacepos[4][1][2][2]=1;reactreplacepos[4][1][2][3]=2;
    reactreplaceval[4][1][2][1]=4;reactreplaceval[4][1][2][2]=6;reactreplaceval[4][1][2][3]=7;
    
    reactrate[4][1][8]=ro2hyd;
    reactinter[4][1][8]=58;
    reactreplacenum[4][1][8]=3;
    reactreplacepos[4][1][8][1]=0;reactreplacepos[4][1][8][2]=1;reactreplacepos[4][1][8][3]=2;
    reactreplaceval[4][1][8][1]=4;reactreplaceval[4][1][8][2]=7;reactreplaceval[4][1][8][3]=6;
    
    
    
    reactrate[4][1][3]=ro2des;
    reactinter[4][1][3]=8;
    reactreplacenum[4][1][3]=3;
    reactreplacepos[4][1][3][1]=0;reactreplacepos[4][1][3][2]=1;reactreplacepos[4][1][3][3]=2;
    reactreplaceval[4][1][3][1]=0;reactreplaceval[4][1][3][2]=0;reactreplaceval[4][1][3][3]=0;
    reactfreenum[4][1][3]=1;
    reactfreename[4][1][3][1]=1;
    reactfreeval[4][1][3][1]=1.0;
    
    
    reactrate[4][1][4]=ro2hyd2;
    reactinter[4][1][4]=580;
    reactlimitnum[4][1][4]=1;
    reactlimitpos[4][1][4][1]=15;
    reactlimitval[4][1][4][1]=0;
    reactreplacenum[4][1][4]=4;
    reactreplacepos[4][1][4][1]=15;reactreplacepos[4][1][4][2]=0;reactreplacepos[4][1][4][3]=1;reactreplacepos[4][1][4][4]=2;
    reactreplaceval[4][1][4][1]=3; reactreplaceval[4][1][4][2]=4;reactreplaceval[4][1][4][3]=6;reactreplaceval[4][1][4][4]=7;
    
    reactrate[4][1][5]=ro2hyd2;
    reactinter[4][1][5]=580;
    reactlimitnum[4][1][5]=1;
    reactlimitpos[4][1][5][1]=16;
    reactlimitval[4][1][5][1]=0;
    reactreplacenum[4][1][5]=4;
    reactreplacepos[4][1][5][1]=16;reactreplacepos[4][1][5][2]=0;reactreplacepos[4][1][5][3]=1;reactreplacepos[4][1][5][4]=2;
    reactreplaceval[4][1][5][1]=3; reactreplaceval[4][1][5][2]=4;reactreplaceval[4][1][5][3]=6;reactreplaceval[4][1][5][4]=7;
    
    reactrate[4][1][6]=ro2hyd2;
    reactinter[4][1][6]=580;
    reactlimitnum[4][1][6]=1;
    reactlimitpos[4][1][6][1]=17;
    reactlimitval[4][1][6][1]=0;
    reactreplacenum[4][1][6]=4;
    reactreplacepos[4][1][6][1]=17;reactreplacepos[4][1][6][2]=0;reactreplacepos[4][1][6][3]=1;reactreplacepos[4][1][6][4]=2;
    reactreplaceval[4][1][6][1]=3; reactreplaceval[4][1][6][2]=4;reactreplaceval[4][1][6][3]=7;reactreplaceval[4][1][6][4]=6;
    
    reactrate[4][1][7]=ro2hyd2;
    reactinter[4][1][7]=580;
    reactlimitnum[4][1][7]=1;
    reactlimitpos[4][1][7][1]=18;
    reactlimitval[4][1][7][1]=0;
    reactreplacenum[4][1][7]=4;
    reactreplacepos[4][1][7][1]=18;reactreplacepos[4][1][7][2]=0;reactreplacepos[4][1][7][3]=1;reactreplacepos[4][1][7][4]=2;
    reactreplaceval[4][1][7][1]=3; reactreplaceval[4][1][7][2]=4;reactreplaceval[4][1][7][3]=7;reactreplaceval[4][1][7][4]=6;
    
    
    reactnum[4][4]=14;
    
    reactrate[4][4][1]=rho2dis;
    reactinter[4][4][1]=136;
    reactlimitnum[4][4][1]=1;
    reactlimitpos[4][4][1][1]=1;
    reactlimitval[4][4][1][1]=7;
    reactreplacenum[4][4][1]=3;
    reactreplacepos[4][4][1][1]=0;reactreplacepos[4][4][1][2]=1;reactreplacepos[4][4][1][3]=2;
    reactreplaceval[4][4][1][1]=0;reactreplaceval[4][4][1][2]=2;reactreplaceval[4][4][1][3]=3;
    
    reactrate[4][4][2]=rho2dis;
    reactinter[4][4][2]=136;
    reactlimitnum[4][4][2]=1;
    reactlimitpos[4][4][2][1]=1;
    reactlimitval[4][4][2][1]=6;
    reactreplacenum[4][4][2]=3;
    reactreplacepos[4][4][2][1]=0;reactreplacepos[4][4][2][2]=1;reactreplacepos[4][4][2][3]=2;
    reactreplaceval[4][4][2][1]=0;reactreplaceval[4][4][2][2]=3;reactreplaceval[4][4][2][3]=2;
    
    reactrate[4][4][3]=0;
    reactreplacenum[4][4][3]=1;
    reactreplacepos[4][4][3][1]=0;
    reactreplaceval[4][4][3][1]=7;
    
    reactrate[4][4][4]=rho2des;
    reactinter[4][4][4]=6;
    reactreplacenum[4][4][4]=3;
    reactreplacepos[4][4][4][1]=0;reactreplacepos[4][4][4][2]=1;reactreplacepos[4][4][4][3]=2;
    reactreplaceval[4][4][4][1]=0;reactreplaceval[4][4][4][2]=0;reactreplaceval[4][4][4][3]=0;
    reactfreenum[4][4][4]=1;
    reactfreename[4][4][4][1]=3;
    reactfreeval[4][4][4][1]=1.0;
    
    reactrate[4][4][5]=rho2hyd2;
    reactinter[4][4][5]=116;
    reactlimitnum[4][4][5]=2;
    reactlimitpos[4][4][5][1]=15;reactlimitpos[4][4][5][2]=1;
    reactlimitval[4][4][5][1]=0;reactlimitval[4][4][5][2]=7;
    reactreplacenum[4][4][5]=4;
    reactreplacepos[4][4][5][1]=15;reactreplacepos[4][4][5][2]=0;reactreplacepos[4][4][5][3]=1;reactreplacepos[4][4][5][4]=2;
    reactreplaceval[4][4][5][1]=3; reactreplaceval[4][4][5][2]=0;reactreplaceval[4][4][5][3]=3;reactreplaceval[4][4][5][4]=3;
    
    reactrate[4][4][6]=rho2hyd2;
    reactinter[4][4][6]=116;
    reactlimitnum[4][4][6]=2;
    reactlimitpos[4][4][6][1]=16;reactlimitpos[4][4][6][2]=1;
    reactlimitval[4][4][6][1]=0;reactlimitval[4][4][6][2]=7;
    reactreplacenum[4][4][6]=4;
    reactreplacepos[4][4][6][1]=16;reactreplacepos[4][4][6][2]=0;reactreplacepos[4][4][6][3]=1;reactreplacepos[4][4][6][4]=2;
    reactreplaceval[4][4][6][1]=3; reactreplaceval[4][4][6][2]=0;reactreplaceval[4][4][6][3]=3;reactreplaceval[4][4][6][4]=3;
    
    reactrate[4][4][7]=rho2hyd2;
    reactinter[4][4][7]=116;
    reactlimitnum[4][4][7]=2;
    reactlimitpos[4][4][7][1]=17;reactlimitpos[4][4][7][2]=2;
    reactlimitval[4][4][7][1]=0;reactlimitval[4][4][7][2]=7;
    reactreplacenum[4][4][7]=4;
    reactreplacepos[4][4][7][1]=17;reactreplacepos[4][4][7][2]=0;reactreplacepos[4][4][7][3]=1;reactreplacepos[4][4][7][4]=2;
    reactreplaceval[4][4][7][1]=3; reactreplaceval[4][4][7][2]=0;reactreplaceval[4][4][7][3]=3;reactreplaceval[4][4][7][4]=3;
    
    reactrate[4][4][8]=rho2hyd2;
    reactinter[4][4][8]=116;
    reactlimitnum[4][4][8]=2;
    reactlimitpos[4][4][8][1]=18;reactlimitpos[4][4][8][2]=2;
    reactlimitval[4][4][8][1]=0;reactlimitval[4][4][8][2]=7;
    reactreplacenum[4][4][8]=4;
    reactreplacepos[4][4][8][1]=18;reactreplacepos[4][4][8][2]=0;reactreplacepos[4][4][8][3]=1;reactreplacepos[4][4][8][4]=2;
    reactreplaceval[4][4][8][1]=3; reactreplaceval[4][4][8][2]=0;reactreplaceval[4][4][8][3]=3;reactreplaceval[4][4][8][4]=3;
    
    reactrate[4][4][9]=rho2hor;
    reactinter[4][4][9]=76;
    reactlimitnum[4][4][9]=2;
    reactlimitpos[4][4][9][1]=15;reactlimitpos[4][4][9][2]=1;
    reactlimitval[4][4][9][1]=3;reactlimitval[4][4][9][2]=6;
    reactreplacenum[4][4][9]=4;
    reactreplacepos[4][4][9][1]=15;reactreplacepos[4][4][9][2]=0;reactreplacepos[4][4][9][3]=1;reactreplacepos[4][4][9][4]=2;
    reactreplaceval[4][4][9][1]=0;reactreplaceval[4][4][9][2]=1;reactreplaceval[4][4][9][3]=5;reactreplaceval[4][4][9][4]=5;
    
    
    reactrate[4][4][10]=rho2hor;
    reactinter[4][4][10]=76;
    reactlimitnum[4][4][10]=2;
    reactlimitpos[4][4][10][1]=16;reactlimitpos[4][4][10][2]=1;
    reactlimitval[4][4][10][1]=3;reactlimitval[4][4][10][2]=6;
    reactreplacenum[4][4][10]=4;
    reactreplacepos[4][4][10][1]=16;reactreplacepos[4][4][10][2]=0;reactreplacepos[4][4][10][3]=1;reactreplacepos[4][4][10][4]=2;
    reactreplaceval[4][4][10][1]=0;reactreplaceval[4][4][10][2]=1;reactreplaceval[4][4][10][3]=5;reactreplaceval[4][4][10][4]=5;
    
    reactrate[4][4][11]=rho2hor;
    reactinter[4][4][11]=76;
    reactlimitnum[4][4][11]=2;
    reactlimitpos[4][4][11][1]=17;reactlimitpos[4][4][11][2]=2;
    reactlimitval[4][4][11][1]=3;reactlimitval[4][4][11][2]=6;
    reactreplacenum[4][4][11]=4;
    reactreplacepos[4][4][11][1]=17;reactreplacepos[4][4][11][2]=0;reactreplacepos[4][4][11][3]=1;reactreplacepos[4][4][11][4]=2;
    reactreplaceval[4][4][11][1]=0;reactreplaceval[4][4][11][2]=1;reactreplaceval[4][4][11][3]=5;reactreplaceval[4][4][11][4]=5;
    
    reactrate[4][4][12]=rho2hor;
    reactinter[4][4][12]=76;
    reactlimitnum[4][4][12]=2;
    reactlimitpos[4][4][12][1]=18;reactlimitpos[4][4][12][2]=2;
    reactlimitval[4][4][12][1]=3;reactlimitval[4][4][12][2]=6;
    reactreplacenum[4][4][12]=4;
    reactreplacepos[4][4][12][1]=18;reactreplacepos[4][4][12][2]=0;reactreplacepos[4][4][12][3]=1;reactreplacepos[4][4][12][4]=2;
    reactreplaceval[4][4][12][1]=0;reactreplaceval[4][4][12][2]=1;reactreplaceval[4][4][12][3]=5;reactreplaceval[4][4][12][4]=5;
    
    reactrate[4][4][13]=rho2dis2;
    reactinter[4][4][13]=36;
    reactlimitnum[4][4][13]=1;
    reactlimitpos[4][4][13][1]=1;
    reactlimitval[4][4][13][1]=6;
    reactreplacenum[4][4][13]=3;
    reactreplacepos[4][4][13][1]=0;reactreplacepos[4][4][13][2]=1;reactreplacepos[4][4][13][3]=2;
    reactreplaceval[4][4][13][1]=0;reactreplaceval[4][4][13][2]=0;reactreplaceval[4][4][13][3]=2;
    
    reactrate[4][4][14]=rho2dis2;
    reactinter[4][4][14]=36;
    reactlimitnum[4][4][14]=1;
    reactlimitpos[4][4][14][1]=2;
    reactlimitval[4][4][14][1]=6;
    reactreplacenum[4][4][14]=3;
    reactreplacepos[4][4][14][1]=0;reactreplacepos[4][4][14][2]=2;reactreplacepos[4][4][14][3]=1;
    reactreplaceval[4][4][14][1]=0;reactreplaceval[4][4][14][2]=0;reactreplaceval[4][4][14][3]=2;
    
    
    reactnum[1][2]=7;
    
    reactrate[1][2][1]=rohyd;
    reactinter[1][2][1]=140;
    reactlimitnum[1][2][1]=6;
    reactlimitpos[1][2][1][1]=1;reactlimitpos[1][2][1][2]=2;reactlimitpos[1][2][1][3]=3;reactlimitpos[1][2][1][4]=4;reactlimitpos[1][2][1][5]=5;reactlimitpos[1][2][1][6]=6;
    reactlimitval[1][2][1][1]=-3;reactlimitval[1][2][1][2]=-3;reactlimitval[1][2][1][3]=-3;reactlimitval[1][2][1][4]=-3;reactlimitval[1][2][1][5]=-3;reactlimitval[1][2][1][6]=-3;
    reactreplacenum[1][2][1]=1;
    reactreplacepos[1][2][1][1]=0;
    reactreplaceval[1][2][1][1]=3;
    
    reactrate[1][2][2]=rohyd2;
    reactinter[1][2][2]=14;
    reactlimitnum[1][2][2]=5;
    reactlimitpos[1][2][2][1]=1;reactlimitpos[1][2][2][2]=3;reactlimitpos[1][2][2][3]=5;reactlimitpos[1][2][2][4]=32;reactlimitpos[1][2][2][5]=42;
    reactlimitval[1][2][2][1]=0;reactlimitval[1][2][2][2]=0;reactlimitval[1][2][2][3]=0;reactlimitval[1][2][2][4]=0;reactlimitval[1][2][2][5]=0;
    reactreplacenum[1][2][2]=2;
    reactreplacepos[1][2][2][1]=0;reactreplacepos[1][2][2][2]=1;
    reactreplaceval[1][2][2][1]=3;reactreplaceval[1][2][2][2]=3;
    
    reactrate[1][2][3]=rohyd2;
    reactinter[1][2][3]=14;
    reactlimitnum[1][2][3]=5;
    reactlimitpos[1][2][3][1]=2;reactlimitpos[1][2][3][2]=4;reactlimitpos[1][2][3][3]=6;reactlimitpos[1][2][3][4]=32;reactlimitpos[1][2][3][5]=34;
    reactlimitval[1][2][3][1]=0;reactlimitval[1][2][3][2]=0;reactlimitval[1][2][3][3]=0;reactlimitval[1][2][3][4]=0;reactlimitval[1][2][3][5]=0;
    reactreplacenum[1][2][3]=2;
    reactreplacepos[1][2][3][1]=0;reactreplacepos[1][2][3][2]=2;
    reactreplaceval[1][2][3][1]=3;reactreplaceval[1][2][3][2]=3;
    
    reactrate[1][2][4]=rohyd2;
    reactinter[1][2][4]=14;
    reactlimitnum[1][2][4]=5;
    reactlimitpos[1][2][4][1]=3;reactlimitpos[1][2][4][2]=5;reactlimitpos[1][2][4][3]=1;reactlimitpos[1][2][4][4]=34;reactlimitpos[1][2][4][5]=36;
    reactlimitval[1][2][4][1]=0;reactlimitval[1][2][4][2]=0;reactlimitval[1][2][4][3]=0;reactlimitval[1][2][4][4]=0;reactlimitval[1][2][4][5]=0;
    reactreplacenum[1][2][4]=2;
    reactreplacepos[1][2][4][1]=0;reactreplacepos[1][2][4][2]=3;
    reactreplaceval[1][2][4][1]=3;reactreplaceval[1][2][4][2]=3;
    
    reactrate[1][2][5]=rohyd2;
    reactinter[1][2][5]=14;
    reactlimitnum[1][2][5]=5;
    reactlimitpos[1][2][5][1]=4;reactlimitpos[1][2][5][2]=6;reactlimitpos[1][2][5][3]=2;reactlimitpos[1][2][5][4]=36;reactlimitpos[1][2][5][5]=38;
    reactlimitval[1][2][5][1]=0;reactlimitval[1][2][5][2]=0;reactlimitval[1][2][5][3]=0;reactlimitval[1][2][5][4]=0;reactlimitval[1][2][5][5]=0;
    reactreplacenum[1][2][5]=2;
    reactreplacepos[1][2][5][1]=0;reactreplacepos[1][2][5][2]=4;
    reactreplaceval[1][2][5][1]=3;reactreplaceval[1][2][5][2]=3;
    
    reactrate[1][2][6]=rohyd2;
    reactinter[1][2][6]=14;
    reactlimitnum[1][2][6]=5;
    reactlimitpos[1][2][6][1]=5;reactlimitpos[1][2][6][2]=1;reactlimitpos[1][2][6][3]=3;reactlimitpos[1][2][6][4]=38;reactlimitpos[1][2][6][5]=40;
    reactlimitval[1][2][6][1]=0;reactlimitval[1][2][6][2]=0;reactlimitval[1][2][6][3]=0;reactlimitval[1][2][6][4]=0;reactlimitval[1][2][6][5]=0;
    reactreplacenum[1][2][6]=2;
    reactreplacepos[1][2][6][1]=0;reactreplacepos[1][2][6][2]=5;
    reactreplaceval[1][2][6][1]=3;reactreplaceval[1][2][6][2]=3;
    
    reactrate[1][2][7]=rohyd2;
    reactinter[1][2][7]=14;
    reactlimitnum[1][2][7]=5;
    reactlimitpos[1][2][7][1]=6;reactlimitpos[1][2][7][2]=2;reactlimitpos[1][2][7][3]=4;reactlimitpos[1][2][7][4]=40;reactlimitpos[1][2][7][3]=42;
    reactlimitval[1][2][7][1]=0;reactlimitval[1][2][7][2]=0;reactlimitval[1][2][7][3]=0;reactlimitval[1][2][7][4]=0;reactlimitval[1][2][7][5]=0;
    reactreplacenum[1][2][7]=2;
    reactreplacepos[1][2][7][1]=0;reactreplacepos[1][2][7][2]=6;
    reactreplaceval[1][2][7][1]=3;reactreplaceval[1][2][7][2]=3;
    
    
    reactnum[1][3]=20;
    
    for(i=1;i<=6;i++)
    {
        j=i+7;
        k=(i+2)%6;
        if (k==0) k=6;
        h=2*k+29;
        jl=k-1;if (jl==0) jl=6;
        jm=k+1;if (jm==7) jm=1;
        jo=h-1;if (jo==30) jo=42;
        jp=h+1;if (jp==43) jp=31;
        reactrate[1][3][j]=rHtrans;
        reactlimitnum[1][3][j]=7;
        reactlimitpos[1][3][j][1]=i;reactlimitpos[1][3][j][2]=k;reactlimitpos[1][3][j][3]=h;reactlimitpos[1][3][j][4]=jl;reactlimitpos[1][3][j][5]=jm;reactlimitpos[1][3][j][6]=jo;reactlimitpos[1][3][j][7]=jp;
        reactlimitval[1][3][j][1]=3;reactlimitval[1][3][j][2]=0;reactlimitval[1][3][j][3]=-3;reactlimitval[1][3][j][4]=-3;reactlimitval[1][3][j][5]=-3;reactlimitval[1][3][j][6]=-3;reactlimitval[1][3][j][7]=-3;
        reactreplacenum[1][3][j]=2;
        reactreplacepos[1][3][j][1]=0;reactreplacepos[1][3][j][2]=k;
        reactreplaceval[1][3][j][1]=0;reactreplaceval[1][3][j][2]=3;
    };
    for(i=1;i<=6;i++)
    {
        j=i+13;
        k=(i+4)%6;
        if (k==0) k=6;
        h=2*k+29;
        jl=k-1;if (jl==0) jl=6;
        jm=k+1;if (jm==7) jm=1;
        jo=h-1;if (jo==30) jo=42;
        jp=h+1;if (jp==43) jp=31;
        reactrate[1][3][j]=rHtrans;
        reactlimitnum[1][3][j]=7;
        reactlimitpos[1][3][j][1]=i;reactlimitpos[1][3][j][2]=k;reactlimitpos[1][3][j][3]=h;reactlimitpos[1][3][j][4]=jl;reactlimitpos[1][3][j][5]=jm;reactlimitpos[1][3][j][6]=jo;reactlimitpos[1][3][j][7]=jp;
        reactlimitval[1][3][j][1]=3;reactlimitval[1][3][j][2]=0;reactlimitval[1][3][j][3]=-3;reactlimitval[1][3][j][4]=-3;reactlimitval[1][3][j][5]=-3;reactlimitval[1][3][j][6]=-3;reactlimitval[1][3][j][7]=-3;
        reactreplacenum[1][3][j]=2;
        reactreplacepos[1][3][j][1]=0;reactreplacepos[1][3][j][2]=k;
        reactreplaceval[1][3][j][1]=0;reactreplaceval[1][3][j][2]=3;
    };
    reactrate[1][3][1]=rhodes;
    reactinter[1][3][1]=2;
    reactreplacenum[1][3][1]=1;
    reactreplacepos[1][3][1][1]=0;
    reactreplaceval[1][3][1][1]=0;
    reactfreenum[1][3][1]=1;
    reactfreename[1][3][1][1]=2;
    reactfreeval[1][3][1][1]=1.0;
    
    reactrate[1][3][2]=rohohr;
    reactinter[1][3][2]=33;
    reactlimitnum[1][3][2]=1;
    reactlimitpos[1][3][2][1]=1;
    reactlimitval[1][3][2][1]=3;
    reactreplacenum[1][3][2]=2;
    reactreplacepos[1][3][2][1]=0;reactreplacepos[1][3][2][2]=1;
    reactreplaceval[1][3][2][1]=0;reactreplaceval[1][3][2][2]=2;
    
    reactrate[1][3][3]=rohohr;
    reactinter[1][3][3]=33;
    reactlimitnum[1][3][3]=1;
    reactlimitpos[1][3][3][1]=2;
    reactlimitval[1][3][3][1]=3;
    reactreplacenum[1][3][3]=2;
    reactreplacepos[1][3][3][1]=0;reactreplacepos[1][3][3][2]=2;
    reactreplaceval[1][3][3][1]=0;reactreplaceval[1][3][3][2]=2;
    
    
    reactrate[1][3][4]=rohohr;
    reactinter[1][3][4]=33;
    reactlimitnum[1][3][4]=1;
    reactlimitpos[1][3][4][1]=3;
    reactlimitval[1][3][4][1]=3;
    reactreplacenum[1][3][4]=2;
    reactreplacepos[1][3][4][1]=0;reactreplacepos[1][3][4][2]=3;
    reactreplaceval[1][3][4][1]=0;reactreplaceval[1][3][4][2]=2;
    
    
    reactrate[1][3][5]=rohohr;
    reactinter[1][3][5]=33;
    reactlimitnum[1][3][5]=1;
    reactlimitpos[1][3][5][1]=4;
    reactlimitval[1][3][5][1]=3;
    reactreplacenum[1][3][5]=2;
    reactreplacepos[1][3][5][1]=0;reactreplacepos[1][3][5][2]=4;
    reactreplaceval[1][3][5][1]=0;reactreplaceval[1][3][5][2]=2;
    
    reactrate[1][3][6]=rohohr;
    reactinter[1][3][6]=33;
    reactlimitnum[1][3][6]=1;
    reactlimitpos[1][3][6][1]=5;
    reactlimitval[1][3][6][1]=3;
    reactreplacenum[1][3][6]=2;
    reactreplacepos[1][3][6][1]=0;reactreplacepos[1][3][6][2]=5;
    reactreplaceval[1][3][6][1]=0;reactreplaceval[1][3][6][2]=2;
    
    reactrate[1][3][7]=rohohr;
    reactinter[1][3][7]=33;
    reactlimitnum[1][3][7]=1;
    reactlimitpos[1][3][7][1]=6;
    reactlimitval[1][3][7][1]=3;
    reactreplacenum[1][3][7]=2;
    reactreplacepos[1][3][7][1]=0;reactreplacepos[1][3][7][2]=6;
    reactreplaceval[1][3][7][1]=0;reactreplaceval[1][3][7][2]=2;
    
    
    reactrate[1][3][20]=rOHoxidation;
    reactinter[1][3][20]=3300;
    reactreplacenum[1][3][20]=1;
    reactreplacepos[1][3][20][1]=0;
    reactreplaceval[1][3][20][1]=2;
    
    reactnum[1][0]=1;
    
    reactrate[1][0][1]=rhoads;
    reactinter[1][0][1]=1;
    /*reactlimitnum[1][0][1]=18;
     for (i=1;i<=6;i++)
     {
     reactlimitpos[1][0][1][i]=i;
     reactlimitval[1][0][1][i]=-3;
     };
     for (i=1;i<=12;i++)
     {
     reactlimitpos[1][0][1][i+6]=i+30;
     reactlimitval[1][0][1][i+6]=-3;
     };*/
    
    reactreplacenum[1][0][1]=1;
    reactreplacepos[1][0][1][1]=0;
    reactreplaceval[1][0][1][1]=3;
    reactfreenum[1][0][1]=1;
    reactfreename[1][0][1][1]=2;
    reactfreeval[1][0][1][1]=-1.0;
    
};                                                                              // End original data
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Editable//////////////////////////
void firstrand()
{srand(time(NULL));};
////////////////////////////////////////////////////////////////////////////////
double getrand()
{
    double r1,r2,r3,rr,ra;
    r1=RAND_MAX;
    r2=RAND_MAX;
    r3=RAND_MAX;
    
    while(r1==RAND_MAX) r1=rand()*1.0;
    while(r2==RAND_MAX) r2=rand()*1.0;
    while(r3==RAND_MAX) r3=rand()*1.0;
    
    rr=r1*RAND_MAX*RAND_MAX+r2*RAND_MAX+r3;
    ra=rr/RAND_MAX/RAND_MAX/RAND_MAX;
    return ra;
};
////////////////////////////////////////////////////////////////////////////////
void sitecover()                                                                //Cover the sites with different elements
{int i,j,k,h;
    for(i=1;i<=totalsite;i++) (site[i]=0);
    for (i=1;i<=xdim;i++)
        for (j=1;j<=ydim;j++)
        {
            k=i-j;
            if (k<0) k=0-k;
            if ((k%3)==0)
            {
                
                h=position[i][j][1];
                site[h]=10;
            };
        };
    
    
    //while
};
////////////////////////////////////////////////////////////////////////////////
int chk(int st, int rcttp)
{
    int i,sttp,cmpnt;
    sttp=sitetype[st];
    cmpnt=site[st];
    for (i=1;i<=reactlimitnum[sttp][cmpnt][rcttp];i++)
    {if ((reactlimitval[sttp][cmpnt][rcttp][i]>=0)&&(site[sn[st][reactlimitpos[sttp][cmpnt][rcttp][i]]]!=reactlimitval[sttp][cmpnt][rcttp][i])) return (0);
        if ((reactlimitval[sttp][cmpnt][rcttp][i]<0)&&(site[sn[st][reactlimitpos[sttp][cmpnt][rcttp][i]]]==(0-reactlimitval[sttp][cmpnt][rcttp][i]))) return (0);};
    return (1);
};
////////////////////////////////////////////////////////////////////////////////
double adjust_site_inter(int posR, int code_site)
{
    int i;
    double sum_OH,c_OH,aFactor,t_rrate,t_ea;
    sum_OH=0;
    
    if (code_site==2)//OH des
    {
        for (i=0;i<=42;i++) if (site[sn[posR][i]]==3) sum_OH=sum_OH+1;
        c_OH=sum_OH/19;
        xOH=c_OH;
        
        deltOHbinding=1.4595*xOH;
        
        t_ea=((dGOHdes-deltOHbinding+1*po)*eFactor+eIntercept)*96486.9;
        
        if (t_ea<(eIntercept*96486.9)) t_ea=eIntercept*96486.9;
        t_rrate=ke*exp(((-1)*t_ea)/(8.314*temp));
        aFactor=t_rrate/rhodes;
        
        return(aFactor);
    };
    
    if (code_site==7)//O2 ads
    {
        
        return(1);
    };
    
    if (code_site==238)//O2 dis
    {
        for (i=0;i<=56;i++) if (site[sn[posR][i]]==3) sum_OH=sum_OH+1;
        c_OH=sum_OH/24;
        xOH=c_OH;
        Eao2dis=3.195*xOH+ActEO2dis;
        if (Eao2dis<0) Eao2dis=0;
        aFactor=exp(-96486.9*(Eao2dis)/(8.314*temp));
        return(aFactor);
    };
    if (code_site==580)
    {
        for (i=0;i<=56;i++) if (site[sn[posR][i]]==3) sum_OH=sum_OH+1;
        c_OH=sum_OH/24;
        xOH=c_OH;
        Eao2hyd2=ActEO2pro+switch_Surf;
        if (xOH>0.11111) Eao2hyd2=9999999999999;
        aFactor=exp(-96486.9*(Eao2hyd2)/(8.314*temp));
        return(aFactor);
    };
    
    if (code_site==58)
    {
        for (i=0;i<=56;i++) if (site[sn[posR][i]]==3) sum_OH=sum_OH+1;
        c_OH=sum_OH/24;
        xOH=c_OH;
        deltOOHbinding=0.7928*xOH;
        deltO2binding=-8.5306*xOH*xOH*xOH-8.522*xOH*xOH+6.2763*xOH;
        t_ea=((dGO2pro+deltOOHbinding-deltO2binding+1*po)*eFactor+eIntercept+switch_Elec)*96486.9;
        if (t_ea<(eIntercept*96486.9)) t_ea=eIntercept*96486.9;
        t_rrate=ke*exp(((-1)*t_ea)/(8.314*temp));
        aFactor=t_rrate/ro2hyd;
        return(aFactor);
    };
    
    if (code_site==8)
    {
        for (i=0;i<=56;i++) if (site[sn[posR][i]]==3) sum_OH=sum_OH+1;
        c_OH=sum_OH/24;
        xOH=c_OH;
        deltO2binding=-8.5306*xOH*xOH*xOH-8.522*xOH*xOH+6.2763*xOH;
        t_ea=96486.9*((-1*dGO2ads)-(1*deltO2binding));
        t_rrate=ksf*exp(((-1)*t_ea)/(8.314*temp));
        aFactor=t_rrate/ro2des;
        return(aFactor);
    };
    
    if (code_site==136)
    {
        for (i=0;i<=56;i++) if (site[sn[posR][i]]==3) sum_OH=sum_OH+1;
        c_OH=sum_OH/24;
        xOH=c_OH;
        Eaho2dis=4.59*xOH+ActEO2Hdis;
        aFactor=exp(-96486.9*(Eaho2dis)/(8.314*temp));
        return(aFactor);
    };
    if (code_site==6)
    {
        for (i=0;i<=56;i++) if (site[sn[posR][i]]==3) sum_OH=sum_OH+1;
        c_OH=sum_OH/24;
        xOH=c_OH;
        deltOOHbinding=0.7928*xOH;
        t_ea=((dGOOHdes-deltOOHbinding+1*po)*eFactor+eIntercept+switch_Elec)*96486.9;
        if (t_ea<(eIntercept*96486.9)) t_ea=eIntercept*96486.9;
        t_rrate=ke*exp(((-1)*t_ea)/(8.314*temp));
        aFactor=t_rrate/rho2des;
        return(aFactor);
    };
    
    if (code_site==116)
    {
        for (i=0;i<=56;i++) if (site[sn[posR][i]]==3) sum_OH=sum_OH+1;
        c_OH=sum_OH/24;
        xOH=c_OH;
        Eaho2hyd2=ActEO2Hpro+switch_Surf+99999999999;
        aFactor=exp(-96486.9*(Eaho2hyd2)/(8.314*temp));
        return(aFactor);
    };
    if (code_site==76)
    {
        for (i=0;i<=56;i++) if (site[sn[posR][i]]==3) sum_OH=sum_OH+1;
        if (sum_OH>=1) sum_OH--;
        c_OH=sum_OH/24;
        xOH=c_OH;
        Eaho2hor=-1.227*xOH+ActHO2HOr+switch_Surf;
        if (Eaho2hor<0) Eaho2hor=0;
        aFactor=exp(-96486.9*(Eaho2hor)/(8.314*temp));
        return(aFactor);
    };
    if (code_site==36)
    {
        for (i=0;i<=56;i++) if (site[sn[posR][i]]==3) sum_OH=sum_OH+1;
        c_OH=sum_OH/24;
        xOH=c_OH;
        deltOOHbinding=0.7928*xOH;
        deltObinding=1.991*xOH;
        t_ea=((dGOOHdis+dGOHdes+deltObinding-deltOOHbinding+1*po)*eFactor+eIntercept+switch_Elec)*96486.9;//+999999999999999999999999999999;
        if (t_ea<(eIntercept*96486.9)) t_ea=eIntercept*96486.9;
        t_rrate=ke*exp(((-1)*t_ea)/(8.314*temp));
        aFactor=t_rrate/rho2dis2;
        return(aFactor);
    };
    if (code_site==140)
    {
        for (i=0;i<=42;i++) if (site[sn[posR][i]]==3) sum_OH=sum_OH+1;
        c_OH=sum_OH/19;
        xOH=c_OH;
        deltObinding=1.991*xOH;
        deltOHbinding=1.4595*xOH;
        
        t_ea=((dGOpro-deltObinding+deltOHbinding+1*po)*eFactor+eIntercept+switch_Elec)*96486.9;
        
        if (t_ea<(eIntercept*96486.9)) t_ea=eIntercept*96486.9;
        t_rrate=ke*exp(((-1)*t_ea)/(8.314*temp));
        aFactor=t_rrate/rohyd;
        return(aFactor);
    };
    if (code_site==14)
    {
        for (i=0;i<=42;i++) if (site[sn[posR][i]]==3) sum_OH=sum_OH;
        c_OH=sum_OH/19;
        xOH=c_OH;
        Eaohyd=1.953*xOH+ActEOpro+switch_Surf;
        if (Eaohyd<0) Eaohyd=0;
        aFactor=exp(-96486.9*Eaohyd/(8.314*temp));
        //aFactor=0;
        return(aFactor);
    };
    if (code_site==33)
    {
        for (i=0;i<=42;i++) if (site[sn[posR][i]]==3) sum_OH=sum_OH+1;
        if (sum_OH>2) sum_OH=sum_OH-2;
        c_OH=sum_OH/19;
        xOH=c_OH;
        Eaohohr=-2.574*xOH+Actohohr+switch_Surf;
        if (Eaohohr<0) Eaohohr=0;
        aFactor=exp(-96486.9*Eaohohr/(8.314*temp));
        //aFactor=0;
        return(aFactor);
    };
    if (code_site==3300)
    {
        for (i=0;i<=42;i++) if (site[sn[posR][i]]==3) sum_OH=sum_OH+1;
        c_OH=sum_OH/19;
        xOH=c_OH;
        deltOHbinding=1.4595*xOH;
        deltObinding=1.991*xOH;
        t_ea=((-dGOpro+deltObinding-deltOHbinding-1*po)*eFactor+eIntercept+switch_Elec)*96486.9;
        if (t_ea<(eIntercept*96486.9)) t_ea=eIntercept*96486.9;
        t_rrate=ke/550*exp(((-1)*t_ea)/(8.314*temp));
        aFactor=t_rrate/rOHoxidation;
        return(aFactor);
    };
    
    if (code_site==2)
    {
        for (i=0;i<=42;i++) if (site[sn[posR][i]]==3) sum_OH=sum_OH+1;
        if (sum_OH>1) sum_OH=sum_OH-1;
        c_OH=sum_OH/19;
        xOH=c_OH;
        deltOHbinding=1.4595*xOH;
        t_ea=((dGOHdes-deltOHbinding+1*po)*eFactor+eIntercept)*96486.9;
        if (t_ea<(eIntercept*96486.9)) t_ea=eIntercept*96486.9;
        t_rrate=ke*exp(((-1)*t_ea)/(8.314*temp));
        aFactor=t_rrate/rhodes;
        
        // if (((step%10000)==0)&(posR>1000)&(posR<1200)) {printf("%f",xOH);printf("%s","XXXXXX");};
        
        
        
        return(aFactor);
    };
    if (code_site==1)
    {
        for (i=0;i<=42;i++) if (site[sn[posR][i]]==3) sum_OH=sum_OH+1;
        c_OH=sum_OH/19;
        xOH=c_OH;
        deltOHbinding=1.4595*xOH;
        t_ea=((-1*dGOHdes+deltOHbinding-1*po)*eFactor+eIntercept)*96486.9;
        if (t_ea<(eIntercept*96486.9)) t_ea=eIntercept*96486.9;
        t_rrate=ke/550*exp(((-1)*t_ea)/(8.314*temp));
        aFactor=t_rrate/rhoads;
        
        return(aFactor);
    };
    
    return(1);
};
////////////////////////////////////////////////////////////////////////////////
double calculaterate(int ss)
{
    int i,cmpnt,sttp,crctnum,acode;
    double  tmrate,temp_rate;
    sttp=sitetype[ss];
    cmpnt=site[ss];
    crctnum=0;
    tmrate=0;
    for (i=1;i<=reactnum[sttp][cmpnt];i++)
        if (chk(ss,i)==1)
        {crctnum++;
            effectreactname[ss][crctnum]=i;
            temp_rate=reactrate[sttp][cmpnt][i]*adjust_site_inter(ss,reactinter[sttp][cmpnt][i]);
            effectreactrate[ss][crctnum]=temp_rate;
            tmrate=tmrate+temp_rate;};
    effectreactnum[ss]=crctnum;
    return(tmrate);
};
////////////////////////////////////////////////////////////////////////////////
void binarytree()
{int i,j;                                                                       //Initialization of tree
    i=1;
    j=2;
    while (j<totalsite)                                                        //Tree level
    {i=i+1;
        j=j*2;};
    treelevel=i;
    treetop=j-1;                                                               //Tree top
    totaltreesite=totalsite+treetop;
    for(i=1;i<=totalsite;i++)
        treeposition[i]=i+treetop;                                             //Sites' position in tree
    for(i=1;i<=totaltreesite;i++)
        treerate[i]=0.0;
    for(i=1;i<=totalsite;i++)
    {treerate[treeposition[i]]=calculaterate(i);
    };
    for(i=totaltreesite;i>=2;i--)
    {
        treerate[i/2]=treerate[i/2]+treerate[i];
        
    };
    
};
////////////////////////////////////////////////////////////////////////////////
void selectsite()
{
    double rrate;
    int i,j,k;
    k=1;
    rrate=treerate[k]*getrand();
    for (i=1;i<=treelevel;i++)
    {j=k;
        if (rrate<=treerate[2*j]){k=k*2;};
        if (rrate>treerate[2*j]){k=k*2+1;rrate=rrate-treerate[2*j];};
        
    };
    selectedtreesite=k;
    selectedsite=k-treetop;
    if (selectedsite>totalsite) selectsite();
    if (selectedsite<1) selectsite();
    if (treerate[selectedtreesite]==0) selectsite();
};
////////////////////////////////////////////////////////////////////////////////
void changefree(int fsttp, int fcmpnt, int frctnum)
{
    int i;
    for (i=1;i<=reactfreenum[fsttp][fcmpnt][frctnum];i++)
        freequantity[reactfreename[fsttp][fcmpnt][frctnum][i]]=freequantity[reactfreename[fsttp][fcmpnt][frctnum][i]]+reactfreeval[fsttp][fcmpnt][frctnum][i];
    
};
////////////////////////////////////////////////////////////////////////////////
void refresh( int rbst)
{int i,rtp,rnn,tst;
    rtp=sitetype[rbst];
    rnn=refreshneighbornum[rtp];
    for (i=1;i<=rnn;i++)
    {tst=sn[rbst][refreshneighbor[rtp][i]];
        if (changedchk[tst]!=step)
        {changednum++;
            treerate[treeposition[tst]]=calculaterate(tst);
            changedsite[changednum]=treeposition[tst];
            changedchk[tst]=step;
        };
        
        
    };
};
////////////////////////////////////////////////////////////////////////////////
void changestcmpnt(int ssttp, int scmpnt, int srctnum)
{
    int i;
    for (i=1;i<=reactreplacenum[ssttp][scmpnt][srctnum];i++)
        site[sn[selectedsite][reactreplacepos[ssttp][scmpnt][srctnum][i]]]=reactreplaceval[ssttp][scmpnt][srctnum][i];
    
    changednum=0;
    
    for (i=1;i<=reactreplacenum[ssttp][scmpnt][srctnum];i++)
        refresh(sn[selectedsite][reactreplacepos[ssttp][scmpnt][srctnum][i]]);
};
////////////////////////////////////////////////////////////////////////////////
void updatebinarytree()
{
    int i,j,k;
    
    for (i=1;i<=changednum;i++)
        for (j=1;j<=treelevel;j++)
        {k=changedsite[i]/2;
            treerate[k]=treerate[2*k]+treerate[2*k+1];
            changedsite[i]=k;
        };
    
    
    
};
////////////////////////////////////////////////////////////////////////////////
void sitereact()
{   int cmpnt,sttp,erctnum,rctnum;
    
    double randrate,tsrate,temprate;
    sttp=sitetype[selectedsite];
    cmpnt=site[selectedsite];
    tsrate=treerate[selectedtreesite];
    randrate=tsrate*getrand();
    temprate=0.0;
    erctnum=0;
    while ((temprate<=randrate)&&(erctnum<effectreactnum[selectedsite]))
    {erctnum++;
        temprate=temprate+effectreactrate[selectedsite][erctnum];
    };
    rctnum=effectreactname[selectedsite][erctnum];
    //if (step==1000) {printf("%d",sttp);printf("%s","   ");printf("%d",cmpnt);printf("%s","   ");printf("%d",rctnum);printf("%s","   ");};
    reactsum[sttp][cmpnt][rctnum]++;
    changefree(sttp,cmpnt,rctnum);
    changestcmpnt(sttp,cmpnt,rctnum);
    updatebinarytree();
    
};
//////////////////////////////d//////////////////////////////////////////////////


int main()
{
    int i,te,ipo,NoIt,MaxofIt,It;
    double tt,ic,rsO2dis,rsO2hyd;
    double aric[101];
    double aric2[101];
    double aric3[101];
    double aric4[101];
    double aric5[101];
    double aric6[101];
    double aric7[101];
    double aric8[101];
    double aric9[101];
    double aric10[101];
    double aric11[101];
    double aric12[101];
    double aric13[101];
    double aric14[101];
    double aric15[101];
    double aric16[101];
    double aric17[101];
    double aric18[101];
    double tsr;
    firstrand();
    MaxofIt=5;
    for (i=1;i<=100;i++) aric[i]=0;
    for (i=1;i<=100;i++) aric2[i]=0;
    for (i=1;i<=100;i++) aric3[i]=0;
    for (i=1;i<=100;i++) aric4[i]=0;
    for (i=1;i<=100;i++) aric5[i]=0;
    for (i=1;i<=100;i++) aric6[i]=0;
    for (i=1;i<=100;i++) aric7[i]=0;
    for (i=1;i<=100;i++) aric8[i]=0;
    for (i=1;i<=100;i++) aric9[i]=0;
    for (i=1;i<=100;i++) aric10[i]=0;
    for (i=1;i<=100;i++) aric11[i]=0;
    for (i=1;i<=100;i++) aric12[i]=0;
    for (i=1;i<=100;i++) aric13[i]=0;
    for (i=1;i<=100;i++) aric14[i]=0;
    for (i=1;i<=100;i++) aric15[i]=0;
    for (i=1;i<=100;i++) aric16[i]=0;
    for (i=1;i<=100;i++) aric17[i]=0;
    for (i=1;i<=100;i++) aric18[i]=0;
    for(NoIt=1;NoIt<=MaxofIt;NoIt++)
    {
        
        ppo=-0.64;
        
        initialization();
        sitecover();
        It=0;
        for (ipo=1;ipo<=54;ipo++)
        {
            It++;
            initialization();                                                       //This initialization is for Structrue
            originaldata();                                                         //Input the original data
            binarytree();
            tt=0;
            step=1;
            while (step<500000)
            {
                selectsite();
                tt=tt+1/treerate[1];
                sitereact();
                step++;
            };
            
            for(i=0;i<=coverspecies;i++)
                componentcoverage[i]=0;
            for(i=1;i<=totalsite;i++)
                componentcoverage[site[i]]++;
            te=reactsum[4][0][1]-reactsum[4][1][3];
            ic=4*te/tt;
            //printf("%f",ic/220000);printf("%s","   ");
            //printf("%3f",tt);printf("%s","   ");
            printf("%d",reactsum[4][4][4]);printf("%s","   ");
            //rsO2dis=reactsum[4][1][1];
            //rsO2hyd=reactsum[4][1][2]+reactsum[4][1][8]+reactsum[4][1][4]+reactsum[4][1][5]+reactsum[4][1][6]+reactsum[4][1][7]-reactsum[4][4][9]-reactsum[4][4][10]-reactsum[4][4][11]-reactsum[4][4][12];
            printf("%f",componentcoverage[2]);printf("%s","   ");
            printf("%f",componentcoverage[4]);printf("%s","   ");
            te=reactsum[4][1][2]+reactsum[4][1][8]+reactsum[4][4][3]+reactsum[4][4][4]+reactsum[4][4][13]+reactsum[4][4][14]+reactsum[1][2][1]+reactsum[1][3][1]-reactsum[1][0][1]-reactsum[1][3][20];
            ic=te/tt;
            //tt=1;
            printf("%f",ppo+0.83);printf("%s","   ");
            printf("%f",tt);printf("%s","   ");
            printf("%f\n",ic/(-233965));
            tsr=reactsum[4][1][1];// O2dis
            aric[It]=aric[It]+(tsr/tt);//aric[It]=componentcoverage[3];
            tsr=reactsum[4][1][2]+reactsum[4][1][8];//O2 pro Ele
            aric2[It]=aric2[It]+(tsr/tt);//aric2[It]=componentcoverage[2];
            tsr=reactsum[4][1][4]+reactsum[4][1][5]+reactsum[4][1][6]+reactsum[4][1][7];//O2 pro Suf
            aric3[It]=aric3[It]+(tsr/tt);//aric3[It]+ic/(-233965);
            tsr=reactsum[4][4][1]+reactsum[4][4][2];//O2H dis suf
            aric4[It]=aric4[It]+(tsr/tt);
            tsr=reactsum[4][4][4];//O2H des
            aric5[It]=aric5[It]+(tsr/tt);
            tsr=reactsum[4][4][5]+reactsum[4][4][6]+reactsum[4][4][7]+reactsum[4][4][8];//O2H pro surf
            aric6[It]=aric6[It]+(tsr/tt);
            tsr=reactsum[4][4][9]+reactsum[4][4][10]+reactsum[4][4][11]+reactsum[4][4][12];//O2H + OH-->H2O+O2 surf
            aric7[It]=aric7[It]+(tsr/tt);
            tsr=0; for (i=8;i<=19;i++) tsr=tsr+reactsum[1][3][i];//H trans
            aric8[It]=aric8[It]+(tsr/tt);
            tsr=reactsum[1][2][1];//*O-->*OH ele
            aric9[It]=aric9[It]+(tsr/tt);
            tsr=reactsum[1][2][2]+reactsum[1][2][3]+reactsum[1][2][4]+reactsum[1][2][5]+reactsum[1][2][6]+reactsum[1][2][7];//O-->OH sur
            aric10[It]=aric10[It]+(tsr/tt);
            tsr=reactsum[4][4][13]+reactsum[4][4][14];//O2H-->*O + OH-
            aric11[It]=aric11[It]+(tsr/tt);
            tsr=reactsum[1][3][1];//OH des
            aric12[It]=aric12[It]+(tsr/tt);
            tsr=reactsum[1][0][1];//OH ads
            aric13[It]=aric13[It]+(tsr/tt);
            tsr=reactsum[1][3][2]+reactsum[1][3][3]+reactsum[1][3][4]+reactsum[1][3][5]+reactsum[1][3][6]+reactsum[1][3][7];//OH_OHR
            aric14[It]=aric14[It]+(tsr/tt);
            tsr=reactsum[1][3][20];//OH oxidation
            aric15[It]=aric15[It]+(tsr/tt);
            aric16[It]=aric16[It]+componentcoverage[3];
            aric17[It]=aric17[It]+componentcoverage[2];
            aric18[It]=aric18[It]+ic/(-233965);
            ppo=ppo+0.02;
        };
        
    };
    
    for(i=1;i<=100;i++) aric[i]=aric[i]/MaxofIt;
    for(i=1;i<=100;i++) aric2[i]=aric2[i]/MaxofIt;
    for(i=1;i<=100;i++) aric3[i]=aric3[i]/MaxofIt;
    for(i=1;i<=100;i++) aric4[i]=aric4[i]/MaxofIt;
    for(i=1;i<=100;i++) aric5[i]=aric5[i]/MaxofIt;
    for(i=1;i<=100;i++) aric6[i]=aric6[i]/MaxofIt;
    for(i=1;i<=100;i++) aric7[i]=aric7[i]/MaxofIt;
    for(i=1;i<=100;i++) aric8[i]=aric8[i]/MaxofIt;
    for(i=1;i<=100;i++) aric9[i]=aric9[i]/MaxofIt;
    for(i=1;i<=100;i++) aric10[i]=aric10[i]/MaxofIt;
    for(i=1;i<=100;i++) aric11[i]=aric11[i]/MaxofIt;
    for(i=1;i<=100;i++) aric12[i]=aric12[i]/MaxofIt;
    for(i=1;i<=100;i++) aric13[i]=aric13[i]/MaxofIt;
    for(i=1;i<=100;i++) aric14[i]=aric14[i]/MaxofIt;
    for(i=1;i<=100;i++) aric15[i]=aric15[i]/MaxofIt;
    for(i=1;i<=100;i++) aric16[i]=aric16[i]/MaxofIt;
    for(i=1;i<=100;i++) aric17[i]=aric17[i]/MaxofIt;
    for(i=1;i<=100;i++) aric18[i]=aric18[i]/MaxofIt;
    
    ofstream fout;
    fout.open("output.txt");
    if (fout.is_open())
    {
        for (i=1;i<=54;i++)
        {
            fout<<aric[i];fout<<" ";
            fout<<aric2[i];fout<<" ";
            fout<<aric3[i];fout<<" ";
            fout<<aric4[i];fout<<" ";
            fout<<aric5[i];fout<<" ";
            fout<<aric6[i];fout<<" ";
            fout<<aric7[i];fout<<" ";
            fout<<aric8[i];fout<<" ";
            fout<<aric9[i];fout<<" ";
            fout<<aric10[i];fout<<" ";
            fout<<aric11[i];fout<<" ";
            fout<<aric12[i];fout<<" ";
            fout<<aric13[i];fout<<" ";
            fout<<aric14[i];fout<<" ";
            fout<<aric15[i];fout<<" ";
            fout<<aric16[i];fout<<" ";
            fout<<aric17[i];fout<<" ";
            fout<<aric18[i];fout<<" ";
            fout<<"\n";
        };
        
        
        
        
        
        
        for (i=1;i<=6534;i++)
            if (site[i]==3) {fout<<coordinatesx[i]+coordinatesy[i]*0.5;fout<<" ";fout<<coordinatesy[i];fout<<"\n";};
        
        
        for (i=1;i<=6534;i++)
            if ((site[i]==0)&&(sitetype[i]==1)) {fout<<coordinatesx[i]+coordinatesy[i]*0.5;fout<<" H2O ";fout<<coordinatesy[i];fout<<"\n";};
        
        
        for (i=1;i<=6534;i++)
            if (site[i]==10) {fout<<coordinatesx[i]+coordinatesy[i]*0.5;fout<<" Pt Pt ";fout<<coordinatesy[i];fout<<"\n";};
        
        
        for (i=1;i<=6534;i++)
            if (site[i]==5) {fout<<coordinatesx[i]+coordinatesy[i]*0.5;fout<<" O2 O2 O2 ";fout<<coordinatesy[i];fout<<"\n";};
        
        
        for (i=1;i<=6534;i++)
            if (site[i]==2) {fout<<coordinatesx[i]+coordinatesy[i]*0.5;fout<<" O O O O ";fout<<coordinatesy[i];fout<<"\n";};
        
        for (i=1;i<=6534;i++)
            if ((site[i]==6)||(site[i]==7)) {fout<<coordinatesx[i]+coordinatesy[i]*0.5;fout<<" O2H O2H O2H O2H O2H ";fout<<coordinatesy[i];fout<<"\n";};
        fout<<"\n";
        
        
        
        fout<<flush;
        fout.close();
    };
    
    return 0;
}

