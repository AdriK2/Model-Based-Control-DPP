

double V =Vinit;
float T1=T0+273;
float B= -0.28*Voc/100;
float a1= q/(Ns*k*T1*a);
 // saturation current
float Is1=(Isc+Isc*Rs/Rp-(Voc+B*(T1-298))/Rp)/(exp(q*(Voc+B*(T1-298))/(Ns*k*T1*a))-exp(q*Isc*Rs/(Ns*k*T1*a)));


//int cnt = 0;  // Initialize the counter
//int cntlim = 10; // Set the counter limit

//do{ 
//for(cnt =0; cnt<=cntlim; cnt++) {
static int cnt = 0;
static int cntlim = 50;

//y = cnt;

if (cnt <= cntlim) {
            
  for(i=0;i<n;i++){
double f= Im1 + V*( -Is1*a1*exp(a1*(V+Im1*Rs))-1/Rsh)/(Rs/Rsh+1+a1*Is1*Rs*exp(a1*(V+Im1*Rs)) );
double f1= 2*(-Is1*a1*exp(a1*(V+Im1*Rs))-1/Rsh)/(Rs/Rsh+1+a1*Is1*Rs*exp(a1*(V+Im1*Rs))) + V*(   ((Rs/Rsh+1+a1*Is1*Rs*exp(a1*(V+Im1*Rs)))*(-Is1*a1*a1*exp(a1*(V+Im1*Rs))-a1*a1*Is1*Rs*
               (-Is1*a1*exp(a1*(V+Im1*Rs))-1/Rsh)/(Rs/Rsh+1+a1*Is1*Rs*exp(a1*(V+Im1*Rs)))*exp(a1*(V+Im1*Rs))) ) - ( (-a1*Is1*exp(a1*(V+Im1*Rs))-1/Rsh)) *( a1*a1*Is1*Rs*exp(a1*(V+Im1*Rs))+a1*a1*Is1*Rs*Rs*
               (-Is1*a1*exp(a1*(V+Im1*Rs))-1/Rsh)/(Rs/Rsh+1+a1*Is1*Rs*exp(a1*(V+Im1*Rs)))*exp(a1*(V+Im1*Rs)) )    
         / ( (Rs/Rsh+1+a1*Is1*Rs*exp(a1*(V+Im1*Rs)))*(Rs/Rsh+1+a1*Is1*Rs*exp(a1*(V+Im1*Rs))))  );

    V1= V-(f/f1);
        // 收斂條件
        if(  fabs(V1-V)<tol )
        {
           break;
        }
     
    V=V1;   // 更新 V
//        y2 = cnt; 
      // y2 = V; 
  } 

    cnt++;     // 外層計數器遞增
   
    //y3=cnt; 
    //y=V;
    Vout0=V;
} 
else {  

Vout=Vout0;  // Final Output    // 計算完成，輸出結果
    //cnt = 0;  // 回到0，繼續循環
  
     //y5=cnt;
     // y4=V;
}    
   //  y1 = V;


  


