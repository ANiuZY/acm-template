//三角形
//已知三边平方，计算夹角
double calR(double a2,double b2,double c2){
	return acos((a2+b2-c2)/(2*sqrt(a2*b2)));
}

//已知三边，计算底边高
double calH(double a,double b,double c){
    double r=calR(a*a,b*b,c*c);
    return a*b*sin(r)/c;
}

//点线
//二维坐标旋转
void calC(double& x,double& y,double r){
    double tx=t,ty=y;
    x=cos(r)*tx-sin(r)*ty;
    y=sin(r)*tx+cos(r)*ty;
}
