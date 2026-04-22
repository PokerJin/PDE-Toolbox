clc,clear;
%设定求解区域
xyt_range=[0,1,0,1,0,0.5];
%设定差分网格
Nxyt=[50,50,100];
%设定方程参数
equ_para.d=@(x,y,t)1;
equ_para.c=@(x,y,t)1;
equ_para.a=@(x,y,t)0;
equ_para.f=@(x, y, t) 2*pi^2*exp(-t).*sin(pi*x).*sin(pi*y);
%设定初始条件
u0=@(x,y)sin(pi*x) .* sin(pi*y);
%设定边界条件
bc_type=[1,1,1,1];
bc_para.h_left=@(x,y,t)1;
bc_para.r_left=@(x,y,t)0;
bc_para.h_right=@(x,y,t)1;
bc_para.r_right=@(x,y,t)0;
bc_para.h_bottom=@(x,y,t)1;
bc_para.r_bottom=@(x,y,t)0;
bc_para.h_top=@(x,y,t)1;
bc_para.r_top=@(x,y,t)0;
%设定求解方法
method='CN';
solve_parabolic(xyt_range,Nxyt, equ_para,u0,bc_type,bc_para,method)
