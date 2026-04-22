function solve_mulhyp(xyt_range,Nmmtp,equ_para,init,bc_type,bc_para,shape)
if bc_type(1)==1
    bc_para.g_left=@(x,y,t)0;
    bc_para.q_left=@(x,y,t)0;
else
    bc_para.h_left=@(x,y,t)0;
    bc_para.r_left=@(x,y,t)0;
end
if bc_type(2)==1
    bc_para.g_right=@(x,y,t)0;
    bc_para.q_right=@(x,y,t)0;
else
    bc_para.h_right=@(x,y,t)0;
    bc_para.r_right=@(x,y,t)0;
end
if bc_type(3)==1
    bc_para.g_bottom=@(x,y,t)0;
    bc_para.q_bottom=@(x,y,t)0;
else
    bc_para.h_bottom=@(x,y,t)0;
    bc_para.r_bottom=@(x,y,t)0;
end
if bc_type(4)==1
    bc_para.g_top=@(x,y,t)0;
    bc_para.q_top=@(x,y,t)0;
else
    bc_para.h_top=@(x,y,t)0;
    bc_para.r_top=@(x,y,t)0;
end
switch shape
    case 'circle'
        solve_mulhyp_circle(xyt_range,Nmmtp,equ_para,init,bc_type,bc_para)
    case 'square'
        solve_mulhyp_square(xyt_range,Nmmtp,equ_para,init,bc_type,bc_para)
    otherwise
        disp('无该形状');
end
end