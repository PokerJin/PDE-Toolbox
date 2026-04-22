function solve_mulell_dir(xy_range,Nxy,equ_para,bc_type,bc_para,method)
if bc_type(1)==1
    bc_para.g_left=@(x,y)0;
    bc_para.q_left=@(x,y)0;
else
    bc_para.h_left=@(x,y)0;
    bc_para.r_left=@(x,y)0;
end
if bc_type(2)==1
    bc_para.g_right=@(x,y)0;
    bc_para.q_right=@(x,y)0;
else
    bc_para.h_right=@(x,y)0;
    bc_para.r_right=@(x,y)0;
end
if bc_type(3)==1
    bc_para.g_bottom=@(x,y)0;
    bc_para.q_bottom=@(x,y)0;
else
    bc_para.h_bottom=@(x,y)0;
    bc_para.r_bottom=@(x,y)0;
end
if bc_type(4)==1
    bc_para.g_top=@(x,y)0;
    bc_para.q_top=@(x,y)0;
else
    bc_para.h_top=@(x,y)0;
    bc_para.r_top=@(x,y)0;
end
switch method
    case 'fivepoint'
        solve_mulell_dirfive(xy_range,Nxy,equ_para,bc_type,bc_para)
    case 'ninepoint'
        solve_mulell_dirnine(xy_range,Nxy,equ_para,bc_type,bc_para)
    otherwise
        disp('输入错误');
end


end