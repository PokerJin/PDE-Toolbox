# PDE-Toolbox

面向材料计算的偏微分方程数值仿真工具箱。

PDE-Toolbox 面向材料学、数学、力学及相关工程计算场景，提供一套轻量化、可交互、可扩展的 PDE 数值仿真工具，重点支持复合材料、多孔介质等非均匀材料中的热传导、波传播及多尺度响应分析。

## 项目背景

偏微分方程广泛出现在科学计算与工程仿真中，例如传热、力学响应、波传播、流体动力学以及材料性能预测。传统工业仿真软件功能强大，但往往依赖国外商业软件生态；国产科学计算平台仍需要更多面向具体工程问题的高层工具箱支持。

本工具箱基于北太天元/类 MATLAB 数值计算环境开发，围绕经典 PDE 与多尺度 PDE 的建模、求解和可视化流程进行封装，目标是降低算法实现门槛，让用户能够通过参数配置和函数调用快速完成数值实验。

## 主要功能

| 模块 | 支持内容 | 说明 |
| --- | --- | --- |
| 经典椭圆型方程 | 五点差分、九点差分 | 支持二维区域上的 Dirichlet/Neumann 边界条件配置 |
| 经典抛物型方程 | 向前差分、向后差分、Crank-Nicolson、ADI | 支持时间范围、初值条件、方程系数与边界条件配置 |
| 经典双曲型方程 | 蛙跳格式 | 支持位移初值、速度初值和时空网格配置 |
| 多尺度椭圆型方程 | 参考解、高阶多尺度解 | 支持圆形/方形夹杂相、材料参数与尺度比设置 |
| 多尺度抛物型方程 | 参考解、高阶多尺度解 | 面向含快速振荡系数的瞬态热传导类问题 |
| 多尺度双曲型方程 | 参考解、高阶多尺度解 | 面向非均匀材料中的波传播类问题 |
| 可视化 | 二维等值线、三维数值解图像 | 便于观察解的空间分布与多尺度振荡特征 |

## 技术特点

- 集成有限差分法，覆盖椭圆型、抛物型、双曲型 PDE 的基础数值求解。
- 集成高阶渐近均匀化方法，用于求解具有快速振荡系数的二维多尺度 PDE。
- 同时提供直接精细数值模拟参考解和高阶多尺度近似解，便于精度对比。
- 支持自定义方程系数、源项、求解区域、网格规模、边界条件和初始条件。
- 针对复合材料和多孔介质等非均匀材料，可设置夹杂相、基体相、几何形状和尺度比。
- 示例材料中的测试结果显示，典型算例的 L2 相对误差可控制在 5% 以下，能够满足一般教学演示与数值实验需求。

## 运行环境

推荐环境：

- Windows 10/11 64 位
- 北太天元 v4.2.1 或兼容 MATLAB 语法的数值计算环境
- MATLAB R2024b / App Designer 环境可用于安装和运行 `.mlappinstall` 应用包

说明书中的开发与测试环境包括：

- CPU: Intel i5-12450H
- 内存: 16 GB
- MATLAB R2023b/R2024b
- MATLAB App Designer

## 仓库结构

```text
PDE-Toolbox/
├── Code/
│   ├── 椭圆型偏微分方程/
│   ├── 抛物型偏微分方程/
│   ├── 双曲型偏微分方程/
│   ├── 多尺度椭圆型_参考解/
│   ├── 多尺度椭圆型_高阶多尺度/
│   ├── 多尺度抛物型_参考解/
│   ├── 多尺度抛物型_高阶多尺度/
│   ├── 多尺度双曲型_参考解/
│   └── 多尺度双曲型_高阶多尺度/
├── Help/
│   ├── 工具箱说明文档.pdf
│   ├── 经典偏微分方程/
│   └── 多尺度偏微分方程/
├── 软件安装包/
│   └── Launcher.mlappinstall
└── README.md
```

每个求解模块目录中通常包含：

- `solve_*.m`: 主函数入口
- `solve_*_*.m`: 具体格式或具体几何形状的子函数
- `Example1.m`: 最小可运行示例

## 快速开始

```bash
git clone https://github.com/PokerJin/PDE-Toolbox.git
cd PDE-Toolbox
```

### 方式一：安装图形界面工具箱

1. 下载或克隆本仓库。
2. 在 MATLAB/App Designer 环境中打开 `Release/PDE-Toolbox-v2.0.mlappinstall`。
3. 按提示完成安装。
4. 启动工具箱后，在模型向导中选择 PDE 类型，配置方程参数、网格、边界条件和初值条件。
5. 点击“求解”生成二维/三维仿真图像。

### 方式二：直接运行函数代码

将 `Code/` 加入运行路径，然后进入相应模块目录运行示例：

```matlab
addpath(genpath('Code'));

cd('Code/椭圆型偏微分方程');
Example1
```

也可以直接调用主函数，例如求解二维椭圆型方程：

```matlab
clc; clear;

xy_range = [0, 1, 0, 1];
Nxy = [41, 41];

equ_para.c = @(x, y) 1;
equ_para.a = @(x, y) 0;
equ_para.f = @(x, y) 1;

bc_type = [1, 1, 1, 1];  % 1: Dirichlet, 2: Neumann
bc_para.h_left = @(x, y) 1;
bc_para.r_left = @(x, y) 0;
bc_para.h_right = @(x, y) 1;
bc_para.r_right = @(x, y) 0;
bc_para.h_bottom = @(x, y) 1;
bc_para.r_bottom = @(x, y) 0;
bc_para.h_top = @(x, y) 1;
bc_para.r_top = @(x, y) 0;

method = 'fivepoint';
solve_elliptic(xy_range, Nxy, equ_para, bc_type, bc_para, method);
```

## 多尺度求解示例

下面示例调用高阶多尺度方法求解二维多尺度椭圆型方程：

```matlab
clc; clear;

xy_range = [0, 1, 0, 1];
Nxy = [100, 50, 200];  % 微观尺度网格、宏观尺度网格、绘图网格

equ_para.c_in = 1;
equ_para.c_ma = 100;
equ_para.epsilon = 0.1;
equ_para.f = @(x, y) 10;
equ_para.len = 0.03;

bc_type = [1, 1, 1, 1];
bc_para.h_left = @(x, y) 1;
bc_para.r_left = @(x, y) 0;
bc_para.h_right = @(x, y) 1;
bc_para.r_right = @(x, y) 0;
bc_para.h_bottom = @(x, y) 1;
bc_para.r_bottom = @(x, y) 0;
bc_para.h_top = @(x, y) 1;
bc_para.r_top = @(x, y) 0;

shape = 'circle';  % circle 或 square
solve_mulell(xy_range, Nxy, equ_para, bc_type, bc_para, shape);
```

## 函数入口

| 方程类型 | 主函数 | 常用方法/参数 |
| --- | --- | --- |
| 经典椭圆型 PDE | `solve_elliptic` | `fivepoint`, `ninepoint` |
| 经典抛物型 PDE | `solve_parabolic` | `forword`, `backword`, `CN`, `adi` |
| 经典双曲型 PDE | `solve_hyperbolic` | `leapfrog` |
| 多尺度椭圆型参考解 | `solve_mulell_dir` | `fivepoint`, `ninepoint` |
| 多尺度抛物型参考解 | `solve_mulpar_dir` | CN 格式 |
| 多尺度双曲型参考解 | `solve_mulhyp_dir` | 蛙跳格式 |
| 多尺度椭圆型高阶解 | `solve_mulell` | `circle`, `square` |
| 多尺度抛物型高阶解 | `solve_mulpar` | `circle`, `square` |
| 多尺度双曲型高阶解 | `solve_mulhyp` | `circle`, `square` |

通用参数约定：

- `xy_range`: 空间求解区域，格式为 `[x_min, x_max, y_min, y_max]`
- `xyt_range`: 时空求解区域，格式为 `[x_min, x_max, y_min, y_max, t_min, t_max]`
- `Nxy` / `Nxyt`: 空间或时空网格点数
- `equ_para`: 方程系数、材料参数、源项和尺度比
- `bc_type`: 四条边界的边界条件类型，`1` 表示 Dirichlet，`2` 表示 Neumann
- `bc_para`: 边界条件参数
- `u0`: 抛物型方程初值
- `init.u0` / `init.v0`: 双曲型方程位移初值与速度初值
- `shape`: 多尺度夹杂相几何形状，支持 `circle` 和 `square`

## 文档

更多说明见仓库中的帮助文档：

- `Help/工具箱说明文档.pdf`
- `Help/经典偏微分方程/函数帮助文档.pdf`
- `Help/经典偏微分方程/工具箱使用示例.pdf`
- `Help/多尺度偏微分方程/函数帮助文档.pdf`
- `Help/多尺度偏微分方程/工具箱使用示例.pdf`

## 适用场景

- 偏微分方程数值方法教学与演示
- 经典 PDE 有限差分格式验证
- 复合材料与多孔介质中的热传导、波传播数值实验
- 多尺度 PDE 的直接精细求解与高阶均匀化近似对比
- 工程仿真平台原型开发与课程设计

## 许可说明

本仓库暂未附带正式开源许可证。若需将代码用于课程、科研或工程项目，请先确认项目组授权方式。
