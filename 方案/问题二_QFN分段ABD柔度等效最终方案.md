# 问题二：QFN 三层分段占域—标量 ABD 柔度等效方案

## 1. 问题概述与问题一承接

问题二要求估计 QFN 封装沿中心至角点对角线方向的等效拉伸杨氏模量、等效弯曲杨氏模量和等效热膨胀系数。题面中的 QFN 本体为三层结构：第一层为纯环氧树脂，第二层为环氧树脂与芯片的分区组合，第三层为环氧树脂、中央铜焊盘及局部焊料几何的分区组合。焊料仅在材料参数上按铜近似，不能因此删去或合并其几何。

本方案直接继承问题一的参数定义和力学口径：

| 问题一内容 | 问题二中的继承与扩展 |
|---|---|
| 沿封装对角线取单位宽度窄条 | 仍沿 QFN 中心至角点的半对角线建模 |
| 用截面轴向刚度 \(S_0\) 定义 \(E_{d,t}=S_0/h\) | 将单一截面推广为分段局部 \(A_r\)，经路径柔度组装得到 \(A_{\mathrm{eq}}\) |
| 用弹性中性轴和弯曲刚度 \(J\) 定义 \(E_{d,b}=12J/h^3\) | 由等效 \(A_{\mathrm{eq}},B_{\mathrm{eq}},D_{\mathrm{eq}}\) 求中性轴及 \(D_n\) |
| 用刚度加权热响应定义 \(\alpha_d=T_0/S_0\) | 将局部热合力、热矩与 ABD 同步组装，得到 \(\alpha_{\mathrm D}^{\mathrm{QFN}}\) |
| 厚度方向各材料共同承载 | 每个局部截面内仍按刚度积分并联承载 |

问题二相对问题一的结构性扩展是：芯片、铜盘和焊料只占据部分对角线路径，因此局部截面随 \(s\) 改变。模型必须先计算各段截面 ABD，再沿 \(s\) 按柔度串联，不能把夹杂刚度铺满整条路径。

最终主结果为

$$
\boxed{E_{\mathrm D,t}^{\mathrm{QFN}}},\qquad
\boxed{E_{\mathrm D,b}^{\mathrm{QFN}}},\qquad
\boxed{\alpha_{\mathrm D}^{\mathrm{QFN}}}.
$$

其中拉伸模量和弯曲模量必须分别报告，不能用一个“等效杨氏模量”替代二者。

---

## 2. 题面证据链与建模边界

### 2.1 几何证据链

本模型严格采用以下题面链条：

$$
\text{QFN 分三层}
\rightarrow
\text{第一层纯环氧}
\rightarrow
\text{第二层为环氧与芯片分区}
\rightarrow
\text{第三层为环氧、铜焊盘和局部焊料分区}
\rightarrow
\text{求 QFN 本体等效参数}.
$$

因此，主模型的等效域是厚度为 \(H_5\) 的三层 QFN 本体。PCB 不构成第四层，焊料也不构成覆盖全长的独立层。

### 2.2 基本假设

1. 各材料在研究温区内服从小变形线弹性热弹性本构，参数取常数，界面理想粘结。
2. 沿对角线截取单位宽度窄条；因题目未给泊松比，采用一个方向上的标量 ABD，而不虚构完整二维本构矩阵。
3. 温度变化准静态且截面内均匀；本模型只估计整体等效参数，不预测局部峰值应力、蠕变、脱层或疲劳寿命。
4. 芯片、铜盘和焊料在其几何域内替换环氧树脂，不与环氧在同一空间重复计量。
5. 焊料几何由 \(LEN_1,H_1\) 独立描述，但按题意取

   $$
   E_{\mathrm{sol}}=E_{\mathrm c},\qquad
   \alpha_{\mathrm{sol}}=\alpha_{\mathrm c}.
   $$

### 2.3 \(H_2\) 至 \(H_5\) 的条件解释

不能把 \(H_5=H_2+H_3+H_4\) 当作题面已明示的硬约束。应先按题图复现各界面位置，再选用相应分支。

取 QFN 本体底面为 \(z_0=0\)。若题图确认 \(H_3,H_4\) 分别就是第二、第三功能层的净厚度，可取

$$
t_3=H_4,\qquad t_2=H_3,\qquad
t_1=H_5-H_3-H_4>0,
\tag{1}
$$

并以 \(H_2=t_1\) 作为“\(H_2\) 是顶部纯环氧净厚度”这一解释的核对条件。只有核对成立时，才可进一步写成 \(H_5=H_2+H_3+H_4\)。

若题图显示 \(H_2\) 的量测范围不同，或 \(H_1\) 改变了第三层界面位置，则直接采用题图读取的真实界面坐标 \(z_0,z_1,z_2,z_3\)，只保留

$$
z_3-z_0=H_5,\qquad t_j=z_j-z_{j-1}>0.
\tag{2}
$$

焊料高度 \(H_1\) 只用于确定第三层内局部焊料子域的竖向边界，不能未经题图确认就加入总厚度，也不能用它替代铜盘厚度 \(H_4\)。

---

## 3. 修正后的符号表

为避免与问题一的材料编号冲突，问题二统一采用材料语义下标，题给编号只在“来源”栏中保留。尤其应注意：问题一的 \(E_1,E_2,E_3\) 分别表示 PCB、焊球和 BGA；问题二题图中的同号参数对应环氧、芯片和铜，不能跨问题直接沿用编号含义。

### 3.1 题给几何与材料参数

| 符号 | 含义 | 单位 | 来源或与问题一的关系 |
|---|---|---:|---|
| \(s\) | 从 QFN 中心指向角点的对角线路径坐标 | m | 继承问题一的对角线坐标思想 |
| \(L_{\mathrm D}=LD_2\) | QFN 中心至角点的半对角线长度 | m | 题图 \(LD_2\)，且 \(0\le s\le L_{\mathrm D}\) |
| \(LD_1\) | 芯片区域沿该路径的边界位置 | m | 题图芯片尺寸 |
| \(LD_3\) | 中央铜焊盘沿该路径的边界位置 | m | 题图铜盘尺寸 |
| \(LEN_1\) | 焊料几何的题给长度 | m | 题图；不预设其就是对角线投影长度 |
| \(H_1\) | 局部焊料几何高度 | m | 题图；保留独立几何身份 |
| \(H_2\) | 题图标注的环氧相关厚度 | m | 仅在确认量测范围后参与厚度闭合 |
| \(H_3\) | 芯片厚度 | m | 题图 |
| \(H_4\) | 中央铜焊盘厚度 | m | 题图；不可与 \(H_1\) 混同 |
| \(H_5\) | QFN 本体总厚度 | m | 题图；三个主结果的统一参考厚度 |
| \(E_{\mathrm m},\alpha_{\mathrm m}\) | 环氧树脂的杨氏模量、CTE | Pa，K\(^{-1}\) | 问题二题给 \(E_1,CTE_1\) |
| \(E_{\mathrm{die}},\alpha_{\mathrm{die}}\) | 芯片的杨氏模量、CTE | Pa，K\(^{-1}\) | 问题二题给 \(E_2,CTE_2\) |
| \(E_{\mathrm c},\alpha_{\mathrm c}\) | 铜的杨氏模量、CTE | Pa，K\(^{-1}\) | 问题二题给 \(E_3,CTE_3\) |
| \(E_{\mathrm{sol}},\alpha_{\mathrm{sol}}\) | 焊料采用的杨氏模量、CTE | Pa，K\(^{-1}\) | 题意简化为铜参数，不改变焊料几何 |
| \(E_{\mathrm{PCB}},\alpha_{\mathrm{PCB}}\) | PCB 的杨氏模量、CTE | Pa，K\(^{-1}\) | 问题二题给 \(E_4,CTE_4\)；不进入 QFN 本体 ABD |

### 3.2 几何占域与分段符号

| 符号 | 含义 | 单位 | 来源或作用 |
|---|---|---:|---|
| \(z\) | QFN 厚度坐标 | m | 继承问题一的截面坐标 |
| \(z_*\) | 所有区段共用的参考线，取 \(H_5/2\) | m | 保证各段广义变形可直接组装 |
| \(z_0,z_1,z_2,z_3\) | 三层的四个实际界面坐标 | m | 由题图复现，且 \(z_3-z_0=H_5\) |
| \(t_1,t_2,t_3\) | 第一、第二、第三功能层净厚度 | m | \(t_j=z_j-z_{j-1}\)；条件分支见式（1）—（2） |
| \(\Omega_{\mathrm Q}\) | 三层 QFN 本体等效域 | — | 主模型外边界 |
| \(\Omega_{\mathrm{die}},\Omega_{\mathrm c},\Omega_{\mathrm S}\) | 芯片、铜盘、焊料的真实子域 | — | 由题图及 \(LD_1,LD_3,LEN_1,H_1\) 定位 |
| \(s_{\mathrm S}^-,s_{\mathrm S}^+\) | 焊料域与对角线的起、止交点 | m | 由真实几何投影确定 |
| \(L_{\mathrm S}\) | 焊料沿对角线的实际投影长度 | m | \(L_{\mathrm S}=s_{\mathrm S}^+-s_{\mathrm S}^-\)；仅在题图确认时令其等于 \(LEN_1\) |
| \(\chi_j(s,z)\) | 材料 \(j\) 的指示函数 | 0 或 1 | 保留真实占域，不另造经验覆盖率 |
| \(\mathcal B\) | 对角线路径完整断点集合 | m | 由所有材料边界组成 |
| \(b_r,I_r,\ell_r\) | 排序后断点、第 \(r\) 段及其长度 | m | \(I_r=(b_{r-1},b_r)\)，\(\ell_r=b_r-b_{r-1}\) |
| \(r,R\) | 路径区段编号、区段总数 | — | \(r=1,2,\ldots,R\) |

### 3.3 局部 ABD、等效量与问题一接口

| 符号 | 含义 | 单位 | 与问题一的关系 |
|---|---|---:|---|
| \(A_r\) | 第 \(r\) 段单位宽度轴向刚度 | N/m | 对应问题一单截面的 \(S_0\) |
| \(B_r\) | 第 \(r\) 段拉伸—弯曲耦合刚度 | N | 对应参考线下一阶刚度矩 |
| \(D_r\) | 第 \(r\) 段关于 \(z_*\) 的弯曲刚度 | N·m | 移至中性轴后对应问题一的 \(J\) |
| \(n_r^T,m_r^T\) | 第 \(r\) 段单位温升热合力、热矩系数 | N/(m·K)，N/K | \(n_r^T\) 是问题一 \(T_0\) 的局部推广 |
| \(\mathbf K_r,\mathbf C_r\) | 第 \(r\) 段刚度矩阵及柔度矩阵 | 混合量纲 | \(\mathbf C_r=\mathbf K_r^{-1}\) |
| \(\varepsilon_{0r},\kappa_r\) | 第 \(r\) 段参考线应变、曲率 | 1，m\(^{-1}\) | 局部广义变形 |
| \(N,M\) | 单位宽度轴力、关于 \(z_*\) 的弯矩 | N/m，N | 沿无分布载荷的各段连续传递 |
| \(\mathbf p\) | 广义内力向量 \([N,M]^{\mathsf T}\) | 混合量纲 | 用于各段串联组装 |
| \(\Delta T\) | 相对参考温度的温度变化 | K | 继承问题一热弹性工况 |
| \(\bar{\mathbf C},\bar{\mathbf g}\) | 路径平均柔度、单位温升平均自由变形向量 | 混合量纲 | 将各段串联为整体 |
| \(\bar{\mathbf e}\) | 路径平均广义变形 \([\bar\varepsilon_0,\bar\kappa]^{\mathsf T}\) | 1，m\(^{-1}\) | 整体本构状态量 |
| \(A_{\mathrm{eq}},B_{\mathrm{eq}},D_{\mathrm{eq}}\) | 等效 ABD 的三个标量分量 | N/m，N，N·m | 问题一单截面刚度的路径推广 |
| \(N_{\mathrm{eq}}^T,M_{\mathrm{eq}}^T\) | 等效单位温升热合力、热矩系数 | N/(m·K)，N/K | 问题一热刚度量的路径推广 |
| \(k_{11},k_{12},k_{22}\) | \(\mathbf K_{\mathrm{eq}}\) 的三个矩阵元 | N/m，N，N·m | 分别等于 \(A_{\mathrm{eq}},B_{\mathrm{eq}},D_{\mathrm{eq}}\)，仅用于集中呈现最终结果 |
| \(z_n\) | 等效弹性中性轴绝对坐标 | m | 继承问题一 \(z_n=S_1/S_0\) 的物理定义 |
| \(D_n\) | 关于等效弹性中性轴的弯曲刚度 | N·m | 对应问题一的 \(J\) |
| \(\bar\varepsilon_n,M_n\) | 等效中性轴处平均应变、关于该轴的弯矩 | 1，N | 用于三个单位工况 |
| \(E_{\mathrm D,t}^{\mathrm{QFN}}\) | QFN 对角线等效拉伸杨氏模量 | Pa | 对应问题一 \(E_{d,t}\) |
| \(E_{\mathrm D,b}^{\mathrm{QFN}}\) | QFN 对角线等效弯曲杨氏模量 | Pa | 对应问题一 \(E_{d,b}\) |
| \(\alpha_{\mathrm D}^{\mathrm{QFN}}\) | QFN 中性轴处自由等效 CTE | K\(^{-1}\) | 对应问题一 \(\alpha_d\) |
| \(\kappa_T\) | 单位温升产生的自由热曲率 | m\(^{-1}\)K\(^{-1}\) | 辅助量，不替代三项主结果 |

本文用 \(r\) 表示路径区段，避免继续使用问题一的 \(q_i\)；问题一中 \(q_i\) 已表示折算模量，若再用作区段编号会造成符号重载。

---

## 4. 三层几何、指示函数与完整分段

### 4.1 三层主模型

令 \(\Omega_1,\Omega_2,\Omega_3\) 分别表示第一、第二、第三功能层：

$$
\Omega_1:\ \text{仅环氧},
$$

$$
\Omega_2:\ \Omega_{\mathrm{die}}\ \text{内为芯片，其余为环氧},
$$

$$
\Omega_3:\ \Omega_{\mathrm c}\ \text{内为铜盘，}\Omega_{\mathrm S}\ \text{内为焊料，其余为环氧}.
$$

焊料域 \(\Omega_{\mathrm S}\) 必须由 \(LEN_1,H_1\) 与题图共同定位，铜盘域仍由 \(LD_3,H_4\) 定位。即使二者采用相同材料参数，也不能删除其中任一几何域。

### 4.2 材料指示函数

定义

$$
\chi_{\mathrm{die}}=\mathbf 1_{\Omega_{\mathrm{die}}},\qquad
\chi_{\mathrm c}=\mathbf 1_{\Omega_{\mathrm c}},\qquad
\chi_{\mathrm S}=\mathbf 1_{\Omega_{\mathrm S}}.
\tag{3}
$$

为防止铜盘与焊料边界接触或投影重合时被重复计量，定义金属占域的并集指示函数

$$
\chi_{\mathrm{met}}=\max(\chi_{\mathrm c},\chi_{\mathrm S}),
\tag{4}
$$

环氧指示函数为

$$
\chi_{\mathrm m}=\chi_{\mathrm Q}-\chi_{\mathrm{die}}-\chi_{\mathrm{met}}.
\tag{5}
$$

于是局部材料场写成

$$
\boxed{
E(s,z)=E_{\mathrm m}\chi_{\mathrm m}
+E_{\mathrm{die}}\chi_{\mathrm{die}}
+E_{\mathrm c}\chi_{\mathrm{met}}},
\tag{6}
$$

$$
\boxed{
(E\alpha)(s,z)=E_{\mathrm m}\alpha_{\mathrm m}\chi_{\mathrm m}
+E_{\mathrm{die}}\alpha_{\mathrm{die}}\chi_{\mathrm{die}}
+E_{\mathrm c}\alpha_{\mathrm c}\chi_{\mathrm{met}}}.
\tag{7}
$$

式（6）—（7）保留“材料指示函数”的核心思想；不再引入没有宽度数据、也未进入最终输出的有限宽窄带占域率 \(f_j\)。

### 4.3 完整分段点集合

沿中心到角点取 \(0\le s\le L_{\mathrm D}=LD_2\)。必须保留完整断点集合

$$
\boxed{
\mathcal B=
\{0,LD_1,LD_3,s_{\mathrm S}^-,s_{\mathrm S}^+,LD_2\}.}
\tag{8}
$$

将集合中的点排序、去重并剔除不在 \([0,LD_2]\) 内的无效值，记为

$$
0=b_0<b_1<\cdots<b_R=LD_2,\qquad
I_r=(b_{r-1},b_r),\qquad
\ell_r=b_r-b_{r-1}.
\tag{9}
$$

该写法不预设 \(LD_1\) 与 \(LD_3\) 的大小关系。每个 \(I_r\) 内材料组合不变，因此每段只需建立一个局部截面模型。

---

## 5. 各段局部一维标量 ABD 模型

所有区段统一取

$$
z_*=\frac{H_5}{2}.
\tag{10}
$$

第 \(r\) 段的对角线应变为

$$
\varepsilon_r(z)=\varepsilon_{0r}+(z-z_*)\kappa_r.
\tag{11}
$$

利用式（6）—（7）在该段沿厚度积分，定义

$$
\boxed{
A_r=\int_{z_0}^{z_3}E_r(z)\,dz,\quad
B_r=\int_{z_0}^{z_3}E_r(z)(z-z_*)\,dz,\quad
D_r=\int_{z_0}^{z_3}E_r(z)(z-z_*)^2\,dz,}
\tag{12}
$$

$$
\boxed{
n_r^T=\int_{z_0}^{z_3}E_r(z)\alpha_r(z)\,dz,\quad
m_r^T=\int_{z_0}^{z_3}E_r(z)\alpha_r(z)(z-z_*)\,dz.}
\tag{13}
$$

由于每段内部均为分片常材料，式（12）—（13）可直接按各材料子域的上下界求和，无须再保留未使用的二维或有限宽积分公式。

记

$$
\mathbf K_r=
\begin{bmatrix}A_r&B_r\\B_r&D_r\end{bmatrix},
\qquad
\mathbf t_r=
\begin{bmatrix}n_r^T\\m_r^T\end{bmatrix}.
\tag{14}
$$

局部热弹性本构为

$$
\boxed{
\begin{bmatrix}N\\M\end{bmatrix}
=\mathbf K_r
\begin{bmatrix}\varepsilon_{0r}\\\kappa_r\end{bmatrix}
-\mathbf t_r\Delta T.}
\tag{15}
$$

这里的 \(2\times2\) ABD 是问题一变换截面模型的直接推广：\(A_r\) 描述拉伸，\(D_r\) 描述弯曲，\(B_r\) 保留非对称截面产生的拉伸—弯曲耦合。

---

## 6. 沿对角线的柔度串联组装

基础模型假设对角线窄条上无分布轴力和分布弯矩，因此广义内力

$$
\mathbf p=\begin{bmatrix}N\\M\end{bmatrix}
$$

在各段连续，而各段应变和曲率可以不同。由式（15）可得

$$
\begin{bmatrix}\varepsilon_{0r}\\\kappa_r\end{bmatrix}
=\mathbf K_r^{-1}(\mathbf p+\mathbf t_r\Delta T).
\tag{16}
$$

各段沿长度是串联关系，应平均柔度而不是直接平均刚度：

$$
\boxed{
\bar{\mathbf C}
=\frac{1}{LD_2}\sum_{r=1}^{R}\ell_r\mathbf K_r^{-1},}
\tag{17}
$$

$$
\boxed{
\bar{\mathbf g}
=\frac{1}{LD_2}\sum_{r=1}^{R}
\ell_r\mathbf K_r^{-1}\mathbf t_r.}
\tag{18}
$$

定义整体等效刚度与热载荷向量

$$
\boxed{
\mathbf K_{\mathrm{eq}}=\bar{\mathbf C}^{-1}
=\begin{bmatrix}
A_{\mathrm{eq}}&B_{\mathrm{eq}}\\
B_{\mathrm{eq}}&D_{\mathrm{eq}}
\end{bmatrix},}
\tag{19}
$$

$$
\boxed{
\mathbf t_{\mathrm{eq}}
=\mathbf K_{\mathrm{eq}}\bar{\mathbf g}
=\begin{bmatrix}N_{\mathrm{eq}}^T\\M_{\mathrm{eq}}^T\end{bmatrix}.}
\tag{20}
$$

整体本构为

$$
\boxed{
\mathbf p
=\mathbf K_{\mathrm{eq}}\bar{\mathbf e}
-\mathbf t_{\mathrm{eq}}\Delta T.}
\tag{21}
$$

式（17）—（20）是问题二的核心新增环节：截面内部通过 ABD 处理厚度并联与耦合，截面之间通过柔度平均处理路径串联。

---

## 7. 中性轴、三个单位工况与主结果

### 7.1 移至等效弹性中性轴

等效中性轴及关于中性轴的弯曲刚度为

$$
\boxed{
z_n=z_*+\frac{B_{\mathrm{eq}}}{A_{\mathrm{eq}}},
\qquad
D_n=D_{\mathrm{eq}}-\frac{B_{\mathrm{eq}}^2}{A_{\mathrm{eq}}}.}
\tag{22}
$$

等效热矩移至中性轴后为

$$
\boxed{
M_{n,\mathrm{eq}}^T
=M_{\mathrm{eq}}^T
-\frac{B_{\mathrm{eq}}}{A_{\mathrm{eq}}}N_{\mathrm{eq}}^T.}
\tag{23}
$$

在中性轴变量 \((\bar\varepsilon_n,\bar\kappa)\) 下，整体本构化为

$$
\boxed{
\begin{bmatrix}N\\M_n\end{bmatrix}
=
\begin{bmatrix}A_{\mathrm{eq}}&0\\0&D_n\end{bmatrix}
\begin{bmatrix}\bar\varepsilon_n\\\bar\kappa\end{bmatrix}
-
\begin{bmatrix}N_{\mathrm{eq}}^T\\M_{n,\mathrm{eq}}^T\end{bmatrix}
\Delta T.}
\tag{24}
$$

### 7.2 单位拉伸工况及拉伸模量结果

令

$$
\bar\varepsilon_n=1,\qquad \bar\kappa=0,\qquad \Delta T=0.
$$

由式（24）有 \(N=A_{\mathrm{eq}}\)。将该 QFN 窄条与厚度同为 \(H_5\) 的均质窄条在相同中性轴应变下等效：

$$
N=E_{\mathrm D,t}^{\mathrm{QFN}}H_5\bar\varepsilon_n.
$$

代入 \(\bar\varepsilon_n=1\) 和 \(N=A_{\mathrm{eq}}\)，得到

$$
\boxed{
E_{\mathrm D,t}^{\mathrm{QFN}}
=\frac{A_{\mathrm{eq}}}{H_5}.}
\tag{25}
$$

该工况规定拉力通过等效中性轴，延续问题一“纯拉伸不附加弯矩”的定义。

### 7.3 单位弯曲工况及弯曲模量结果

令

$$
\bar\varepsilon_n=0,\qquad \bar\kappa=1,\qquad \Delta T=0.
$$

由式（24）有 \(M_n=D_n\)。对同厚度均质窄条，关于中性轴的弯矩为

$$
M_n=E_{\mathrm D,b}^{\mathrm{QFN}}
\left(\frac{H_5^3}{12}\right)\bar\kappa.
$$

代入 \(\bar\kappa=1\) 和 \(M_n=D_n\)，得到

$$
\boxed{
E_{\mathrm D,b}^{\mathrm{QFN}}
=\frac{12D_n}{H_5^3}
=\frac{12}{H_5^3}
\left(D_{\mathrm{eq}}-
\frac{B_{\mathrm{eq}}^2}{A_{\mathrm{eq}}}\right).}
\tag{26}
$$

### 7.4 单位温升工况及等效 CTE 结果

令

$$
\Delta T=1,\qquad N=0,\qquad M_n=0.
$$

由式（24）的轴力方程

$$
0=A_{\mathrm{eq}}\bar\varepsilon_n-N_{\mathrm{eq}}^T\Delta T
$$

可得中性轴处自由热应变与温升之比：

$$
\boxed{
\alpha_{\mathrm D}^{\mathrm{QFN}}
=\bar\varepsilon_n
=\frac{N_{\mathrm{eq}}^T}{A_{\mathrm{eq}}},}
\tag{27}
$$

并可同时得到辅助量

$$
\boxed{
\kappa_T
=\bar\kappa
=\frac{M_{n,\mathrm{eq}}^T}{D_n}.}
\tag{28}
$$

式（25）—（27）是论文主表应报告的三个结果；式（28）只用于说明非对称 QFN 的自由热翘曲，不将其升级为题目要求之外的主目标。

### 7.5 可直接交付的三项结果

为使三项参数不再分散在推导中，令

$$
\mathbf K_{\mathrm{eq}}
=\bar{\mathbf C}^{-1}
=\begin{bmatrix}k_{11}&k_{12}\\k_{12}&k_{22}\end{bmatrix},
\qquad
\mathbf t_{\mathrm{eq}}=\mathbf K_{\mathrm{eq}}\bar{\mathbf g}
=\begin{bmatrix}N_{\mathrm{eq}}^T\\M_{\mathrm{eq}}^T\end{bmatrix}.
\tag{29}
$$

其中 \(\bar{\mathbf C}\) 和 \(\bar{\mathbf g}\) 已由每段题图几何、材料参数通过式（12）—（18）唯一确定。于是三项最终可交付的显式符号结果为

$$
\boxed{
E_{\mathrm D,t}^{\mathrm{QFN}}=\frac{k_{11}}{H_5},\qquad
E_{\mathrm D,b}^{\mathrm{QFN}}
=\frac{12}{H_5^3}\left(k_{22}-\frac{k_{12}^2}{k_{11}}\right),\qquad
\alpha_{\mathrm D}^{\mathrm{QFN}}=\frac{N_{\mathrm{eq}}^T}{k_{11}}.}
\tag{30}
$$

| 交付参数 | 可直接报告的最终表达式 | 单位 | 定义工况 | 数值结果栏 |
|---|---|---:|---|---|
| \(E_{\mathrm D,t}^{\mathrm{QFN}}\) | \(k_{11}/H_5\) | Pa | \(\bar\varepsilon_n=1,\bar\kappa=0,\Delta T=0\) | 待将已确认题图参数代入 |
| \(E_{\mathrm D,b}^{\mathrm{QFN}}\) | \(\dfrac{12}{H_5^3}\left(k_{22}-\dfrac{k_{12}^2}{k_{11}}\right)\) | Pa | \(\bar\varepsilon_n=0,\bar\kappa=1,\Delta T=0\) | 待将已确认题图参数代入 |
| \(\alpha_{\mathrm D}^{\mathrm{QFN}}\) | \(N_{\mathrm{eq}}^T/k_{11}\) | K\(^{-1}\) | \(\Delta T=1,N=0,M_n=0\) | 待将已确认题图参数代入 |

数值交付时，应将表中最后一列改为“数值 + 单位”，并在表注注明采用的 \(H_2\) 厚度解释，以及 \(LEN_1\) 到 \([s_{\mathrm S}^-,s_{\mathrm S}^+]\) 的真实几何投影规则。这样读者能同时看到最终答案、物理工况和结果来源。

---

## 8. PCB 参数为何不直接进入 QFN 本体 ABD

QFN 本体 ABD 的积分域是 \(\Omega_{\mathrm Q}\)，其外边界和参考厚度均由 QFN 的 \(H_5\) 决定。PCB 位于该域之外。若把 \(E_{\mathrm{PCB}},\alpha_{\mathrm{PCB}}\) 直接加入式（12）—（13），就会同时改变总厚度、中性轴和弯曲惯性矩，所得结果将成为“QFN—焊料—PCB 装配体”的等效参数，而不再是题目要求的 QFN 封装参数。

因此：

- \(E_{\mathrm{PCB}}\) 不进入问题二主模型的 \(\mathbf K_r\)；只有另建装配状态约束模型且已知 PCB 厚度、宽度和边界条件时才可使用。
- \(\alpha_{\mathrm{PCB}}\) 可在求得 QFN 本体参数后作为热失配参照：

  $$
  \Delta\alpha
  =\alpha_{\mathrm D}^{\mathrm{QFN}}-\alpha_{\mathrm{PCB}}.
  \tag{31}
  $$

式（31）是问题二结果与后续板级连接分析的接口，不属于 QFN 本体 ABD 的组成部分。

---

## 9. 必做的解析退化与自洽性检查

以下检查不需要新增实验数据，且都直接检验主公式。

### 9.1 芯片几何退化

- 当 \(LD_1\to0\) 时，芯片占域消失，第二层相应区域退化为环氧。
- 当 \(LD_1\to LD_2\) 时，芯片沿整条对角线路径铺满第二层。

两种极限都应由断点集合与指示函数自动实现，不得另改主公式。

### 9.2 铜盘几何退化

- 当 \(LD_3\to0\) 时，中央铜盘贡献消失。
- 当 \(LD_3\to LD_2\) 时，铜盘沿整条路径铺满第三层的相应厚度区域。

### 9.3 焊料几何退化

- 当 \(s_{\mathrm S}^+-s_{\mathrm S}^-\to0\) 或 \(H_1\to0\) 时，焊料贡献应连续消失。
- 当焊料域铺满第三层中由题图指定的可占区域时，\(\chi_{\mathrm S}\) 应趋于该区域的指示函数；由于 \(E_{\mathrm{sol}}=E_{\mathrm c}\)、\(\alpha_{\mathrm{sol}}=\alpha_{\mathrm c}\)，该区域应等价为完整铜参数区。

“铺满”只改变 \(\Omega_{\mathrm S}\)，不能把 \(LEN_1/H_1\) 擅自改成覆盖全厚或全长的铜层，也不能抹去中央铜盘 \(\Omega_{\mathrm c}\) 的独立几何身份。

### 9.4 均质与无分段退化

若全部材料具有相同 \(E,\alpha\)，则必须得到

$$
\boxed{
E_{\mathrm D,t}^{\mathrm{QFN}}
=E_{\mathrm D,b}^{\mathrm{QFN}}=E,\qquad
\alpha_{\mathrm D}^{\mathrm{QFN}}=\alpha,\qquad
\kappa_T=0.}
\tag{32}
$$

若所有区段的 \(\mathbf K_r,\mathbf t_r\) 完全相同，则

$$
\bar{\mathbf C}=\mathbf K_r^{-1},\qquad
\mathbf K_{\mathrm{eq}}=\mathbf K_r,\qquad
\mathbf t_{\mathrm{eq}}=\mathbf t_r.
\tag{33}
$$

此时式（25）—（27）分别退化为问题一的

$$
E_{d,t}=\frac{S_0}{h},\qquad
E_{d,b}=\frac{12J}{h^3},\qquad
\alpha_d=\frac{T_0}{S_0},
$$

这构成问题一至问题二最直接的承接验证。

### 9.5 正定性与量纲

每段都应满足

$$
A_r>0,\qquad D_r>0,\qquad A_rD_r-B_r^2>0,
\tag{34}
$$

并检查 \(A\) 为 N/m、\(B\) 为 N、\(D\) 为 N·m，最终两个模量均为 Pa、CTE 为 K\(^{-1}\)。

---

## 10. 闭合的基础建模流程

1. 根据题图确定三层界面、芯片域、铜盘域与焊料域，并按第 2.3 节核对厚度解释。
2. 由式（8）构造完整断点集合，在每个区段确定三层的实际材料组成。
3. 用式（12）—（15）得到每段局部 \(A_r,B_r,D_r,n_r^T,m_r^T\)。
4. 用式（17）—（20）完成路径柔度组装，得到等效 ABD 与等效热载荷。
5. 用式（22）—（24）移至等效中性轴，通过三个单位工况定义并输出式（25）—（27）。
6. 依次执行第 9 节的消失、铺满、均质、无分段、正定性和量纲检查。

至此，题图输入、模型中间量、最终输出和解析验证形成完整闭环。几何与材料数据齐全时可直接代入；若 \(H_2\)、\(LEN_1\) 的量测方向或焊料竖向位置尚未由题图确认，则保留相应几何符号，不虚构数值关系。

---

## 11. 本版明确删除的内容

以下内容不进入问题二基础方案：

1. **灵敏度分析公式**：题目未给参数波动范围，本阶段也不产生灵敏度结果，保留只会形成无数据支撑的空章节。
2. **有限元误差公式与网格流程**：当前没有有限元模型、网格数据或对照结果，不能预先书写貌似完整但无法兑现的误差指标。
3. **有限宽窄带占域率 \(f_j\)**：没有窄带宽度及横向几何数据，且该量不进入最终三项输出，故删除。
4. **Mori–Tanaka、自洽均匀化及经验修正系数**：芯片和铜盘是单个中央大占域，真实边界已由 \(LD_1,LD_3\) 给出，无需再假设随机小夹杂；加入这类模型既会丢失路径边界，也不会增加可验证结果。
5. **未使用的 Voigt 上界、多套装配边界和伪代码**：本方案只保留一套与问题一同口径、可闭合的本体参数定义。装配约束可在后续问题需要时以式（31）为接口另建模型。

---

## 12. 完成度与适用范围

本方案已完成问题二基础数学模型所需的四个闭环：三层几何闭环、局部 ABD 本构闭环、路径柔度组装闭环以及三个单位工况的参数定义闭环；同时给出了不依赖额外数据的解析退化检查。其完成度属于“可形成论文模型建立部分、可在数据齐备后直接计算”的符号模型。

尚未完成且本阶段不应伪装成已完成的内容只有数值代入；三项结果的定义、推导和交付表已由第 7.5 节集中给出。数值栏必须依赖题图的最终尺寸口径和材料数值填写。本方案也不声称能够给出角点峰值应力、焊料疲劳寿命或 PCB 约束反力；这些均超出问题二 QFN 本体等效参数估计的基础范围。

## 参考依据

[1] A. T. Nettles. *Basic Mechanics of Laminated Composite Plates*. NASA Reference Publication 1351, 1994.

[2] J. N. Reddy. *Mechanics of Laminated Composite Plates and Shells: Theory and Analysis*, 2nd ed. CRC Press, 2004.

[3] R. A. Schapery. “Thermal Expansion Coefficients of Composite Materials Based on Energy Principles.” *Journal of Composite Materials*, 2(3), 380–404, 1968.
