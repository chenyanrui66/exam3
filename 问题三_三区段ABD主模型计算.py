# -*- coding: utf-8 -*-
"""
问题三 实际BGA封装角点热弹性参数估计
=====================================
三区段 ABD 柔度等效主模型（角点窗口融合版）完整计算代码。

模型结构（对应《问题三_方案B_角点窗口融合版.md》正文）：
    第 1 步  桶形焊球等柔度面积        —— 式(3-7)~(3-11)
    第 2 步  离散连接器 → 连续层接口   —— 式(3-12)~(3-14)
    第 3 步  437/441 阵列坐标生成      —— 式(3-2)~(3-4)
    第 4 步  角点窗口段局部折算模量    —— 式(3-17)~(3-20)
    第 5 步  三(两)区段 ABD 柔度组装  —— 式(3-21)~(3-28)
    第 6 步  结果输出：主结果、437/441 对比、窗口宽度敏感性

单位制约定：力 N，长度 mm，温度 K。
    弹性模量输入用 GPa，程序内统一换算为 N/mm^2（1 GPa = 1e3 N/mm^2）；
    热膨胀系数输入用 ppm/K，程序内统一换算为 1/K（1 ppm/K = 1e-6 /K）。
    于是 A 刚度单位 N/mm，B 刚度单位 N，D 刚度单位 N·mm，
    输出等效模量时再除 1e3 换回 GPa，等效热膨胀系数乘 1e6 换回 ppm/K。

依赖：numpy、scipy（仅用于焊盘热膨胀系数的反演校准；若直接取
      alpha_pad = 17.0 ppm/K，可删除该校准步骤，scipy 亦非必需）。
运行：python 问题三_三区段ABD主模型计算.py
"""

import numpy as np

# ======================================================================
# 0. 输入参数
# ======================================================================

# ---------- 0.1 题目给定数据（A题问题3） ----------
# PCB 层
E_PCB   = 28.6e3      # PCB 杨氏模量 [N/mm^2]（题给 28.6 GPa）
A_PCB   = 17.0e-6     # PCB 热膨胀系数 [1/K]（题给 17.0 ppm/C）
T_PCB   = 1.6         # PCB 厚度 [mm]

# BGA 基板（laminate）层
E_LAM   = 73.3e3      # 基板杨氏模量 [N/mm^2]（题给 73.3 GPa）
A_LAM   = 21.0e-6     # 基板热膨胀系数 [1/K]（题给 21.0 ppm/C）
T_LAM   = 0.8         # 基板厚度 [mm]

# 覆膜（overmold）层：芯片上封装材料参数与覆膜相同，合并为满铺层
E_OM    = 11.7e3      # 覆膜杨氏模量 [N/mm^2]（题给 11.7 GPa）
A_OM    = 23.0e-6     # 覆膜热膨胀系数 [1/K]（题给 23.0 ppm/C）
T_OM    = 1.17        # 覆膜总厚度 [mm]（0.20 芯片结 + 0.97 芯片上封装）

# 芯片结（die）：位于封装正中心，替换覆膜层底部 0.2 mm
E_DIE   = 130.0e3     # 芯片杨氏模量 [N/mm^2]（题给 130 GPa）
A_DIE   = 2.6e-6      # 芯片热膨胀系数 [1/K]（题给 2.6 ppm/C）
T_DIE   = 0.2         # 芯片结厚度 [mm]
L_DIE   = 4.5         # 芯片结边长 [mm]

# 封装与焊球几何
L_BGA   = 26.0        # BGA 封装边长 [mm]
N_BALL  = 437         # 实际焊球总数（题给硬约束）
H_CONN  = 0.50        # 焊球连接体总高 Total_h [mm]
H_PAD_L = 0.03        # 下焊盘厚度 low_pad_h [mm]（题给 030 um）
H_PAD_U = 0.02        # 上焊盘厚度 high_pad_h [mm]（题给 020 um）
R_PAD   = 0.25        # 焊盘半径 [mm]（pad 直径 500 um）
R_MAX   = 0.35        # 焊球最大半径 [mm]（球径 700 um）

# ---------- 0.2 可替换示例参数（非题给，见正文假设 7） ----------
P_PITCH = 1.0         # 焊球节距 [mm]（图3.2对应的23格外框示意值）
E_SOLD  = 50.0e3      # 焊料杨氏模量 [N/mm^2]（示例 50 GPa）
NU_SOLD = 0.35        # 焊料泊松比（示例）
A_SOLD  = 22.0e-6     # 焊料热膨胀系数 [1/K]（示例 22 ppm/K）
E_PAD   = 110.0e3     # 焊盘（铜）杨氏模量 [N/mm^2]（示例）
NU_PAD  = 0.34        # 焊盘（铜）泊松比（示例）
A_PAD   = 17.0e-6     # 焊盘（铜）热膨胀系数 [1/K]（示例；
                      # 该值同时使两区段模型复现原文 19.28973 ppm/K，
                      # 见第 6 节校准函数）

# 派生几何量
LD      = L_BGA / np.sqrt(2)      # 中心到角点的路径长度 L_D [mm]
S_DIE   = L_DIE / np.sqrt(2)      # 芯片段终点对角线坐标 s_die [mm]
H_BAR   = H_CONN - H_PAD_L - H_PAD_U   # 桶形焊料本体高度 0.45 mm
H_TOT   = T_PCB + H_CONN + T_LAM + T_OM  # 组合截面总厚度 4.07 mm
Z_STAR  = H_TOT / 2               # ABD 参考面（取截面中面）

# 焊料与焊盘剪切模量（各向同性换算 G = E / 2(1+nu)）
G_SOLD  = E_SOLD / (2 * (1 + NU_SOLD))
G_PAD   = E_PAD / (2 * (1 + NU_PAD))


# ======================================================================
# 1. 桶形焊球等柔度面积（式 3-7 ~ 3-11）
# ======================================================================
# 物理图像：沿高度把焊球切成薄圆片，同一剪切力依次穿过所有薄片，
# 各薄片柔度 1/(G·A(z)) 串联求和——所以是"柔度的积分"而非"面积的平均"。
# 母线用满足两端半径 r_pad、中点最大半径 r_max 的抛物线重构。

def barrel_area(n_gauss=96):
    """返回桶体等柔度面积 A_bar [mm^2] 与积分 ∫dz/A(z) [1/mm]。"""
    # 96 点 Gauss-Legendre 求积节点与权重（标准区间 [-1,1]）
    xi, wi = np.polynomial.legendre.leggauss(n_gauss)
    # 映射到 [0, H_BAR]
    z  = H_BAR / 2 * (1 + xi)
    # 抛物线母线：z=0 与 z=H_BAR 处 r=r_pad，z=H_BAR/2 处 r=r_max
    r  = R_PAD + 4 * (R_MAX - R_PAD) * (z / H_BAR) * (1 - z / H_BAR)
    A  = np.pi * r**2
    # ∫ dz/A(z)
    integral = np.sum(wi * (H_BAR / 2) / A)
    # 等柔度面积：把变截面桶体折算成"同高、等柔度"的等截面柱
    A_bar = H_BAR / integral
    return A_bar, integral


def connector_stiffness():
    """
    单连接器（下焊盘 + 桶形焊料 + 上焊盘）串联总柔度与剪切刚度。
    返回:
        C_G   总剪切柔度 [mm/N]
        k_G   剪切刚度 [N/mm]（数值上 1 N/mm = 1e3 N/m，注意换算）
        A_Geq 折算到总高 H_CONN 上的等柔度面积 [mm^2]
    """
    A_bar, integral = barrel_area()
    A_pad = np.pi * R_PAD**2
    # 三段柔度串联：下焊盘 + 焊料桶体 + 上焊盘
    C_G = (H_PAD_L / (G_PAD * A_pad)
           + integral / G_SOLD
           + H_PAD_U / (G_PAD * A_pad))
    k_G   = 1.0 / C_G
    A_Geq = H_CONN / (G_SOLD * C_G)   # 式(3-10)
    return C_G, k_G, A_Geq


# ======================================================================
# 2. 离散连接器 → 连续连接层接口（式 3-12 ~ 3-14）
# ======================================================================

def alpha_conn(a_pad=A_PAD):
    """连接层热膨胀系数：三段材料按厚度做自由串联伸长等效（式 3-14）。"""
    return (A_SOLD * H_BAR + a_pad * (H_PAD_L + H_PAD_U)) / H_CONN


def E_conn_global(n_ball):
    """
    全局均匀化连接层模量（式 3-12）：
    让连续层与 N 个离散连接器具有相同的平均剪切柔度。
        G_c = N·k_G·H / L^2 ，  E_c = 2(1+nu_s)·G_c
    """
    G_c = n_ball * K_G * H_CONN / L_BGA**2
    return 2 * (1 + NU_SOLD) * G_c


# ======================================================================
# 3. 437 / 441 阵列坐标（式 3-2 ~ 3-4）
# ======================================================================
# 23×23 外框 - 13×13 中部空区 + 9×9 中心阵列，再去四个外角球。
# 验证：23^2 - 13^2 + 9^2 - 4 = 529 - 169 + 81 - 4 = 437（精确吻合题给球数）

def ball_set(n_ball):
    """
    返回焊球格点列表 [(i,j), ...]，实际坐标为 (i·p, j·p) [mm]。
    n_ball = 441：保留四个外角球；n_ball = 437：去掉四个外角球。
    """
    S = [(i, j) for i in range(-11, 12) for j in range(-11, 12)
         if max(abs(i), abs(j)) <= 4 or max(abs(i), abs(j)) >= 7]
    if n_ball == 441:
        return S
    corners = [(-11, -11), (-11, 11), (11, -11), (11, 11)]
    return [b for b in S if b not in corners]


# ======================================================================
# 4. 角点窗口段：局部折算模量（式 3-17 ~ 3-20）
# ======================================================================
# 窗口取角三角形 Omega_w = {s >= L_D - w}（斜边 w·sqrt2，面积恰为 w^2），
# 与路径坐标 s 天然对齐，其在路径上的投影就是第三区段 (L_D-w, L_D]。
# 折算口径与全局式(3-12)完全一致，仅把"全球数/全面积"换成"窗口球数/窗口面积"：
#   窗口密度 = 全局密度  ⇒  E_win = E_global（退化为两区段模型，一致性条件）

def window_count(balls, w, p=P_PITCH):
    """统计落在角三角形窗口内的焊球数 N_c（式 3-18）。"""
    return sum(1 for (i, j) in balls if (i + j) * p / np.sqrt(2) >= LD - w)


def E_conn_window(n_c, w):
    """
    窗口连接层局部折算模量 [N/mm^2]（式 3-19）。
    N_c = 0 时窗口段连接层缺席（零刚度、零热载荷），返回 0。
    """
    if n_c == 0:
        return 0.0
    G_c = n_c * K_G * H_CONN / w**2
    return 2 * (1 + NU_SOLD) * G_c


def corner_diagnostics(balls, p=P_PITCH):
    """角点窗口诊断量：最近球距离 d_min [mm]（纯几何先验指标）。"""
    c = L_BGA / 2   # 角点坐标 (13, 13)
    return min(np.hypot(c - i * p, c - j * p) for (i, j) in balls)


# ======================================================================
# 5. 区段 ABD 刚度与路径柔度组装（式 3-21 ~ 3-28）
# ======================================================================

def make_layers(E_c, with_die, a_c):
    """
    构造一个区段的叠层列表 [(z0, z1, E, alpha), ...]，z 自 PCB 底面起算。
    with_die=True ：芯片段——覆膜底部 0.2 mm 替换为芯片材料；
    with_die=False：非芯片段——覆膜 1.17 mm 满铺。
    E_c = 0 表示连接层缺席（437 阵列的角点窗口段）。
    """
    layers = [(0.0, T_PCB, E_PCB, A_PCB)]                     # PCB 层
    if E_c > 0:
        layers.append((T_PCB, T_PCB + H_CONN, E_c, a_c))      # 连接层
    z2 = T_PCB + H_CONN
    layers.append((z2, z2 + T_LAM, E_LAM, A_LAM))             # 基板层
    z3 = z2 + T_LAM
    if with_die:
        layers.append((z3, z3 + T_DIE, E_DIE, A_DIE))         # 芯片结
        layers.append((z3 + T_DIE, z3 + T_OM, E_OM, A_OM))    # 芯片上覆膜
    else:
        layers.append((z3, z3 + T_OM, E_OM, A_OM))            # 完整覆膜
    return layers


def segment_ABD(layers):
    """
    对单区段沿厚度积分，得局部 ABD 刚度矩阵 K 与单位温升热载荷向量 t（式 3-22）。
    参考面取 z_* = h/2。D 中包含各层自身惯性项 t^3/12 与杠杆项 t·(zc-z*)^2。
    """
    A = B = D = n_th = m_th = 0.0
    for z0, z1, E, al in layers:
        t  = z1 - z0
        zc = 0.5 * (z0 + z1)
        A    += E * t
        B    += E * t * (zc - Z_STAR)
        D    += E * (t**3 / 12 + t * (zc - Z_STAR)**2)
        n_th += E * al * t
        m_th += E * al * t * (zc - Z_STAR)
    K = np.array([[A, B], [B, D]])
    t_vec = np.array([n_th, m_th])
    return K, t_vec


def assemble(segments):
    """
    串联路径柔度组装（式 3-25 ~ 3-28）：
    各段传递相同广义内力 (N,M)，故对柔度 K^{-1} 做长度加权平均，再求逆。
    segments: [(ell_r, K_r, t_r), ...]
    返回字典：路径等效刚度、热载荷、中性轴位置与三个最终等效参数。
    """
    L_tot = sum(s[0] for s in segments)
    # 柔度与"单位温升自由应变"的长度加权平均
    C_bar = sum(s[0] * np.linalg.inv(s[1]) for s in segments) / L_tot
    g_bar = sum(s[0] * (np.linalg.inv(s[1]) @ s[2]) for s in segments) / L_tot
    K_eq  = np.linalg.inv(C_bar)
    t_eq  = K_eq @ g_bar

    A_eq, B_eq, D_eq = K_eq[0, 0], K_eq[0, 1], K_eq[1, 1]
    # 参考面移到等效弹性中性轴（式 3-27），消去拉伸-弯曲耦合
    z_n = Z_STAR + B_eq / A_eq
    D_n = D_eq - B_eq**2 / A_eq

    return {
        'A_eq': A_eq, 'B_eq': B_eq, 'D_eq': D_eq, 'D_n': D_n,
        'z_n': z_n, 'N_th': t_eq[0], 'M_th': t_eq[1],
        # 三个最终等效参数（式 3-28），换回 GPa 与 ppm/K
        'E_dt':  A_eq / H_TOT / 1e3,            # 等效拉伸杨氏模量 [GPa]
        'E_db':  12 * D_n / H_TOT**3 / 1e3,     # 等效弯曲杨氏模量 [GPa]
        'alpha': t_eq[0] / A_eq * 1e6,          # 等效热膨胀系数 [ppm/K]
    }


# ======================================================================
# 6. 两种路径划分方案
# ======================================================================

def two_segment(n_ball, a_c):
    """原方案：两区段（芯片段 + 非芯片段），连接层全程全局均匀化。"""
    E_c = E_conn_global(n_ball)
    K1, t1 = segment_ABD(make_layers(E_c, True,  a_c))
    K2, t2 = segment_ABD(make_layers(E_c, False, a_c))
    return assemble([(S_DIE, K1, t1), (LD - S_DIE, K2, t2)])


def three_segment(n_ball, w, a_c, p=P_PITCH):
    """
    本版方案：三区段 = 芯片段 + 中段 + 角点窗口段。
    窗口段连接层改用局部折算模量 E_conn_window。
    返回 (结果字典, 窗口球数 N_c, 窗口模量 [GPa])。
    """
    balls = ball_set(n_ball)
    E_c   = E_conn_global(n_ball)
    N_c   = window_count(balls, w, p)
    E_win = E_conn_window(N_c, w)

    K1, t1 = segment_ABD(make_layers(E_c,   True,  a_c))   # I1 芯片段
    K2, t2 = segment_ABD(make_layers(E_c,   False, a_c))   # I2 中段
    K3, t3 = segment_ABD(make_layers(E_win, False, a_c))   # I3 角点窗口段

    r = assemble([(S_DIE, K1, t1),
                  (LD - w - S_DIE, K2, t2),
                  (w, K3, t3)])
    return r, N_c, E_win / 1e3


def calibrate_pad_alpha(target=19.28973):
    """
    反演焊盘热膨胀系数，使两区段模型复现原文 alpha_d = 19.28973 ppm/K。
    校准结果恰为 17.0 ppm/K——说明原文即取该值，输入闭合。
    若不需校准，直接令 A_PAD = 17.0e-6 即可跳过本函数。
    """
    from scipy.optimize import brentq
    f = lambda ap: two_segment(437, alpha_conn(ap))['alpha'] - target
    return brentq(f, 10e-6, 25e-6)


# ======================================================================
# 7. 主程序
# ======================================================================

if __name__ == '__main__':

    # ---------- 第 1 步：单球力学标定 ----------
    A_bar, _ = barrel_area()
    C_G, K_G, A_Geq = connector_stiffness()
    # 把 k_G 提升为全局变量，供接口函数使用（N/mm）
    globals()['K_G'] = K_G
    A_C = alpha_conn()

    print('=' * 64)
    print('第 1 步  桶形焊球等柔度修正（式 3-8 ~ 3-11）')
    print(f'  桶体等柔度面积 A_bar   = {A_bar:.6f} mm^2   [基准复现值 0.305935]')
    print(f'  连接器等柔度面积 A_Geq = {A_Geq:.6f} mm^2   [基准复现值 0.315299]')
    print(f'  单球剪切刚度 k_G       = {K_G * 1e-3:.5f} MN/m  [基准复现值 11.67775]')
    print(f'  连接层热膨胀系数 a_c   = {A_C * 1e6:.2f} ppm/K')

    # ---------- 第 2 步：连接层全局折算 ----------
    print('=' * 64)
    print('第 2 步  离散连接器 → 连续层接口（式 3-12）')
    print(f'  E_c(437) = {E_conn_global(437) / 1e3:.5f} GPa   [基准复现值 10.19126]')
    print(f'  E_c(441) = {E_conn_global(441) / 1e3:.5f} GPa   [基准复现值 10.28454]')

    # ---------- 第 3 步：阵列自检 ----------
    s437, s441 = ball_set(437), ball_set(441)
    print('=' * 64)
    print('第 3 步  阵列拓扑自检（式 3-4：23^2-13^2+9^2-4=437）')
    print(f'  |S_437| = {len(s437)}，|S_441| = {len(s441)}')
    print(f'  角点诊断 d_min：437 = {corner_diagnostics(s437):.4f} mm，'
          f'441 = {corner_diagnostics(s441):.4f} mm')

    # ---------- 第 4 步：原两区段模型复现 ----------
    r2_437, r2_441 = two_segment(437, A_C), two_segment(441, A_C)
    print('=' * 64)
    print('第 4 步  原两区段模型（对照组，复现原文结果）')
    for tag, r in [('437', r2_437), ('441', r2_441)]:
        print(f'  {tag}: E_dt = {r["E_dt"]:.5f} GPa, '
              f'E_db = {r["E_db"]:.5f} GPa, '
              f'alpha = {r["alpha"]:.5f} ppm/K')
    print('  [基准复现值 437: 31.02254 / 22.85975 / 19.28973]')

    # ---------- 第 5 步：三区段主模型（基准 w = 3p） ----------
    W_BASE = 3.0 * P_PITCH
    r3_437, nc437, ew437 = three_segment(437, W_BASE, A_C)
    r3_441, nc441, ew441 = three_segment(441, W_BASE, A_C)
    print('=' * 64)
    print(f'第 5 步  三区段主模型（角点窗口 w = 3p = {W_BASE:.0f} mm）')
    for tag, r, nc, ew in [('437', r3_437, nc437, ew437),
                           ('441', r3_441, nc441, ew441)]:
        print(f'  {tag}: N_c = {nc}, E_c^win = {ew:.3f} GPa | '
              f'E_dt = {r["E_dt"]:.5f} GPa, '
              f'E_db = {r["E_db"]:.5f} GPa, '
              f'alpha = {r["alpha"]:.5f} ppm/K')

    pct = lambda a, b: (a - b) / b * 100
    print('-' * 64)
    print('  437 相对 441（三区段）:  '
          f'E_dt {pct(r3_437["E_dt"], r3_441["E_dt"]):+.4f}%, '
          f'E_db {pct(r3_437["E_db"], r3_441["E_db"]):+.4f}%, '
          f'alpha {pct(r3_437["alpha"], r3_441["alpha"]):+.4f}%')
    print('  437 相对 441（两区段）:  '
          f'E_dt {pct(r2_437["E_dt"], r2_441["E_dt"]):+.5f}%, '
          f'E_db {pct(r2_437["E_db"], r2_441["E_db"]):+.5f}%, '
          f'alpha {pct(r2_437["alpha"], r2_441["alpha"]):+.5f}%')
    print('  三区段相对两区段（437）: '
          f'E_dt {pct(r3_437["E_dt"], r2_437["E_dt"]):+.4f}%, '
          f'E_db {pct(r3_437["E_db"], r2_437["E_db"]):+.4f}%, '
          f'alpha {pct(r3_437["alpha"], r2_437["alpha"]):+.4f}%')

    # ---------- 第 6 步：窗口宽度敏感性 ----------
    print('=' * 64)
    print('第 6 步  窗口宽度敏感性扫描 w = 2p / 3p / 4p')
    print(f'  {"w":>4} {"阵列":>4} {"N_c":>4} {"E_win/GPa":>10} '
          f'{"E_dt/GPa":>10} {"E_db/GPa":>10} {"alpha/ppmK":>11}')
    for w in (2.0 * P_PITCH, 3.0 * P_PITCH, 4.0 * P_PITCH):
        for n in (437, 441):
            r, nc, ew = three_segment(n, w, A_C)
            print(f'  {w:>4.0f} {n:>6} {nc:>4} {ew:>10.3f} '
                  f'{r["E_dt"]:>10.5f} {r["E_db"]:>10.5f} {r["alpha"]:>11.5f}')

    # ---------- 可选：焊盘热膨胀系数校准 ----------
    print('=' * 64)
    try:
        ap_cal = calibrate_pad_alpha()
        print(f'可选校准：复现原文 19.28973 ppm/K 所需 alpha_pad '
              f'= {ap_cal * 1e6:.2f} ppm/K（恰为整数 17，输入闭合）')
    except ImportError:
        print('可选校准：未安装 scipy，跳过（直接取 alpha_pad = 17 ppm/K 即可）')
