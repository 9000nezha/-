/**
 * 枪炮内弹道学数值解法程序
 * 基于《枪炮内弹道学》第2章第2.4节内容实现
 * 
 * 功能：
 * 1. 使用四阶龙格-库塔法求解量纲为1的内弹道方程组
 * 2. 计算并定位特殊点：最大压力点、燃烧分裂点、燃烧结束点、炮口点
 * 3. 生成完整的计算结果表格（表2-8格式）
 * 4. 生成p-t、p-l、v-t、v-l曲线数据
 * 
 * 作者：9000nezha
 * 日期：2026年
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>

using namespace std;

// ============================================================================
// 常量定义
// ============================================================================
const double GOLDEN_RATIO = 0.618033988749895;  // 黄金分割比例
const double EPSILON = 1e-10;                   // 精度判别阈值

// ============================================================================
// 数据结构定义
// ============================================================================

// 输入参数结构体
struct InputParams {
    // 火炮构造及弹丸诸元
    double S;      // 炮膛截面积 (m^2)
    double V0;     // 药室容积 (m^3)
    double lg;     // 炮膛全长 (m)
    double m;      // 弹丸质量 (kg)
    
    // 装药条件
    double f;      // 火药力 (J/kg)
    double omega;  // 装药量 (kg)
    double alpha;  // 余容 (m^3/kg)
    double rho_p;  // 火药密度 (kg/m^3)
    double theta;  // 火药特征量
    double u1;     // 燃速系数 (m/(s·Pa^n))
    double n;      // 燃速指数
    double e1;     // 火药弧厚的一半 (m)
    double chi;    // 火药形状特征量
    double lambda; // 火药形状特征量
    double mu;     // 火药形状特征量
    double chi_s;  // 火药形状特征量
    double lambda_s; // 火药形状特征量
    
    // 起始条件
    double p0;     // 起始压力 (Pa)
    
    // 计算常数
    double phi1;   // 次要功系数
    double lambda2; // 系数
    
    // 计算条件
    double h;      // 步长
};

// 计算常量结构体
struct CalculatedConstants {
    double phi;    // 次要功系数
    double Delta;  // 装填密度
    double l0;     // 药室长度
    double vj;     // 特征速度
    double B;      // 系数
    double lg_bar; // 量纲为1的炮膛全长
    double Zk;     // 火药燃尽时相对燃烧厚度
    double psi_s;  // 分裂点相对燃烧厚度
    double xi_s;   // 分裂点相对燃烧厚度
};

// 状态变量结构体（量纲为1）
struct State {
    double t_bar;  // 量纲为1的时间
    double l_bar;  // 量纲为1的行程
    double v_bar;  // 量纲为1的速度
    double p_bar;  // 量纲为1的压力
    double psi;    // 相对燃烧厚度
    double Z;      // 相对燃烧厚度
};

// 结果点结构体
struct ResultPoint {
    double l_dm;   // 行程 (dm)
    double t_ms;   // 时间 (ms)
    double v;      // 速度 (m/s)
    double p;      // 平均压力 (MPa)
    double p_d;    // 膛底压力 (MPa)
    double p_t;    // 膛口压力 (MPa)
    double psi;    // 相对燃烧厚度
    double Z;      // 相对燃烧厚度
    string label;  // 标签
};

// 全局变量
InputParams input;
CalculatedConstants consts;
vector<ResultPoint> allResults;    // 所有计算点
vector<ResultPoint> specialPoints; // 特殊点

// ============================================================================
// 辅助函数
// ============================================================================

// 计算psi(Z) - 相对燃烧厚度函数（公式2-71）
double calculatePsi(double Z) {
    if (Z < 1.0) {
        // Z < 1: psi = chi*Z*(1 + lambda*Z + mu*Z^2)
        return input.chi * Z * (1.0 + input.lambda * Z + input.mu * Z * Z);
    } else if (Z < consts.Zk) {
        // 1 <= Z < Zk: psi = chi_s * Z/Zk * (1 + lambda_s * Z/Zk)
        double Z_ratio = Z / consts.Zk;
        return input.chi_s * Z_ratio * (1.0 + input.lambda_s * Z_ratio);
    } else {
        // Z >= Zk: psi = 1
        return 1.0;
    }
}

// 计算dpsi/dZ
double calculateDpsiDZ(double Z) {
    if (Z < 1.0) {
        // dpsi/dZ = chi * (1 + 2*lambda*Z + 3*mu*Z^2)
        return input.chi * (1.0 + 2.0 * input.lambda * Z + 3.0 * input.mu * Z * Z);
    } else if (Z < consts.Zk) {
        // dpsi/dZ = chi_s/Zk * (1 + 2*lambda_s*Z/Zk)
        return input.chi_s / consts.Zk * (1.0 + 2.0 * input.lambda_s * Z / consts.Zk);
    } else {
        return 0.0;
    }
}

// 计算l_psi_bar - 自由容积修正系数（公式2-72下方）
double calculateLpsiBar(double psi) {
    double term1 = 1.0 - 1.0 / consts.Delta / input.rho_p;
    double term2 = -(input.alpha - 1.0 / input.rho_p) * psi;
    return term1 + term2;
}

// 量纲转换函数
double toL(double l_bar) { return l_bar * consts.l0; }           // m
double toT(double t_bar) { return t_bar * consts.l0 / consts.vj; } // s
double toV(double v_bar) { return v_bar * consts.vj; }           // m/s
double toP(double p_bar) { return p_bar * input.f * consts.Delta; } // Pa

// 计算膛底压力和膛口压力
double calculatePd(double p, double v, double l) {
    if (l < EPSILON) return p;
    double term = consts.phi * input.m * v * v / (2.0 * input.S * p * l * 1e6); // p in MPa
    return p * (1.0 + term);
}

double calculatePt(double p, double v, double l) {
    if (l < EPSILON) return p;
    double term = consts.phi * input.m * v * v / (2.0 * input.S * p * l * 1e6);
    return p * (1.0 - term);
}

// ============================================================================
// 右端函数计算（公式2-72）
// ============================================================================

// 计算dZ/dt_bar
double dZdt(double Z, double p_bar) {
    if (Z >= consts.Zk) return 0.0;
    double term = sqrt(input.theta / (2.0 * consts.B));
    return term * pow(p_bar, input.n);
}

// 计算dl_bar/dt_bar
double dldt(double v_bar) {
    return v_bar;
}

// 计算dv_bar/dt_bar
double dvdt(double p_bar) {
    return input.theta / 2.0 * p_bar;
}

// 计算dp_bar/dt_bar（公式2-72）
double dpdt(double Z, double p_bar, double v_bar, double l_bar, double psi) {
    double l_psi_bar = calculateLpsiBar(psi);
    double denominator = (l_bar + l_psi_bar) * v_bar;
    
    if (fabs(denominator) < EPSILON) return 0.0;
    
    double dZ = dZdt(Z, p_bar);
    double dpsi_dt = calculateDpsiDZ(Z) * dZ;
    
    // dp̄/dt̄ = [l₀/(vj·(l̄+l̄_ψ)·v̄)] * [1+Δ(α-1/ρ_p)·p̄]·(dψ/dt̄) - [(1+θ)/(l̄+l̄_ψ)]·p̄·v̄
    double term1 = consts.l0 / (consts.vj * denominator);
    double term2 = 1.0 + consts.Delta * (input.alpha - 1.0 / input.rho_p) * p_bar;
    double first_part = term1 * term2 * dpsi_dt;
    
    double second_part = (1.0 + input.theta) / (l_bar + l_psi_bar) * p_bar * v_bar;
    
    return first_part - second_part;
}

// ============================================================================
// 四阶龙格-库塔法（公式2-73）
// ============================================================================

void rungeKutta4(State& state, double h) {
    double Z = state.Z;
    double l_bar = state.l_bar;
    double v_bar = state.v_bar;
    double p_bar = state.p_bar;
    double psi = calculatePsi(Z);
    
    // K1
    double k1_Z = dZdt(Z, p_bar);
    double k1_l = dldt(v_bar);
    double k1_v = dvdt(p_bar);
    double k1_p = dpdt(Z, p_bar, v_bar, l_bar, psi);
    
    // K2
    double Z2 = Z + h * k1_Z / 2.0;
    double l2 = l_bar + h * k1_l / 2.0;
    double v2 = v_bar + h * k1_v / 2.0;
    double p2 = p_bar + h * k1_p / 2.0;
    double psi2 = calculatePsi(Z2);
    
    double k2_Z = dZdt(Z2, p2);
    double k2_l = dldt(v2);
    double k2_v = dvdt(p2);
    double k2_p = dpdt(Z2, p2, v2, l2, psi2);
    
    // K3
    double Z3 = Z + h * k2_Z / 2.0;
    double l3 = l_bar + h * k2_l / 2.0;
    double v3 = v_bar + h * k2_v / 2.0;
    double p3 = p_bar + h * k2_p / 2.0;
    double psi3 = calculatePsi(Z3);
    
    double k3_Z = dZdt(Z3, p3);
    double k3_l = dldt(v3);
    double k3_v = dvdt(p3);
    double k3_p = dpdt(Z3, p3, v3, l3, psi3);
    
    // K4
    double Z4 = Z + h * k3_Z;
    double l4 = l_bar + h * k3_l;
    double v4 = v_bar + h * k3_v;
    double p4 = p_bar + h * k3_p;
    double psi4 = calculatePsi(Z4);
    
    double k4_Z = dZdt(Z4, p4);
    double k4_l = dldt(v4);
    double k4_v = dvdt(p4);
    double k4_p = dpdt(Z4, p4, v4, l4, psi4);
    
    // 更新状态
    state.Z += h / 6.0 * (k1_Z + 2.0 * k2_Z + 2.0 * k3_Z + k4_Z);
    state.l_bar += h / 6.0 * (k1_l + 2.0 * k2_l + 2.0 * k3_l + k4_l);
    state.v_bar += h / 6.0 * (k1_v + 2.0 * k2_v + 2.0 * k3_v + k4_v);
    state.p_bar += h / 6.0 * (k1_p + 2.0 * k2_p + 2.0 * k3_p + k4_p);
    state.t_bar += h;
    state.psi = calculatePsi(state.Z);
}

// ============================================================================
// 从给定状态积分指定的时间步长
// ============================================================================

State integrateForTime(State initial, double delta_t, double h) {
    State state = initial;
    double target_t = state.t_bar + delta_t;
    int max_steps = 100000;
    int steps = 0;
    
    while (state.t_bar < target_t && steps < max_steps) {
        double remaining = target_t - state.t_bar;
        double step = min(h, remaining);
        if (step < EPSILON) break;
        rungeKutta4(state, step);
        steps++;
    }
    return state;
}

// ============================================================================
// 黄金分割法求最大压力点（图2-7）
// ============================================================================

ResultPoint findMaxPressurePoint(State state_before, State state_after, double h) {
    double ta = state_before.t_bar;
    double tb = state_after.t_bar;
    double h_current = tb - ta;
    
    if (h_current < EPSILON) {
        ResultPoint rp;
        double l = toL(state_before.l_bar);
        double v = toV(state_before.v_bar);
        double p = toP(state_before.p_bar) / 1e6;
        rp.l_dm = l * 10;  // convert to dm
        rp.t_ms = toT(state_before.t_bar) * 1000;
        rp.v = v;
        rp.p = p;
        rp.p_d = calculatePd(p, v, l);
        rp.p_t = calculatePt(p, v, l);
        rp.psi = state_before.psi;
        rp.Z = state_before.Z;
        rp.label = "最大压力点";
        return rp;
    }
    
    double dt1 = (1.0 - GOLDEN_RATIO) * h_current;
    double dt2 = GOLDEN_RATIO * h_current;
    
    State s1 = integrateForTime(state_before, dt1, h);
    State s2 = integrateForTime(state_before, dt2, h);
    
    double p1 = s1.p_bar;
    double p2 = s2.p_bar;
    
    int max_iter = 50;
    int iter = 0;
    
    while (h_current > 1e-8 && iter < max_iter) {
        if (p1 > p2) {
            tb = ta + dt2;
            dt2 = dt1;
            p2 = p1;
            s2 = s1;
            h_current = tb - ta;
            dt1 = (1.0 - GOLDEN_RATIO) * h_current;
            s1 = integrateForTime(state_before, dt1, h);
            p1 = s1.p_bar;
        } else {
            ta = ta + dt1;
            dt1 = dt2;
            p1 = p2;
            s1 = s2;
            h_current = tb - ta;
            dt2 = GOLDEN_RATIO * h_current;
            s2 = integrateForTime(state_before, dt2, h);
            p2 = s2.p_bar;
        }
        iter++;
    }
    
    State state_max = (p1 > p2) ? s1 : s2;
    
    ResultPoint rp;
    double l = toL(state_max.l_bar);
    double v = toV(state_max.v_bar);
    double p = toP(state_max.p_bar) / 1e6;
    rp.l_dm = l * 10;
    rp.t_ms = toT(state_max.t_bar) * 1000;
    rp.v = v;
    rp.p = p;
    rp.p_d = calculatePd(p, v, l);
    rp.p_t = calculatePt(p, v, l);
    rp.psi = state_max.psi;
    rp.Z = state_max.Z;
    rp.label = "最大压力点";
    return rp;
}

// ============================================================================
// 输入输出函数
// ============================================================================

bool readInput(const string& filename) {
    ifstream fin(filename);
    if (!fin.is_open()) {
        cerr << "错误：无法打开输入文件 " << filename << endl;
        return false;
    }
    
    fin >> input.S >> input.V0 >> input.lg >> input.m;
    fin >> input.f >> input.omega >> input.alpha >> input.rho_p >> input.theta;
    fin >> input.u1 >> input.n >> input.e1;
    fin >> input.chi >> input.lambda >> input.mu;
    fin >> input.chi_s >> input.lambda_s;
    fin >> input.p0;
    fin >> input.phi1 >> input.lambda2;
    fin >> input.h;
    
    fin.close();
    return true;
}

void printInput() {
    cout << "\n========== 输入参数 ==========" << endl;
    cout << "\n【火炮构造及弹丸诸元】" << endl;
    cout << "炮膛截面积 S = " << input.S << " m^2" << endl;
    cout << "药室容积 V0 = " << input.V0 << " m^3" << endl;
    cout << "炮膛全长 lg = " << input.lg << " m" << endl;
    cout << "弹丸质量 m = " << input.m << " kg" << endl;
    
    cout << "\n【装药条件】" << endl;
    cout << "火药力 f = " << input.f << " J/kg" << endl;
    cout << "装药量 omega = " << input.omega << " kg" << endl;
    cout << "余容 alpha = " << input.alpha << " m^3/kg" << endl;
    cout << "火药密度 rho_p = " << input.rho_p << " kg/m^3" << endl;
    cout << "火药特征量 theta = " << input.theta << endl;
    cout << "燃速系数 u1 = " << input.u1 << " m/(s·Pa^n)" << endl;
    cout << "燃速指数 n = " << input.n << endl;
    cout << "火药弧厚的一半 e1 = " << input.e1 << " m" << endl;
    cout << "火药形状特征量 chi = " << input.chi << endl;
    cout << "火药形状特征量 lambda = " << input.lambda << endl;
    cout << "火药形状特征量 mu = " << input.mu << endl;
    cout << "火药形状特征量 chi_s = " << input.chi_s << endl;
    cout << "火药形状特征量 lambda_s = " << input.lambda_s << endl;
    
    cout << "\n【起始条件】" << endl;
    cout << "起始压力 p0 = " << input.p0 << " Pa" << endl;
    
    cout << "\n【计算常数】" << endl;
    cout << "次要功系数 phi1 = " << input.phi1 << endl;
    cout << "系数 lambda2 = " << input.lambda2 << endl;
    
    cout << "\n【计算条件】" << endl;
    cout << "步长 h = " << input.h << endl;
}

void printConstants() {
    cout << "\n========== 计算常量 ==========" << endl;
    cout << "次要功系数 phi = " << consts.phi << endl;
    cout << "装填密度 Delta = " << consts.Delta << endl;
    cout << "药室长度 l0 = " << consts.l0 << " m" << endl;
    cout << "特征速度 vj = " << consts.vj << " m/s" << endl;
    cout << "系数 B = " << consts.B << endl;
    cout << "量纲为1的炮膛全长 lg_bar = " << consts.lg_bar << endl;
    cout << "火药燃尽时相对燃烧厚度 Zk = " << consts.Zk << endl;
}

// 创建结果点
ResultPoint createResultPoint(const State& state, const string& label) {
    ResultPoint rp;
    double l = toL(state.l_bar);
    double v = toV(state.v_bar);
    double p = toP(state.p_bar) / 1e6;  // Convert to MPa
    rp.l_dm = l * 10;  // convert to dm
    rp.t_ms = toT(state.t_bar) * 1000;  // convert to ms
    rp.v = v;
    rp.p = p;
    rp.p_d = calculatePd(p, v, l);
    rp.p_t = calculatePt(p, v, l);
    rp.psi = state.psi;
    rp.Z = state.Z;
    rp.label = label;
    return rp;
}

void saveResultsCSV(const string& filename) {
    ofstream fout(filename);
    if (!fout.is_open()) {
        cerr << "警告：无法保存结果到 " << filename << endl;
        return;
    }
    
    fout << "l_dm,t_ms,v_m_s,p_MPa,pd_MPa,pt_MPa,psi,Z,label" << endl;
    
    for (const auto& rp : allResults) {
        fout << fixed << setprecision(4) << rp.l_dm << ","
             << fixed << setprecision(4) << rp.t_ms << ","
             << fixed << setprecision(2) << rp.v << ","
             << fixed << setprecision(2) << rp.p << ","
             << fixed << setprecision(2) << rp.p_d << ","
             << fixed << setprecision(2) << rp.p_t << ","
             << fixed << setprecision(6) << rp.psi << ","
             << fixed << setprecision(6) << rp.Z << ","
             << rp.label << endl;
    }
    
    fout.close();
    cout << "\n结果已保存到: " << filename << endl;
}

// 保存曲线数据
void saveCurveData(const string& prefix) {
    // p-t 曲线数据
    ofstream fpt(prefix + "_p_t.csv");
    fpt << "t_ms,p_MPa" << endl;
    for (const auto& rp : allResults) {
        fpt << rp.t_ms << "," << rp.p << endl;
    }
    fpt.close();
    
    // p-l 曲线数据
    ofstream fpl(prefix + "_p_l.csv");
    fpl << "l_dm,p_MPa" << endl;
    for (const auto& rp : allResults) {
        fpl << rp.l_dm << "," << rp.p << endl;
    }
    fpl.close();
    
    // v-t 曲线数据
    ofstream fvt(prefix + "_v_t.csv");
    fvt << "t_ms,v_m_s" << endl;
    for (const auto& rp : allResults) {
        fvt << rp.t_ms << "," << rp.v << endl;
    }
    fvt.close();
    
    // v-l 曲线数据
    ofstream fvl(prefix + "_v_l.csv");
    fvl << "l_dm,v_m_s" << endl;
    for (const auto& rp : allResults) {
        fvl << rp.l_dm << "," << rp.v << endl;
    }
    fvl.close();
    
    cout << "曲线数据已保存到: " << prefix << "_*.csv" << endl;
}

// ============================================================================
// 常量计算（公式2-71下方）
// ============================================================================

void calculateConstants() {
    // phi = phi1 + lambda2 * omega / m
    consts.phi = input.phi1 + input.lambda2 * input.omega / input.m;
    
    // Delta = omega / V0
    consts.Delta = input.omega / input.V0;
    
    // l0 = V0 / S
    consts.l0 = input.V0 / input.S;
    
    // vj = sqrt(2 * f * omega / (theta * phi * m))
    consts.vj = sqrt(2.0 * input.f * input.omega / (input.theta * consts.phi * input.m));
    
    // B = S^2 * e1^2 / (f * omega * phi * m * u1^2) * (f * Delta)^(2*(1-n))
    double term1 = input.S * input.S * input.e1 * input.e1 / 
                   (input.f * input.omega * consts.phi * input.m * input.u1 * input.u1);
    double term2 = pow(input.f * consts.Delta, 2.0 * (1.0 - input.n));
    consts.B = term1 * term2;
    
    // lg_bar = lg / l0
    consts.lg_bar = input.lg / consts.l0;
    
    // Zk = 1.5 (根据火药形状确定)
    consts.Zk = 1.5;
    
    // psi_s = chi * (1 + lambda + mu)
    consts.psi_s = input.chi * (1.0 + input.lambda + input.mu);
    
    // xi_s = 1 / Zk
    consts.xi_s = 1.0 / consts.Zk;
}

// ============================================================================
// 初值计算（公式2-71下方）
// ============================================================================

State calculateInitialValues() {
    State state;
    
    state.t_bar = 0.0;
    state.l_bar = 0.0;
    state.v_bar = 0.0;
    state.p_bar = input.p0 / (input.f * consts.Delta);
    
    // psi0 = (1/Delta - 1/rho_p) / (f/p0 + (alpha - 1/rho_p))
    double numerator = 1.0 / consts.Delta - 1.0 / input.rho_p;
    double denominator = input.f / input.p0 + (input.alpha - 1.0 / input.rho_p);
    state.psi = numerator / denominator;
    
    // Z0 = (sqrt(1 + 4*lambda*psi0/chi) - 1) / (2*lambda)
    if (input.lambda > EPSILON) {
        double term = sqrt(1.0 + 4.0 * input.lambda * state.psi / input.chi);
        state.Z = (term - 1.0) / (2.0 * input.lambda);
    } else {
        state.Z = state.psi / input.chi;
    }
    
    return state;
}

// ============================================================================
// 主弹道循环计算
// ============================================================================

void ballisticCalculation() {
    State state = calculateInitialValues();
    double h = input.h;
    
    cout << "\n========== 开始弹道计算 ==========" << endl;
    
    // 保存初始点
    allResults.push_back(createResultPoint(state, "起始点"));
    
    // 状态标志
    bool maxPressureFound = false;
    bool splitPointFound = false;
    bool burnoutPointFound = false;
    bool muzzlePointFound = false;
    
    double p_bar_prev = state.p_bar;
    
    int max_steps = 500000;
    int step = 0;
    int output_counter = 0;
    
    // 主循环
    while (step < max_steps && !muzzlePointFound) {
        State state_prev = state;
        
        // 龙格-库塔法积分一步
        rungeKutta4(state, h);
        
        // 每一定步数保存一个点（用于绘制曲线）
        output_counter++;
        if (output_counter >= 50) {  // 每50步保存一个点
            allResults.push_back(createResultPoint(state, ""));
            output_counter = 0;
        }
        
        // 检查最大压力点（必须在Z<1时）
        if (!maxPressureFound && state.Z < 1.0 && state.p_bar < p_bar_prev) {
            ResultPoint rp_max = findMaxPressurePoint(state_prev, state, h);
            allResults.push_back(rp_max);
            specialPoints.push_back(rp_max);
            maxPressureFound = true;
            cout << "已找到最大压力点: p = " << fixed << setprecision(1) << rp_max.p << " MPa" << endl;
        }
        p_bar_prev = state.p_bar;
        
        // 检查燃烧分裂点 (Z = 1)
        if (!splitPointFound && state.Z >= 1.0 && state_prev.Z < 1.0) {
            double ratio = (1.0 - state_prev.Z) / (state.Z - state_prev.Z);
            State state_split;
            state_split.t_bar = state_prev.t_bar + ratio * (state.t_bar - state_prev.t_bar);
            state_split.l_bar = state_prev.l_bar + ratio * (state.l_bar - state_prev.l_bar);
            state_split.v_bar = state_prev.v_bar + ratio * (state.v_bar - state_prev.v_bar);
            state_split.p_bar = state_prev.p_bar + ratio * (state.p_bar - state_prev.p_bar);
            state_split.psi = state_prev.psi + ratio * (state.psi - state_prev.psi);
            state_split.Z = 1.0;
            
            ResultPoint rp_split = createResultPoint(state_split, "燃烧分裂点");
            allResults.push_back(rp_split);
            specialPoints.push_back(rp_split);
            splitPointFound = true;
            cout << "已找到燃烧分裂点: Z = 1.0" << endl;
        }
        
        // 检查燃烧结束点 (Z = Zk)
        if (!burnoutPointFound && state.Z >= consts.Zk && state_prev.Z < consts.Zk) {
            double ratio = (consts.Zk - state_prev.Z) / (state.Z - state_prev.Z);
            State state_burnout;
            state_burnout.t_bar = state_prev.t_bar + ratio * (state.t_bar - state_prev.t_bar);
            state_burnout.l_bar = state_prev.l_bar + ratio * (state.l_bar - state_prev.l_bar);
            state_burnout.v_bar = state_prev.v_bar + ratio * (state.v_bar - state_prev.v_bar);
            state_burnout.p_bar = state_prev.p_bar + ratio * (state.p_bar - state_prev.p_bar);
            state_burnout.psi = state_prev.psi + ratio * (state.psi - state_prev.psi);
            state_burnout.Z = consts.Zk;
            
            ResultPoint rp_burnout = createResultPoint(state_burnout, "燃烧结束点");
            allResults.push_back(rp_burnout);
            specialPoints.push_back(rp_burnout);
            burnoutPointFound = true;
            cout << "已找到燃烧结束点: Z = " << consts.Zk << endl;
        }
        
        // 检查炮口点 (l_bar = lg_bar)
        if (state.l_bar >= consts.lg_bar) {
            double ratio = (consts.lg_bar - state_prev.l_bar) / (state.l_bar - state_prev.l_bar);
            State state_muzzle;
            state_muzzle.t_bar = state_prev.t_bar + ratio * (state.t_bar - state_prev.t_bar);
            state_muzzle.l_bar = consts.lg_bar;
            state_muzzle.v_bar = state_prev.v_bar + ratio * (state.v_bar - state_prev.v_bar);
            state_muzzle.p_bar = state_prev.p_bar + ratio * (state.p_bar - state_prev.p_bar);
            state_muzzle.psi = state_prev.psi + ratio * (state.psi - state_prev.psi);
            state_muzzle.Z = state_prev.Z + ratio * (state.Z - state_prev.Z);
            
            ResultPoint rp_muzzle = createResultPoint(state_muzzle, "炮口点");
            allResults.push_back(rp_muzzle);
            specialPoints.push_back(rp_muzzle);
            muzzlePointFound = true;
            cout << "已找到炮口点: v = " << fixed << setprecision(1) << rp_muzzle.v << " m/s" << endl;
        }
        
        step++;
    }
    
    if (step >= max_steps) {
        cout << "警告：达到最大迭代次数，计算可能未收敛" << endl;
    }
    
    cout << "弹道计算完成，共迭代 " << step << " 步" << endl;
}

// ============================================================================
// 主函数
// ============================================================================

int main(int argc, char* argv[]) {
    cout << "========================================" << endl;
    cout << "    枪炮内弹道学数值解法程序" << endl;
    cout << "    基于四阶龙格-库塔法" << endl;
    cout << "========================================" << endl;
    
    string inputFile = "input_clean.txt";
    if (argc > 1) {
        inputFile = argv[1];
    }
    
    cout << "\n读取输入文件: " << inputFile << endl;
    
    if (!readInput(inputFile)) {
        cerr << "程序退出" << endl;
        return 1;
    }
    
    printInput();
    calculateConstants();
    printConstants();
    
    ballisticCalculation();
    
    // 保存结果
    saveResultsCSV("./output/results_table.csv");
    saveCurveData("./output/curve");
    
    cout << "\n程序执行完毕" << endl;
    
    return 0;
}