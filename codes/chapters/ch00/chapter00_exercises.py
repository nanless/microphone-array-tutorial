"""Run E00-01..E00-16 mathematical foundation exercises with pure NumPy; never write files or access network.

Usage:
    .venv/bin/python -m codes.chapters.ch00.chapter00_exercises
    .venv/bin/python codes/chapters/ch00/chapter00_exercises.py --list
    .venv/bin/python codes/chapters/ch00/chapter00_exercises.py --exercise E00-11 --json
"""
from __future__ import annotations

if __name__ == "__main__" and not __package__:
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

import argparse
import json
import math
from typing import Callable
import numpy as np

from codes.chapters.ch00.core.chapter00_math import (
    bessel_legendre_and_spherical_harmonics_example,
    dft_parseval_and_zero_padded_conv_example,
    hungarian_bipartite_assignment_example,
    nnls_kkt_and_mm_reweighted_l1_example,
    orthogonal_projection_and_four_subspaces_example,
    phasor_delay_and_decibel_budget_example,
    proper_complex_gaussian_and_crlb_example,
    rayleigh_quotient_and_gev_diagonalization_example,
    so3_rotation_matrix_composition_example,
    spherical_and_fresnel_geometry_example,
    stft_cola_window_reconstruction_example,
    svd_condition_and_tikhonov_loading_example,
    vec_kronecker_and_khatri_rao_example,
    wiener_orthogonality_and_kalman_joseph_example,
    wirtinger_lagrange_mvdr_and_lcmv_example,
    woodbury_and_schur_complement_example,
)

EXERCISE_IDS = tuple(f"E00-{i:02d}" for i in range(1, 17))

EXERCISE_TITLES: dict[str, str] = {
    "E00-01": "球面/AES69 坐标变换、空间方向余弦与二阶菲涅耳球面波程差近似",
    "E00-02": "三维正交旋转矩阵 SO(3) 复合、逆变换保范性与刚体坐标映射",
    "E00-03": "复相量合成、微秒级时延相位斜率与幅值/功率分贝换算",
    "E00-04": "DFT 酉性帕塞瓦尔能量守恒、循环卷积混叠与补零线性卷积精确等价",
    "E00-05": "汉宁窗与开方汉宁分析—合成窗的 COLA 重叠相加完美重建检验",
    "E00-06": "复列空间正交投影算子、幂等厄米性与四基本子空间正交分解",
    "E00-07": "奇异值分解二范数条件数、病态逆噪声放大与 Tikhonov 对角加载正则化",
    "E00-08": "Sherman–Morrison–Woodbury 矩阵求逆引理、秩一协方差递推与分块 Schur 补行列式",
    "E00-09": "厄米矩阵瑞利商极值、广义特征值问题（GEV）与白化合同联合对角化",
    "E00-10": "列主序向量化 vec(AXB)、Kronecker 积谱性质与 Khatri–Rao 积对角协方差恒等式",
    "E00-11": "Wirtinger 共轭偏导数验证与复线性等式约束拉格朗日乘子法闭式求解 MVDR/LCMV",
    "E00-12": "非负最小二乘（NNLS）KKT 互补松弛条件检验与主化—最小化（MM）替代上界单调下降",
    "E00-13": "圆对称复高斯（Proper vs Improper）伪协方差判据、实频点例外与标量时延 CRLB",
    "E00-14": "LMMSE 维纳—霍夫正交原理验证与线性高斯卡尔曼预测—新息—Joseph 协方差更新",
    "E00-15": "置换不变代价矩阵、贪心局部指派次优反例与匈牙利二分图全局最优指派",
    "E00-16": "柱/球贝塞尔函数递推、奇次 Legendre 多项式区间正交解耦与球谐函数加法定理检验",
}

EXERCISE_DISPATCH: dict[str, Callable[[], dict]] = {
    "E00-01": spherical_and_fresnel_geometry_example,
    "E00-02": so3_rotation_matrix_composition_example,
    "E00-03": phasor_delay_and_decibel_budget_example,
    "E00-04": dft_parseval_and_zero_padded_conv_example,
    "E00-05": stft_cola_window_reconstruction_example,
    "E00-06": orthogonal_projection_and_four_subspaces_example,
    "E00-07": svd_condition_and_tikhonov_loading_example,
    "E00-08": woodbury_and_schur_complement_example,
    "E00-09": rayleigh_quotient_and_gev_diagonalization_example,
    "E00-10": vec_kronecker_and_khatri_rao_example,
    "E00-11": wirtinger_lagrange_mvdr_and_lcmv_example,
    "E00-12": nnls_kkt_and_mm_reweighted_l1_example,
    "E00-13": proper_complex_gaussian_and_crlb_example,
    "E00-14": wiener_orthogonality_and_kalman_joseph_example,
    "E00-15": hungarian_bipartite_assignment_example,
    "E00-16": bessel_legendre_and_spherical_harmonics_example,
}


def _json(value):
    if isinstance(value, np.ndarray):
        if np.iscomplexobj(value):
            stacked = np.stack([value.real, value.imag], axis=-1)
            return _json(stacked.tolist())
        return [_json(v) for v in value.tolist()]
    if isinstance(value, (complex, np.complexfloating)):
        r_val, i_val = float(value.real), float(value.imag)
        if not (math.isfinite(r_val) and math.isfinite(i_val)):
            raise ValueError("non-finite complex scalar encountered during JSON serialization")
        return [r_val, i_val]
    if isinstance(value, np.generic):
        return _json(value.item())
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("non-finite float encountered during JSON serialization")
        return value
    if isinstance(value, dict):
        return {str(key): _json(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json(item) for item in value]
    return value


def run_exercise(exercise_id: str) -> dict:
    if exercise_id not in EXERCISE_DISPATCH:
        raise ValueError(f"unknown Chapter 0 exercise ID: {exercise_id!r}")
    return _json(EXERCISE_DISPATCH[exercise_id]())


def run_experiments() -> dict:
    if tuple(EXERCISE_DISPATCH) != EXERCISE_IDS or tuple(EXERCISE_TITLES) != EXERCISE_IDS:
        raise AssertionError("Chapter 0 exercise registry differs from E00-01..E00-16")
    exercises = {ex_id: EXERCISE_DISPATCH[ex_id]() for ex_id in EXERCISE_IDS}
    return _json(
        {
            "schema_version": 1,
            "scope": "Chapter 0 deterministic pure-NumPy mathematical foundations exercises (E00-01..E00-16)",
            "complex_encoding": "complex scalar -> [real,imag], complex arrays -> trailing [real,imag] axis",
            "exercises": exercises,
        }
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--list", action="store_true", help="list available E00-01..E00-16 exercise IDs")
    parser.add_argument("--exercise", choices=EXERCISE_IDS, help="run a single exercise by ID")
    parser.add_argument("--json", action="store_true", help="emit strict finite UTF-8 JSON output")
    args = parser.parse_args(argv)

    if args.list:
        if args.json:
            payload = {
                "schema_version": 1,
                "exercise_ids": list(EXERCISE_IDS),
                "titles": EXERCISE_TITLES,
            }
            print(json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False))
        else:
            for ex_id in EXERCISE_IDS:
                print(f"{ex_id}: {EXERCISE_TITLES[ex_id]}")
        return 0

    if args.exercise:
        payload = {
            "schema_version": 1,
            "scope": "Chapter 0 deterministic pure-NumPy mathematical foundations exercises (E00-01..E00-16)",
            "complex_encoding": "complex scalar -> [real,imag], complex arrays -> trailing [real,imag] axis",
            "exercises": {args.exercise: run_exercise(args.exercise)},
        }
    else:
        payload = run_experiments()

    print(json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
