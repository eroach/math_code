from sympy import symbols, expand, Matrix, Rational, Integer, collect, zeros, simplify
from sympy.ntheory import isprime
import time
import warnings
import sys
import os
from itertools import product
from functools import lru_cache

# 忽略警告信息，保持输出整洁
warnings.filterwarnings("ignore")

# ===================== 全局配置 =====================
MIN_PRIME = 5  # 最小质数≥5，确保n-3≥2
MAX_PRIME = 1000  # 扩展到1000以内的质数
MAX_ENUM_SOLUTIONS = None  # 输出所有解，不限制数量
PRINT_PROGRESS_INTERVAL = 100  # 进度打印间隔
BASE_SAVE_DIR = "prime_results_1000"  # 基础保存目录
REGULAR_DIR = os.path.join(BASE_SAVE_DIR, "regular_primes")  # 正则素数目录
IRREGULAR_DIR = os.path.join(BASE_SAVE_DIR, "irregular_primes")  # 非正则素数目录
m = symbols('m', integer=True)  # 全局符号变量

# 1000以内的非正则素数列表
IRREGULAR_PRIMES_1000 = [
    37, 59, 67, 101, 103, 131, 149, 151, 157, 163, 167, 173, 179, 181, 191,
    193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271,
    277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367,
    373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457,
    461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563,
    569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647,
    653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751,
    757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857,
    859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967,
    971, 977, 983, 991, 997
]


# ===================== 正则素数判定函数 =====================
def is_regular_prime(p):
    """判断质数p是否为正则素数（基于预定义的非正则素数列表）"""
    if not isprime(p):
        return False
    if p == 2 or p == 3:
        return True
    return p not in IRREGULAR_PRIMES_1000


def get_prime_type_info(p):
    """获取质数类型信息"""
    if is_regular_prime(p):
        return "regular", "✅ 正则素数 (Regular Prime)", REGULAR_DIR
    else:
        return "irregular", "❌ 非正则素数 (Irregular Prime)", IRREGULAR_DIR


# ===================== 工具函数：解过滤 =====================
def is_trivial_solution(mj_sol):
    """判断是否为平凡解（所有m_j都等于0）"""
    if not mj_sol:
        return False
    return all(value == 0 for value in mj_sol.values())


def is_non_trivial_solution(mj_sol):
    """判断是否为非平凡解（至少有一个m_j不等于0）"""
    return not is_trivial_solution(mj_sol)


def filter_non_trivial_solutions(all_solutions):
    """
    过滤掉所有m_j都等于0的平凡解，只保留非平凡解
    返回：(非平凡解列表, 过滤掉的平凡解数量)
    """
    # 先去重（防止重复生成相同解）
    unique_solutions = []
    seen = set()
    for sol in all_solutions:
        sol_tuple = tuple(sorted(sol.items()))
        if sol_tuple not in seen:
            seen.add(sol_tuple)
            unique_solutions.append(sol)

    trivial_solutions = [sol for sol in unique_solutions if is_trivial_solution(sol)]
    trivial_count = len(trivial_solutions)
    non_trivial_solutions = [sol for sol in unique_solutions if is_non_trivial_solution(sol)]

    print(f"🔍 解过滤：总共找到{len(unique_solutions)}个唯一解")
    print(f"   - 移除{trivial_count}个平凡解（所有m_j=0）")
    print(f"   - 保留{len(non_trivial_solutions)}个非平凡解（至少一个m_j≠0）")

    return non_trivial_solutions, trivial_count


# ===================== 初始化保存目录 =====================
def init_save_directories():
    """初始化正则/非正则素数保存目录"""
    if not os.path.exists(BASE_SAVE_DIR):
        os.makedirs(BASE_SAVE_DIR)
        print(f"📁 创建基础保存目录：{BASE_SAVE_DIR}")

    if not os.path.exists(REGULAR_DIR):
        os.makedirs(REGULAR_DIR)
        print(f"📁 创建正则素数目录：{REGULAR_DIR}")
    else:
        print(f"📁 正则素数目录已存在：{REGULAR_DIR}")

    if not os.path.exists(IRREGULAR_DIR):
        os.makedirs(IRREGULAR_DIR)
        print(f"📁 创建非正则素数目录：{IRREGULAR_DIR}")
    else:
        print(f"📁 非正则素数目录已存在：{IRREGULAR_DIR}")


# ===================== 工具函数：质数生成 =====================
def generate_primes(min_p, max_p):
    """生成[min_p, max_p]内的所有奇质数（≥5）"""
    primes = []
    print(f"正在生成{min_p}~{max_p}范围内的质数列表...")

    for num in range(min_p, max_p + 1):
        if isprime(num) and num >= 5:
            primes.append(num)

    # 分类统计
    regular_primes = [p for p in primes if is_regular_prime(p)]
    irregular_primes_in_list = [p for p in primes if p in IRREGULAR_PRIMES_1000]

    print(f"生成{min_p}~{max_p}范围内的奇质数列表：共{len(primes)}个")
    print(f"  - 正则素数：{len(regular_primes)}个")
    print(f"  - 非正则素数：{len(irregular_primes_in_list)}个")
    print(f"前10个质数：{primes[:10]}...")
    print(f"前10个非正则素数：{irregular_primes_in_list[:10]}...")

    return primes


# ===================== 工具函数：格式化输出 =====================
def print_separator(title):
    """打印分隔符，仅保留关键标题"""
    print("\n" + "=" * 80)
    print(f"【{title}】")
    print("=" * 80)


def print_matrix(mat, name="矩阵"):
    """简化矩阵打印，仅显示维度和关键行"""
    if len(mat) == 0:
        print(f"\n{name}：空矩阵")
        return
    if len(mat) <= 5:
        print(f"\n{name}（{len(mat)}行×{len(mat[0])}列）：")
        for i, row in enumerate(mat):
            print(f"  第{i + 1}行：{row}")
    else:
        print(f"\n{name}：{len(mat)}行×{len(mat[0])}列")


# ===================== 模运算工具函数 =====================
def mod_pow_custom(base, exp, mod):
    """自定义快速模幂运算"""
    if mod == 1:
        return 0
    result = 1
    base = base % mod
    while exp > 0:
        if exp % 2 == 1:
            result = (result * base) % mod
        exp = exp // 2
        base = (base * base) % mod
    return result


def mod_pow_builtin(base, exp, mod):
    """使用Python内置pow函数的快速模幂（更高效）"""
    if mod == 1:
        return 0
    return pow(base, exp, mod)


def compute_a_value(mj_dict, prime):
    """计算a = ∏(k^m_k) mod n²"""
    prime_square = prime * prime
    var_count = (prime - 1) // 2
    max_k = 1 + var_count
    a = 1
    for k in range(2, max_k + 1):
        mj_key = f"m{k}"
        mj_val = mj_dict.get(mj_key, 0)
        if mj_val == 0:
            continue
        term = mod_pow_custom(k, mj_val, prime_square)
        a = (a * term) % prime_square
    print(f"  a = ∏(k^m_k) mod {prime}² = {a}")
    return a


def compute_a_n_minus_1_double_verify(mj_dict, prime):
    """双验证计算a^(n-1) mod n²"""
    prime_square = prime * prime
    exp = prime - 1

    a_mod_nsq = compute_a_value(mj_dict, prime)
    a_exp_custom = mod_pow_custom(a_mod_nsq, exp, prime_square)
    a_exp_builtin = mod_pow_builtin(a_mod_nsq, exp, prime_square)

    is_consistent = (a_exp_custom == a_exp_builtin)
    print(f"  a^{exp} mod {prime}² = {a_exp_custom}")
    if not is_consistent:
        print(f"  ❌ 警告：内置计算结果={a_exp_builtin}（不一致）")

    return {
        "a_mod_nsq": a_mod_nsq,
        "a_exp_custom": a_exp_custom,
        "a_exp_builtin": a_exp_builtin,
        "is_consistent": is_consistent,
        "quotient_custom": a_exp_custom // prime if prime != 0 else 0,
        "remainder_custom": a_exp_custom % prime if prime != 0 else 0,
        "quotient_builtin": a_exp_builtin // prime if prime != 0 else 0,
        "remainder_builtin": a_exp_builtin % prime if prime != 0 else 0
    }


# ===================== 核心工具函数：符号表达式模运算 =====================
def mod_symbolic_expr(expr, prime):
    """对带符号m的表达式进行模运算处理"""
    if prime == 1:
        return Integer(0)

    expr_expanded = expand(expr, deep=True, rational=True)
    expr_collected = collect(expr_expanded, m)

    coeffs = expr_collected.as_coefficients_dict()
    new_coeffs = {}

    for var, coeff in coeffs.items():
        if isinstance(coeff, Rational):
            numer = int(coeff.numerator) % prime
            denom = int(coeff.denominator) % prime

            if denom == 0:
                new_coeff = Integer(0)
            else:
                inv_denom = pow(denom, prime - 2, prime)
                new_coeff = Integer((numer * inv_denom) % prime)
        elif isinstance(coeff, (int, Integer)):
            new_coeff = Integer(int(coeff) % prime)
        else:
            new_coeff = Integer(int(coeff) % prime)

        if new_coeff < 0:
            new_coeff += prime

        new_coeffs[var] = new_coeff

    new_expr = sum([coeff * var for var, coeff in new_coeffs.items()])
    return new_expr


# ===================== P_k计算（核心逻辑） =====================
def factorial_mod(n_val, prime):
    """计算阶乘并即时模n，避免大数溢出"""
    if n_val < 0 or prime == 1:
        return 0
    if n_val == 0 or n_val == 1:
        return 1
    res = 1
    for i in range(2, n_val + 1):
        res = (res * i) % prime
        if res == 0:  # 提前终止，阶乘包含prime因子
            break
    return res


def e_k_optimized(k, prime):
    """优化版e(k)计算，减少输出"""
    if k < 1 or prime == 1:
        return Integer(0)

    # 计算分子：(m-1)(m-2)...(m-k)
    numer = Integer(1)
    for i in range(1, k + 1):
        numer = numer * (m - Integer(i))

    # 计算分母：(k+1)!
    denom_fact = factorial_mod(k + 1, prime)
    if denom_fact == 0:
        ek = Integer(0)
    else:
        # 计算分母逆元
        inv_denom = pow(denom_fact, prime - 2, prime)
        ek = numer * Integer(inv_denom)

    ek_mod = mod_symbolic_expr(ek, prime)
    return ek_mod


def compute_pk_optimized(prime):
    """计算每个质数prime专属的P_k(m)表达式，简化输出"""
    if prime < 5:
        print(f"❌ 质数{prime}小于5，跳过")
        return [], 0

    max_pk_index = prime - 3
    print_separator(f"处理质数 {prime} - 步骤1：计算Pk表达式")
    pk_list = [None] * (max_pk_index + 1)

    # 计算P1
    pk1 = (m - 1) / 2
    pk1_mod = mod_symbolic_expr(pk1, prime)
    pk_list[1] = pk1_mod

    # 递推计算Pk（k≥2）
    for k in range(2, max_pk_index + 1):
        pk = Integer(0)

        for i in range(1, k):
            if k - i > max_pk_index or pk_list[k - i] is None:
                continue

            sign = Integer((-1) ** (i + 1))
            ei = e_k_optimized(i, prime)
            p_ki = pk_list[k - i]

            term = sign * ei * p_ki
            term_mod = mod_symbolic_expr(term, prime)

            pk += term_mod
            pk = mod_symbolic_expr(pk, prime)

        # 最终项
        sign_final = Integer((-1) ** (k + 1))
        ek = e_k_optimized(k, prime)
        final_term = sign_final * Integer(k) * ek
        final_term_mod = mod_symbolic_expr(final_term, prime)

        pk += final_term_mod
        pk = mod_symbolic_expr(pk, prime)
        pk_list[k] = pk

    print(f"✅ P1~P{max_pk_index} 表达式计算完成")
    return pk_list, max_pk_index


def build_mod_n_matrix_optimized(pk_list, max_pk_index, prime):
    """构建系数矩阵，简化输出"""
    if prime < 5 or max_pk_index < 1:
        return [], []

    var_count = (prime - 1) // 2
    M_VALUES = list(range(2, 2 + var_count))

    print_separator(f"处理质数 {prime} - 步骤2：构建系数矩阵")
    print(f"变量数：{var_count} | 代入m值：{M_VALUES}")

    A_list = []
    # 目标PK列表
    target_pk = [1, 2] + list(range(4, max_pk_index + 1, 2))
    target_pk = target_pk[:var_count]

    for row_idx, k in enumerate(target_pk):
        if k >= len(pk_list) or pk_list[k] is None:
            row_values = [0] * var_count
        else:
            pk = pk_list[k]
            row_values = []

            for col_idx, mj in enumerate(M_VALUES):
                val = pk.subs(m, Integer(mj))
                val = simplify(val)

                # 统一转为模prime整数
                if isinstance(val, Rational):
                    numer = int(val.numerator) % prime
                    denom = int(val.denominator) % prime
                    inv_denom = pow(denom, prime - 2, prime) if (denom % prime != 0 and prime > 1) else 0
                    val_mod = (numer * inv_denom) % prime
                else:
                    val_mod = int(val) % prime

                val_mod = val_mod if val_mod >= 0 else val_mod + prime
                row_values.append(val_mod)

        A_list.append(row_values)

    print(f"✅ {var_count}×{var_count} 系数矩阵构建完成")
    return A_list, M_VALUES


# ===================== 矩阵求解 =====================
def compute_determinant_mod_n(A_list, prime):
    """计算行列式，简化输出"""
    if not A_list or len(A_list) != len(A_list[0]) or prime == 1:
        print(f"❌ 无效矩阵或质数，无法计算行列式")
        return 1, None

    var_count = len(A_list)
    print_separator(f"处理质数 {prime} - 步骤3：行列式分析")

    try:
        A_sympy = Matrix(A_list)
        det = A_sympy.det()
        det_mod = det % prime

        print(f"行列式模{prime} = {det_mod}")
        if det_mod == 0:
            print(f"结论：存在非平凡解")
        else:
            print(f"结论：仅存在唯一平凡解（全零解）")
        return det_mod, A_sympy
    except Exception as e:
        print(f"❌ 行列式计算错误：{e}")
        return None, None


def mod_p_rref(matrix, prime):
    """有限域𝔽ₚ上的RREF计算，简化输出"""
    if not matrix or prime == 1:
        print(f"❌ 无效矩阵或质数，无法计算RREF")
        return [], [], 0

    var_count = len(matrix)
    print_separator(f"处理质数 {prime} - 步骤4：RREF计算")

    mat = [row.copy() for row in matrix]
    n_rows = len(mat)
    n_cols = len(mat[0]) if n_rows > 0 else 0
    rank = 0
    pivots = []

    for col in range(n_cols):
        # 找主元行
        pivot_row = None
        for r in range(rank, n_rows):
            if mat[r][col] % prime != 0:
                pivot_row = r
                break

        if pivot_row is None:
            continue

        # 交换行
        if pivot_row != rank:
            mat[rank], mat[pivot_row] = mat[pivot_row], mat[rank]

        # 主元归一化
        pivot_val = mat[rank][col]
        if pivot_val == 0:
            continue
        inv_pivot = pow(pivot_val, prime - 2, prime)
        for c in range(col, n_cols):
            mat[rank][c] = (mat[rank][c] * inv_pivot) % prime

        # 消去其他行
        for r in range(n_rows):
            if r != rank and mat[r][col] % prime != 0:
                factor = mat[r][col]
                for c in range(col, n_cols):
                    mat[r][c] = (mat[r][c] - factor * mat[rank][c]) % prime

        pivots.append(col)
        rank += 1

    # 最终模处理
    for r in range(n_rows):
        for c in range(n_cols):
            mat[r][c] = mat[r][c] % prime

    free_vars_idx = [i for i in range(n_cols) if i not in pivots]
    print(f"✅ RREF计算完成 | 秩={rank} | 自由变量数={len(free_vars_idx)}")
    return mat, pivots, rank


def solve_mod_n_equations(A_list, det_mod, prime, M_VALUES):
    """
    求解模n同余方程组，输出所有解
    核心修复：行列式≠0时，仅返回1个平凡解，而非错误的2个
    """
    if det_mod is None or prime == 1:
        trivial_sol = {f"m{j}": 0 for j in M_VALUES}
        return [trivial_sol], 1  # 1个平凡解

    var_count = len(M_VALUES)
    mj_max = prime - 1
    print_separator(f"处理质数 {prime} - 步骤5：求解同余方程")

    # 情况1：行列式≠0 → 仅存在唯一平凡解（全零解）
    if det_mod != 0:
        trivial_sol = {f"m{j}": 0 for j in M_VALUES}
        all_solutions = [trivial_sol]
        print(f"✅ 仅存在唯一平凡解：{trivial_sol}")

        # 过滤平凡解
        non_trivial_solutions, trivial_count = filter_non_trivial_solutions(all_solutions)
        print(f"✅ 方程求解完成 | 总有效解数：{len(non_trivial_solutions)}（已过滤{trivial_count}个平凡解）")

        if non_trivial_solutions:
            print(f"   所有非平凡解：")
            for i, sol in enumerate(non_trivial_solutions):
                print(f"     解#{i + 1}：{sol}")
        else:
            print(f"   无有效非平凡解")

        return non_trivial_solutions, trivial_count

    # 情况2：行列式=0 → 存在非平凡解
    rref_mat, pivots, real_rank = mod_p_rref(A_list, prime)
    free_vars_idx = [i for i in range(var_count) if i not in pivots]
    n_free = len(free_vars_idx)

    # 无自由变量（理论上不会走到这里，因为det_mod=0）
    if n_free == 0:
        trivial_sol = {f"m{j}": 0 for j in M_VALUES}
        all_solutions = [trivial_sol]
        print(f"✅ 仅存在平凡解：{trivial_sol}")

        non_trivial_solutions, trivial_count = filter_non_trivial_solutions(all_solutions)
        print(f"✅ 方程求解完成 | 总有效解数：{len(non_trivial_solutions)}（已过滤{trivial_count}个平凡解）")

        if non_trivial_solutions:
            print(f"   所有非平凡解：")
            for i, sol in enumerate(non_trivial_solutions):
                print(f"     解#{i + 1}：{sol}")
        else:
            print(f"   无有效非平凡解")

        return non_trivial_solutions, trivial_count

    # 枚举自由变量
    free_var_ranges = [range(0, mj_max + 1) for _ in free_vars_idx]
    total_comb = (mj_max + 1) ** n_free
    print(f"开始枚举自由变量（总组合数：{total_comb}）")

    all_solutions = []
    enum_count = 0

    for idx, free_vals in enumerate(product(*free_var_ranges)):
        # 进度打印
        if idx % PRINT_PROGRESS_INTERVAL == 0:
            print(f"  进度：{idx}/{total_comb} | 已找到有效解：{enum_count}")

        # 初始化变量值
        var_vals = [0] * var_count
        for i, free_idx in enumerate(free_vars_idx):
            if i >= len(free_vals) or free_idx >= len(var_vals):
                continue
            var_vals[free_idx] = free_vals[i]

        # 计算主变量值
        for r in range(real_rank):
            if r >= len(rref_mat):
                continue
            pivot_col = pivots[r] if r < len(pivots) else -1
            if pivot_col < 0 or pivot_col >= var_count:
                continue

            row = rref_mat[r]
            # 计算右侧值
            rhs = 0
            for c in range(var_count):
                if c != pivot_col and c < len(row):
                    rhs = (rhs - row[c] * var_vals[c]) % prime

            # 主元系数应为1（已归一化）
            var_vals[pivot_col] = rhs % prime

        # 验证取值范围
        if not all(0 <= v <= mj_max for v in var_vals):
            continue

        # 二次验证：代入原矩阵
        sol_valid = True
        for row_idx in range(var_count):
            if row_idx >= len(A_list):
                continue
            dot_product = 0
            for c in range(var_count):
                if c >= len(A_list[row_idx]) or c >= len(var_vals):
                    continue
                dot_product = (dot_product + A_list[row_idx][c] * var_vals[c]) % prime
            if dot_product != 0:
                sol_valid = False
                break

        if sol_valid:
            enum_count += 1
            mj_sol = {f"m{M_VALUES[c]}": var_vals[c] for c in range(var_count) if c < len(M_VALUES)}
            all_solutions.append(mj_sol)

    # 过滤平凡解
    non_trivial_solutions, trivial_count = filter_non_trivial_solutions(all_solutions)

    print(f"✅ 方程求解完成 | 总有效解数：{len(non_trivial_solutions)}（已过滤{trivial_count}个平凡解）")
    if non_trivial_solutions:
        print(f"   所有非平凡解：")
        for i, sol in enumerate(non_trivial_solutions):
            print(f"     解#{i + 1}：{sol}")
    else:
        print(f"   无有效非平凡解")

    return non_trivial_solutions, trivial_count


# ===================== 多项式验证 =====================
@lru_cache(maxsize=1024)
def build_factor_poly(j, mod_n, max_deg):
    """构建因子多项式P_j(x)，修复后的正确版本"""
    if j <= 0 or mod_n == 1 or max_deg <= 0:
        return tuple([0] * max_deg)

    # 正确的多项式递推公式：(1 - x)^j = 1 - C(j,1)x + C(j,2)x² - ... + (-1)^j x^j
    poly = [0] * (max_deg + 1)
    poly[0] = 1  # 初始为1

    for step in range(j):
        new_poly = [0] * (max_deg + 1)
        # (1 - x) * 当前多项式
        for d in range(max_deg):
            if poly[d] == 0:
                continue
            new_poly[d] = (new_poly[d] + poly[d]) % mod_n  # 1 * x^d 项
            # 关键修复：减法转加法，确保符号正确
            sub_term = (mod_n - poly[d] % mod_n) % mod_n
            new_poly[d + 1] = (new_poly[d + 1] + sub_term) % mod_n  # -x * x^d 项
        poly = new_poly

    # 分子：1 - (1 - x)^j
    numerator = [0] * (max_deg + 1)
    numerator[0] = 1
    for d in range(max_deg + 1):
        numerator[d] = (numerator[d] - poly[d]) % mod_n

    # 最终多项式（去掉常数项，返回x^1到x^max_deg的系数）
    result = [0] * max_deg
    for d in range(1, max_deg + 1):
        result[d - 1] = numerator[d] % mod_n

    return tuple(result)


def poly_mult(p1, p2, mod_n, max_deg):
    """多项式乘法（模运算）"""
    result = [0] * (max_deg + 1)
    for d1 in range(len(p1)):
        if d1 > max_deg or p1[d1] == 0:
            continue
        for d2 in range(len(p2)):
            if d2 > max_deg or p2[d2] == 0:
                continue
            d = d1 + d2
            if d > max_deg:
                continue
            result[d] = (result[d] + p1[d1] * p2[d2]) % mod_n
    return result


def poly_pow(poly, exp, mod_n, max_deg):
    """多项式快速幂（模运算）"""
    if exp == 0:
        result = [0] * (max_deg + 1)
        result[0] = 1
        return result

    result = [0] * (max_deg + 1)
    result[0] = 1
    current = list(poly)

    while exp > 0:
        if exp % 2 == 1:
            result = poly_mult(result, current, mod_n, max_deg)
        current = poly_mult(current, current, mod_n, max_deg)
        exp = exp // 2

    return result


def verify_polynomial_only(mj_dict, prime, M_VALUES, print_detail=True):
    """验证多项式条件，简化输出"""
    if prime == 1:
        if print_detail:
            print(f"❌ 质数为1，无法验证")
        return False, "失败：无效质数"

    try:
        mj = {int(k.replace("m", "")): v for k, v in mj_dict.items()}
        mod_n = prime
        max_deg = prime - 1
        final_poly = [0] * (max_deg + 1)
        final_poly[0] = 1

        # 逐个因子计算乘积
        for j in M_VALUES:
            mj_val = mj.get(j, 0)
            if mj_val == 0:
                continue

            base_poly = build_factor_poly(j, mod_n, max_deg)
            pow_poly = poly_pow(base_poly, mj_val, mod_n, max_deg)
            final_poly = poly_mult(final_poly, pow_poly, mod_n, max_deg)

        # 检查条件：除常数项外所有系数≡0 mod n
        invalid_terms = [d for d in range(1, max_deg + 1) if final_poly[d] % mod_n != 0]
        if invalid_terms:
            if print_detail:
                print(f"❌ 多项式验证失败：b{invalid_terms} ≠ 0 mod{prime}")
            return False, f"失败：非零项={invalid_terms}"

        if print_detail:
            print(f"✅ 多项式验证通过")
        return True, "通过"

    except Exception as e:
        error_msg = f"异常：{str(e)}"
        if print_detail:
            print(f"❌ 多项式验证异常：{error_msg}")
        return False, error_msg


# ===================== 整合验证 =====================
def verify_mj_with_a_exp_calc(mj_solutions, prime, M_VALUES):
    """验证所有m_j解并计算a^(n-1) mod n²"""
    print_separator(f"处理质数 {prime} - 步骤6：验证与计算")
    valid_results = []
    total_solutions = len(mj_solutions)

    print(f"待验证非平凡解数：{total_solutions}")

    # 验证所有解（不分批，完整输出）
    for global_idx, mj_sol in enumerate(mj_solutions, 1):
        print(f"\n----- 验证解#{global_idx} -----")
        # 多项式验证
        is_poly_valid, poly_msg = verify_polynomial_only(mj_sol, prime, M_VALUES)

        if is_poly_valid:
            print(f"  解#{global_idx}：{mj_sol} → ✅ 通过验证")
            # 计算a^(n-1) mod n²
            a_exp_result = compute_a_n_minus_1_double_verify(mj_sol, prime)

            valid_results.append({
                "mj": mj_sol,
                "poly_verify": poly_msg,
                "a_mod_nsq": a_exp_result["a_mod_nsq"],
                "a_exp_custom": a_exp_result["a_exp_custom"],
                "a_exp_builtin": a_exp_result["a_exp_builtin"],
                "is_consistent": a_exp_result["is_consistent"]
            })
        else:
            print(f"  解#{global_idx}：{mj_sol} → ❌ 验证失败：{poly_msg}")

    print(f"\n验证总结：")
    print(f"  总非平凡解数：{total_solutions}")
    print(f"  通过多项式验证：{len(valid_results)}")

    return valid_results


# ===================== 结果保存 =====================
def save_prime_results(prime, pk_list, A_list, det_mod, mj_solutions, trivial_count, valid_results, elapsed, M_VALUES,
                       max_pk_index):
    """
    分目录保存结果：
    1. 正则素数保存到 regular_primes 目录
    2. 非正则素数保存到 irregular_primes 目录
    3. 完整保存所有有效解，不截断
    4. 正确统计过滤后的非平凡解数
    """
    # 获取质数类型信息
    prime_type, prime_type_desc, save_dir = get_prime_type_info(prime)
    is_regular = (prime_type == "regular")

    # 构建文件名
    filename = os.path.join(save_dir, f"prime_{prime}_results.txt")
    prime_square = prime * prime

    try:
        with open(filename, "w", encoding="utf-8") as f:
            # 基础信息（包含正则素数标识）
            f.write(f"质数 n = {prime} 完整验证结果\n")
            f.write(f"质数类型：{prime_type_desc}\n")
            f.write(f"生成时间：{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}\n")
            f.write(f"PK范围=P1~P{max_pk_index} | 变量数={(prime - 1) // 2} | 耗时={elapsed:.2f}秒\n")
            f.write(f"⚠️  已过滤{trivial_count}个平凡解（所有m_j=0），仅保留非平凡解\n")
            f.write("=" * 80 + "\n\n")

            # 1. 质数基本信息
            f.write("【1. 质数基本信息】\n")
            f.write(f"质数：{prime}\n")
            f.write(f"是否正则素数：{is_regular} ({prime_type})\n")
            f.write(f"所属范围：1000以内的奇质数\n")
            f.write(f"最小验证质数限制：{MIN_PRIME}\n")
            f.write(f"平凡解过滤：已移除{trivial_count}个平凡解（所有m_j=0）\n")
            f.write("\n")

            # 2. 矩阵信息
            f.write("【2. 系数矩阵信息】\n")
            f.write(f"矩阵维度：{len(A_list)}行×{len(A_list[0]) if A_list else 0}列\n")
            f.write(f"行列式模{prime}：{det_mod if det_mod is not None else 'N/A'}\n")
            f.write(f"矩阵内容：\n")
            for i, row in enumerate(A_list):
                f.write(f"  第{i + 1}行：{row}\n")
            f.write("\n")

            # 3. 所有非平凡解（已过滤平凡解）
            f.write("【3. 所有非平凡解（已过滤平凡解）】\n")
            total_before_filter = len(mj_solutions) + trivial_count
            f.write(f"总解数（过滤前）：{total_before_filter}\n")
            f.write(f"过滤掉的平凡解数：{trivial_count}\n")
            f.write(f"非平凡解总数：{len(mj_solutions)}\n")

            if mj_solutions:
                f.write("所有非平凡解列表（完整，未截断）：\n")
                for i, sol in enumerate(mj_solutions, 1):
                    f.write(f"  解#{i}：{sol}\n")
                f.write("\n")
            else:
                f.write("  无有效非平凡解\n")
            f.write("\n")

            # 4. 多项式验证结果
            f.write("【4. 多项式验证结果】\n")
            f.write(f"通过验证的有效解数：{len(valid_results)} / {len(mj_solutions)}\n")
            if valid_results:
                for i, res in enumerate(valid_results, 1):
                    f.write(f"\n===== 有效解#{i} =====\n")
                    f.write(f"m_j：{res['mj']}\n")
                    f.write(f"多项式验证：{res['poly_verify']}\n")
                    f.write(f"a mod {prime}² = {res['a_mod_nsq']}\n")
                    f.write(f"a^{prime - 1} mod {prime}² = {res['a_exp_custom']}\n")
                    f.write(f"内置函数计算结果：{res['a_exp_builtin']}\n")
                    f.write(f"结果一致性：{'✅ 一致' if res['is_consistent'] else '❌ 不一致'}\n")
            else:
                f.write("  无通过验证的有效解\n")
            f.write("\n")

            # 5. 性能统计
            f.write("【5. 性能统计】\n")
            f.write(f"- 总计算耗时：{elapsed:.2f}秒\n")
            f.write(f"- 变量范围：m{min(M_VALUES) if M_VALUES else '无'}~m{max(M_VALUES) if M_VALUES else '无'}\n")

            # 计算自由变量数
            if A_list and len(A_list) > 0:
                try:
                    mat = Matrix(A_list)
                    rank = mat.rank()
                    free_vars = len(A_list) - rank
                except:
                    free_vars = 0
            else:
                free_vars = 0
            f.write(f"- 自由变量数：{free_vars}\n")

            # 解枚举总组合数
            if det_mod != 0:
                enum_count = 1
            else:
                enum_count = (prime) ** free_vars if free_vars > 0 else 1
            f.write(f"- 解枚举总组合数：{enum_count}\n")
            f.write(f"- 过滤前总解数：{len(mj_solutions) + trivial_count}\n")
            f.write(f"- 过滤掉的平凡解数：{trivial_count}\n")
            f.write(f"- 过滤后非平凡解数：{len(mj_solutions)}\n")

        print(f"\n📄 质数{prime}({prime_type})的完整结果已保存到：{filename}")
        return True, filename

    except Exception as e:
        print(f"❌ 保存质数{prime}结果失败：{str(e)}")
        # 尝试保存到基础目录作为备用
        backup_filename = os.path.join(BASE_SAVE_DIR, f"prime_{prime}_{prime_type}_backup.txt")
        try:
            with open(backup_filename, "w", encoding="utf-8") as f:
                f.write(f"质数 {prime} ({prime_type}) 结果（备份）\n")
                f.write(f"过滤前总解数：{len(mj_solutions) + trivial_count}\n")
                f.write(f"过滤掉的平凡解数：{trivial_count}\n")
                f.write(f"非平凡解数：{len(mj_solutions)}\n")
                f.write(f"有效解数：{len(valid_results)}\n")
                if mj_solutions:
                    f.write(f"所有非平凡解：\n")
                    for i, sol in enumerate(mj_solutions, 1):
                        f.write(f"  解#{i}：{sol}\n")
            print(f"📄 已保存备份文件：{backup_filename}")
            return True, backup_filename
        except:
            print(f"❌ 备份文件保存也失败！")
            return False, ""


# ===================== 单质数处理主流程 =====================
def process_single_prime(prime, progress_idx, total_primes):
    """处理单个质数的完整流程，带进度显示"""
    print(f"\n" + "*" * 80)
    print(f"开始处理质数 {prime}（{progress_idx}/{total_primes}）")
    prime_type, _, _ = get_prime_type_info(prime)
    print(f"质数类型：{prime_type}")
    print(f"*" * 80)

    if prime < 5:
        print(f"\n❌ 跳过质数{prime}（需≥5）")
        save_prime_results(prime, [], [], None, [], 0, [], 0.0, [], 0)
        return None

    start_time = time.time()
    pk_list = []
    A_list = []
    det_mod = None
    mj_solutions = []
    trivial_count = 0
    valid_results = []
    M_VALUES = []
    max_pk_index = 0

    try:
        # 1. 计算Pk表达式
        pk_list, max_pk_index = compute_pk_optimized(prime)
        if not pk_list:
            raise Exception("Pk表达式计算失败")

        # 2. 构建模n矩阵
        A_list, M_VALUES = build_mod_n_matrix_optimized(pk_list, max_pk_index, prime)
        if not A_list or not M_VALUES:
            raise Exception("系数矩阵构建失败")

        # 3. 计算行列式
        det_mod, _ = compute_determinant_mod_n(A_list, prime)
        if det_mod is None:
            raise Exception("行列式计算失败")

        # 4. 求解同余方程
        mj_solutions, trivial_count = solve_mod_n_equations(A_list, det_mod, prime, M_VALUES)

        # 5. 验证解并计算a^(n-1) mod n²
        valid_results = verify_mj_with_a_exp_calc(mj_solutions, prime, M_VALUES)

        # 6. 保存结果
        elapsed = time.time() - start_time
        save_success, filename = save_prime_results(
            prime, pk_list, A_list, det_mod, mj_solutions, trivial_count,
            valid_results, elapsed, M_VALUES, max_pk_index
        )

        # 输出处理完成信息
        print_separator(f"质数 {prime} 处理完成")
        total_before_filter = len(mj_solutions) + trivial_count
        print(f"核心统计：")
        print(f"  - PK范围：P1~P{max_pk_index} | 变量数：{(prime - 1) // 2}")
        print(f"  - 行列式模{prime}：{det_mod}")
        print(f"  - 过滤前总解数：{total_before_filter}")
        print(f"  - 过滤掉的平凡解数：{trivial_count}")
        print(f"  - 剩余非平凡解数：{len(mj_solutions)}")
        print(f"  - 通过验证的有效解数：{len(valid_results)}")
        print(f"  - 总耗时：{elapsed:.2f}秒")
        print(f"  - 结果文件：{filename}")

        return {
            "prime": prime,
            "is_regular": (prime_type == "regular"),
            "max_pk_index": max_pk_index,
            "var_count": len(M_VALUES),
            "var_range": f"m{min(M_VALUES) if M_VALUES else '无'}~m{max(M_VALUES) if M_VALUES else '无'}",
            "det_mod": det_mod,
            "total_solutions_before_filter": total_before_filter,
            "trivial_solutions_filtered": trivial_count,
            "non_trivial_solutions": len(mj_solutions),
            "valid_solutions": len(valid_results),
            "elapsed_time": elapsed,
            "save_success": save_success,
            "filename": filename
        }

    except Exception as e:
        print(f"❌ 处理质数 {prime} 时出错：{str(e)}")
        import traceback
        traceback.print_exc()

        # 保存错误信息
        elapsed = time.time() - start_time
        save_prime_results(
            prime, pk_list, A_list, det_mod, mj_solutions, trivial_count,
            valid_results, elapsed, M_VALUES, max_pk_index
        )
        return None


# ===================== 批量处理主程序 =====================
def main():
    """主程序：批量处理1000以内的奇质数"""
    print("=" * 80)
    print("1000以内奇质数批量验证系统（v2.0 - 修复统计错误）")
    print("=" * 80)

    # 初始化目录
    init_save_directories()

    # 生成质数列表
    primes = generate_primes(MIN_PRIME, MAX_PRIME)
    if not primes:
        print("❌ 未找到符合条件的质数（≥5）")
        return

    # 批量处理质数
    batch_results = []
    total_start = time.time()
    total_primes = len(primes)

    print(f"\n开始批量处理 {total_primes} 个质数...")
    print(f"📁 正则素数结果保存到：{REGULAR_DIR}")
    print(f"📁 非正则素数结果保存到：{IRREGULAR_DIR}")

    # 可选：测试模式（仅处理前N个质数）
    test_mode = False
    test_count = 5
    primes_to_process = primes[:test_count] if test_mode else primes

    # 遍历处理每个质数
    for idx, prime in enumerate(primes_to_process):
        progress_idx = idx + 1
        result = process_single_prime(prime, progress_idx, len(primes_to_process))
        if result:
            batch_results.append(result)

        # 进度统计
        if progress_idx % 5 == 0 or progress_idx == len(primes_to_process):
            elapsed = time.time() - total_start
            avg_time = elapsed / progress_idx
            remaining_time = avg_time * (len(primes_to_process) - progress_idx)
            print(f"\n📊 批量进度：{progress_idx}/{len(primes_to_process)} 个质数已处理")
            print(f"   已耗时：{elapsed:.2f}秒 | 预计剩余：{remaining_time:.2f}秒")

            # 汇总统计
            if batch_results:
                total_before = sum([r['total_solutions_before_filter'] for r in batch_results])
                total_trivial = sum([r['trivial_solutions_filtered'] for r in batch_results])
                total_non_trivial = sum([r['non_trivial_solutions'] for r in batch_results])
                total_valid = sum([r['valid_solutions'] for r in batch_results])
                print(
                    f"   累计统计：过滤前={total_before} | 过滤平凡解={total_trivial} | 非平凡解={total_non_trivial} | 有效解={total_valid}")

    # 生成汇总报告
    total_elapsed = time.time() - total_start
    generate_batch_summary(batch_results, total_elapsed, primes_to_process)

    # 最终统计
    print_separator("批量处理完成")
    print(f"总耗时：{total_elapsed:.2f}秒")
    print(f"处理质数数：{len(batch_results)} / {len(primes_to_process)}")

    if batch_results:
        regular_count = len([r for r in batch_results if r['is_regular']])
        irregular_count = len([r for r in batch_results if not r['is_regular']])
        total_before = sum([r['total_solutions_before_filter'] for r in batch_results])
        total_trivial = sum([r['trivial_solutions_filtered'] for r in batch_results])
        total_non_trivial = sum([r['non_trivial_solutions'] for r in batch_results])
        total_valid = sum([r['valid_solutions'] for r in batch_results])

        print(f"\n📈 最终统计：")
        print(f"  - 正则素数：{regular_count} 个")
        print(f"  - 非正则素数：{irregular_count} 个")
        print(f"  - 过滤前总解数：{total_before}")
        print(f"  - 过滤掉的平凡解数：{total_trivial}")
        print(f"  - 剩余非平凡解数：{total_non_trivial}")
        print(f"  - 通过验证的有效解数：{total_valid}")

    print(f"\n📁 结果目录：{BASE_SAVE_DIR}")
    print(f"📄 汇总报告：{os.path.join(BASE_SAVE_DIR, 'batch_summary_1000.txt')}")


def generate_batch_summary(batch_results, total_elapsed, all_primes):
    """生成批量汇总报告，包含正确的过滤统计"""
    summary_filename = os.path.join(BASE_SAVE_DIR, "batch_summary_1000.txt")

    # 统计正则/非正则素数
    regular_primes = [r for r in batch_results if r['is_regular']]
    irregular_primes = [r for r in batch_results if not r['is_regular']]

    # 统计解的总数
    total_before_filter = sum([r['total_solutions_before_filter'] for r in batch_results]) if batch_results else 0
    total_trivial_filtered = sum([r['trivial_solutions_filtered'] for r in batch_results]) if batch_results else 0
    total_non_trivial = sum([r['non_trivial_solutions'] for r in batch_results]) if batch_results else 0
    total_valid = sum([r['valid_solutions'] for r in batch_results]) if batch_results else 0

    try:
        with open(summary_filename, "w", encoding="utf-8") as f:
            f.write("1000以内奇质数批量验证汇总报告（v2.0 - 修复统计错误）\n")
            f.write(f"生成时间：{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}\n")
            f.write(f"质数范围：{MIN_PRIME}~{MAX_PRIME} | 处理数量：{len(batch_results)}/{len(all_primes)}\n")
            f.write(f"总运行时间：{total_elapsed:.2f}秒\n")
            f.write("=" * 80 + "\n\n")

            # 整体统计
            f.write("【整体统计】\n")
            f.write(f"处理的质数总数：{len(batch_results)}\n")
            f.write(f"正则素数数量：{len(regular_primes)} / {len(batch_results)}\n")
            f.write(f"非正则素数数量：{len(irregular_primes)} / {len(batch_results)}\n")
            f.write(f"正则素数比例：{len(regular_primes) / len(batch_results) * 100:.1f}% (若分母不为0)\n")
            f.write(f"总计过滤前解数：{total_before_filter} 个\n")
            f.write(f"总计过滤掉的平凡解数：{total_trivial_filtered} 个\n")
            f.write(f"总计剩余非平凡解数：{total_non_trivial} 个\n")
            f.write(f"总计通过验证的有效解数：{total_valid} 个\n")
            if total_non_trivial > 0:
                f.write(f"有效解比例：{total_valid / total_non_trivial * 100:.1f}% (非平凡解中)\n")
            else:
                f.write(f"有效解比例：0% (无ufer非平凡解)\n")
            f.write("\n")

            # 质数列表
            f.write("【质数列表】\n")
            f.write(f"正则素数列表：{[r['prime'] for r in regular_primes]}\n")
            f.write(f"非正则素数列表：{[r['prime'] for r in irregular_primes]}\n")
            f.write("\n")

            # 详细统计
            f.write("【各质数详细统计】\n")
            f.write(
                f"{'质数':<6} {'类型':<8} {'变量数':<6} {'行列式':<8} {'过滤前':<8} {'平凡解':<8} {'非平凡解':<8} {'有效解':<8} {'耗时(秒)':<8}\n")
            f.write("-" * 80 + "\n")

            for res in batch_results:
                prime_type = "正则" if res['is_regular'] else "非正则"
                f.write(f"{res['prime']:<6} {prime_type:<8} {res['var_count']:<6} {res['det_mod']:<8} "
                        f"{res['total_solutions_before_filter']:<8} {res['trivial_solutions_filtered']:<8} "
                        f"{res['non_trivial_solutions']:<8} {res['valid_solutions']:<8} {res['elapsed_time']:<8.2f}\n")

            f.write("\n【统计修正说明】\n")
            f.write(f"1. 修复了行列式≠0时过滤前总解数统计错误（从2修正为1）\n")
            f.write(f"2. 增加了解的去重逻辑，避免重复统计相同解\n")
            f.write(f"3. 确保过滤前总解数 = 非平凡解数 + 过滤掉的平凡解数\n")
            f.write(f"4. 行列式≠0时，强制返回唯一平凡解（过滤前总解数=1）\n")

        print(f"\n📄 批量汇总报告已保存到：{summary_filename}")

    except Exception as e:
        print(f"❌ 生成汇总报告失败：{str(e)}")


if __name__ == "__main__":
    main()
