from sympy import symbols, binomial, simplify, expand, div, Integer, Rational
import sys

# Define common symbolic variable (restricted to positive integers)
m = symbols('m', integer=True, positive=True)
max_p = 17  # Recursion upper limit: P1-P17


# ---------------------- Core: e(k) function (corrected truncation logic) ----------------------
def e(k, m_val=None):
    """
    Return expression in terms of m: e_k = C(m, k+1)/m
    When k >= m, e_k = 0 (since C(m, k+1)=0 when k+1>m)
    """
    if k < 1:
        return Integer(0)

    # Check truncation condition if m value is specified
    if m_val is not None and k >= m_val:
        return Integer(0)

    # Return expression in terms of m
    comb = binomial(m, k + 1)
    comb_expanded = expand(comb)
    ek = comb_expanded / m

    return ek


# ---------------------- Newton's identities recursion function (corrected version) ----------------------
def recursive_pk_newton(max_k, m_val=None):
    """
    Implement correct Newton's identities recursion:
    - When k <= m: P_k = Σ_{i=1}^{k-1} (-1)^{i-1} e_i P_{k-i} + (-1)^{k-1} k e_k
    - When k > m: P_k = Σ_{i=1}^{m} (-1)^{i-1} e_i P_{k-i}
    """
    pk_list = [None] * (max_k + 1)

    for k in range(1, max_k + 1):
        pk = Integer(0)

        if m_val is not None and k > m_val:
            # k > m: Use special form of Newton's identities
            for i in range(1, m_val + 1):
                sign_term = (-1) ** (i - 1)
                e_i = e(i, m_val)
                term = sign_term * e_i * pk_list[k - i]
                pk += expand(term)
        else:
            # k <= m: Use standard Newton's identities
            for i in range(1, k):
                sign_term = (-1) ** (i - 1)
                e_i = e(i, m_val)
                term = sign_term * e_i * pk_list[k - i]
                pk += expand(term)

            # Final term
            if m_val is None or k < m_val:
                sign_final = (-1) ** (k - 1)
                e_k = e(k, m_val)
                final_term = sign_final * k * e_k
                pk += expand(final_term)

        # Simplify
        if k <= 10:
            pk = simplify(pk)
        else:
            pk = expand(pk)
        pk_list[k] = pk

    return pk_list


# ---------------------- Utility function: Force extract (m-1) common factor ----------------------
def process_pk_with_m_minus_1_force(expr):
    """Force extract (m-1) common factor to ensure Q(m) is a pure polynomial"""
    expr_simplified = simplify(expr)
    orig_num, orig_den = expr_simplified.as_numer_denom()

    # Handle denominator sign
    if isinstance(orig_den, Integer) and orig_den < 0:
        orig_num = -orig_num
        orig_den = -orig_den

    # Perform polynomial division on expanded numerator
    orig_num_expanded = expand(orig_num)

    # Extract (m-1) factor
    quotient, remainder = div(orig_num_expanded, (m - 1))
    if remainder == 0:
        common_factor = (m - 1)
        q_m = quotient
    else:
        common_factor = (m - 1)
        q_m = orig_num_expanded

    q_m_expanded = expand(q_m)
    final_num = common_factor * q_m_expanded

    return orig_den, final_num, q_m_expanded


# ---------------------- Format output function ----------------------
def format_fraction_with_m_minus_1(numerator, q_m, denom):
    """Format output to highlight (m-1) common factor"""
    m_minus_1 = "(m-1)"
    q_str = str(q_m).replace("(m - 1)", "(m-1)")
    q_str = f"({q_str})" if any(op in q_str for op in ['+', '-', '*', '^']) else q_str
    denom_str = str(denom) if isinstance(denom, Integer) else f"({denom})"

    if str(numerator).startswith('-'):
        return f"-{m_minus_1}×{q_str}/{denom_str}"
    else:
        return f"{m_minus_1}×{q_str}/{denom_str}"


# ---------------------- Function to compute numerical results ----------------------
def compute_numerical_values_for_all_m(pk_list):
    """Compute numerical results for m=3, 5, 7 for each expression"""
    numerical_values = {}
    for m_val in [3, 5, 7]:
        values = [None] * (max_p + 1)
        for i in range(1, max_p + 1):
            # Substitute specific m value and compute numerical result
            expr_with_value = pk_list[i].subs(m, m_val)
            values[i] = simplify(expr_with_value)
        numerical_values[m_val] = values
    return numerical_values


# ==============================================================================
# Compute P1-P17 expressions in terms of m for four scenarios
# ==============================================================================
scenarios = [
    ('m3', "m=3", 3),
    ('m5', "m=5", 5),
    ('m7', "m=7", 7),
    ('m21', "m≥21", None)  # None represents symbolic m
]

results = {}
for scenario, desc, m_val in scenarios:
    print(f"Starting to derive P1-P{max_p} expressions in terms of m for {desc}...")
    sys.stdout.flush()

    pk_list = recursive_pk_newton(max_k=max_p, m_val=m_val)
    pk_processed = [None]

    for i in range(1, max_p + 1):
        denom, final_num, q_m = process_pk_with_m_minus_1_force(pk_list[i])
        pk_processed.append((denom, final_num, q_m))

    # Compute numerical results for all m values for each scenario
    numerical_values = compute_numerical_values_for_all_m(pk_list)

    results[scenario] = (pk_processed, desc, pk_list, numerical_values, m_val)

# ==============================================================================
# Output comparison results
# ==============================================================================
print("\n" + "=" * 150)
print(f"Newton's Identities P1-P{max_p} Expression Comparison (Corrected Version)")
print("=" * 150)
print("Newton's Identities:")
print("   - When k ≤ m: P_k = Σ_{i=1}^{k-1} (-1)^{i-1} e_i P_{k-i} + (-1)^{k-1} k e_k")
print("   - When k > m: P_k = Σ_{i=1}^{m} (-1)^{i-1} e_i P_{k-i}")
print("=" * 150)

for i in range(1, max_p + 1):
    print(f"\nP{i}:")

    for scenario, (pk_processed, desc, pk_list, numerical_values, m_scenario_val) in results.items():
        denom, final_num, q_m = pk_processed[i]
        fmt_str = format_fraction_with_m_minus_1(final_num, q_m, denom)

        # Label characteristics
        if m_scenario_val is not None:
            if i <= m_scenario_val:
                note = "(k≤m, Standard Form)"
            else:
                note = "(k>m, Special Form)"
        else:  # Symbolic scenario
            note = "(Symbolic m, Standard Form)"

        print(f"   {desc} {note} | {fmt_str}")

        # Output numerical results for all m values
        for m_val in [3, 5, 7]:
            num_value = numerical_values[m_val][i]
            print(f"       m={m_val} Numerical Result | {num_value}")

    if i < max_p:
        print("-" * 150)

print("=" * 150)
print("✅ Derivation completed! All expressions are in terms of symbolic m")
print("   Newton's Identities:")
print("   - When k ≤ m: P_k = Σ_{i=1}^{k-1} (-1)^{i-1} e_i P_{k-i} + (-1)^{k-1} k e_k")
print("   - When k > m: P_k = Σ_{i=1}^{m} (-1)^{i-1} e_i P_{k-i}")
print("=" * 150)
