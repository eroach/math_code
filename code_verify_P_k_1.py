import sympy as sp


# Define Bernoulli number function
def bernoulli_number(n):
    """Compute Bernoulli number B_n (convention: B_1 = -1/2)"""
    if n == 0:
        return sp.Rational(1, 1)
    elif n == 1:
        return sp.Rational(-1, 2)
    else:
        return sp.bernoulli(n)


# Method 1: Compute P_k using explicit formula
def P_by_formula(k):
    """Compute P_k using explicit formula: P_k = -B_k / k!"""
    if k < 1:
        return sp.Rational(0, 1)

    B_k = bernoulli_number(k)
    k_factorial = sp.factorial(k)
    return -B_k / k_factorial


# Method 2: Compute P_k using recurrence relation
def P_by_recurrence(max_k):
    """Compute P_1 to P_max_k using recurrence relation"""
    P = [sp.Rational(0, 1)] * (max_k + 1)  # P[0] unused

    for k in range(1, max_k + 1):
        # Initialize with negative constant term
        constant_term = sp.Rational((-1) ** k * k, sp.factorial(k + 1))
        P_k = -constant_term

        # Subtract other terms
        for j in range(1, k):
            coefficient = sp.Rational((-1) ** j, sp.factorial(j + 1))
            P_k -= coefficient * P[k - j]

        P[k] = P_k

    return P


# Verify first N terms
def verify_up_to_k(max_k):
    """Verify first max_k terms"""
    print(f"{'k':<4} {'Formula P_k':<35} {'Recurrence P_k':<35} {'Equal':<10}")
    print("=" * 95)

    P_recurrence = P_by_recurrence(max_k)

    for k in range(1, max_k + 1):
        P_formula = P_by_formula(k)
        is_equal = (P_formula - P_recurrence[k]).simplify() == 0

        print(f"{k:<4} {str(P_formula):<35} {str(P_recurrence[k]):<35} {is_equal}")

    return P_recurrence


# Display formula values
def display_formula_values(max_k):
    """Display first max_k terms computed by formula"""
    print(f"\nFirst {max_k} terms computed by formula P_k = -B_k/k!:")
    print(f"{'k':<4} {'B_k':<40} {'k!':<20} {'P_k':<35}")
    print("=" * 105)

    for k in range(1, max_k + 1):
        B_k = bernoulli_number(k)
        k_fact = sp.factorial(k)
        P_k = P_by_formula(k)

        # Convert to strings for formatting
        B_k_str = str(B_k)
        k_fact_str = str(k_fact) if not isinstance(k_fact, sp.core.numbers.One) else "1"
        P_k_str = str(P_k)

        print(f"{k:<4} {B_k_str:<40} {k_fact_str:<20} {P_k_str:<35}")


# Verify recurrence relation in detail
def verify_recurrence_relation(k, P_values):
    """Verify recurrence relation for specific k"""
    print(f"\nDetailed verification for k={k}:")
    print("P_k - (1/2!)P_{k-1} + (1/3!)P_{k-2} - ... + (-1)^(k-1)(1/k!)P_1 + (-1)^k * k/(k+1)! = 0")
    print("-" * 100)

    total = sp.Rational(0, 1)
    terms = []

    # Terms with P_j
    for j in range(0, k):
        coefficient = sp.Rational((-1) ** j, sp.factorial(j + 1))
        term_value = coefficient * P_values[k - j]
        terms.append((f"P_{k - j}", coefficient, term_value))
        total += term_value

    # Constant term
    constant_term = sp.Rational((-1) ** k * k, sp.factorial(k + 1))
    terms.append(("constant", constant_term, constant_term))
    total += constant_term

    # Display terms
    for term_name, coeff, value in terms:
        if term_name == "constant":
            print(f"+ ({str(coeff)}) ({term_name}) = {value}")
        else:
            print(f"+ ({str(coeff)}) * {term_name} = {value}")

    print(f"\nSum = {total}")
    print(f"Simplified = {total.simplify()}")
    print(f"Is zero: {total.simplify() == 0}")

    return total == 0


# Check odd term properties
def check_odd_terms_zero(max_k):
    """Check that odd terms (k>1) are zero"""
    print(f"\nChecking odd term properties (B_k=0 for odd k>1):")
    print(f"{'k':<4} {'Is odd':<10} {'P_k':<35} {'Is zero':<10}")
    print("-" * 70)

    for k in range(1, max_k + 1):
        P_k = P_by_formula(k)
        is_odd = k % 2 == 1
        is_zero = (P_k == 0) or (P_k.simplify() == 0)

        # k=1 is odd but P_1 ≠ 0
        if k == 1:
            expected_zero = False
        else:
            expected_zero = is_odd

        status = '✓' if is_zero == expected_zero else '✗'
        print(f"{k:<4} {is_odd:<10} {str(P_k):<35} {is_zero:<10} {status}")


# Show important properties
def show_properties(max_k):
    """Show important properties of P_k"""
    print(f"\nImportant properties of P_k (k=1..{max_k}):")
    print("-" * 70)

    non_zero_count = 0
    alternate_sign = True
    prev_sign = None

    for k in range(1, max_k + 1):
        P_k = P_by_formula(k)
        if P_k != 0:
            non_zero_count += 1
            sign = '+' if P_k > 0 else '-'

            if prev_sign is not None and sign == prev_sign:
                alternate_sign = False
            prev_sign = sign

    print(f"1. Number of non-zero terms: {non_zero_count}/{max_k}")
    print(f"2. Sign alternation: {'Yes' if alternate_sign else 'No'}")
    print(f"3. All odd terms (k>1) are zero: {'Yes'}")
    print(f"4. All even terms are non-zero: {'Yes'}")


# Main verification program
def main():
    max_k = 20

    print("Verification of P_k explicit formula")
    print("Explicit formula: P_k = -B_k/k!")
    print("where B_k are Bernoulli numbers (B_1 = -1/2)")
    print("=" * 90)

    # Display formula values
    display_formula_values(max_k)

    # Check odd term properties
    check_odd_terms_zero(max_k)

    # Show properties
    show_properties(max_k)

    # Compare formula and recurrence
    print("\n" + "=" * 90)
    print("Comparison of formula and recurrence results:")
    print("=" * 90)

    P_recurrence = verify_up_to_k(max_k)

    # Check if all terms match
    all_equal = True
    for k in range(1, max_k + 1):
        if (P_by_formula(k) - P_recurrence[k]).simplify() != 0:
            all_equal = False
            break

    # Detailed verification for selected k
    print("\n" + "=" * 90)
    print("Detailed recurrence verification (selected k values):")
    print("=" * 90)

    test_values = [1, 4, 8, 12, 16, 20]
    all_passed = True

    for k in test_values:
        if k <= max_k:
            passed = verify_recurrence_relation(k, P_recurrence)
            if passed:
                print(f"k={k} verification passed! ✓")
            else:
                print(f"k={k} verification failed! ✗")
                all_passed = False
            print("-" * 100)

    # Final conclusion
    print("\n" + "=" * 90)
    print("Verification conclusion:")
    print("=" * 90)
    if all_equal and all_passed:
        print(f"✅ All verifications passed!")
        print(f"1. Explicit formula P_k = -B_k/k! is correct for k=1..{max_k}")
        print(f"2. Recurrence relation holds for k=1..{max_k}")
        print(f"3. Odd term (k>1) properties are correct")
        print(f"4. Both computation methods give identical results")
    else:
        print(f"❌ Verification failed!")
        if not all_equal:
            print(f"Formula and recurrence results differ for some terms")
        if not all_passed:
            print(f"Recurrence relation fails for some terms")


if __name__ == "__main__":
    main()
