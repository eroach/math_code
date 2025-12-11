import sympy


def create_matrix(m):
    """Create a (m-1)-order matrix with base numbers 2, 3, ..., m"""
    size = m - 1
    M = sympy.zeros(size, size)
    for i in range(size):  # Row index i, base number a = i+2
        a = i + 2
        for j in range(size):  # Column index j, corresponding column number j+1
            if j == 0:
                M[i, j] = 1
            elif j == 1:
                M[i, j] = a
            else:
                exp1 = 2 * (j - 1)
                exp2 = 2 * (j - 1) + 1
                M[i, j] = a ** exp1 + a ** exp2
    return M


def det_and_factorization(m):
    """Calculate determinant and its prime factorization, then verify divisibility by prime (2m-1)"""
    M = create_matrix(m)
    det = M.det()
    # Prime factorization
    factors = sympy.factorint(det)
    # Format as product string
    factor_str = ' * '.join(f'{p}^{e}' if e > 1 else str(p) for p, e in factors.items())

    # Step 1: Calculate n = 2m - 1 and check if it is a prime number
    n = 2 * m - 1
    is_prime = sympy.isprime(n)

    # Step 2: Verify if n divides the determinant
    if is_prime:
        is_divisible = det % n == 0
        divisible_str = f"{n} (prime) divides the determinant" if is_divisible else f"{n} (prime) does NOT divide the determinant"
    else:
        is_divisible = None
        divisible_str = f"{n} is NOT a prime number, skip divisibility check"

    return det, factor_str, n, is_prime, is_divisible, divisible_str


# Calculate results for m from 3 to 10
results = {}
for m in range(3, 11):
    det, factor_str, n, is_prime, is_divisible, divisible_str = det_and_factorization(m)
    results[m] = (det, factor_str, n, is_prime, is_divisible)
    # Print detailed results
    print(f"m = {m}:")
    print(f"  Matrix order = {m - 1}")
    print(f"  Determinant value = {det}")
    print(f"  Prime factorization = {factor_str}")
    print(f"  Check number: n = 2m-1 = {n}")
    print(f"  Divisibility check: {divisible_str}")
    print("-" * 60)  # Separator for readability
