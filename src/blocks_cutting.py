"""
Cutting blocks.
"""

import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import time
matplotlib.rcParams.update({'font.size': 22})

#===================================================================================================

# Find minimal cuts to extract from block n x m x k (n >= m >= k)
# part of size t.

# Global mem (we can process blocks with n <= 10, m <= 10, k <= 10).
min_cuts_for_extract_part_n = 11
min_cuts_for_extract_part_shape = (min_cuts_for_extract_part_n,
                                   min_cuts_for_extract_part_n,
                                   min_cuts_for_extract_part_n,
                                   min_cuts_for_extract_part_n**3)
min_cuts_for_extract_part_mem = 0 - np.ones(min_cuts_for_extract_part_shape, dtype=int)

#---------------------------------------------------------------------------------------------------

def min_cuts_for_extract_part_1d(n, t):
    """
    Find minimum cuts count to extract from block n
    part of size t.

    Parameters
    ----------
    n : int
        Block side size.
    t : int
        Target part size.

    Returns
    -------
    int
        Minimum cuts count.
    """

    if t > n:
        return math.inf

    # Use memory.
    if min_cuts_for_extract_part_mem[n][1][1][t] >= 0:
        return min_cuts_for_extract_part_mem[n][1][1][t]

    if (t == 0) or (t == n):
        r = 0
    else:
        r = 1

    # Use memory.
    min_cuts_for_extract_part_mem[n][1][1][t] = r

    return r

#---------------------------------------------------------------------------------------------------

def min_cuts_for_extract_part_2d(n, m, t):
    """
    Find minimal cut count to extract from block n x m
    part of size t.

    Parameters
    ----------
    n : int
        Block side size.
    m : int
        Block side size.
    t : int
        Target part size.

    Returns
    -------
    int
        Minimal cuts count.
    """

    s = n * m

    if t > s:
        return math.inf

    # Normalize n >= m, t <= s / 2.
    if not (n >= m):
        n, m = m, n
    if t > s // 2:
        t = s - t

    # Use memory.
    if min_cuts_for_extract_part_mem[n][m][1][t] >= 0:
        return min_cuts_for_extract_part_mem[n][m][1][t]

    if (t == 0) or (t == s):
        r = 0
    elif n == 1:
        r = min_cuts_for_extract_part_1d(m, t)
    elif m == 1:
        r = min_cuts_for_extract_part_1d(n, t)
    else:
        r = math.inf
        for t1 in range(t):
            for n1 in range(1, n):
                r = min(r,
                        min_cuts_for_extract_part_2d(n1, m, t1) \
                        + min_cuts_for_extract_part_2d(n - n1, m, t - t1))
            for m1 in range(1, m):
                r = min(r,
                        min_cuts_for_extract_part_2d(n, m1, t1) \
                        + min_cuts_for_extract_part_2d(n, m - m1, t - t1))
        r = r + 1

    # Use memory.
    min_cuts_for_extract_part_mem[n][m][1][t] = r

    return r

#---------------------------------------------------------------------------------------------------

def min_cuts_for_extract_part_3d(n, m, k, t):
    """
    Find minimal cuts count to extract from block n x m x k
    part of size t.

    Parameters
    ----------
    n : int
        Block side size.
    m : int
        Block side size.
    k : int
        Block side size.
    t : int
        Target part size.

    Returns
    -------
    int
        Minimal cuts count.
    """

    s = n * m * k

    if t > s:
        return math.inf

    # Normalize n >= m >= k, t <= s / 2
    if not ((n >= m) and (m >= k)):
        a = [n, m, k]
        a.sort()
        [k, m, n] = a
    if t > s // 2:
        t = s - t

    # Use memory.
    if min_cuts_for_extract_part_mem[n][m][k][t] >= 0:
        return min_cuts_for_extract_part_mem[n][m][k][t]

    if (t == 0) or (t == s):
        r = 0
    elif n == 1:
        r = min_cuts_for_extract_part_2d(m, k, t)
    elif m == 1:
        r = min_cuts_for_extract_part_2d(n, k, t)
    elif k == 1:
        r = min_cuts_for_extract_part_2d(n, m, t)
    else:
        r = math.inf
        for t1 in range(t):
            for n1 in range(1, n):
                r = min(r,
                        min_cuts_for_extract_part_3d(n1, m, k, t1) \
                        + min_cuts_for_extract_part_3d(n - n1, m, k, t - t1))
            for m1 in range(1, m):
                r = min(r,
                        min_cuts_for_extract_part_3d(n, m1, k, t1) \
                        + min_cuts_for_extract_part_3d(n, m - m1, k, t - t1))
            for k1 in range(1, k):
                r = min(r,
                        min_cuts_for_extract_part_3d(n, m, k1, t1) \
                        + min_cuts_for_extract_part_3d(n, m, k - k1, t - t1))
        r = r + 1

    # Use memory.
    min_cuts_for_extract_part_mem[n][m][k][t] = r

    return r

#===================================================================================================

def test():
    """
    Tests.
    """

    # Minimal cuts count to extract from block n x 1 x 1 part of size t.
    # 1d
    assert math.isinf(min_cuts_for_extract_part_1d(5, 6))
    assert min_cuts_for_extract_part_1d(5, 0) == 0
    assert min_cuts_for_extract_part_1d(5, 5) == 0
    assert min_cuts_for_extract_part_1d(5, 3) == 1
    # 2d
    assert math.isinf(min_cuts_for_extract_part_2d(3, 3, 10))
    assert min_cuts_for_extract_part_2d(3, 3, 0) == 0
    assert min_cuts_for_extract_part_2d(3, 3, 9) == 0
    assert min_cuts_for_extract_part_2d(3, 3, 6) == 1
    assert min_cuts_for_extract_part_2d(3, 3, 4) == 2
    assert min_cuts_for_extract_part_2d(3, 3, 5) == 2
    # 3d
    assert math.isinf(min_cuts_for_extract_part_3d(3, 3, 3, 28))
    assert min_cuts_for_extract_part_3d(3, 3, 3, 0) == 0
    assert min_cuts_for_extract_part_3d(3, 3, 3, 27) == 0
    assert min_cuts_for_extract_part_3d(3, 3, 3, 9) == 1
    assert min_cuts_for_extract_part_3d(5, 5, 5, 50) == 1
    assert min_cuts_for_extract_part_3d(6, 6, 6, 53) == 5
    assert min_cuts_for_extract_part_3d(6, 6, 6, 59) == 5
    assert min_cuts_for_extract_part_3d(6, 6, 6, 89) == 5

#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    test()

    # Draw plot.
    plt.figure(figsize=(10, 6), dpi=100)
    n = 6
    xs = range(n**3 + 1)
    plt.plot(xs, [min_cuts_for_extract_part_3d(n, n, n, x) for x in xs], marker='o', linewidth=2.0)
    plt.xlabel(r'Размер выделяемой части $t$')
    plt.ylabel(r'Минимальное количество разрезов $P_{6 \times 6 \times 6}^{t}$')
    plt.grid(True)
    plt.show()

#===================================================================================================
