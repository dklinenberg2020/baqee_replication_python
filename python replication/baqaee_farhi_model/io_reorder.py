"""
Python translation of replication/baqaee_farhi_model/IO_reorder.m
(the `initial_tariff_index == 1`, i.e. no-initial-tariff, data path used by
the paper's main driver script).

Loads the WIOD 2008 world input-output table (wiott2008.mat), the WIOD
socio-economic accounts (wiodsea2008.mat) and industry-level trade
elasticities (trade_elast_2008.mat), collapses the raw 31-industry WIOD
classification to the 30-sector classification used in the paper (merging
Textiles+Leather, Motor-vehicle-sales+Wholesale-trade, and sectors 8+9;
dropping "Other community/social/personal services" and "Private households
with employed persons" into "Public admin"), and reorganizes the world
economy into `len(keep_c)` individually tracked countries plus a single
aggregate "rest of the world" (ROW) country appended at the end.

All arrays are 0-indexed numpy arrays; comments give the 1-indexed MATLAB
line each block corresponds to. This is a direct, index-by-index port (not
a reformulation) because IO_reorder.m's correctness rests entirely on the
precise WIOD-2008 row/column layout, which is safest to preserve literally.
Validated element-wise against the original MATLAB code (GNU Octave 8.4.0)
for the full 41-country run -- see baqaee_farhi_model/README.md.
"""
import numpy as np
import scipy.io as sio


def _matlab_range(start, step, stop):
    """0-indexed numpy array equivalent to MATLAB's 1-indexed start:step:stop."""
    n = (stop - start) // step + 1
    return start - 1 + step * np.arange(n)


def io_reorder(keep_c, data_dir='.'):
    """
    Parameters
    ----------
    keep_c : 1-indexed (WIOD numbering) array-like of country numbers to
        track individually, ascending, e.g. np.r_[1:35, 36:42] for "all
        countries" (WIOD country 35 is always excluded from keep_c: it is
        handled as the aggregate rest-of-world country instead).

    Returns
    -------
    Omega, beta, alpha_VA, alpha, trade_elast, GDP_weights -- see
    IO_reorder.m's header comment / main_load_data_rev.m for definitions.
    """
    keep_c = np.asarray(keep_c, dtype=int)

    wiott2008 = sio.loadmat(f'{data_dir}/wiott2008.mat')['wiott2008']
    trade_elast_2008 = sio.loadmat(f'{data_dir}/trade_elast_2008.mat')['trade_elast_2008'].ravel()
    wiodsea2008 = sio.loadmat(f'{data_dir}/wiodsea2008.mat')['wiodsea2008']

    N = 31
    J = 41

    # --- Data cleaning (IO_reorder.m lines 30-60) --------------------------
    IO = wiott2008.copy()
    IO[IO < 0] = 0

    Nb_row = N + 4  # 35: row-block size (31 industries + 4 extra rows)
    Nb_col = N + 9  # 40: column-block size (31 industries + 9 extra columns)

    r_N = _matlab_range(N, Nb_row, (J - 1) * Nb_row + N)
    r_N3 = _matlab_range(N + 3, Nb_row, (J - 1) * Nb_row + N + 3)
    r_N4 = _matlab_range(N + 4, Nb_row, (J - 1) * Nb_row + N + 4)
    IO[r_N, :] = IO[r_N, :] + IO[r_N3, :] + IO[r_N4, :]

    c_N = _matlab_range(N, Nb_col, (J - 1) * Nb_col + N)
    c_N3 = _matlab_range(N + 3, Nb_col, (J - 1) * Nb_col + N + 3)
    c_N4 = _matlab_range(N + 4, Nb_col, (J - 1) * Nb_col + N + 4)
    IO[:, c_N] = IO[:, c_N] + IO[:, c_N3] + IO[:, c_N4]

    r4 = _matlab_range(4, Nb_row, (J - 1) * Nb_row + 4)
    r5 = _matlab_range(5, Nb_row, (J - 1) * Nb_row + 5)
    IO[r4, :] = IO[r4, :] + IO[r5, :]
    c4 = _matlab_range(4, Nb_col, (J - 1) * Nb_col + 4)
    c5 = _matlab_range(5, Nb_col, (J - 1) * Nb_col + 5)
    IO[:, c4] = IO[:, c4] + IO[:, c5]

    r19 = _matlab_range(19, Nb_row, (J - 1) * Nb_row + 19)
    r20 = _matlab_range(20, Nb_row, (J - 1) * Nb_row + 20)
    IO[r19, :] = IO[r19, :] + IO[r20, :]
    c19 = _matlab_range(19, Nb_col, (J - 1) * Nb_col + 19)
    c20 = _matlab_range(20, Nb_col, (J - 1) * Nb_col + 20)
    IO[:, c19] = IO[:, c19] + IO[:, c20]

    r9 = _matlab_range(9, Nb_row, (J - 1) * Nb_row + 9)
    r8 = _matlab_range(8, Nb_row, (J - 1) * Nb_row + 8)
    IO[r9, :] = IO[r9, :] + IO[r8, :]
    c9 = _matlab_range(9, Nb_col, (J - 1) * Nb_col + 9)
    c8 = _matlab_range(8, Nb_col, (J - 1) * Nb_col + 8)
    IO[:, c9] = IO[:, c9] + IO[:, c8]

    rows_to_drop = np.concatenate([r_N3, r_N4, r5, r20, r8])
    cols_to_drop = np.concatenate([c_N3, c_N4, c5, c20, c8])
    IO = np.delete(IO, rows_to_drop, axis=0)
    IO = np.delete(IO, cols_to_drop, axis=1)

    IO = IO.T

    N = 30

    # --- Final consumption columns and value added (lines 62-121) ---------
    Nb5 = N + 5  # 35: row-block size after transpose, with 5 extra rows now
    rN1 = _matlab_range(N + 1, Nb5, (J - 1) * Nb5 + N + 1)
    rN2 = _matlab_range(N + 2, Nb5, (J - 1) * Nb5 + N + 2)
    rN3b = _matlab_range(N + 3, Nb5, (J - 1) * Nb5 + N + 3)
    rN4b = _matlab_range(N + 4, Nb5, (J - 1) * Nb5 + N + 4)
    rN5b = _matlab_range(N + 5, Nb5, (J - 1) * Nb5 + N + 5)

    JN = J * N
    H_cons_h = (IO[np.ix_(rN1, np.arange(JN))].T + IO[np.ix_(rN4b, np.arange(JN))].T
                + IO[np.ix_(rN5b, np.arange(JN))].T)
    NP_cons_h = IO[np.ix_(rN2, np.arange(JN))].T
    Gov_cons_h = IO[np.ix_(rN3b, np.arange(JN))].T
    Tot_cons_h = H_cons_h + NP_cons_h + Gov_cons_h

    rows_drop2 = np.concatenate([rN5b, rN4b, rN3b, rN2, rN1])
    IO = np.delete(IO, rows_drop2, axis=0)

    IO[:JN, -1] = IO[:, :JN].sum(axis=0) + Tot_cons_h.sum(axis=1) - IO[:JN, :JN].sum(axis=1)

    VA_all = IO[:N * J, -1].copy()
    alpha = np.maximum(0.0, np.nan_to_num(VA_all))
    Omega = IO[:JN, :JN].copy()
    Omega = np.nan_to_num(Omega, posinf=1.0, neginf=1.0)

    beta = Tot_cons_h[:JN, :].copy()

    sea = np.nan_to_num(wiodsea2008)
    N = 31
    blk = (J - 1) * (N + 4)
    col0 = sea[0:blk, 0]
    col1 = sea[blk:2 * blk, 0]
    col2 = sea[2 * blk:3 * blk, 0]
    col3 = sea[3 * blk:4 * blk, 0]
    col4 = sea[4 * blk:5 * blk, 0]
    VA = np.column_stack([col0, col1 * col3, col1 * col2, col1 * col4])

    Nb4 = N + 4
    VA_new = np.zeros((J * Nb4, 4))
    VA_new[:34 * Nb4, :] = VA[:34 * Nb4, :]
    VA_new[34 * Nb4:35 * Nb4, :] = 250 * np.ones((Nb4, 4))
    VA_new[35 * Nb4:, :] = VA[34 * Nb4:, :]
    VA = VA_new

    v4 = _matlab_range(4, Nb4, (J - 1) * Nb4 + 4)
    v5 = _matlab_range(5, Nb4, (J - 1) * Nb4 + 5)
    VA[v4, :] = VA[v4, :] + VA[v5, :]
    v19 = _matlab_range(19, Nb4, (J - 1) * Nb4 + 19)
    v20 = _matlab_range(20, Nb4, (J - 1) * Nb4 + 20)
    VA[v19, :] = VA[v19, :] + VA[v20, :]
    v9 = _matlab_range(9, Nb4, (J - 1) * Nb4 + 9)
    v8 = _matlab_range(8, Nb4, (J - 1) * Nb4 + 8)
    VA[v9, :] = VA[v9, :] + VA[v8, :]
    vN3 = _matlab_range(N + 3, Nb4, (J - 1) * Nb4 + N + 3)
    vN4 = _matlab_range(N + 4, Nb4, (J - 1) * Nb4 + N + 4)
    VA = np.delete(VA, np.concatenate([vN3, vN4, v5, v20, v8]), axis=0)

    VA = np.maximum(VA, 0)
    alpha_VA = np.maximum(0.0, np.nan_to_num(VA))

    N = 30
    zero_va_rows = alpha_VA.sum(axis=1) == 0
    alpha[zero_va_rows] = 0

    # --- Reorganize into `keep_c` countries + aggregate ROW (lines 123-266) --
    M = 4
    J_new = len(keep_c)

    IO_new = np.zeros((J_new * N, J_new * N))
    beta_new = np.zeros((J_new * N, J_new))
    alpha_VA_new = np.zeros((J_new * N, M))
    alpha_new = np.zeros(J_new * N)

    for a, i in enumerate(keep_c):
        for b, j in enumerate(keep_c):
            IO_new[a * N:(a + 1) * N, b * N:(b + 1) * N] = IO[(i - 1) * N:i * N, (j - 1) * N:j * N]
            beta_new[a * N:(a + 1) * N, b] = beta[(i - 1) * N:i * N, j - 1]
        alpha_VA_new[a * N:(a + 1) * N, :] = alpha_VA[(i - 1) * N:i * N, :]
        alpha_new[a * N:(a + 1) * N] = alpha[(i - 1) * N:i * N]

    keep_mask = np.zeros(J + 1, dtype=bool)
    keep_mask[keep_c] = True  # index 1..J (0 unused) -> True if kept

    IO_hh = IO_new
    IO_ff_mod = IO[:J * N, :J * N].copy()
    for j in range(1, J + 1):
        if keep_mask[j]:
            IO_ff_mod[:, (j - 1) * N:j * N] = 0
            IO_ff_mod[(j - 1) * N:j * N, :] = 0
    IO_ff = np.zeros((N, N))
    for k in range(J):
        for l in range(J):
            IO_ff += IO_ff_mod[k * N:(k + 1) * N, l * N:(l + 1) * N]

    IO_mod_fh = IO[:J * N, :J * N].copy()
    IO_fh = np.zeros((J_new * N, N))
    IO_mod_hf = IO[:J * N, :J * N].copy()
    IO_hf = np.zeros((N, J_new * N))
    for j in range(1, J + 1):
        if keep_mask[j]:
            IO_mod_fh[(j - 1) * N:j * N, (j - 1) * N:j * N] = 0
            IO_mod_fh[:, (j - 1) * N:j * N] = 0
            IO_mod_hf[(j - 1) * N:j * N, (j - 1) * N:j * N] = 0
            IO_mod_hf[(j - 1) * N:j * N, :] = 0
        else:
            IO_mod_fh[(j - 1) * N:j * N, :] = 0
            IO_mod_hf[:, (j - 1) * N:j * N] = 0
    for k, j in enumerate(keep_c):
        for l in range(1, J + 1):
            IO_fh[k * N:(k + 1) * N, :] += IO_mod_fh[(j - 1) * N:j * N, (l - 1) * N:l * N]
            IO_hf[:, k * N:(k + 1) * N] += IO_mod_hf[(l - 1) * N:l * N, (j - 1) * N:j * N]

    beta_hh = beta_new
    beta_fh = np.zeros((N, J_new))
    beta_ff = np.zeros((N, 1))
    beta_hf = np.zeros((J_new * N, 1))
    for j in range(1, J + 1):
        for l, k in enumerate(keep_c):
            if not keep_mask[j]:
                beta_fh[:, l] += beta[(j - 1) * N:j * N, k - 1]
                beta_hf[l * N:(l + 1) * N, 0] += beta[(k - 1) * N:k * N, j - 1]
        if not keep_mask[j]:
            for k in range(1, J + 1):
                if not keep_mask[k]:
                    beta_ff[:, 0] += beta[(j - 1) * N:j * N, k - 1]

    alpha_h = alpha_new
    alpha_f = np.zeros(N)
    alpha_VA_h = alpha_VA_new
    alpha_VA_f = np.zeros((N, M))
    for k in range(1, J + 1):
        if not keep_mask[k]:
            alpha_f += alpha[(k - 1) * N:k * N]
            alpha_VA_f += alpha_VA[(k - 1) * N:k * N, :]

    IO_new = np.block([[IO_hh, IO_fh], [IO_hf, IO_ff]])
    alpha_new = np.concatenate([alpha_h, alpha_f])
    alpha_VA_new = np.vstack([alpha_VA_h, alpha_VA_f])
    beta_new = np.block([[beta_hh, beta_hf], [beta_fh, beta_ff]])

    alpha_VA = np.divide(alpha_VA_new, alpha_VA_new.sum(axis=1, keepdims=True),
                          out=np.zeros_like(alpha_VA_new), where=alpha_VA_new.sum(axis=1, keepdims=True) != 0)

    row_sums = IO_new.sum(axis=1, keepdims=True)
    Omega = np.divide(IO_new, row_sums, out=np.zeros_like(IO_new), where=row_sums != 0)

    col_sums = beta_new.sum(axis=0, keepdims=True)
    beta = np.divide(beta_new, col_sums, out=np.zeros_like(beta_new), where=col_sums != 0)

    alpha = alpha_new / (alpha_new + IO_new.sum(axis=1))
    alpha = np.nan_to_num(alpha)

    GDP_weights = beta_new.sum(axis=0)
    GDP_weights = GDP_weights / GDP_weights.sum()

    Omega = np.nan_to_num(Omega, posinf=1.0, neginf=1.0)
    beta = np.nan_to_num(beta, posinf=1.0, neginf=1.0)

    trade_elast = np.concatenate([trade_elast_2008[:8], trade_elast_2008[9:]])

    # WIOD SEA's medium-skill/high-skill columns are in the opposite order
    # from what main_load_data.py's standard-form construction expects --
    # this is a WIOD-specific data-semantics fix, not a generic model step,
    # so it belongs here (in the WIOD-specific loader) rather than in the
    # source-agnostic main_load_data(). A different source's own loader is
    # responsible for handing its own alpha_VA columns to main_load_data()
    # in whatever order that source's own factor categories are meant to be
    # used in -- there is no universal "correct" factor-column order beyond
    # whatever the source loader itself defines.
    alpha_VA = alpha_VA.copy()
    alpha_VA[:, [2, 3]] = alpha_VA[:, [3, 2]]

    return Omega, beta, alpha_VA, alpha, trade_elast, GDP_weights
