# Plan: full proof of Proposition (Asymptotic Validity of Lower Limit under Sample Splitting)

File: `partialID_Overleaf/Methods/derivations.tex`, inside the existing `proof` environment, after the outline (outline kept unchanged, including its Step 5 wording).

Approved new notation:
- `Z_{k,n} = (\hat L_{k,eval} - L_k)/\sqrt{s_{k,n}}`, k = 1, 2
- `c_{1,n}(\omega) = \omega\sqrt{s_{1,n}/v_n(\omega)}`, `c_{2,n}(\omega) = (1-\omega)\sqrt{s_{2,n}/v_n(\omega)}`; `c_1^2 + c_2^2 = 1`, `c_k \in [0,1]`
- `T_n(\omega)` (known variances), `\hat T_n(\omega)` (estimated variances)
- `\varphi_{k,n}(t)` characteristic function of `Z_{k,n}`, `\varphi(t) = e^{-t^2/2}`
- `\Phi` standard normal CDF

Steps written out:
1. Condition on `D_tr(n)`; reduce to deterministic sequences.
2. Standardise: `T_n(\omega) = c_{1,n} Z_{1,n} + c_{2,n} Z_{2,n}`.
3. CF factorises (independence); add and subtract; modulus <= 1; `c_k \in [0,1]` gives a bound free of omega. State explicitly that it holds for the worst-case sequence.
4. Lévy's continuity theorem along any (incl. worst-case) sequence.
5. Variance ratio as a weighted average; uniform in omega via (A4).
6. Slutsky; event equivalence; continuity of Phi; (A2).
7-8. Already written in full in the outline; refer back.

No citations added; `% TODO: cite` for uniform-on-compacts convergence of CFs.
Verify: compile `compiled.tex`, check log.
