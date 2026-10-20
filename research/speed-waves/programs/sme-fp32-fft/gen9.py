# gen9.py -- v9: drain-overlapped schedule for the 256-point kernel (see fft256.c, Kernel v9).
import sys
exec(open('gen.py').read().split('out = []')[0])   # reuse pa(), pb()
def A_group(g):   # pass A^T of the next transform, tiles 0,1 (data left = X rows, F right)
    L = [f"{{ svfloat32x4_t r4_{g} = svld1_f32_x4(pc, xr + {64*g}); svfloat32x4_t i4_{g} = svld1_f32_x4(pc, xi + {64*g});"]
    body = []
    for m in range(4): body += pb(4*g+m, f"svget4(r4_{g},{m})", f"svget4(i4_{g},{m})", 0, 1)
    return L, body, ["}"]
def B_group(j):   # pass B^T of the current transform, tiles 2,3 (F left, data right = scratch slots 4j..4j+3 = n2 j+4m)
    L = [f"{{ svfloat32x4_t s4_{j} = svld1_f32_x4(pc, S + {64*j}); svfloat32x4_t t4_{j} = svld1_f32_x4(pc, S + 256 + {64*j});"]
    body = []
    for m in range(4): body += pa(j + 4*m, f"svget4(s4_{j},{m})", f"svget4(t4_{j},{m})", 2, 3)
    return L, body, ["}"]
def TW_STR(j):
    L = [f"TW({j})"]
    for m in range(4):
        n2 = j + 4*m
        L += [f"svstr_za({4*n2+0}, S + {16*(4*j+m)});", f"svstr_za({4*n2+1}, S + {256 + 16*(4*j+m)});"]
    return L
def OUT_pieces():
    P = []
    for t, base in ((2, "yo_r"), (3, "yo_i")):
        for c in range(4):
            P.append(f"svst1_f32_x4(pc, {base} + {64*c}, svread_hor_za32_f32_vg4({t}, {4*c}));")
    return P
def interleave(a, b):
    out = []; i = j = 0
    na, nb = len(a), len(b)
    while i < na or j < nb:
        # spread b evenly over a
        if j < nb and (i >= na or j * na <= i * nb): out.append(b[j]); j += 1
        else: out.append(a[i]); i += 1
    return out
def block(groups):   # combine several (open, body, close) groups into straight code with all loads first
    opens = sum((g[0] for g in groups), []); bodies = [g[1] for g in groups]; closes = sum((g[2] for g in groups), [])
    return opens, bodies, closes
lines = []
emit = lines.append
# step 1
emit("#define V9_STEP1 " + " ".join(TW_STR(0) + TW_STR(1)))
# step 2: B groups 0,1 interleaved with TW_STR(2), TW_STR(3)
o, bodies, c = block([B_group(0), B_group(1)])
seq = interleave(bodies[0] + bodies[1], TW_STR(2) + TW_STR(3))
emit("#define V9_STEP2 " + " ".join(o + seq + c))
emit("#define V9_STEP2_LAST " + " ".join(o + bodies[0] + bodies[1] + c + TW_STR(2) + TW_STR(3)))
# step 4: B groups 2,3 interleaved FMOPA by FMOPA with A groups 0,1
oB, bB, cB = block([B_group(2), B_group(3)])
oA, bA, cA = block([A_group(0), A_group(1)])
b1 = bB[0] + bB[1]; a1 = bA[0] + bA[1]
seq = interleave(b1, a1)
emit("#define V9_STEP4 " + " ".join(oB + oA + seq + cA + cB))
emit("#define V9_STEP4_LAST " + " ".join(oB + b1 + cB))
# step 5: A groups 2,3 interleaved with the output stores
oA2, bA2, cA2 = block([A_group(2), A_group(3)])
seq = interleave(bA2[0] + bA2[1], OUT_pieces())
emit("#define V9_STEP5 " + " ".join(oA2 + seq + cA2))
emit("#define V9_STEP5_LAST " + " ".join(OUT_pieces()))
open(sys.argv[1], "w").write("\n".join(lines) + "\n")
