# gen.py -- emits the unrolled pass A / pass B FMOPA sequences (sizeless SVE types cannot live in arrays)
import sys
def pa(n, ar, ai, t0=0, t1=1):
    # left = F column n, right = data row; tiles t0 (re), t1 (im)
    if n in (0, 8):
        return [f"svmopa_za32_f32_m({t0},pg,pg,fr{n},{ar});", f"svmopa_za32_f32_m({t1},pg,pg,fr{n},{ai});"]
    if n < 8:
        return [f"svmopa_za32_f32_m({t0},pg,pg,fr{n},{ar});", f"svmops_za32_f32_m({t0},pg,pg,fi{n},{ai});",
                f"svmopa_za32_f32_m({t1},pg,pg,fr{n},{ai});", f"svmopa_za32_f32_m({t1},pg,pg,fi{n},{ar});"]
    m = 16 - n
    return [f"svmopa_za32_f32_m({t0},pg,pg,fr{m},{ar});", f"svmopa_za32_f32_m({t0},pg,pg,fi{m},{ai});",
            f"svmopa_za32_f32_m({t1},pg,pg,fr{m},{ai});", f"svmops_za32_f32_m({t1},pg,pg,fi{m},{ar});"]
def pb(n, yr, yi, t0=2, t1=3):
    # left = data column, right = F row n
    if n in (0, 8):
        return [f"svmopa_za32_f32_m({t0},pg,pg,{yr},fr{n});", f"svmopa_za32_f32_m({t1},pg,pg,{yi},fr{n});"]
    if n < 8:
        return [f"svmopa_za32_f32_m({t0},pg,pg,{yr},fr{n});", f"svmops_za32_f32_m({t0},pg,pg,{yi},fi{n});",
                f"svmopa_za32_f32_m({t1},pg,pg,{yr},fi{n});", f"svmopa_za32_f32_m({t1},pg,pg,{yi},fr{n});"]
    m = 16 - n
    return [f"svmopa_za32_f32_m({t0},pg,pg,{yr},fr{m});", f"svmopa_za32_f32_m({t0},pg,pg,{yi},fi{m});",
            f"svmops_za32_f32_m({t1},pg,pg,{yr},fi{m});", f"svmopa_za32_f32_m({t1},pg,pg,{yi},fr{m});"]
out = []
def emit(name, lines):
    out.append(f"#define {name} \\")
    for l in lines: out.append("  " + l + " \\")
    out.append("")
# pass A groups: rows n1 = 4g+m from x4 loads r4, i4
for g in range(4):
    L = [f"{{ svfloat32x4_t r4 = svld1_f32_x4(pc, xr + {64*g}); svfloat32x4_t i4 = svld1_f32_x4(pc, xi + {64*g});"]
    for m in range(4):
        L += pa(4*g+m, f"svget4(r4,{m})", f"svget4(i4,{m})")
    L.append("}")
    emit(f"PASSA_G{g}", L)
for g in range(4):
    L = [f"{{ svfloat32x4_t r4 = svread_ver_za32_f32_vg4(0, {4*g}); svfloat32x4_t i4 = svread_ver_za32_f32_vg4(1, {4*g});"]
    for m in range(4):
        L += pb(4*g+m, f"svget4(r4,{m})", f"svget4(i4,{m})")
    L.append("}")
    emit(f"PASSB_G{g}", L)
open(sys.argv[1], "w").write("\n".join(out) + "\n")
# --- v2: interleaved pass B(t) on tiles 2,3 (left operand = Y' columns loaded from scratch sr/si)
#         with pass A(t+1) on tiles 0,1 (right operand = rows of X(t+1))
out2 = []
def emit2(name, lines):
    out2.append(f"#define {name} \\")
    for l in lines: out2.append("  " + l + " \\")
    out2.append("")
for g in range(4):
    L = [f"{{ svfloat32x4_t r4 = svld1_f32_x4(pc, xr + {64*g}); svfloat32x4_t i4 = svld1_f32_x4(pc, xi + {64*g});",
         f"  svfloat32x4_t s4 = svld1_f32_x4(pc, sr + {64*g}); svfloat32x4_t t4 = svld1_f32_x4(pc, si + {64*g});"]
    for m in range(4):
        a = pa(4*g+m, f"svget4(r4,{m})", f"svget4(i4,{m})", 0, 1)
        b = pb(4*g+m, f"svget4(s4,{m})", f"svget4(t4,{m})", 2, 3)
        for x, y in zip(b, a): L += [x, y]
    L.append("}")
    emit2(f"AB_G{g}", L)
for g in range(4):
    L = [f"{{ svfloat32x4_t s4 = svld1_f32_x4(pc, sr + {64*g}); svfloat32x4_t t4 = svld1_f32_x4(pc, si + {64*g});"]
    for m in range(4):
        L += pb(4*g+m, f"svget4(s4,{m})", f"svget4(t4,{m})", 2, 3)
    L.append("}")
    emit2(f"BM_G{g}", L)
open(sys.argv[1].replace('.h', '2.h'), "w").write("\n".join(out2) + "\n")
# --- v3: transposed tiles so that every store is a whole ZA array vector (STR ZA):
#   pass A^T: tile(n2,k1) += X[n1,:] (x) F[:,n1]   (data left, F right = pb pattern) into tiles 0,1
#   pass B^T: tile(k2,k1) += F[:,n2] (x) Y'[:,n2]  (F left, data right = pa pattern) into tiles 2,3
out3 = []
def emit3(name, lines):
    out3.append(f"#define {name} \\")
    for l in lines: out3.append("  " + l + " \\")
    out3.append("")
for g in range(4):
    L = [f"{{ svfloat32x4_t r4 = svld1_f32_x4(pc, xr + {64*g}); svfloat32x4_t i4 = svld1_f32_x4(pc, xi + {64*g});"]
    for m in range(4):
        L += pb(4*g+m, f"svget4(r4,{m})", f"svget4(i4,{m})", 0, 1)
    L.append("}")
    emit3(f"T_A_G{g}", L)
for g in range(4):
    L = [f"{{ svfloat32x4_t r4 = svld1_f32_x4(pc, xr + {64*g}); svfloat32x4_t i4 = svld1_f32_x4(pc, xi + {64*g});",
         f"  svfloat32x4_t s4 = svld1_f32_x4(pc, sr + {64*g}); svfloat32x4_t t4 = svld1_f32_x4(pc, si + {64*g});"]
    for m in range(4):
        a = pb(4*g+m, f"svget4(r4,{m})", f"svget4(i4,{m})", 0, 1)
        b = pa(4*g+m, f"svget4(s4,{m})", f"svget4(t4,{m})", 2, 3)
        for x, y in zip(b, a): L += [x, y]
    L.append("}")
    emit3(f"T_AB_G{g}", L)
for g in range(4):
    L = [f"{{ svfloat32x4_t s4 = svld1_f32_x4(pc, sr + {64*g}); svfloat32x4_t t4 = svld1_f32_x4(pc, si + {64*g});"]
    for m in range(4):
        L += pa(4*g+m, f"svget4(s4,{m})", f"svget4(t4,{m})", 2, 3)
    L.append("}")
    emit3(f"T_B_G{g}", L)
open(sys.argv[1].replace('.h', '3.h'), "w").write("\n".join(out3) + "\n")
