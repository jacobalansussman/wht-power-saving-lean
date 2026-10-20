// as-4: SIMULATION MODEL ONLY of the subset of the Xilinx 7-series DSP48E1 used by gen2.py, written from the
// UG479 semantics [recalled]: A_INPUT = DIRECT, USE_DPORT = FALSE, INMODE = 0 (multiplier takes A2 = A[24:0] and B2),
// 25 x 18 signed multiplier, M sign-extended to 48 bits; ALUMODE = 0000 (P = Z + X + Y + CIN), X = Y = M
// (OPMODE[3:0] = 0101), Z = 0 / PCIN / P / C / PCIN>>>17 / P>>>17 (OPMODE[6:4] = 000 / 001 / 010 / 011 / 101 / 110),
// CARRYIN = 0; AREG, BREG, CREG, MREG, PREG each 0 or 1; CE = 1, RST = 0 assumed.
// Cross-checked against Xilinx's unisim DSP48E1.v (XilinxUnisimLibrary, Apache-2.0) by tb2.py UNISIM=1.
module DSP48E1 #(
  parameter A_INPUT = "DIRECT", B_INPUT = "DIRECT", USE_DPORT = "FALSE", USE_MULT = "MULTIPLY", USE_SIMD = "ONE48",
  parameter USE_PATTERN_DETECT = "NO_PATDET",
  parameter integer AREG = 1, BREG = 1, ACASCREG = 1, BCASCREG = 1, CREG = 1, MREG = 1, PREG = 1, DREG = 1, ADREG = 1,
  parameter integer OPMODEREG = 1, ALUMODEREG = 1, INMODEREG = 1, CARRYINREG = 1, CARRYINSELREG = 1
) (
  input CLK, input [29:0] A, input [17:0] B, input [47:0] C, input [24:0] D, input [47:0] PCIN,
  input [29:0] ACIN, input [17:0] BCIN, input [6:0] OPMODE, input [3:0] ALUMODE, input [4:0] INMODE,
  input [2:0] CARRYINSEL, input CARRYIN, input CARRYCASCIN, input MULTSIGNIN,
  input CEA1, CEA2, CEAD, CEALUMODE, CEB1, CEB2, CEC, CECARRYIN, CECTRL, CED, CEINMODE, CEM, CEP,
  input RSTA, RSTALLCARRYIN, RSTALUMODE, RSTB, RSTC, RSTCTRL, RSTD, RSTINMODE, RSTM, RSTP,
  output [47:0] P, output [47:0] PCOUT
);
  reg [29:0] a_r = 0; reg [17:0] b_r = 0; reg [47:0] c_r = 0; reg [42:0] m_r = 0; reg [47:0] p_r = 0;
  always @(posedge CLK) begin a_r <= A; b_r <= B; c_r <= C; end
  wire [29:0] a_m = AREG ? a_r : A;
  wire [17:0] b_m = BREG ? b_r : B;
  wire [47:0] c_m = CREG ? c_r : C;
  wire signed [42:0] m_c = $signed(a_m[24:0]) * $signed(b_m);
  always @(posedge CLK) m_r <= m_c;
  wire [42:0] m = MREG ? m_r : m_c;
  wire [47:0] m48 = {{5{m[42]}}, m};
  wire [47:0] z = (OPMODE[6:4] == 3'b001) ? PCIN : (OPMODE[6:4] == 3'b011) ? c_m : (OPMODE[6:4] == 3'b010) ? p_r :
                  (OPMODE[6:4] == 3'b101) ? {{17{PCIN[47]}}, PCIN[47:17]} : (OPMODE[6:4] == 3'b110) ? {{17{p_r[47]}}, p_r[47:17]} : 48'd0;
  wire [47:0] alu = m48 + z;
  always @(posedge CLK) p_r <= alu;
  assign P = PREG ? p_r : alu;
  assign PCOUT = P;
endmodule
