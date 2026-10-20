// SIMULATION / PROOF MODEL ONLY: the subset of the Xilinx DSP48E2 used by gen.py (UG579 semantics, [recalled]):
// A_INPUT = DIRECT, INMODE = 0 (multiplier takes A[26:0] and B), ALUMODE = 0000 (P = W + X + Y + Z + CIN),
// X = Y = M (OPMODE[3:0] = 0101), W = 0, Z = 0 / PCIN / C (OPMODE[6:4] = 000 / 001 / 011), CARRYIN = 0,
// AREG, BREG, CREG, MREG, PREG each 0 or 1. CE = 1, RST = 0 assumed. Yosys has no behavioural DSP48E2
// (cells_xtra.v declares it as a black box), so synthesis uses the black box and simulation uses this file.
module DSP48E2 #(
  parameter A_INPUT = "DIRECT", B_INPUT = "DIRECT", AMULTSEL = "A", BMULTSEL = "B", USE_MULT = "MULTIPLY",
  parameter USE_SIMD = "ONE48",
  parameter integer AREG = 1, BREG = 1, ACASCREG = 1, BCASCREG = 1, CREG = 1, MREG = 1, PREG = 1, DREG = 1, ADREG = 1,
  parameter integer OPMODEREG = 1, ALUMODEREG = 1, INMODEREG = 1, CARRYINREG = 1, CARRYINSELREG = 1
) (
  input CLK, input [29:0] A, input [17:0] B, input [47:0] C, input [26:0] D, input [47:0] PCIN,
  input [29:0] ACIN, input [17:0] BCIN, input [8:0] OPMODE, input [3:0] ALUMODE, input [4:0] INMODE,
  input [2:0] CARRYINSEL, input CARRYIN, input CARRYCASCIN, input MULTSIGNIN,
  input CEA1, CEA2, CEAD, CEALUMODE, CEB1, CEB2, CEC, CECARRYIN, CECTRL, CED, CEINMODE, CEM, CEP,
  input RSTA, RSTALLCARRYIN, RSTALUMODE, RSTB, RSTC, RSTCTRL, RSTD, RSTINMODE, RSTM, RSTP,
  output [47:0] P, output [47:0] PCOUT
);
  reg [29:0] a_r = 0; reg [17:0] b_r = 0; reg [47:0] c_r = 0; reg [44:0] m_r = 0; reg [47:0] p_r = 0;
  always @(posedge CLK) begin a_r <= A; b_r <= B; c_r <= C; end
  wire [29:0] a_m = AREG ? a_r : A;
  wire [17:0] b_m = BREG ? b_r : B;
  wire [47:0] c_m = CREG ? c_r : C;
  wire signed [44:0] m_c = $signed(a_m[26:0]) * $signed(b_m);
  always @(posedge CLK) m_r <= m_c;
  wire [44:0] m = MREG ? m_r : m_c;
  wire [47:0] m48 = {{3{m[44]}}, m};
  wire [47:0] z = (OPMODE[6:4] == 3'b001) ? PCIN : (OPMODE[6:4] == 3'b011) ? c_m : (OPMODE[6:4] == 3'b010) ? p_r : 48'd0;
  wire [47:0] alu = m48 + z;
  always @(posedge CLK) p_r <= alu;
  assign P = PREG ? p_r : alu;
  assign PCOUT = P;
endmodule
